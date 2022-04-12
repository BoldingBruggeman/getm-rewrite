import operator
from typing import Union, Optional, Type, List
import itertools
import logging
import datetime
import timeit

import numpy
import numpy.typing

from .constants import *
from . import _pygetm
from . import core
from . import domain
from . import parallel
from . import output
from . import pyfabm
import pygetm.mixing
import pygetm.density
import pygetm.airsea
import pygetm.radiation

BAROTROPIC = BAROTROPIC_2D = 1
BAROTROPIC_3D = 2
FROZEN_DENSITY = 3
BAROCLINIC = 4

# Must match number in .../operators/advection.F90 - do you use them Jorn?
HSIMT = 1
MUSCL = 2
P2_PDM = 3
SPLMAX13 = 4
SUPERBEE = 5
UPSTREAM = 6 

class Boundaries:
    __slots__ = '_tracer', '_type', 'values'
    def __init__(self, tracer: 'Tracer'):
        self._tracer = tracer
        self._type = ZERO_GRADIENT
        self.values: Optional[core.Array] = None

    @property
    def type(self) -> int:
        return self._type

    @type.setter
    def type(self, value: int):
        self._type = value
        self.values = None if self._type == ZERO_GRADIENT else self._tracer.grid.array(name='%s_bdy' % self._tracer.name, z=CENTERS, on_boundary=True, attrs={'_3d_only': True})

    def update(self):
        self._tracer.update_boundary(self._type, self.values)

class Tracer(core.Array):
    __slots__ = 'source', 'source_scale', 'vertical_velocity', 'boundaries', 'river_values', 'river_follow', 'rivers'
    def __init__(self, grid: domain.Grid, data: Optional[numpy.ndarray]=None, source: Optional[core.Array]=None, source_scale: float=1., vertical_velocity: Optional[core.Array]=None, rivers_follow_target_cell: bool=False, **kwargs):
        super().__init__(grid=grid, shape=grid.hn.all_values.shape, **kwargs)
        if data is None:
            data = numpy.full_like(grid.hn.all_values, numpy.nan)
        self.wrap_ndarray(data)
        self.register()
        assert source is None or (source.grid is self.grid and source.z == CENTERS)
        assert vertical_velocity is None or (vertical_velocity.grid is self.grid and vertical_velocity.z == CENTERS)
        self.source = source
        self.source_scale = source_scale
        self.vertical_velocity = vertical_velocity
        self.boundaries = Boundaries(self)
        self.river_values = numpy.zeros((len(grid.domain.rivers),))
        self.river_follow = numpy.full((len(grid.domain.rivers),), rivers_follow_target_cell, dtype=bool)
        self.rivers = {}
        for iriver, river in enumerate(grid.domain.rivers.values()):
            river_tracer = domain.RiverTracer(grid, river.name, self.name, self.river_values[..., iriver], self.river_follow[..., iriver], units=self.units, attrs={'_3d_only': True})
            river._tracers[self.name] = river_tracer
            self.rivers[river.name] = river_tracer

class Simulation(_pygetm.Simulation):
    _momentum_arrays = 'U', 'V', 'fU', 'fV', 'advU', 'advV', 'u1', 'v1', 'bdyu', 'bdyv', 'uk', 'vk', 'ru', 'rru', 'rv', 'rrv', 'pk', 'qk', 'ww', 'advpk', 'advqk', 'Ui', 'Vi', 'SS', 'fpk', 'fqk'
    _pressure_arrays = 'dpdx', 'dpdy', 'idpdx', 'idpdy'
    _sealevel_arrays = 'zbdy',
    _time_arrays = 'timestep', 'macrotimestep', 'split_factor', 'timedelta', 'time', 'istep', 'report'
    _all_fortran_arrays = tuple(['_%s' % name for name in _momentum_arrays + _pressure_arrays + _sealevel_arrays]) + ('uadv', 'vadv', 'uua', 'uva', 'vua', 'vva', 'uua3d', 'uva3d', 'vua3d', 'vva3d')
    __slots__ = _all_fortran_arrays + ('output_manager', 'input_manager', 'fabm_model', '_fabm_interior_diagnostic_arrays', '_fabm_horizontal_diagnostic_arrays', 'fabm_sources_interior', 'fabm_sources_surface', 'fabm_sources_bottom', 'fabm_vertical_velocity', 'fabm_conserved_quantity_totals', '_yearday', 'tracers', 'tracer_totals', 'logger', 'airsea', 'turbulence', 'density', 'temp', 'salt', 'pres', 'rad', 'par', 'par0', 'rho', 'sst', 'sss', 'SS', 'NN', 'u_taus', 'u_taub', 'z0s', 'z0b', 'fwf', 'vertical_diffusion', 'tracer_advection', '_cum_river_height_increase', '_start_time', '_profile', 'radiation') + _time_arrays

    def __init__(self, dom: domain.Domain, runtype: int, advection_scheme: int=4, apply_bottom_friction: bool=True, fabm: Union[bool, str, None]=None, gotm: Union[str, None]=None,
        turbulence: Optional[pygetm.mixing.Turbulence]=None, airsea: Optional[pygetm.airsea.Fluxes]=None, density: Optional[pygetm.density.Density]=None,
        logger: Optional[logging.Logger]=None, log_level: int=logging.INFO, A=0.7, g1=1., g2=15.):

        self.logger = dom.root_logger
        self.output_manager = output.OutputManager(rank=dom.tiling.rank, logger=self.logger.getChild('output_manager'))
        self.input_manager = dom.input_manager
        dom.field_manager = self.output_manager

        self.input_manager.set_logger(self.logger.getChild('input_manager'))

        # Disable bottom friction if physical bottom roughness is 0 everywhere
        if apply_bottom_friction and (numpy.ma.array(dom.z0b_min, mask=dom.mask==0) == 0.).any():
            self.logger.warning('Disabling bottom friction because bottom roughness is 0 in one or more points.')
            apply_bottom_friction = False

        assert not dom.initialized
        _pygetm.Simulation.__init__(self, dom, runtype, apply_bottom_friction)
        self.logger.info('Maximum dt = %.3f s' % dom.maxdt)
        dom.T.hn.fabm_standard_name = 'cell_thickness'
        if dom.T.lon is not None:
            dom.T.lon.fabm_standard_name = 'longitude'
        if dom.T.lat is not None:
            dom.T.lat.fabm_standard_name = 'latitude'

        array_args = {
            'uk': dict(units='m s-1', long_name='velocity in Eastward direction', fill_value=FILL_VALUE),
            'vk': dict(units='m s-1', long_name='velocity in Northward direction', fill_value=FILL_VALUE),
            'ww': dict(units='m s-1', long_name='vertical velocity', fill_value=FILL_VALUE),
            'SS': dict(units='s-2', long_name='shear frequency squared', fill_value=FILL_VALUE)
        }

        for name in Simulation._momentum_arrays:
            setattr(self, '_%s' % name, self.wrap(core.Array(name=name, **array_args.get(name, {})), name.encode('ascii'), source=1))
        for name in Simulation._pressure_arrays:
            setattr(self, '_%s' % name, self.wrap(core.Array(name=name, **array_args.get(name, {})), name.encode('ascii'), source=2))
        for name in Simulation._sealevel_arrays:
            setattr(self, '_%s' % name, self.wrap(core.Array(name=name, **array_args.get(name, {})), name.encode('ascii'), source=3))

        if runtype > BAROTROPIC_2D:
            self.uk.all_values.fill(0.)
            self.vk.all_values.fill(0.)
            self.ww.all_values.fill(0.)

        self.update_depth()
        self._cum_river_height_increase = numpy.zeros((len(self.domain.rivers),))

        self.airsea = airsea or pygetm.airsea.FluxesFromMeteo()
        assert isinstance(self.airsea, pygetm.airsea.Fluxes)
        self.airsea.initialize(self.domain)

        self.fwf = dom.T.array(name='fwf', units='m s-1', long_name='freshwater flux', fill_value=FILL_VALUE)
        self.fwf.fill(0.)
        self.uadv = _pygetm.Advection(dom.U, scheme=advection_scheme)
        self.vadv = _pygetm.Advection(dom.V, scheme=advection_scheme)

        self.uua = dom.UU.array(fill=numpy.nan)
        self.uva = dom.UV.array(fill=numpy.nan)
        self.vua = dom.VU.array(fill=numpy.nan)
        self.vva = dom.VV.array(fill=numpy.nan)

        self.uua3d = dom.UU.array(fill=numpy.nan, z=CENTERS)
        self.uva3d = dom.UV.array(fill=numpy.nan, z=CENTERS)
        self.vua3d = dom.VU.array(fill=numpy.nan, z=CENTERS)
        self.vva3d = dom.VV.array(fill=numpy.nan, z=CENTERS)

        self.tracers: List[Tracer] = []
        self.tracer_totals = []

        self.fabm_model = None

        if runtype > BAROTROPIC_2D:
            self.tracer_advection = _pygetm.Advection(dom.T, scheme=advection_scheme)

            # Turbulence and associated fields
            self.turbulence = turbulence or pygetm.mixing.GOTM(self.domain, nml_path=gotm)
            self.NN = dom.T.array(fill=0., z=INTERFACES, name='NN', units='s-2', long_name='buoyancy frequency squared', fill_value=FILL_VALUE)
            self.u_taus = dom.T.array(fill=0., name='u_taus', units='m s-1', long_name='shear velocity (surface)', fill_value=FILL_VALUE)
            self.u_taub = dom.T.array(fill=0., name='u_taub', units='m s-1', long_name='shear velocity (bottom)', fill_value=FILL_VALUE)
            self.z0s = dom.T.array(fill=0.1, name='z0s', units='m', long_name='hydrodynamic roughness (surface)', fill_value=FILL_VALUE)
            self.z0b = dom.T.array(fill=0.1, name='z0b', units='m', long_name='hydrodynamic roughness (bottom)', fill_value=FILL_VALUE)

            self.vertical_diffusion = _pygetm.VerticalDiffusion(dom.T, cnpar=1.)

            if fabm:
                def fabm_variable_to_array(variable, send_data: bool=False, **kwargs):
                    ar = core.Array(name=variable.output_name, units=variable.units, long_name=variable.long_path, fill_value=variable.missing_value, dtype=self.fabm_model.fabm.dtype, grid=self.domain.T, **kwargs)
                    if send_data:
                        ar.wrap_ndarray(variable.data)
                    ar.register()
                    return ar

                shape = self.domain.T.hn.all_values.shape
                self.fabm_model = pyfabm.Model(fabm if isinstance(fabm, str) else 'fabm.yaml', shape=shape, libname='fabm_c',
                    start=(0, self.domain.halo, self.domain.halo), stop=(shape[0], shape[1] - self.domain.halo, shape[2] - self.domain.halo))
                self.fabm_sources_interior = numpy.zeros_like(self.fabm_model.interior_state)
                self.fabm_sources_surface = numpy.zeros_like(self.fabm_model.surface_state)
                self.fabm_sources_bottom = numpy.zeros_like(self.fabm_model.bottom_state)
                self.fabm_vertical_velocity = numpy.zeros_like(self.fabm_model.interior_state)
                for i, variable in enumerate(self.fabm_model.interior_state_variables):
                    ar_w = core.Array(grid=self.domain.T)
                    ar_w.wrap_ndarray(self.fabm_vertical_velocity[i, ...])
                    self.create_tracer(data=variable.data, vertical_velocity=ar_w, name=variable.output_name, units=variable.units, long_name=variable.long_path, 
                        fill_value=variable.missing_value, rivers_follow_target_cell=variable.no_river_dilution)
                for variable in itertools.chain(self.fabm_model.surface_state_variables, self.fabm_model.bottom_state_variables):
                    ar = fabm_variable_to_array(variable, send_data=True)
                self._fabm_interior_diagnostic_arrays = [fabm_variable_to_array(variable, shape=self.domain.T.hn.shape) for variable in self.fabm_model.interior_diagnostic_variables]
                self._fabm_horizontal_diagnostic_arrays = [fabm_variable_to_array(variable, shape=self.domain.T.H.shape) for variable in self.fabm_model.horizontal_diagnostic_variables]
                self.fabm_model.link_mask(self.domain.T.mask.all_values)
                self.fabm_model.link_cell_thickness(self.domain.T.hn.all_values)

                self.fabm_conserved_quantity_totals = numpy.empty((len(self.fabm_model.conserved_quantities),) + self.domain.T.H.all_values.shape, dtype=self.fabm_sources_interior.dtype)
                for i, variable in enumerate(self.fabm_model.conserved_quantities):
                    ar = core.Array(name=variable.output_name, units=variable.units, long_name=variable.long_name, fill_value=variable.missing_value, dtype=self.fabm_model.fabm.dtype, grid=self.domain.T)
                    ar.wrap_ndarray(self.fabm_conserved_quantity_totals[i, ...])
                    self.tracer_totals.append(ar)

            self.pres = dom.T.array(z=CENTERS, name='pres', units='dbar', long_name='pressure', fabm_standard_name='pressure', fill_value=FILL_VALUE)

        self.sst = dom.T.array(name='sst', units='degrees_Celsius', long_name='sea surface temperature', fill_value=FILL_VALUE)

        if runtype == BAROCLINIC:
            self.temp = self.create_tracer(name='temp', units='degrees_Celsius', long_name='conservative temperature', fabm_standard_name='temperature', fill_value=FILL_VALUE, source=dom.T.array(fill=0., z=CENTERS, units='W m-2'), source_scale=1. / (RHO0 * CP), rivers_follow_target_cell=True)
            self.salt = self.create_tracer(name='salt', units='-', long_name='absolute salinity', fabm_standard_name='practical_salinity', fill_value=FILL_VALUE, source=dom.T.array(fill=0., z=CENTERS))
            self.pres.saved = True
            self.temp.fill(5.)
            self.salt.fill(35.)
            self.rho = dom.T.array(z=CENTERS, name='rho', units='kg m-3', long_name='density', fabm_standard_name='density')
            self.tracer_totals.append(self.salt)
            self.sss = self.salt.isel(-1)

            self.radiation = pygetm.radiation.TwoBand(dom.T)
            self.density = density or pygetm.density.Density()
        else:
            self.sss = None

    def __getitem__(self, key: str) -> core.Array:
        return self.output_manager.fields[key]

    def get_fabm_dependency(self, name: str):
        variable = self.fabm_model.dependencies.find(name)
        if len(variable.shape) == 0:
            return variable
        arr = self.domain.T.array(name=variable.output_name, units=variable.units, long_name=variable.long_path, z=len(variable.shape) == 3)
        variable.link(arr.all_values)
        return arr

    def create_tracer(self, **kwargs):
        """Add a tracer that will be subject to advection and diffusion.
        The optional source array must after multiplication with source_scale have tracer units per time, multiplied by layer thickness."""
        tracer = Tracer(grid=self.domain.T, **kwargs)
        self.tracers.append(tracer)
        return tracer

    def start(self, time: datetime.datetime, timestep: float, split_factor: int=1, report: int=10, save: bool=True, profile: str=False):
        """This should be called after the output configuration is complete (because we need toknow when variables need to be saved),
        and after the FABM model has been provided with all dependencies"""
        self.logger.info('Starting simulation at %s' % time)
        self.timestep = timestep
        self.split_factor = split_factor
        self.macrotimestep = self.timestep * self.split_factor
        self.timedelta = datetime.timedelta(seconds=timestep)
        self.time = time
        self.istep = 0
        self.report = report

        if self.runtype > BAROTROPIC_2D:
            # ensure ho and hn are up to date and identical
            self.start_3d()

        if self.fabm_model:
            # Tell FABM which diagnostics are saved. FABM will allocate and manage memory only for those that are.
            # This MUST be done before calling self.fabm_model.start
            for variable, ar in zip(itertools.chain(self.fabm_model.interior_diagnostic_variables, self.fabm_model.horizontal_diagnostic_variables), itertools.chain(self._fabm_interior_diagnostic_arrays, self._fabm_horizontal_diagnostic_arrays)):
                variable.save = ar.saved

            # Transfer GETM fields with a standard name to FABM
            for field in self.domain.field_manager.fields.values():
                if field.fabm_standard_name:
                    try:
                        variable = self.fabm_model.dependencies.find(field.fabm_standard_name)
                    except KeyError:
                        continue
                    field.saved = True
                    variable.link(field.all_values)

            try:
                self._yearday = self.fabm_model.dependencies.find('number_of_days_since_start_of_the_year')
                self._yearday.value = (self.time - datetime.datetime(self.time.year, 1, 1)).total_seconds() / 86400.
            except KeyError:
                self._yearday = None

            # Start FABM. This verifies whether all dependencies are fulfilled and freezes the set of dsiagsntoics that will be saved.
            assert self.fabm_model.start(), 'FABM failed to start. Likely its configuration is incomplete.'

            # Fill GETM placeholder arrays for all FABM diagnostics that will be computed/saved.
            for variable, ar in zip(itertools.chain(self.fabm_model.interior_diagnostic_variables, self.fabm_model.horizontal_diagnostic_variables), itertools.chain(self._fabm_interior_diagnostic_arrays, self._fabm_horizontal_diagnostic_arrays)):
                if ar.saved:
                    ar.wrap_ndarray(variable.data)

            # Apply mask to state variables
            for variable in itertools.chain(self.fabm_model.interior_state_variables, self.fabm_model.surface_state_variables, self.fabm_model.bottom_state_variables):
                variable.value[..., self.domain.T.mask.all_values == 0] = variable.missing_value

        if self.runtype == BAROCLINIC:
            # ensure density is in sync with initial T & S
            if self.rho.saved:
                self.density.get_density(self.salt, self.temp, p=self.pres, out=self.rho)
            self.density.get_potential_temperature(self.salt.isel(-1), self.temp.isel(-1), self.sst)

        # Update all input fields
        self.domain.input_manager.update(time, include_3d=True)

        # Ensure heat and momentum fluxes have sensible values
        self.airsea(self.time, self.sst, self.sss, calculate_heat_flux=self.runtype == BAROCLINIC)
        if self.runtype == BAROCLINIC:
            assert self.airsea.shf.require_set(self.logger)
            self.radiation(self.airsea.swr)

        # Update elevation at the open boundaries
        self.update_sealevel_boundaries(self.timestep)

        if self.fabm_model:
            self.update_fabm_sources()

        if self.domain.open_boundaries:
            for tracer in self.tracers:
                tracer.boundaries.update()

        self.output_manager.start(save=save)

        self._start_time = timeit.default_timer()

        self._profile = None
        if profile:
            import cProfile
            pr = cProfile.Profile()
            self._profile = (profile, pr)
            pr.enable()

    def advance(self):
        self.time += self.timedelta
        self.istep += 1
        macro_active = self.istep % self.split_factor == 0

        # Update all inputs
        self.domain.input_manager.update(self.time, include_3d=macro_active)

        # Update air-sea fluxes of heat and momentum (T grid for all, U and V grid for x and y stresses respectively)
        self.airsea(self.time, self.sst, self.sss, calculate_heat_flux=macro_active)

        # Update elevation at the open boundaries
        self.update_sealevel_boundaries(self.timestep)

        # Calculate the surface pressure gradient in the U and V points.
        # This requires elevation and surface pressure (both on T grid) to be valid in the halos
        self.airsea.sp.update_halos()
        self.update_surface_pressure_gradient(self.domain.T.z, self.airsea.sp)

        # Update momentum using surface stresses and pressure gradients. Inputs and outputs on U and V grids.
        self.uv_momentum_2d(self.timestep, self.airsea.taux_U, self.airsea.tauy_V, self.dpdx, self.dpdy)

        # Update halos of U and V as the sea level update below requires valid transports within the halos
        self.U.update_halos()
        self.V.update_halos()

        self.update_freshwater_fluxes(self.timestep)

        # Update sea level on T grid, and from that calculate sea level and water depth on all grids
        self.update_sealevel(self.timestep, self.U, self.V, self.fwf)
        self.update_depth()

        if self.report != 0 and self.istep % self.report == 0:
            self.logger.info(self.time)

        if self.runtype > BAROTROPIC_2D and macro_active:
            # Depth-integrated transports have been summed over all microtimesteps. Now average them.
            self.Ui.all_values *= 1. / self.split_factor
            self.Vi.all_values *= 1. / self.split_factor

            # Use previous source terms for biogeochemistry to update tracers
            # This should be done before the tracer concentrations change due to transport or rivers,
            # as the source terms are only valid for the current tracer concentrations.
            if self.fabm_model:
                self.update_fabm(self.macrotimestep)

            # Update of layer thicknesses and tracer concentrations to account for river inflow over the past macrotimestep.
            self.add_rivers_3d()

            # Update 3D elevations and layer thicknesses. New elevation on T grid will match elevation at end of 2D timestep,
            # thicknesses on T grid will match. Elevation and thicknesses on U/V grids will be 1/2 macrotimestep behind. Old
            # elevations and thicknesses will be one macrotimestep behind new elevations aand thicknesses.
            self.start_3d()

            # Update presssure gradient for start of the 3D time step. JB: presumably airsea.sp needs to be at start of the macrotimestep - not currently the case if we update meteo on 2D timestep!
            self.update_surface_pressure_gradient(self.domain.T.zio, self.airsea.sp)

            # JB: at what time should airsea.taux_U and airsea.tauy_V formally be defined? Start of the macrotimestep like dpdx/dpdy? TODO
            self.uvw_momentum_3d(self.macrotimestep, self.airsea.taux_U, self.airsea.tauy_V, self.dpdx, self.dpdy, self.idpdx, self.idpdy, self.turbulence.num)

            if self.runtype == BAROCLINIC:
                # Buoyancy frequency and turbulence (T grid - interfaces)
                self.density.get_buoyancy_frequency(self.salt, self.temp, p=self.pres, out=self.NN)
                self.u_taus.all_values[...] = (self.airsea.taux.all_values**2 + self.airsea.tauy.all_values**2)**0.25 / numpy.sqrt(RHO0)
                self.turbulence(self.macrotimestep, self.u_taus, self.u_taub, self.z0s, self.z0b, self.NN, self.SS)

                # Temperature and salinity (T grid - centers)
                self.radiation(self.airsea.swr)
                self.temp.source.all_values[...] = self.radiation.swr_abs.all_values
                self.temp.source.all_values[-1, ...] += self.airsea.shf.all_values

                # Transport of passive tracers (including biogeochemical ones)
                w_if = None
                for tracer in self.tracers:
                    if tracer.source is not None:
                        tracer.source.all_values *= self.macrotimestep * tracer.source_scale
                    w = self.ww
                    if tracer.vertical_velocity is not None:
                        w_if = tracer.vertical_velocity.interp(w_if or self.domain.T.array(fill=0., z=INTERFACES))
                        w_if.all_values += self.ww.all_values
                        w = w_if
                    self.tracer_advection.apply_3d(self.uk, self.vk, self.ww, self.macrotimestep, tracer, w_var=w)
                    self.vertical_diffusion(self.turbulence.nuh, self.macrotimestep, tracer, ea4=tracer.source)
                    if self.domain.open_boundaries:
                        tracer.boundaries.update()

                # Update density to keep it in sync with T and S.
                # Unlike buoyancy frequency, density does not influence hydrodynamics directly.
                # Therefore, it is computed only if needed for output (or as input for biogeochemistry)
                if self.rho.saved:
                    self.density.get_density(self.salt, self.temp, p=self.pres, out=self.rho)

                # From conservative temperature to in-situ sea surface temperature, needed to compute heat/momentum fluxes at the surface
                self.density.get_potential_temperature(self.salt.isel(-1), self.temp.isel(-1), out=self.sst)

                # Update source terms of biogeochemistry, using the new tracer concentrations
                if self.fabm_model:
                    self.update_fabm_sources()

            self.report_domain_integrals()

            # Reset depth-integrated transports that will be incremented over subsequent 3D timestep.
            self.Ui.all_values.fill(0.)
            self.Vi.all_values.fill(0.)

        self.output_manager.save()
        return macro_active

    def finish(self):
        if self._profile:
            import pstats
            name, pr = self._profile
            pr.disable()
            with open('%s-%03i.prof' % (name, self.domain.tiling.rank), 'w') as f:
                ps = pstats.Stats(pr, stream=f).sort_stats(pstats.SortKey.TIME)
                ps.print_stats()
        self.logger.info('Time spent in main loop: %.3f s' % (timeit.default_timer() - self._start_time,))
        self.output_manager.close()

    def update_freshwater_fluxes(self, timestep: float):
        height_increase_rate = self.domain.rivers.flow * self.domain.rivers.iarea
        self._cum_river_height_increase += height_increase_rate * timestep
        self.fwf.all_values[self.domain.rivers.j, self.domain.rivers.i] = height_increase_rate

    def add_rivers_3d(self):
        """Update layer thicknesses and tracer concentrations to account for river inflow."""
        h = self.domain.T.hn.all_values[:, self.domain.rivers.j, self.domain.rivers.i]
        river_active = numpy.full(h.shape, True)
        # JB TODO: customize river_active by flagging layers where the river does not go with False
        river_depth = (h * river_active).copy(order='F').sum(axis=0)   # ensure pairwise summation also for more than 1 river - needed for reproducibility for different divisions!
        withdrawal = self._cum_river_height_increase < 0.
        h_increase = numpy.where(river_active, h * self._cum_river_height_increase / river_depth, 0.)
        for tracer in self.tracers:
            river_follow = numpy.logical_or(tracer.river_follow, withdrawal)
            if not river_follow.all():
                tracer_old = tracer.all_values[:, self.domain.rivers.j, self.domain.rivers.i]
                river_values = numpy.where(river_follow, tracer_old, tracer.river_values)
                tracer_new = (tracer_old * river_depth + river_values * self._cum_river_height_increase) / (river_depth + self._cum_river_height_increase)
                tracer.all_values[:, self.domain.rivers.j, self.domain.rivers.i] = numpy.where(river_active, tracer_new, tracer_old)
        self.domain.T.hn.all_values[:, self.domain.rivers.j, self.domain.rivers.i] = h + h_increase
        self.domain.T.zin.all_values[self.domain.rivers.j, self.domain.rivers.i] += self._cum_river_height_increase
        self._cum_river_height_increase.fill(0.)

    def update_fabm_sources(self):
        """Update FABM sources, vertical velocities, and diagnostics. This does not update the state variables themselves; that is done by update_fabm"""
        if self._yearday:
            self._yearday.value = (self.time - datetime.datetime(self.time.year, 1, 1)).total_seconds() / 86400.
        self.fabm_model.get_sources(out=(self.fabm_sources_interior, self.fabm_sources_surface, self.fabm_sources_bottom))
        self.fabm_model.get_vertical_movement(self.fabm_vertical_velocity)

    def update_fabm(self, timestep: float):
        """Time-integrate source terms of all FABM state variables (3D pelagic tracers as well as bottom- and surface-attached variables)"""
        self.fabm_sources_interior *= timestep
        self.fabm_sources_surface *= timestep
        self.fabm_sources_bottom *= timestep
        self.fabm_model.interior_state += self.fabm_sources_interior
        self.fabm_model.surface_state += self.fabm_sources_surface
        self.fabm_model.bottom_state += self.fabm_sources_bottom

    def report_domain_integrals(self):
        """Write totals of selected variables over the global domain (those in list self.tracer_totals) to the log."""
        total_volume = (self.domain.T.D * self.domain.T.area).global_sum(where=self.domain.T.mask != 0)
        if total_volume is not None:
            self.logger.info('Integrals over global domain:')
            self.logger.info('  volume: %.15e m3' % total_volume)
        if self.fabm_model:
            self.fabm_model.get_conserved_quantities(out=self.fabm_conserved_quantity_totals)
        for var in self.tracer_totals:
            total = var * var.grid.area
            if total.ndim == 3:
                total.all_values *= var.grid.hn.all_values
            total = total.global_sum(where=var.grid.mask != 0)
            if total is not None:
                self.logger.info('  %s: %.15e %s m3 (per volume: %s %s)' % (var.name, total, var.units, total / total_volume, var.units))

    def uv_momentum_2d(self, timestep: float, tausx: core.Array, tausy: core.Array, dpdx: core.Array, dpdy: core.Array):
        # compute velocities at time=n-1/2
        self.u1.all_values[:, :] = self.U.all_values / self.U.grid.D.all_values
        self.v1.all_values[:, :] = self.V.all_values / self.V.grid.D.all_values

        itimestep = 1. / timestep

        # Advect U using velocities interpolated to its own advection grids
        self.U.interp(self.uua)
        self.V.interp(self.uva)
        self.uua.all_values /= self.domain.UU.D.all_values
        self.uva.all_values /= self.domain.UV.D.all_values
        self.uadv(self.uua, self.uva, timestep, self.u1)
        self.advU.all_values[...] = (self.u1.all_values * self.uadv.D - self.U.all_values) * itimestep

        # Advect V using velocities interpolated to its own advection grids
        self.U.interp(self.vua)
        self.V.interp(self.vva)
        self.vua.all_values /= self.domain.VU.D.all_values
        self.vva.all_values /= self.domain.VV.D.all_values
        self.vadv(self.vua, self.vva, timestep, self.v1)
        self.advV.all_values[...] = (self.v1.all_values * self.vadv.D - self.V.all_values) * itimestep

        # Restore velocity at time=n-1/2
        self.u1.all_values[:, :] = self.U.all_values / self.U.grid.D.all_values
        self.v1.all_values[:, :] = self.V.all_values / self.V.grid.D.all_values

        _pygetm.Simulation.uv_momentum_2d(self, timestep, tausx, tausy, dpdx, dpdy)

    def uvw_momentum_3d(self, timestep: float, tausx: core.Array, tausy: core.Array, dpdx: core.Array, dpdy: core.Array, idpdx: core.Array, idpdy: core.Array, viscosity: core.Array):
        # Do the halo exchange for viscosity, as this needs to be interpolated to the U and V grids. For that, information from the halos is used.
        viscosity.update_halos()

        _pygetm.Simulation.uvw_momentum_3d(self, timestep, tausx, tausy, dpdx, dpdy, idpdx, idpdy, viscosity.interp(self.domain.U), viscosity.interp(self.domain.V))

        itimestep = 1. / timestep

        # Update the halos so that transports (and more importantly, the velocities derived subsequently) are valid there.
        # This is needed to advect quantities defined on the T grid, as this requires velocities at the boundaries of every T cell
        # of the subdomain interior; this includes cells at the very Western and Southern boundary, which for U and V grids lie within the halo
        # Moreover, the velocities will also be interpolated to the advection grids for momentum (UU, UV, VU, VV), which again requires halos values.
        # The vertical velocity is already computed at all points surrounding U and V points, so no further halo exchange of w is needed
        # to support interpolation to U and V grids later on.
        self.pk.update_halos()
        self.qk.update_halos()

        # Compute 3D velocities (m s-1) from 3D transports (m2 s-1) by dividing by layer heights
        numpy.divide(self.pk.all_values, self.U.grid.hn.all_values, where=self.pk.grid.mask.all_values != 0, out=self.uk.all_values)
        numpy.divide(self.qk.all_values, self.V.grid.hn.all_values, where=self.qk.grid.mask.all_values != 0, out=self.vk.all_values)

        # Use updated velocities (uk, vk) to compute shear frequency (SS)
        self.update_shear_frequency(viscosity)

        # Interpolate 3D velocities to advection grids.
        # This needs to be done before uk/vk are changed by the advection operator (apply_3d).
        self.uk.interp(self.uua3d)
        self.vk.interp(self.uva3d)
        self.uk.interp(self.vua3d)
        self.vk.interp(self.vva3d)

        # Advect 3D u velocity from time=n-1/2 to n+1/2 using velocities interpolated to its own advection grids
        # JB the alternative would be to interpolate transports and then divide by (colocated) layer heights, like we do for 2D
        self.uadv.apply_3d(self.uua3d, self.uva3d, self.ww.interp(self.uk.grid), timestep, self.uk, new_h=True, skip_initial_halo_exchange=True)
        self.advpk.all_values[...] = (self.uk.all_values * self.uadv.h - self.pk.all_values) * itimestep

        # Advect 3D v velocity from time=n-1/2 to n+1/2 using velocities interpolated to its own advection grids
        # JB the alternative would be to interpolate transports and then divide by (colocated) layer heights, like we do for 2D
        self.vadv.apply_3d(self.vua3d, self.vva3d, self.ww.interp(self.vk.grid), timestep, self.vk, new_h=True, skip_initial_halo_exchange=True)
        self.advqk.all_values[...] = (self.vk.all_values * self.vadv.h - self.qk.all_values) * itimestep

        # Restore velocity at time=n-1/2
        numpy.divide(self.pk.all_values, self.U.grid.hn.all_values, where=self.pk.grid.mask.all_values != 0, out=self.uk.all_values)
        numpy.divide(self.qk.all_values, self.V.grid.hn.all_values, where=self.qk.grid.mask.all_values != 0, out=self.vk.all_values)

    def start_3d(self):
        """Update surface elevations and layer thicknesses for the 3D time step, starting from elevations at the end of the most recent 2D time step.
        Note: this uses sea level on T grid as computed by the 2D time step. This has to be up to date in the halos too!
        """
        # Store current elevations as previous elevations (on the 3D time step)
        self.domain.T.zio.all_values[...] = self.domain.T.zin.all_values
        self.domain.U.zio.all_values[...] = self.domain.U.zin.all_values
        self.domain.V.zio.all_values[...] = self.domain.V.zin.all_values
        self.domain.X.zio.all_values[...] = self.domain.X.zin.all_values

        # Synchronize new elevations on the 3D time step to those of the 2D time step that has just completed.
        self.domain.T.zin.all_values[...] = self.domain.T.z.all_values

        # Compute elevations on U, V, X grids.
        # Note that this must be at time=n+1/2, whereas elevations on the T grid is now at time=n+1.
        zi_T_half = 0.5 * (self.domain.T.zio + self.domain.T.zin)
        zi_T_half.interp(self.domain.U.zin)
        zi_T_half.interp(self.domain.V.zin)
        zi_T_half.interp(self.domain.X.zin)

        # Clip newly inferred elevations to ensure the minimum depth is respected.
        self.domain.U.zin.all_values.clip(min=-self.domain.U.H.all_values + self.domain.Dmin, out=self.domain.U.zin.all_values)
        self.domain.V.zin.all_values.clip(min=-self.domain.V.H.all_values + self.domain.Dmin, out=self.domain.V.zin.all_values)
        self.domain.X.zin.all_values.clip(min=-self.domain.X.H.all_values + self.domain.Dmin, out=self.domain.X.zin.all_values)

        # Halo exchange for elevation on U, V grids, needed because the very last points in the halos
        # (x=-1 for U, y=-1 for V) are not valid after interpolating from the T grid above.
        # These elevations are needed to later compute velocities from transports
        # (by dividing by layer thicknesses, which are computed from elevation)
        # These velocities will be advected, and therefore need to be valid througout the halos.
        # We do not need to halo-exchange elevation on the X grid, since that need to be be valid
        # at the innermost halo point only, which is ensured by T.zin exchange.
        self.domain.U.zin.update_halos(parallel.LEFT_RIGHT)
        self.domain.V.zin.update_halos(parallel.TOP_BOTTOM)

        # Update layer thicknesses (hn) using bathymetry H and new elevations zin (on the 3D timestep)
        # This routine also sets ho to the previous value of hn
        # All points (interior and halo) are processed, so ho and hn will be valid in halos after this completes,
        # provided zin was valid in the halo beforehand.
        self.domain.do_vertical()

        # Update thicknesses on advection grids. These must be at time=n+1/2
        # That's already the case for the X grid, but for the T grid (now at t=n+1) we explicitly compute thicknesses at time=n+1/2.
        # Note that UU.hn and VV.hn will miss the x=-1 and y=-1 strips, respectively (the last strip of values within their halos);
        # fortunately these values are not needed for advection.
        h_half = self.domain.T.ho.all_values + self.domain.T.hn.all_values
        self.domain.UU.hn.all_values[:, :, :-1] = h_half[:, :, 1:]
        self.domain.VV.hn.all_values[:, :-1, :] = h_half[:, 1:, :]
        self.domain.UV.hn.all_values[:, :, :] = self.domain.VU.hn.all_values[:, :, :] = self.domain.X.hn.all_values[:, 1:, 1:]

        if self.pres.saved:
            # Update pressure (dbar) at layer centers, assuming it is equal to depth in m
            _pygetm.thickness2center_depth(self.domain.T.mask, self.domain.T.hn, self.pres)

        if self.domain.zc_bdy.saved:
            # Update vertical coordinate at open boundary, used to interpolate inputs on z grid to dynamic model depths
            self.domain.zc_bdy.all_values[...] = self.domain.T.zc.all_values[:, self.domain.bdy_j, self.domain.bdy_i].T

    def update_depth(self):
        """Use surface elevation on T grid to update elevations on U,V,X grids and subsequently update total water depth D on all grids."""
        # Halo exchange for sea level on T grid
        self.domain.T.z.update_halos()

        # Compute sea level on U, V, X grids.
        # Note that this must be at time=n+1/2, whereas sea level on T grid is now at time=n+1.
        z_T_half = 0.5 * (self.domain.T.zo + self.domain.T.z)
        z_T_half.interp(self.domain.U.z)
        z_T_half.interp(self.domain.V.z)
        z_T_half.interp(self.domain.X.z)
        self.domain.U.z.all_values.clip(min=-self.domain.U.H.all_values + self.domain.Dmin, out=self.domain.U.z.all_values)
        self.domain.V.z.all_values.clip(min=-self.domain.V.H.all_values + self.domain.Dmin, out=self.domain.V.z.all_values)
        self.domain.X.z.all_values.clip(min=-self.domain.X.H.all_values + self.domain.Dmin, out=self.domain.X.z.all_values)
        z_T_half.all_values.clip(min=-self.domain.T.H.all_values + self.domain.Dmin, out=z_T_half.all_values)

        # Halo exchange for elevation on U, V grids, needed because the very last points in the halos
        # (x=-1 for U, y=-1 for V) are not valid after interpolating from the T grid above.
        # These elevations are needed to later compute velocities from transports
        # (by dividing by layer thicknesses, which are computed from elevation)
        # These velocities will be advected, and therefore need to be valid througout the halos.
        # We do not need to halo-exchange elevation on the X grid, since that need to be be valid
        # at the innermost halo point only, which is ensured by T.zin exchange.
        self.domain.U.z.update_halos(parallel.LEFT_RIGHT)
        self.domain.V.z.update_halos(parallel.TOP_BOTTOM)

        # Update total water depth D on T, U, V, X grids
        # This also processes the halos; no further halo exchange needed.
        self.domain.update_depths()

        # Update column depth on advection grids. These must be at time=n+1/2.
        # That's already the case for the X grid, but for the T grid we explicitly compute and use D at time=n+1/2.
        D_T_half = self.domain.T.H.all_values + z_T_half.all_values
        self.domain.UU.D.all_values[:, :-1] = D_T_half[:, 1:]
        self.domain.VV.D.all_values[:-1, :] = D_T_half[1:, :]
        self.domain.UV.D.all_values[:, :] = self.domain.VU.D.all_values[:, :] = self.domain.X.D.all_values[1:, 1:]

        self.u1.all_values[:, :] = self.U.all_values / self.U.grid.D.all_values
        self.v1.all_values[:, :] = self.V.all_values / self.V.grid.D.all_values

    @property
    def Ekin(self, rho0: float=RHO0):
        dom = self.domain
        U = self.U.interp(dom.T)
        V = self.V.interp(dom.T)
        vel2_D2 = U**2 + V**2
        return 0.5 * rho0 * dom.T.area * vel2_D2 / dom.T.D

for membername in Simulation._all_fortran_arrays:
    setattr(Simulation, membername[1:], property(operator.attrgetter(membername)))
