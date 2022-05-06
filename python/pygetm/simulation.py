import operator
from typing import Union, Optional, Type, List
import itertools
import logging
import datetime
import timeit

import numpy
import numpy.typing
import cftime
import enum

from .constants import *
from . import _pygetm
from . import core
from . import domain
from . import parallel
from . import output
from . import operators
from . import pyfabm
import pygetm.mixing
import pygetm.density
import pygetm.airsea
import pygetm.radiation

BAROTROPIC = BAROTROPIC_2D = 1
BAROTROPIC_3D = 2
FROZEN_DENSITY = 3
BAROCLINIC = 4

class InternalPressure(enum.IntEnum):
    BLUMBERG_MELLOR = 1
#    BLUMBERG_MELLOR_LIN=2
#    Z_INTERPOL=3
#    SONG_WRIGHT=4
#    CHU_FAN=5
#    SHCHEPETKIN_MCWILLIAMS=6
#    STELLING_VANKESTER=7


class OpenBoundaries:
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
    __slots__ = 'source', 'source_scale', 'vertical_velocity', 'open_boundaries', 'river_values', 'river_follow', 'rivers'
    def __init__(self, grid: domain.Grid, data: Optional[numpy.ndarray]=None, source: Optional[core.Array]=None, source_scale: float=1., vertical_velocity: Optional[core.Array]=None, rivers_follow_target_cell: bool=False, **kwargs):
        kwargs.setdefault('attrs', {}).update(_part_of_state=True, _3d_only=True)
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
        self.open_boundaries = OpenBoundaries(self)
        self.river_values = numpy.zeros((len(grid.domain.rivers),))
        self.river_follow = numpy.full((len(grid.domain.rivers),), rivers_follow_target_cell, dtype=bool)
        self.rivers = {}
        for iriver, river in enumerate(grid.domain.rivers.values()):
            river_tracer = domain.RiverTracer(grid, river.name, self.name, self.river_values[..., iriver], self.river_follow[..., iriver], units=self.units, attrs={'_3d_only': True})
            river._tracers[self.name] = river_tracer
            self.rivers[river.name] = river_tracer

class Simulation(_pygetm.Simulation):
    _momentum_arrays = 'U', 'V', 'fU', 'fV', 'advU', 'advV', 'diffu1', 'diffv1', 'u1', 'v1', 'uk', 'vk', 'ru', 'rru', 'rv', 'rrv', 'pk', 'qk', 'ww', 'advpk', 'advqk', 'diffpk', 'diffqk', 'Ui', 'Vi', 'SS', 'fpk', 'fqk', 'ustar2_s', 'ustar2_b', 'SxB', 'SyB', 'SxA', 'SyA', 'SxD', 'SyD'
    _pressure_arrays = 'dpdx', 'dpdy', 'idpdx', 'idpdy'
    _sealevel_arrays = ()
    _time_arrays = 'timestep', 'macrotimestep', 'split_factor', 'timedelta', 'time', 'istep', 'report'
    _all_fortran_arrays = tuple(['_%s' % name for name in _momentum_arrays + _pressure_arrays + _sealevel_arrays]) + ('uadv', 'vadv', 'uua', 'uva', 'vua', 'vva', 'uua3d', 'uva3d', 'vua3d', 'vva3d')
    __slots__ = _all_fortran_arrays + ('output_manager', 'input_manager', 'fabm_model', '_fabm_interior_diagnostic_arrays', '_fabm_horizontal_diagnostic_arrays', 'fabm_sources_interior', 'fabm_sources_surface', 'fabm_sources_bottom', 'fabm_vertical_velocity', 'fabm_conserved_quantity_totals', '_yearday', 'tracers', 'tracer_totals', 'logger', 'airsea', 'turbulence', 'density', 'buoy', 'temp', 'salt', 'pres', 'rad', 'par', 'par0', 'rho', 'sst', 'sss', 'SS', 'NN', 'ustar_s', 'ustar_b', 'taub', 'z0s', 'z0b', 'fwf', 'vertical_diffusion', 'tracer_advection', '_cum_river_height_increase', '_start_time', '_profile', 'radiation', '_ufirst', '_u3dfirst', 'diffuse_momentum', 'apply_bottom_friction') + _time_arrays

    def __init__(self, dom: domain.Domain, runtype: int, advection_scheme: operators.AdvectionScheme=operators.AdvectionScheme.HSIMT, apply_bottom_friction: bool=True, fabm: Union[bool, str, None]=None, gotm: Union[str, None]=None,
        turbulence: Optional[pygetm.mixing.Turbulence]=None, airsea: Optional[pygetm.airsea.Fluxes]=None, density: Optional[pygetm.density.Density]=None,
        logger: Optional[logging.Logger]=None, log_level: int=logging.INFO, A=0.7, g1=1., g2=15., internal_pressure_method=0, Am: float=0.):

        self.logger = dom.root_logger
        self.logger.setLevel(log_level)
        self.output_manager = output.OutputManager(rank=dom.tiling.rank, logger=self.logger.getChild('output_manager'))
        self.input_manager = dom.input_manager
        dom.field_manager = self.output_manager

        self.input_manager.set_logger(self.logger.getChild('input_manager'))

        # Disable bottom friction if physical bottom roughness is 0 everywhere
        if apply_bottom_friction and (numpy.ma.array(dom.z0b_min, mask=dom.mask==0) == 0.).any():
            self.logger.warning('Disabling bottom friction because bottom roughness is 0 in one or more points.')
            apply_bottom_friction = False
        self.apply_bottom_friction = apply_bottom_friction
        self.diffuse_momentum = Am > 0.
        if not self.diffuse_momentum:
            self.logger.info('Diffusion of momemntum if off because Am is 0')

        assert not dom.initialized
        super().__init__(dom, runtype, internal_pressure_method=internal_pressure_method, Am0=Am)
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

        self.domain.open_boundaries.z = self.wrap(self.domain.open_boundaries.z, b'zbdy', source=3)
        self.domain.open_boundaries.u = self.wrap(self.domain.open_boundaries.u, b'bdyu', source=1)
        self.domain.open_boundaries.v = self.wrap(self.domain.open_boundaries.v, b'bdyv', source=1)
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
            self.SS.fill(0.)

        # Derive old and new elevations, water depths and thicknesses from current surface elevation on T grid
        self.domain.update_depth_and_thicknesses()
        self.domain.update_depth_and_thicknesses()

        self._cum_river_height_increase = numpy.zeros((len(self.domain.rivers),))

        self.airsea = airsea or pygetm.airsea.FluxesFromMeteo()
        assert isinstance(self.airsea, pygetm.airsea.Fluxes)
        self.airsea.initialize(self.domain)

        self.fwf = dom.T.array(name='fwf', units='m s-1', long_name='freshwater flux', fill_value=FILL_VALUE)
        self.fwf.fill(0.)
        self.uadv = operators.Advection(dom.U, scheme=advection_scheme)
        self.vadv = operators.Advection(dom.V, scheme=advection_scheme)

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

        self._ufirst = False
        self._u3dfirst = False

        self.fabm_model = None

        if runtype > BAROTROPIC_2D:
            self.tracer_advection = operators.Advection(dom.T, scheme=advection_scheme)

            # Turbulence and associated fields
            self.turbulence = turbulence or pygetm.mixing.GOTM(nml_path=gotm)
            self.turbulence.initialize(self.domain.T)
            self.NN = dom.T.array(fill=0., z=INTERFACES, name='NN', units='s-2', long_name='buoyancy frequency squared', fill_value=FILL_VALUE)
            self.ustar_s = dom.T.array(fill=0., name='ustar_s', units='m s-1', long_name='shear velocity (surface)', fill_value=FILL_VALUE)
            self.ustar_b = dom.T.array(fill=0., name='ustar_b', units='m s-1', long_name='shear velocity (bottom)', fill_value=FILL_VALUE)
            self.z0s = dom.T.array(fill=0.1, name='z0s', units='m', long_name='hydrodynamic roughness (surface)', fill_value=FILL_VALUE)
            self.taub = dom.T.array(fill=0., name='taub', units='Pa', long_name='bottom shear stress', fill_value=FILL_VALUE, fabm_standard_name='bottom_stress')

            self.vertical_diffusion = operators.VerticalDiffusion(dom.T, cnpar=1.)

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

            self.pres = dom.depth

        self.sst = dom.T.array(name='sst', units='degrees_Celsius', long_name='sea surface temperature', fill_value=FILL_VALUE)

        if runtype == BAROCLINIC:
            self.temp = self.create_tracer(name='temp', units='degrees_Celsius', long_name='conservative temperature', fabm_standard_name='temperature', fill_value=FILL_VALUE, source=dom.T.array(fill=0., z=CENTERS, units='W m-2'), source_scale=1. / (RHO0 * CP), rivers_follow_target_cell=True)
            self.salt = self.create_tracer(name='salt', units='-', long_name='absolute salinity', fabm_standard_name='practical_salinity', fill_value=FILL_VALUE, source=dom.T.array(fill=0., z=CENTERS))
            self.pres.saved = True
            self.temp.fill(5.)
            self.salt.fill(35.)
            self.rho = dom.T.array(z=CENTERS, name='rho', units='kg m-3', long_name='density', fabm_standard_name='density')
            self.buoy = dom.T.array(z=CENTERS, name='buoy', units='m s-2', long_name='buoyancy')
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

    def start(self, time: Union[cftime.datetime, datetime.datetime], timestep: float, split_factor: int=1, report: int=10, save: bool=True, profile: str=False):
        """This should be called after the output configuration is complete (because we need toknow when variables need to be saved),
        and after the FABM model has been provided with all dependencies"""
        if isinstance(time, datetime.datetime):
            time = cftime.datetime(time.year, time.month, time.day, time.hour, time.minute, time.second, time.microsecond)
        self.logger.info('Starting simulation at %s' % time)
        self.timestep = timestep
        self.split_factor = split_factor
        self.macrotimestep = self.timestep * self.split_factor
        self.timedelta = datetime.timedelta(seconds=timestep)
        self.time = time
        self.istep = 0
        self.report = report

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
                self._yearday.value = (self.time - cftime.datetime(self.time.year, 1, 1)).total_seconds() / 86400.
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

        # Update inputs and forcing variables based on the current time and state
        self.update_forcing(macro_active=True)

        if self.apply_bottom_friction:
            self.bottom_friction_3d()

        # Start output manager
        self.output_manager.start(self.istep, self.time, save=save)

        # Record true start time for performance analysis
        self._start_time = timeit.default_timer()

        # Start profiing if requested
        self._profile = None
        if profile:
            import cProfile
            pr = cProfile.Profile()
            self._profile = (profile, pr)
            pr.enable()

    def advance(self):
        # Update momentum from time=-1/2 to +1/2, using surface stresses and pressure gradients defined at time=0
        # Inputs and outputs on U and V grids. As depth-averaged velocities u1 and v1 are already in sync
        # with transports U and V and water depths D, we skip their computation by passing u1_ready=True
        self.update_2d_momentum(self.timestep, self.airsea.taux_U, self.airsea.tauy_V, self.dpdx, self.dpdy, u1_ready=True)

        # Update surface elevation on T grid from time=0 to time=1 using transports U and V at time=1/2 and
        # freshwater fluxes at time=0. Exchange elevation halos so that depths and thicknesses can be computed
        # everywhere without further halo exchange
        self.update_sealevel(self.timestep, self.U, self.V, self.fwf)
        self.domain.T.z.update_halos()

        # Track cumulative increase in elevation due to river inflow over the current macrotimestep
        self._cum_river_height_increase += self.domain.rivers.flow * self.domain.rivers.iarea * self.timestep

        # Update the time
        self.time += self.timedelta
        self.istep += 1
        macro_active = self.istep % self.split_factor == 0
        if self.report != 0 and self.istep % self.report == 0:
            self.logger.info(self.time)

        if self.runtype > BAROTROPIC_2D and macro_active:
            # Depth-integrated transports have been summed over all microtimesteps. Now average them.
            self.Ui.all_values *= 1. / self.split_factor
            self.Vi.all_values *= 1. / self.split_factor

            # Use previous source terms for biogeochemistry (valid for the start of the current macrotimestep)
            # to update tracers. This should be done before the tracer concentrations change due to transport
            # or rivers, as the source terms are only valid for the current tracer concentrations.
            if self.fabm_model:
                self.update_fabm(self.macrotimestep)

            # Update layer thicknesses and tracer concentrations to account for river inflow
            # between start and end of the current macrotimestep.
            self.add_rivers_3d()

            # Update 3D elevations and layer thicknesses. New elevation on T grid will match elevation at end of 2D timestep,
            # thicknesses on T grid will match. Elevation and thicknesses on U/V grids will be 1/2 macrotimestep behind. Old
            # elevations zio and thicknesses ho will be one macrotimestep behind new elevations zin and thicknesses hn.
            self.domain.update_depth_and_thicknesses()

            # Update presssure gradient for start of the 3D time step
            # Note that we use previously-recorded elevations and surface pressure at the start of the 3D timestep
            self.update_surface_pressure_gradient(self.domain.T.zio, self.airsea.spo)

            # Update momentum from time=-1/2 to 1/2 of the macrotimestep, using forcing defined at time=0
            # For this purpose, surface stresses at the end of the previous macrotimestep were saved (taux_Uo, tauy_Vo)
            # Pressure gradients dpdx and dpdy have just been updated to match the start of the current macrotimestep
            # Internal pressure idpdx and idpdy were calculated at the end of the previous macrotimestep and are therefore ready as-is.
            self.update_3d_momentum(self.macrotimestep, self.airsea.taux_Uo, self.airsea.tauy_Vo, self.dpdx, self.dpdy, self.idpdx, self.idpdy, self.turbulence.num)

            if self.runtype == BAROCLINIC:
                # Update total stresses (x and y combined) and calculate the friction velocities (m s-1)
                # This is for turbulence (GOTM), so all on the T grid
                # Note that surface stress is currently defined at the start of the MICROtimestep
                # (only one microtimestep before the end of the current macrotimestep),
                # whereas bottom stress is at 1/2 of the macrotimestep, since it is computed from velocities that have just been updated.
                self.update_stresses(self.airsea.taux, self.airsea.tauy)
                numpy.sqrt(self.ustar2_s.all_values, out=self.ustar_s.all_values)
                numpy.sqrt(self.ustar2_b.all_values, out=self.ustar_b.all_values)
                self.taub.all_values[...] = self.ustar2_b.all_values * RHO0

                # Update turbulent quantities (T grid - interfaces) from time=0 to time=1 (macrotimestep),
                # using surface/buoyancy-related forcing at time=0, and velocity-related forcing at time=1/2
                #self.domain.T.z0b.all_values[1:, 1:] = 0.5 * (numpy.maximum(self.domain.U.z0b.all_values[1:, 1:], self.domain.U.z0b.all_values[1:, :-1]) + numpy.maximum(self.domain.V.z0b.all_values[:-1, 1:], self.domain.V.z0b.all_values[1:, :-1]))
                self.turbulence(self.macrotimestep, self.ustar_s, self.ustar_b, self.z0s, self.domain.T.z0b, self.NN, self.SS)

                # Temperature sources (T grid), defined at the start of the macrotimestep
                self.temp.source.all_values[...] = self.radiation.swr_abs.all_values
                self.temp.source.all_values[-1, ...] += self.airsea.shf.all_values

                # Advection of passive tracers (including biogeochemical ones)from time=0 to time=1 (macrotimestep),
                # using velocities defined at time=1/2
                w_res = []
                for tracer in self.tracers:
                    w = self.ww
                    if tracer.vertical_velocity is not None and tracer.vertical_velocity.all_values.any():
                        w = self.domain.T.array(fill=0., z=INTERFACES)
                        tracer.vertical_velocity.interp(w)
                        w.all_values += self.ww.all_values
                    w_res.append(w)
                self.tracer_advection.apply_3d_batch(self.uk, self.vk, self.ww, self.macrotimestep, self.tracers, w_vars=w_res)

                # Diffusion of passive tracers (including biogeochemical ones)
                # This simultaneously time-integrates source terms, if specified - but BGC sources have already been handled separately.
                for tracer in self.tracers:
                    if tracer.source is not None:
                        tracer.source.all_values *= self.macrotimestep * tracer.source_scale
                    self.vertical_diffusion(self.turbulence.nuh, self.macrotimestep, tracer, ea4=tracer.source)

            self.report_domain_integrals()

            # Reset depth-integrated transports that will be incremented over subsequent 3D timestep.
            self.Ui.all_values.fill(0.)
            self.Vi.all_values.fill(0.)

        # Update all inputs and fluxes that will drive the next state update
        self.update_forcing(macro_active)

        self.output_manager.save(self.timestep * self.istep, self.istep, self.time)

        return macro_active

    def update_forcing(self, macro_active: bool):
        # Update all inputs
        self.domain.input_manager.update(self.time, include_3d=macro_active)

        if self.runtype == BAROCLINIC and macro_active:
            # Update tracer values at open boundaries. This must be done after input_manager.update,
            # but before diagnostics/forcing variables derived from the tracers are calculated
            if self.domain.open_boundaries:
                for tracer in self.tracers:
                    tracer.open_boundaries.update()

            # Update density, buoyancy and internal pressure to keep them in sync with T and S.
            self.density.get_density(self.salt, self.temp, p=self.pres, out=self.rho)
            self.rho.update_halos(parallel.Neighbor.TOP_AND_RIGHT)   # valid rho around all U and V needed for internal pressure; not yet valid becasue T&S were not valid in halos when rho was calculated
            self.buoy.all_values[...] = (-GRAVITY / RHO0) * (self.rho.all_values - RHO0)
            self.update_internal_pressure_gradient(self.buoy, self.SxB, self.SyB)

            # From conservative temperature to in-situ sea surface temperature,
            # needed to compute heat/momentum fluxes at the surface
            self.density.get_potential_temperature(self.salt.isel(-1), self.temp.isel(-1), out=self.sst)

            # Calculate squared buoyancy frequency NN (T grid, interfaces between layers)
            self.density.get_buoyancy_frequency(self.salt, self.temp, p=self.pres, out=self.NN)

        # Update surface elevation z on U, V, X grids and water depth D on all grids
        # This is based on old and new elevation (T grid) for the microtimestep.
        # Thus, for grids lagging 1/2 a timestep behind (U, V, X grids), the elevations
        # and water depths will be representative for 1/2 a MICROtimestep ago.
        self.domain.update_depth()

        # Calculate 2D (depth-averaged) velocities using updated water depths
        # Note that as long as transports U and V are 0 in masked points=, u1 and v1 will be too,
        # provided D contains only non-zero finite values
        numpy.divide(self.U.all_values, self.U.grid.D.all_values, out=self.u1.all_values)
        numpy.divide(self.V.all_values, self.V.grid.D.all_values, out=self.v1.all_values)

        # Calculate bottom friction using updated depth-averaged velocities u1 and v1
        # Warning: this uses velocities u1 and v1 at masked points, which therefore need to be kept at 0
        if self.apply_bottom_friction:
            self.bottom_friction_2d()

        # Update freshwater fluxes (TODO: add precipitation)
        self.fwf.all_values[self.domain.rivers.j, self.domain.rivers.i] = self.domain.rivers.flow * self.domain.rivers.iarea

        # Update air-sea fluxes of heat and momentum (T grid for all, U and V grid for x and y stresses respectively)
        self.airsea(self.time, self.sst, self.sss, calculate_heat_flux=macro_active and self.runtype == BAROCLINIC)

        # Update elevation at the open boundaries. This must be done before update_surface_pressure_gradient
        self.update_sealevel_boundaries(self.timestep)

        # Calculate the surface pressure gradient in the U and V points.
        # Note: this requires elevation and surface air pressure (both on T grid) to be valid in the halos,
        # which is guaranteed for elevation (halo exchange happens just after update), and for air pressure
        # if it is managed by the input manager (e.g. read from file)
        self.airsea.sp.update_halos(parallel.Neighbor.TOP_AND_RIGHT)
        self.update_surface_pressure_gradient(self.domain.T.z, self.airsea.sp)

        if self.runtype == BAROCLINIC and macro_active:
            # Update radiation. This must come after the airsea update, which is responsible for calculating swr
            self.radiation(self.airsea.swr)

            # Update source terms of biogeochemistry, using the new tracer concentrations
            # Do this last because FABM could depend on any of the variables computed before
            if self.fabm_model:
                self.update_fabm_sources()

            # Save forcing variables for the next baroclinic update
            self.airsea.spo.all_values[...] = self.airsea.sp.all_values
            self.airsea.taux_Uo.all_values[...] = self.airsea.taux_U.all_values
            self.airsea.tauy_Vo.all_values[...] = self.airsea.tauy_V.all_values

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
            tracer.update_halos_start(parallel.Neighbor.LEFT_AND_RIGHT)   # to prepare for advection
        self.domain.T.hn.all_values[:, self.domain.rivers.j, self.domain.rivers.i] = h + h_increase
        self.domain.T.zin.all_values[self.domain.rivers.j, self.domain.rivers.i] += self._cum_river_height_increase
        self._cum_river_height_increase.fill(0.)

    def update_fabm_sources(self):
        """Update FABM sources, vertical velocities, and diagnostics. This does not update the state variables themselves; that is done by update_fabm"""
        if self._yearday:
            self._yearday.value = (self.time - cftime.datetime(self.time.year, 1, 1)).total_seconds() / 86400.
        self.fabm_model.get_sources(out=(self.fabm_sources_interior, self.fabm_sources_surface, self.fabm_sources_bottom))
        self.fabm_model.get_vertical_movement(self.fabm_vertical_velocity)

    def update_fabm(self, timestep: float, repair: bool=True):
        """Time-integrate source terms of all FABM state variables (3D pelagic tracers as well as bottom- and surface-attached variables)"""
        self.fabm_sources_interior *= timestep
        self.fabm_sources_surface *= timestep
        self.fabm_sources_bottom *= timestep
        self.fabm_model.interior_state += self.fabm_sources_interior
        self.fabm_model.surface_state += self.fabm_sources_surface
        self.fabm_model.bottom_state += self.fabm_sources_bottom
        self.fabm_model.check_state(repair=repair)

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

    def advect_2d_momentum(self, U: core.Array, V: core.Array, timestep: float, advU: core.Array, advV: core.Array, diffu1: core.Array, diffv1: core.Array, u1_ready: bool=False):
        # Advect and optionally diffuse depth-averaged u and v velocity (u1, v1),
        # using velocities reconstructed from transports interpolated to their own advection grids
        # Then reconstruct the resulting change in transports U and V and return this (advU, advV)

        if not u1_ready:
            numpy.divide(U.all_values, U.grid.D.all_values, out=self.u1.all_values)
            numpy.divide(V.all_values, V.grid.D.all_values, out=self.v1.all_values)

        if self.diffuse_momentum:
            # Compute velocity diffusion contribution to transport sources.
            # This uses depth-averaged velocities u1 and v1, which therefore have to be up to date
            # Water depths should be in sync with velocities, which means they should lag 1/2 a timestep behind the tracer/T grid
            self.momentum_diffusion_driver(self.domain.D_T_half, self.domain.X.D, self.u1,self.v1, diffu1, diffv1)

        itimestep = 1. / timestep

        # Advection of u velocity (u1)
        U.interp(self.uua)
        V.interp(self.uva)
        self.uua.all_values /= self.domain.UU.D.all_values
        self.uva.all_values /= self.domain.UV.D.all_values
        self.uadv(self.uua, self.uva, timestep, self.u1, skip_initial_halo_exchange=True)
        advU.all_values[...] = (self.u1.all_values * self.uadv.D - U.all_values) * itimestep

        # Advection of v velocity (v1)
        self.U.interp(self.vua)
        self.V.interp(self.vva)
        self.vua.all_values /= self.domain.VU.D.all_values
        self.vva.all_values /= self.domain.VV.D.all_values
        self.vadv(self.vua, self.vva, timestep, self.v1, skip_initial_halo_exchange=True)
        advV.all_values[...] = (self.v1.all_values * self.vadv.D - V.all_values) * itimestep

    def update_2d_momentum(self, timestep: float, tausx: core.Array, tausy: core.Array, dpdx: core.Array, dpdy: core.Array, u1_ready: bool=False):
        """Update depth-integrated transports (U, V) and depth-averaged velocities (u1, v1).
        This will also update their halos. Note that velocities u1 and v1 are assumed to be in sync
        with current transports and water depths (u1=U/U.grid.D and v1=V/V.grid.D)"""

        # Advect depth-averaged u and v velocity (u1, v1) from t-1/2 to t+1/2,
        # store the resulting change in transports U and V (advU, advV),
        # to be taken into account by transport update further down.
        # Since self.u1=self.U/self.U.grid.D (and same for v1/V), we state that with u1_ready=True
        self.advect_2d_momentum(self.U, self.V, timestep, self.advU, self.advV, self.diffu1, self.diffv1, u1_ready=u1_ready)

        # Update 2D transports from t-1/2 to t+1/2
        if self._ufirst:
            self.u_2d(timestep, tausx, dpdx)
            self.U.update_halos()
            self.coriolis_fu()
            self.v_2d(timestep, tausy, dpdy)
            self.V.update_halos()
            self.coriolis_fv()
        else:
            self.v_2d(timestep, tausy, dpdy)
            self.V.update_halos()
            self.coriolis_fv()
            self.u_2d(timestep, tausx, dpdx)
            self.U.update_halos()
            self.coriolis_fu()
        self._ufirst = not self._ufirst

        self.Ui.all_values += self.U.all_values
        self.Vi.all_values += self.V.all_values

    def update_3d_momentum(self, timestep: float, tausx: core.Array, tausy: core.Array, dpdx: core.Array, dpdy: core.Array, idpdx: core.Array, idpdy: core.Array, viscosity: core.Array):
        """Update depth-explicit transports (pk, qk) and velocities (uk, vk). This will also update their halos."""
        # Do the halo exchange for viscosity, as this needs to be interpolated to the U and V grids.
        # For that, information from the halos is used.
        viscosity.update_halos(parallel.Neighbor.TOP_AND_RIGHT)

        # Update horizontal transports. Also update the halos so that transports (and more importantly, the velocities
        # derived subsequently) are valid there. Information from these halos is needed for many reasons:
        # - the Coriolis update requires horizontal velocities at the four points surrounding each U/V point
        # - to advect the horizontal velocities themselves, for which they need to be valid in the halos in the direction of transport
        # - to advect quantities defined on the T grid, as this requires horizontal velocities at the boundaries of every T cell
        #   of the subdomain interior; this includes cells at the very Western and Southern boundary,
        #   which for U and V grids lie within the halo
        # - to allow interpolation of horizontal velocities to the advection grids for momentum (UU, UV, VU, VV),
        #   which again requires halos values
        # - to calculate vertical velocities, which requires horizontal transports at the four interfaces around every T point
        if self._u3dfirst:
            self.pk_3d(timestep, tausx, dpdx, idpdx, viscosity.interp(self.domain.U))
            self.pk.update_halos()
            self.coriolis_fpk()
            self.qk_3d(timestep, tausy, dpdy, idpdy, viscosity.interp(self.domain.V))
            self.qk.update_halos()
            self.coriolis_fqk()
        else:
            self.qk_3d(timestep, tausy, dpdy, idpdy, viscosity.interp(self.domain.V))
            self.qk.update_halos()
            self.coriolis_fqk()
            self.pk_3d(timestep, tausx, dpdx, idpdx, viscosity.interp(self.domain.U))
            self.pk.update_halos()
            self.coriolis_fpk()
        self._u3dfirst = not self._u3dfirst

        # Infer vertical velocity from horizontal transports and desired layer height change.
        # This is done at all points surrounding U and V points, so no further halo exchange of w is needed
        # to support interpolation to U and V grids later on. This does require that transports are up to date in halos.
        self.w_3d(timestep)

        itimestep = 1. / timestep

        # Compute 3D velocities (m s-1) from 3D transports (m2 s-1) by dividing by layer heights
        # Both velocities and U/V thicknesses are now at time 1/2
        numpy.divide(self.pk.all_values, self.U.grid.hn.all_values, where=self.pk.grid.mask.all_values != 0, out=self.uk.all_values)
        numpy.divide(self.qk.all_values, self.V.grid.hn.all_values, where=self.qk.grid.mask.all_values != 0, out=self.vk.all_values)

        # Use updated velocities (uk, vk) to compute shear frequency (SS) at T points (interior only, not in halos)
        self.update_shear_frequency(viscosity)

        # Calculate bottom friction from updated velocities (and syncronized layer thicknesses hn)
        # This needs to be done before derived quantities such as bottom stress are calculated
        if self.apply_bottom_friction:
            self.bottom_friction_3d()

        # Interpolate 3D velocities to advection grids.
        # This needs to be done before uk/vk are changed by the advection operator (apply_3d).
        self.uk.interp(self.uua3d)
        self.vk.interp(self.uva3d)
        self.uk.interp(self.vua3d)
        self.vk.interp(self.vva3d)

        # Advect 3D u velocity from time=1/2 to 1 1/2 using velocities interpolated to its own advection grids
        # JB the alternative would be to interpolate transports and then divide by (colocated) layer heights, like we do for 2D
        self.uadv.apply_3d(self.uua3d, self.uva3d, self.ww.interp(self.uk.grid), timestep, self.uk, new_h=True, skip_initial_halo_exchange=True)
        self.advpk.all_values[...] = (self.uk.all_values * self.uadv.h - self.pk.all_values) * itimestep

        # Advect 3D v velocity from time=1/2 to 1 1/2 using velocities interpolated to its own advection grids
        # JB the alternative would be to interpolate transports and then divide by (colocated) layer heights, like we do for 2D
        self.vadv.apply_3d(self.vua3d, self.vva3d, self.ww.interp(self.vk.grid), timestep, self.vk, new_h=True, skip_initial_halo_exchange=True)
        self.advqk.all_values[...] = (self.vk.all_values * self.vadv.h - self.qk.all_values) * itimestep

        # Restore velocity at time=1/2 (the final value at the end of the current timestep)
        numpy.divide(self.pk.all_values, self.U.grid.hn.all_values, where=self.pk.grid.mask.all_values != 0, out=self.uk.all_values)
        numpy.divide(self.qk.all_values, self.V.grid.hn.all_values, where=self.qk.grid.mask.all_values != 0, out=self.vk.all_values)

        if self.diffuse_momentum:
            # Note that thicknesses should be in sync with velocities uk and vk
            # This means they should lag 1/2 a timestep behind the T grid (already the case for X, but for T we use 1/2(ho+hn))
            self.momentum_diffusion_driver(self.domain.h_T_half, self.domain.X.hn, self.uk, self.vk, self.diffpk, self.diffqk)

        # Compute slow (3D) advection contribution to 2D advection.
        # This is done by comparing the previously calculated depth-integrated 3D transport
        # (between centers of the current and next macrotime step) with the newly calculated
        # depth-integrated transport based on accumulated 2D transports (accumulated over the
        # current macrotimestep, and thus representative for its center).
        self.advect_2d_momentum(self.Ui, self.Vi, timestep, self.SxA, self.SyA, self.SxD, self.SyD)
        self.SxA.all_values[...] = self.advpk.all_values.sum(axis=0) - self.SxA.all_values
        self.SyA.all_values[...] = self.advqk.all_values.sum(axis=0) - self.SyA.all_values
        self.SxD.all_values[...] = self.diffpk.all_values.sum(axis=0) - self.SxD.all_values
        self.SyD.all_values[...] = self.diffqk.all_values.sum(axis=0) - self.SyD.all_values

    @property
    def Ekin(self, rho0: float=RHO0):
        dom = self.domain
        U = self.U.interp(dom.T)
        V = self.V.interp(dom.T)
        vel2_D2 = U**2 + V**2
        return 0.5 * rho0 * dom.T.area * vel2_D2 / dom.T.D

for membername in Simulation._all_fortran_arrays:
    setattr(Simulation, membername[1:], property(operator.attrgetter(membername)))
