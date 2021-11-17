import operator
from typing import Union, Optional
import itertools
import logging
import datetime

import numpy

from . import constants
from . import _pygetm
from . import core
from . import domain
from . import output
from . import pyfabm
from . import mixing
import pygetm.airsea

class Simulation(_pygetm.Simulation):
    _momentum_arrays = 'U', 'V', 'fU', 'fV', 'advU', 'advV', 'u1', 'v1', 'bdyu', 'bdyv'
    _pressure_arrays = 'dpdx', 'dpdy'
    _sealevel_arrays = 'zbdy',
    _time_arrays = 'timestep', 'macrotimestep', 'split_factor', 'timedelta', 'time', 'istep', 'report'
    _all_fortran_arrays = tuple(['_%s' % name for name in _momentum_arrays + _pressure_arrays + _sealevel_arrays]) + ('uadv', 'vadv', 'uua', 'uva', 'vua', 'vva')
    __slots__ = _all_fortran_arrays + ('output_manager', 'input_manager', 'fabm_model', '_fabm_interior_diagnostic_arrays', '_fabm_horizontal_diagnostic_arrays', 'fabm_sources_interior', 'fabm_sources_surface', 'fabm_sources_bottom', 'tracers', 'logger', 'airsea', 'turbulence', 'temp', 'salt', 'sst', 'temp_source', 'salt_source', 'shf', 'SS', 'NN', 'u_taus', 'u_taub', 'z0s', 'z0b', 'vertical_diffusion') + _time_arrays

    def __init__(self, dom: domain.Domain, runtype: int, advection_scheme: int=4, apply_bottom_friction: bool=True, fabm: Union[bool, str, None]=None, gotm: Union[str, None]=None, turbulence: Optional[mixing.Turbulence]=None, airsea: Optional[pygetm.airsea.Fluxes]=None, logger: Optional[logging.Logger]=None, log_level: int=logging.INFO):
        if logger is None:
            logging.basicConfig(level=log_level)
            self.logger = logging.getLogger()

        self.output_manager = output.OutputManager(rank=dom.tiling.rank, logger=self.logger.getChild('output_manager'))
        self.input_manager = dom.input_manager
        dom.field_manager = self.output_manager

        self.input_manager.set_logger(self.logger.getChild('input_manager'))
        dom.logger = self.logger.getChild('domain')

        assert not dom.initialized
        _pygetm.Simulation.__init__(self, dom, runtype, apply_bottom_friction)

        for name in Simulation._momentum_arrays:
            setattr(self, '_%s' % name, self.wrap(core.Array(name=name), name.encode('ascii'), source=1))
        for name in Simulation._pressure_arrays:
            setattr(self, '_%s' % name, self.wrap(core.Array(name=name), name.encode('ascii'), source=2))
        for name in Simulation._sealevel_arrays:
            setattr(self, '_%s' % name, self.wrap(core.Array(name=name), name.encode('ascii'), source=3))

        self.update_depth()

        self.airsea = airsea or pygetm.airsea.FluxesFromMeteo(self.domain)

        self.uadv = _pygetm.Advection(dom.U, scheme=advection_scheme)
        self.vadv = _pygetm.Advection(dom.V, scheme=advection_scheme)

        self.uua = dom.UU.array(fill=numpy.nan)
        self.uva = dom.UV.array(fill=numpy.nan)
        self.vua = dom.VU.array(fill=numpy.nan)
        self.vva = dom.VV.array(fill=numpy.nan)

        self.tracers = []

        self.fabm_model = None

        if runtype == 4:
            self.temp = dom.T.array(fill=5., is_3d=True, name='temp', units='degrees_Celsius', long_name='conservative temperature', fabm_standard_name='temperature')
            self.salt = dom.T.array(fill=35., is_3d=True, name='salt', units='-', long_name='absolute salinity', fabm_standard_name='practical_salinity')
            self.sst = self.temp.isel(z=-1)
            self.temp_source = dom.T.array(fill=0., is_3d=True, units='W m-2')
            self.salt_source = dom.T.array(fill=0., is_3d=True)
            self.shf = self.temp_source.isel(z=-1, name='shf', long_name='surface heat flux')

            self.SS = dom.W.array(fill=0., is_3d=True, name='SS', units='s-2', long_name='shear frequency squared')
            self.NN = dom.W.array(fill=0., is_3d=True, name='NN', units='s-2', long_name='buoyancy frequency squared')
            self.u_taus = dom.T.array(fill=0., name='u_taus', units='m s-1', long_name='surface shear velocity')
            self.u_taub = dom.T.array(fill=0., name='u_taub', units='m s-1', long_name='bottom shear velocity')
            self.z0s = dom.T.array(fill=0.1, name='z0s', units='m', long_name='hydrodynamic surface roughness')
            self.z0b = dom.T.array(fill=0.1, name='z0b', units='m', long_name='hydrodynamic bottom roughness')

            self.turbulence = turbulence or mixing.GOTM(self.domain, nml_path=gotm)

            self.vertical_diffusion = pygetm.VerticalDiffusion(dom.T, cnpar=1.)

            # ensure ho and hn are up to date and identical
            dom.do_vertical()

            if fabm:
                def fabm_variable_to_array(variable, send_data: bool=False, **kwargs):
                    ar = core.Array(name=variable.output_name, units=variable.units, long_name=variable.long_name, fill_value=variable.missing_value, dtype=self.fabm_model.fabm.dtype, grid=self.domain.T, **kwargs)
                    if send_data:
                        ar.wrap_ndarray(variable.data)
                    ar.register()
                    return ar

                self.fabm_model = pyfabm.Model(fabm if isinstance(fabm, str) else 'fabm.yaml', shape=self.domain.T.hn.all_values.shape, libname='fabm_c')
                for variable in itertools.chain(self.fabm_model.interior_state_variables, self.fabm_model.surface_state_variables, self.fabm_model.bottom_state_variables):
                    ar = fabm_variable_to_array(variable, send_data=True)
                    if ar.ndim == 3:
                        self.add_tracer(ar)
                self._fabm_interior_diagnostic_arrays = [fabm_variable_to_array(variable, shape=self.domain.T.hn.shape) for variable in self.fabm_model.interior_diagnostic_variables]
                self._fabm_horizontal_diagnostic_arrays = [fabm_variable_to_array(variable, shape=self.domain.T.H.shape) for variable in self.fabm_model.horizontal_diagnostic_variables]
                self.fabm_model.link_mask(self.domain.T.mask.all_values)
                self.fabm_model.link_cell_thickness(self.domain.T.hn.all_values)

                self.fabm_sources_interior = numpy.empty_like(self.fabm_model.interior_state)
                self.fabm_sources_surface = numpy.empty_like(self.fabm_model.surface_state)
                self.fabm_sources_bottom = numpy.empty_like(self.fabm_model.bottom_state)
        else:
            self.sst = self.airsea.t2m

    def get_fabm_dependency(self, name):
        variable = self.fabm_model.dependencies.find(name)
        arr = self.domain.T.array(name=variable.output_name, units=variable.units, long_name=variable.long_name, is_3d=len(variable.shape) == 3)
        variable.link(arr.all_values)
        return arr

    def add_tracer(self, array: core.Array, source: Optional[core.Array]=None):
        assert array.grid is self.domain.T
        assert source is None or source.grid is self.domain.T
        self.tracers.append((array, source))

    def start(self, time: datetime.datetime, timestep, split_factor=1, report=10, save: bool=True):
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
                    variable.link(field.all_values)

            # Start FABM. This verifies whether all dependencies are fulfilled and freezes the set of dsiagsntoics that will be saved.
            assert self.fabm_model.start(), 'FABM failed to start. Likely its configuration is incomplete.'

            # Fill GETM placeholder arrays for all FABM diagnostics that will be computed/saved.
            for variable, ar in zip(itertools.chain(self.fabm_model.interior_diagnostic_variables, self.fabm_model.horizontal_diagnostic_variables), itertools.chain(self._fabm_interior_diagnostic_arrays, self._fabm_horizontal_diagnostic_arrays)):
                if ar.saved:
                    ar.wrap_ndarray(variable.data)

            # Apply mask to state variables
            for variable in itertools.chain(self.fabm_model.interior_state_variables, self.fabm_model.surface_state_variables, self.fabm_model.bottom_state_variables):
                variable.value[..., self.domain.T.mask.all_values == 0] = variable.missing_value

            self.update_fabm_sources()

        self.domain.input_manager.update(time)
        self.output_manager.start()

    def advance(self):
        self.time += self.timedelta

        # Update all inputs (todo: separate 2D and 3D inputs, as the latter are needed much less frequently)
        self.domain.input_manager.update(self.time)

        # Update air-sea fluxes of heat and momentum (T grid for all, U and V grid for x and y stresses respectively)
        self.airsea(self.sst)

        # Update elevation at the open boundaries
        self.update_sealevel_boundaries(self.timestep)

        # Calculate the surface pressure gradient (T grid)
        # This requires elevation and surface pressure (both on T grid) to be valid in the halos
        self.airsea.sp.update_halos()
        self.update_surface_pressure_gradient(self.domain.T.z, self.airsea.sp)

        # Update momentum using surface stresses and pressure gradients
        self.uv_momentum_2d(self.timestep, self.airsea.taux_U, self.airsea.tauy_V, self.dpdx, self.dpdy)

        # Update halos of U and V as the sea level update below will go into the halos
        self.U.update_halos()
        self.V.update_halos()

        # Update sea level on T grid, and from that calculate sea level and water depth on all grids
        self.update_sealevel(self.timestep, self.U, self.V)
        self.update_depth()

        self.istep += 1
        if self.report != 0 and self.istep % self.report == 0:
            self.logger.info(self.time)

        if self.runtype == 4 and self.istep % self.split_factor == 0:
            # Buoyancy frequency and turbulence
            pygetm.mixing.get_buoyancy_frequency(self.salt, self.temp, out=self.NN)
            self.u_taus.all_values[...] = (self.airsea.taux.all_values**2 + self.airsea.tauy.all_values**2)**0.25 / numpy.sqrt(constants.rho0)
            self.turbulence(self.macrotimestep, self.u_taus, self.u_taub, self.z0s, self.z0b, self.NN, self.SS)

            # Temperature and salinity
            self.shf.all_values[...] = self.airsea.qe.all_values + self.airsea.qh.all_values
            self.vertical_diffusion(self.turbulence.nuh, self.macrotimestep, self.temp, ea4=self.temp_source * (self.macrotimestep / (constants.rho0 * constants.cp)))
            self.vertical_diffusion(self.turbulence.nuh, self.macrotimestep, self.salt, ea4=self.salt_source)

            # Transport of passive tracers (including biogeochemical ones)
            for array, src in self.tracers:
                self.vertical_diffusion(self.turbulence.nuh, self.macrotimestep, array, ea4=src)

            # Time-integrate source terms of biogeochemistry
            if self.fabm_model:
                self.update_fabm(self.macrotimestep)
            self.update_fabm_sources()

        self.output_manager.save()

    def finish(self):
        self.output_manager.close()

    def update_fabm_sources(self):
        self.fabm_model.get_sources(out=(self.fabm_sources_interior, self.fabm_sources_surface, self.fabm_sources_bottom))

    def update_fabm(self, timestep: float):
        self.fabm_model.interior_state += self.fabm_sources_interior * timestep
        self.fabm_model.surface_state += self.fabm_sources_surface * timestep
        self.fabm_model.bottom_state += self.fabm_sources_bottom * timestep

    def uv_momentum_2d(self, timestep: float, tausx: core.Array, tausy: core.Array, dpdx: core.Array, dpdy: core.Array):
        # compute velocities at time=n-1/2
        self.u1.all_values[:, :] = self.U.all_values / self.U.grid.D.all_values
        self.v1.all_values[:, :] = self.V.all_values / self.V.grid.D.all_values

        itimestep = 1. / timestep

        # Advect U using velocities interpolated to its own advection grids
        self.U.interp(self.uua)
        self.V.interp(self.uva)
        self.uua.all_values[...] /= self.domain.UU.D.all_values
        self.uva.all_values[...] /= self.domain.UV.D.all_values
        self.uadv(self.uua, self.uva, timestep, self.u1)
        self.advU.all_values[...] = (self.u1.all_values * self.uadv.D - self.U.all_values) * itimestep

        # Advect V using velocities interpolated to its own advection grids
        self.U.interp(self.vua)
        self.V.interp(self.vva)
        self.vua.all_values[...] /= self.domain.VU.D.all_values
        self.vva.all_values[...] /= self.domain.VV.D.all_values
        self.vadv(self.vua, self.vva, timestep, self.v1)
        self.advV.all_values[...] = (self.v1.all_values * self.vadv.D - self.V.all_values) * itimestep

        # Restore velocity at time=n-1/2
        self.u1.all_values[:, :] = self.U.all_values / self.U.grid.D.all_values
        self.v1.all_values[:, :] = self.V.all_values / self.V.grid.D.all_values

        _pygetm.Simulation.uv_momentum_2d(self, timestep, tausx, tausy, dpdx, dpdy)

    def update_depth(self):
        # Halo exchange for sea level on T grid
        self.domain.T.z.update_halos()

        # Compute sea level on U, V, X grids.
        # Note that this must be at time=n+1/2, whereas sea level on T grid is now at time=n+1.
        z_T_half = 0.5 * (self.domain.T.zo + self.domain.T.z)
        z_T_half.interp(self.domain.U.z)
        z_T_half.interp(self.domain.V.z)
        z_T_half.interp(self.domain.X.z)

        # Halo exchange for sea level on U, V, X grids
        self.domain.U.z.update_halos()
        self.domain.V.z.update_halos()
        self.domain.X.z.update_halos()

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
    def Ekin(self, rho0: float=1025.):
        dom = self.domain
        U = self.U.interp(dom.T)
        V = self.V.interp(dom.T)
        vel2_D2 = U**2 + V**2
        return 0.5 * rho0 * dom.T.area * vel2_D2 / dom.T.D

for membername in Simulation._all_fortran_arrays:
    setattr(Simulation, membername[1:], property(operator.attrgetter(membername)))
