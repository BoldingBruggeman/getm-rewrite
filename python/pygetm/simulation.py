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

import xarray

from .constants import *
from . import _pygetm
from . import core
from . import domain
from . import parallel
from . import output
from . import operators
from . import tracer
import pygetm.input
import pygetm.mixing
import pygetm.density
import pygetm.airsea
import pygetm.radiation
import pygetm.fabm

BAROTROPIC = BAROTROPIC_2D = 1
BAROTROPIC_3D = 2
FROZEN_DENSITY = 3
BAROCLINIC = 4

class InternalPressure(enum.IntEnum):
    OFF = 0                    #: internal pressure is disabled
    BLUMBERG_MELLOR = 1        #: Blumberg and Mellor
#    BLUMBERG_MELLOR_LIN=2
#    Z_INTERPOL=3
#    SONG_WRIGHT=4
#    CHU_FAN=5
    SHCHEPETKIN_MCWILLIAMS=6   #: Shchepetkin and McWilliams (2003)
#    STELLING_VANKESTER=7

class Simulation(_pygetm.Simulation):
    _momentum_arrays = 'U', 'V', 'fU', 'fV', 'advU', 'advV', 'diffu1', 'diffv1', 'u1', 'v1', 'uk', 'vk', 'ru', 'rru', 'rv', 'rrv', 'pk', 'qk', 'ww', 'advpk', 'advqk', 'diffpk', 'diffqk', 'Ui', 'Vi', 'SS', 'fpk', 'fqk', 'ustar2_s', 'ustar2_b', 'SxB', 'SyB', 'SxA', 'SyA', 'SxD', 'SyD', 'SxF', 'SyF'
    _pressure_arrays = 'dpdx', 'dpdy', 'idpdx', 'idpdy'
    _sealevel_arrays = ()
    _time_arrays = 'timestep', 'macrotimestep', 'split_factor', 'timedelta', 'time', 'istep', 'report', 'report_totals', 'default_time_reference'
    _all_fortran_arrays = tuple(['_%s' % name for name in _momentum_arrays + _pressure_arrays + _sealevel_arrays]) + ('uadv', 'vadv', 'uua', 'uva', 'vua', 'vva', 'uua3d', 'uva3d', 'vua3d', 'vva3d')
    __slots__ = _all_fortran_arrays + ('output_manager', 'input_manager', 'fabm', '_yearday', 'tracers', 'tracer_totals', 'logger', 'airsea', 'turbulence', 'density', 'buoy', 'temp', 'salt', 'pres', 'rad', 'par', 'par0', 'rho', 'sst', 'sss', 'NN', 'ustar_s', 'ustar_b', 'taub', 'z0s', 'z0b', 'fwf', '_cum_river_height_increase', '_start_time', '_profile', 'radiation', '_ufirst', '_u3dfirst', 'diffuse_momentum', 'apply_bottom_friction', 'Ui_tmp', 'Vi_tmp', '_initialized_variables') + _time_arrays

    _array_args = {
        'U': dict(units='m2 s-1', long_name='depth-integrated transport in Eastward direction', fill_value=FILL_VALUE, attrs={'_part_of_state': True, '_mask_output': True}),
        'V': dict(units='m2 s-1', long_name='depth-integrated transport in Northward direction', fill_value=FILL_VALUE, attrs={'_part_of_state': True, '_mask_output': True}),
        'Ui': dict(units='m2 s-1', fill_value=FILL_VALUE, attrs={'_part_of_state': True, '_mask_output': True}),
        'Vi': dict(units='m2 s-1', fill_value=FILL_VALUE, attrs={'_part_of_state': True, '_mask_output': True}),
        'u1': dict(units='m s-1', long_name='depth-averaged velocity in Eastward direction', fill_value=FILL_VALUE, attrs={'_mask_output': True}),
        'v1': dict(units='m s-1', long_name='depth-averaged velocity in Northward direction', fill_value=FILL_VALUE, attrs={'_mask_output': True}),
        'pk': dict(units='m2 s-1', long_name='layer-integrated transport in Eastward direction', fill_value=FILL_VALUE, attrs={'_part_of_state': True, '_mask_output': True}),
        'qk': dict(units='m2 s-1', long_name='layer-integrated transport in Northward direction', fill_value=FILL_VALUE, attrs={'_part_of_state': True, '_mask_output': True}),
        'uk': dict(units='m s-1', long_name='velocity in Eastward direction', fill_value=FILL_VALUE, attrs={'_mask_output': True}),
        'vk': dict(units='m s-1', long_name='velocity in Northward direction', fill_value=FILL_VALUE, attrs={'_mask_output': True}),
        'ww': dict(units='m s-1', long_name='vertical velocity', fill_value=FILL_VALUE),
        'SS': dict(units='s-2', long_name='shear frequency squared', fill_value=FILL_VALUE)
    }

    def __init__(self, dom: domain.Domain, runtype: int, advection_scheme: operators.AdvectionScheme=operators.AdvectionScheme.HSIMT, apply_bottom_friction: bool=True, fabm: Union[pygetm.fabm.FABM, bool, str, None]=None, gotm: Union[str, None]=None,
        turbulence: Optional[pygetm.mixing.Turbulence]=None, airsea: Optional[pygetm.airsea.Fluxes]=None, density: Optional[pygetm.density.Density]=None,
        logger: Optional[logging.Logger]=None, log_level: int=logging.INFO, internal_pressure_method: InternalPressure=InternalPressure.OFF, Am: float=0.):

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
            self.logger.info('Diffusion of momentum if off because Am is 0')

        assert not dom._initialized
        super().__init__(dom, runtype, internal_pressure_method=internal_pressure_method, Am0=Am)
        dom.T.hn.fabm_standard_name = 'cell_thickness'
        if dom.T.lon is not None:
            dom.T.lon.fabm_standard_name = 'longitude'
        if dom.T.lat is not None:
            dom.T.lat.fabm_standard_name = 'latitude'
        dom.T.z.attrs['_part_of_state'] = True
        dom.T.zo.attrs['_part_of_state'] = True
        if self.runtype > BAROTROPIC_2D:
            dom.T.zio.attrs['_part_of_state'] = True
            dom.T.zin.attrs['_part_of_state'] = True
            dom.T.ho.attrs['_part_of_state'] = True    # ho cannot be computed from zio, because rivers modify ho-from-zio before it is stored
        if self.runtype == BAROTROPIC_2D:
            dom.U.z0b.attrs['_part_of_state'] = True
            dom.V.z0b.attrs['_part_of_state'] = True

        self.domain.open_boundaries.z = self.wrap(self.domain.open_boundaries.z, b'zbdy', source=3)
        self.domain.open_boundaries.u = self.wrap(self.domain.open_boundaries.u, b'bdyu', source=1)
        self.domain.open_boundaries.v = self.wrap(self.domain.open_boundaries.v, b'bdyv', source=1)
        for name in Simulation._momentum_arrays:
            setattr(self, '_%s' % name, self.wrap(core.Array(name=name, **Simulation._array_args.get(name, {})), name.encode('ascii'), source=1))
        for name in Simulation._pressure_arrays:
            setattr(self, '_%s' % name, self.wrap(core.Array(name=name, **Simulation._array_args.get(name, {})), name.encode('ascii'), source=2))
        for name in Simulation._sealevel_arrays:
            setattr(self, '_%s' % name, self.wrap(core.Array(name=name, **Simulation._array_args.get(name, {})), name.encode('ascii'), source=3))

        self.U.all_values.fill(0.)
        self.V.all_values.fill(0.)
        self.u1.all_values.fill(0.)
        self.v1.all_values.fill(0.)
        self.Ui.all_values.fill(0.)
        self.Vi.all_values.fill(0.)
        self.Ui_tmp = numpy.zeros_like(self.Ui.all_values)
        self.Vi_tmp = numpy.zeros_like(self.Vi.all_values)
        if runtype > BAROTROPIC_2D:
            self.pk.all_values.fill(0.)
            self.qk.all_values.fill(0.)
            self.uk.all_values.fill(0.)
            self.vk.all_values.fill(0.)
            self.ww.all_values.fill(0.)
            self.SS.fill(0.)             # for surface/bottom interfaces, which are not updated

        self._cum_river_height_increase = numpy.zeros((len(self.domain.rivers),))

        #: Provider of air-water fluxes of heat and momentum. This must inherit from :class:`pygetm.airsea.Fluxes` and should be provided as argument airsea to :class:`Simulation`.
        self.airsea = airsea or pygetm.airsea.FluxesFromMeteo()
        assert isinstance(self.airsea, pygetm.airsea.Fluxes), 'airsea argument should be of type pygetm.airsea.Fluxes, but is %s' % type(self.airsea)
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

        #: Collection of tracers that are to be transported. Optionally they can have sources, open boundary conditions and riverine concentrations set.
        self.tracers: tracer.TracerCollection = tracer.TracerCollection(self.domain.T, advection_scheme=advection_scheme)

        #: List of variables for which the domain-integrated total needs to be reported. These can be depth-integrated (2D) or depth-explicit (3D).
        self.tracer_totals: List[core.Array] = []

        self._ufirst = False   #: Whether to start the depth-integrated (2D) momentum update with u (as opposed to v)
        self._u3dfirst = False #: Whether to start the depth-explicit (3D) momentum update with u (as opposed to v)

        self.fabm = None

        if runtype > BAROTROPIC_2D:
            #: Provider of turbulent viscosity and diffusivity. This must inherit from :class:`pygetm.mixing.Turbulence` and should be provided as argument turbulence to :class:`Simulation`.
            self.turbulence = turbulence or pygetm.mixing.GOTM(nml_path=gotm)
            self.turbulence.initialize(self.domain.T)
            self.NN = dom.T.array(fill=0., z=INTERFACES, name='NN', units='s-2', long_name='buoyancy frequency squared', fill_value=FILL_VALUE)
            self.ustar_s = dom.T.array(fill=0., name='ustar_s', units='m s-1', long_name='shear velocity (surface)', fill_value=FILL_VALUE)
            self.ustar_b = dom.T.array(fill=0., name='ustar_b', units='m s-1', long_name='shear velocity (bottom)', fill_value=FILL_VALUE)
            self.z0s = dom.T.array(fill=0.1, name='z0s', units='m', long_name='hydrodynamic roughness (surface)', fill_value=FILL_VALUE)
            self.taub = dom.T.array(fill=0., name='taub', units='Pa', long_name='bottom shear stress', fill_value=FILL_VALUE, fabm_standard_name='bottom_stress')

            if fabm:
                if not isinstance(fabm, pygetm.fabm.FABM):
                    fabm = pygetm.fabm.FABM(fabm if isinstance(fabm, str) else 'fabm.yaml')
                self.fabm = fabm
                self.fabm.initialize(self.domain.T, self.tracers, self.tracer_totals, self.logger.getChild('FABM'))

            self.pres = dom.depth
            self.pres.fabm_standard_name = 'pressure'

        self.sst = dom.T.array(name='sst', units='degrees_Celsius', long_name='sea surface temperature', fill_value=FILL_VALUE)

        if runtype == BAROCLINIC:
            self.radiation = pygetm.radiation.TwoBand(dom.T)
            self.temp = self.tracers.add(name='temp', units='degrees_Celsius', long_name='conservative temperature', fabm_standard_name='temperature', fill_value=FILL_VALUE,
                source=self.radiation.swr_abs, surface_flux=self.airsea.shf, source_scale=1. / (RHO0 * CP), rivers_follow_target_cell=True, precipitation_follows_target_cell=True)
            self.salt = self.tracers.add(name='salt', units='-', long_name='absolute salinity', fabm_standard_name='practical_salinity', fill_value=FILL_VALUE)
            self.pres.saved = True
            self.temp.fill(5.)
            self.salt.fill(35.)
            self.rho = dom.T.array(z=CENTERS, name='rho', units='kg m-3', long_name='density', fabm_standard_name='density')
            self.buoy = dom.T.array(z=CENTERS, name='buoy', units='m s-2', long_name='buoyancy')
            self.tracer_totals.append(self.salt)
            self.sss = self.salt.isel(-1)

            self.density = density or pygetm.density.Density()
        else:
            self.sss = None

        # Derive old and new elevations, water depths and thicknesses from current surface elevation on T grid
        # This must be done after self.pres.saved is set
        self.domain.update_depth(_3d=runtype > BAROTROPIC_2D)
        self.domain.update_depth(_3d=runtype > BAROTROPIC_2D)

        self.default_time_reference: Optional[cftime.datetime] = None
        self._initialized_variables = set()

    def __getitem__(self, key: str) -> core.Array:
        return self.output_manager.fields[key]

    def load_restart(self, path: str, time: Optional[cftime.datetime]=None, **kwargs):
        """Load the model state from a restart file."""
        kwargs.setdefault('decode_times', True)
        kwargs['use_cftime'] = True
        with xarray.open_dataset(path, **kwargs) as ds:
            timevar = ds['zt'].getm.time
            time_coord = timevar.values
            self.default_time_reference = cftime.num2date(0., units=timevar.encoding['units'], calendar=timevar.encoding['calendar'])
            assert time_coord.size == 1, 'Currently a restart file should contain a single time point only'
        with xarray.open_dataset(path, **kwargs) as ds:
            for name, field in self.output_manager.fields.items():
                if field.attrs.get('_part_of_state', False):
                    if name not in ds:
                        raise Exception('Field %s is part of state but not found in %s' % path)
                    field.set(ds[name], on_grid=pygetm.input.OnGrid.ALL, mask=True)
                    self._initialized_variables.add(name)
        if self.runtype > BAROTROPIC_2D:
            # Restore elevation from before open boundary condition was applied
            self.domain.T.z.all_values[...] = self.domain.T.zin.all_values
        return time_coord[0]

    def start(self, time: Union[cftime.datetime, datetime.datetime], timestep: float, split_factor: int=1, report: Union[int, datetime.timedelta]=10, report_totals: Union[int, datetime.timedelta]=datetime.timedelta(days=1), save: bool=True, profile: Optional[str]=None):
        """Start a simulation by configuring the time, zeroing velocities, updating diagnostics to match the start time, and optionally saving output.

        This should be called after the output configuration is complete (because we need to know when variables need to be saved),
        and after the FABM model has been provided with all dependencies.
        
        Args:
            time (:class:`cftime.datetime`): start time
            timestep: micro time step (s) used for 2D barotropic processes
            split_factor: number of microtimesteps per macrotimestep
            report: time interval or number of microtimesteps between reporting of the current time, used as indicator of simulation progress
            report_totals: time interval or number of microtimesteps between reporting of integrals over the global domain
            save: whether to save the model state and diagnostics at the very start of the simulation
            profile: base name for the file to write profiling results to. The proces rank and extension ``.prof`` will be appended,
                so that the final name becomes ``<profile>-<rank>.prof``. If the argument is not provided, profiling is disabled.
        """
        if isinstance(time, datetime.datetime):
            time = cftime.datetime(time.year, time.month, time.day, time.hour, time.minute, time.second, time.microsecond)
        self.logger.info('Starting simulation at %s' % time)
        self.timestep = timestep
        self.split_factor = split_factor
        self.macrotimestep = self.timestep * self.split_factor
        self.timedelta = datetime.timedelta(seconds=timestep)
        self.time = time
        self.istep = 0
        if isinstance(report, datetime.timedelta):
            report = int(round(report.total_seconds() / self.timestep))
        self.report = report
        if isinstance(report_totals, datetime.timedelta):
            report_totals = int(round(report_totals.total_seconds() / self.timestep))
        self.report_totals = report_totals

        # Ensure transports and velocities are 0 in masked points
        # NB velocities will be computed from transports, but only in unmasked points, so zeroing them here is needed.
        zero_masked = [self.U, self.V, self.u1, self.v1]
        if self.runtype > BAROTROPIC_2D:
            zero_masked += [self.Ui, self.Vi, self.pk, self.qk, self.uk, self.vk, self.ww]
        for array in zero_masked:
            array.all_values[..., array.grid.mask.all_values == 0] = 0.

        if self.fabm:
            self.fabm.start(self.time)

        # First (out of two) 2D depth update based on old elevations zo
        z_backup = self.domain.T.z.all_values.copy()
        self.domain.T.z.all_values[...] = self.domain.T.zo.all_values
        self.domain.update_depth(_3d=False)

        if self.runtype > BAROTROPIC_2D:
            zin_backup = self.domain.T.zin.all_values.copy()
            ho_T_backup = self.domain.T.ho.all_values.copy()

            # First (out of two) 3D depth/thickness update based on zio.
            # This serves to generate T.ho when T.zio is set, but T.ho is not available.
            # Since we do not have the preceding (2 time steps before start) zi/h, we explicitly set them (here: T.zio/T.ho)
            # to NaN to make it easier to detect algorithms depending on them.
            # As a result of that, all new metrics on the U, V, X grids will be NaN too!
            self.domain.T.z.all_values[...] = self.domain.T.zio.all_values  # to become T.zin when update_depth is called
            self.domain.T.zio.fill(numpy.nan)
            self.domain.T.ho.fill(numpy.nan)
            self.domain.update_depth(_3d=True)

            # Second 3D depth/thickness update based on zin.
            # Override T.ho with user-provided value if available, since this may incorporate river inflow impacts that
            # our previously calculated ho cannot account for.
            # New metrics for U, V, X grids will be calculated from valid old and new metrics on T grid; therefore they will be valid too.
            # However, old metrics (ho/zio) for U, V, X grids will still be NaN and should not be used.
            self.domain.T.z.all_values[...] = zin_backup  # to become T.zin when update_depth is called
            if 'hot' in self._initialized_variables:
                self.domain.T.hn.all_values[...] = ho_T_backup
            self.domain.update_depth(_3d=True)  # this moves our zin backup into zin, and at the same time moves the current zin (originally zio) to zio
            self.update_3d_momentum_diagnostics(self.macrotimestep, self.turbulence.num)

        # Update all forcing, which includes the final 2D depth update based on (original) z
        self.domain.T.z.all_values[...] = z_backup
        self.update_forcing(macro_active=self.runtype > BAROTROPIC_2D)

        # Start output manager
        self.output_manager.start(self.istep, self.time, save=save, default_time_reference=self.default_time_reference)

        # Record true start time for performance analysis
        self._start_time = timeit.default_timer()

        # Start profiling if requested
        self._profile = None
        if profile:
            import cProfile
            pr = cProfile.Profile()
            self._profile = (profile, pr)
            pr.enable()

    def advance(self):
        """Advance the model state by one microtimestep.
        If this completes the current macrotimestep, the part of the state associated with that timestep will be advanced too.
        """

        # Update transports U and V from time=-1/2 to +1/2, using surface stresses and pressure gradients defined at time=0
        # Inputs and outputs are on U and V grids. Stresses and pressure gradients have already been updated by the call to
        # update_forcing at the end of the previous time step.
        self.advance_2d_momentum(self.timestep, self.airsea.taux_U, self.airsea.tauy_V, self.dpdx, self.dpdy)

        # Update surface elevation on T grid from time=0 to time=1 using transports U and V at time=1/2 and
        # freshwater fluxes at time=0. This also updates halos so that depths and thicknesses can be computed
        # everywhere without further halo exchange
        self.advance_surface_elevation(self.timestep, self.U, self.V, self.fwf)

        # Track cumulative increase in elevation due to river inflow over the current macrotimestep
        self._cum_river_height_increase += self.domain.rivers.flow * self.domain.rivers.iarea * self.timestep

        # Update the time
        self.time += self.timedelta
        self.istep += 1
        macro_active = self.istep % self.split_factor == 0
        if self.report != 0 and self.istep % self.report == 0:
            self.logger.info(self.time)

        if self.runtype > BAROTROPIC_2D and macro_active:
            # Use previous source terms for biogeochemistry (valid for the start of the current macrotimestep)
            # to update tracers. This should be done before the tracer concentrations change due to transport
            # or rivers, as the source terms are only valid for the current tracer concentrations.
            if self.fabm:
                self.fabm.advance(self.macrotimestep)

            # Update layer thicknesses and tracer concentrations to account for precipitation, evaporation and river inflow
            # between start and end of the current macrotimestep.
            self.add_freshwater_inputs(self.macrotimestep)

            # Update 3D elevations and layer thicknesses. New elevation zin on T grid will match elevation at end of the microtimestep,
            # thicknesses on T grid will match. Elevation and thicknesses on U/V grids will be 1/2 MACROtimestep behind, as they
            # are calculated by averaging zio and zin. Old elevations zio and thicknesses ho will store the previous values of zin and hn,
            # and thus are one macrotimestep behind new elevations zin and thicknesses hn.
            self.domain.update_depth(_3d=True)

            # Update presssure gradient for start of the 3D time step
            # Note that we use previously-recorded elevations and surface pressure at the start of the current macrotimestep
            self.update_surface_pressure_gradient(self.domain.T.zio, self.airsea.spo)

            # Update momentum from time=-1/2 to 1/2 of the macrotimestep, using forcing defined at time=0
            # For this purpose, surface stresses at the end of the previous macrotimestep were saved (taux_Uo, tauy_Vo)
            # Pressure gradients dpdx and dpdy have just been updated to match the start of the current macrotimestep
            # Internal pressure idpdx and idpdy were calculated at the end of the previous macrotimestep and are therefore ready as-is.
            self.advance_3d_momentum(self.macrotimestep, self.split_factor, self.airsea.taux_Uo, self.airsea.tauy_Vo, self.dpdx, self.dpdy, self.idpdx, self.idpdy, self.turbulence.num)

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
                self.turbulence.advance(self.macrotimestep, self.ustar_s, self.ustar_b, self.z0s, self.domain.T.z0b, self.NN, self.SS)

                # Advect and diffuse tracers. Source terms are optionally handled too, as part of the diffusion update.
                self.tracers.advance(self.macrotimestep, self.uk, self.vk, self.ww, self.turbulence.nuh)

        if self.report_totals != 0 and self.istep % self.report_totals == 0:
            self.report_domain_integrals()

        # Update all inputs and fluxes that will drive the next state update
        self.update_forcing(macro_active, skip_2d_coriolis=True, update_z0b=self.runtype==BAROTROPIC_2D)

        self.output_manager.save(self.timestep * self.istep, self.istep, self.time)

        return macro_active

    def update_forcing(self, macro_active: bool, skip_2d_coriolis: bool=False, update_z0b: bool=False):
        """Update all inputs and fluxes that will drive the next state update.

        Args:
            macro_active: update all quantities associated with the macrotimestep
            skip_2d_coriolis: whether to skip the update of depth-integrated Coriolis terms,
                typically because they have already been recalculated as part of the transport update
        """
        # Update all inputs.
        self.domain.input_manager.update(self.time, include_3d=macro_active)

        if self.runtype == BAROCLINIC and macro_active:
            # Update tracer values at open boundaries. This must be done after input_manager.update,
            # but before diagnostics/forcing variables derived from the tracers are calculated
            if self.domain.open_boundaries:
                for tracer in self.tracers:
                    tracer.open_boundaries.update()

            # Update density, buoyancy and internal pressure to keep them in sync with T and S.
            self.density.get_density(self.salt, self.temp, p=self.pres, out=self.rho)
            self.rho.update_halos(parallel.Neighbor.LEFT_AND_RIGHT_AND_TOP_AND_BOTTOM)   # valid rho around all U and V needed for internal pressure; not yet valid because T&S were not valid in halos when rho was calculated. Note BM needs only right/top, SMcW needs left/right/top/bottom
            self.buoy.all_values[...] = (-GRAVITY / RHO0) * (self.rho.all_values - RHO0)
            self.update_internal_pressure_gradient(self.buoy, self.SxB, self.SyB)

            # From conservative temperature to in-situ sea surface temperature,
            # needed to compute heat/momentum fluxes at the surface
            self.density.get_potential_temperature(self.sss, self.temp.isel(-1), out=self.sst)

            # Calculate squared buoyancy frequency NN (T grid, interfaces between layers)
            self.density.get_buoyancy_frequency(self.salt, self.temp, p=self.pres, out=self.NN)

        # Update surface elevation z on U, V, X grids and water depth D on all grids
        # This is based on old and new elevation (T grid) for the microtimestep.
        # Thus, for grids lagging 1/2 a timestep behind (U, V, X grids), the elevations
        # and water depths will be representative for 1/2 a MICROtimestep ago.
        # Note that T grid elevations at the open boundary have not yet been updated,
        # so the derived elevations and water depths calculated here will not take those into account.
        # This is intentional: it ensures that water depths on the U and V grids are in sync with the already-updated transports,
        # so that velocities can be calculated correctly.
        # The call to update_surface_elevation_boundaries is made later.
        self.domain.update_depth()

        # Calculate advection and diffusion tendencies of transports, bottom friction and, if needed, Coriolis terms
        self.update_2d_momentum_diagnostics(self.timestep, skip_coriolis=skip_2d_coriolis, update_z0b=update_z0b)

        # Update air-sea fluxes of heat and momentum (T grid for all, U and V grid for x and y stresses respectively)
        # Note SST is the true in-situ/potential temperature. SSS currently is absolute salinity - not practical salinity.
        self.airsea(self.time, self.sst, self.sss, calculate_heat_flux=macro_active and self.runtype == BAROCLINIC)

        # Update depth-integrated freshwater fluxes (precipitation, evaporation, condensation, rivers)
        self.fwf.all_values[...] = self.airsea.pe.all_values
        self.fwf.all_values[self.domain.rivers.j, self.domain.rivers.i] += self.domain.rivers.flow * self.domain.rivers.iarea

        # Update elevation at the open boundaries. This must be done before update_surface_pressure_gradient
        self.update_surface_elevation_boundaries(self.timestep)

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
            if self.fabm:
                self.fabm.update_sources(self.time)

            # Save forcing variables for the next baroclinic update
            self.airsea.spo.all_values[...] = self.airsea.sp.all_values
            self.airsea.taux_Uo.all_values[...] = self.airsea.taux_U.all_values
            self.airsea.tauy_Vo.all_values[...] = self.airsea.tauy_V.all_values

    def finish(self):
        """Clean-up after simulation: save profiling result (if any), write output where appropriate (restarts), and close output files"""
        if self._profile:
            import pstats
            name, pr = self._profile
            pr.disable()
            profile_path = '%s-%03i.prof' % (name, self.domain.tiling.rank)
            self.logger.info('Writing profiling report to %s' % (profile_path,))
            with open(profile_path, 'w') as f:
                ps = pstats.Stats(pr, stream=f).sort_stats(pstats.SortKey.TIME)
                ps.print_stats()
        self.logger.info('Time spent in main loop: %.3f s' % (timeit.default_timer() - self._start_time,))
        self.output_manager.close(self.timestep * self.istep, self.time)

    def add_freshwater_inputs(self, timestep: float):
        """Update layer thicknesses and tracer concentrations to account for precipitation, evaporation and river inflow."""

        # Precipitation and evaporation (surface layer only)
        # First update halos for the net freshwater flux, as we need to ensure that the layer heights updated as a result remain valid in the halos.
        self.airsea.pe.update_halos()
        h_increase_pe = self.airsea.pe.all_values * timestep
        h = self.domain.T.hn.all_values[-1, :, :]
        h_new = h + h_increase_pe
        for tracer in self.tracers:
            if not tracer.precipitation_follows_target_cell:
                tracer.all_values[-1, :, :] *= h / h_new
        h[:, :] = h_new
        self.domain.T.zin.all_values += h_increase_pe

        # Rivers, potentially entering anywhere in the water column
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

    def report_domain_integrals(self):
        """Write totals of selected variables over the global domain (those in :attr:`tracer_totals`) to the log."""
        total_volume = (self.domain.T.D * self.domain.T.area).global_sum(where=self.domain.T.mask != 0)
        if total_volume is not None:
            self.logger.info('Integrals over global domain:')
            self.logger.info('  volume: %.15e m3' % total_volume)
        if self.fabm:
            self.fabm.update_totals()
        for var in self.tracer_totals:
            total = var * var.grid.area
            if total.ndim == 3:
                total.all_values *= var.grid.hn.all_values
            total = total.global_sum(where=var.grid.mask != 0)
            if total is not None:
                self.logger.info('  %s: %.15e %s m3 (per volume: %s %s)' % (var.name, total, var.units, total / total_volume, var.units))

    def advance_surface_elevation(self, timestep: float, U: core.Array, V: core.Array, fwf: core.Array):
        """Advance surface elevation (T grid only)

        Args:
            timestep: time step (s)
            U: depth-integrated velocity in x direction (m2 s-1)
            V: depth-integrated velocity in y direction (m2 s-1)
            fwf: freshwater flux (m s-1)

        This also updates the surface elevation halos.
        This method does `not` update elevation on the U, V, X grids, nor water depths, layer thicknesses or vertical coordinates.
        This is done by :meth:`~pygetm.domain.Domain.update_depth` instead.
        """
        super().advance_surface_elevation(timestep, U, V, fwf)
        self.domain.T.z.update_halos()

    def transport_2d_momentum(self, U: core.Array, V: core.Array, timestep: float, advU: core.Array, advV: core.Array, diffU: core.Array, diffV: core.Array, update_z0b: bool):
        """Advect and optionally diffuse depth-integrated transports in x and y direction (arguments ``U`` and ``V``).
        From these, first the depth-averaged velocities are calculated and stored in :attr:`u1` and :attr:`v1`.
        This routine also updates bottom friction :attr:`ru` and  :attr:`rv`.

        Args:
            U: depth-integrated velocity (m2 s-1) in x direction
            V: depth-integrated velocity (m2 s-1) in y direction
            timestep: time step (s) to calculate advection over
            advU: array for storing the change in transport ``U`` due to advection (m2 s-2)
            advV: array for storing the change in transport ``V`` due to advection (m2 s-2)
            diffU: array for storing the change in transport ``U`` due to diffusion (m2 s-2)
            diffV: array for storing the change in transport ``V`` due to diffusion (m2 s-2)
            update_z0b: whether to iteratively update hydrodynamic bottom roughness
        """
        numpy.divide(U.all_values, U.grid.D.all_values, out=self.u1.all_values)
        numpy.divide(V.all_values, V.grid.D.all_values, out=self.v1.all_values)

        if self.diffuse_momentum:
            # Compute velocity diffusion contribution to transport sources.
            # This uses depth-averaged velocities u1 and v1, which therefore have to be up to date
            # Water depths should be in sync with velocities, which means they should lag 1/2 a timestep behind the tracer/T grid
            self.momentum_diffusion_driver(self.domain.D_T_half, self.domain.X.D, self.u1, self.v1, diffU, diffV)

        # Calculate bottom friction (ru and rv) using updated depth-averaged velocities u1 and v1
        # Warning: this uses velocities u1 and v1 at masked points, which therefore need to be kept at 0
        if self.apply_bottom_friction:
            self.bottom_friction_2d(update_z0b)

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

        # Restore depth-averaged velocities as they need to be valid on exit
        numpy.divide(U.all_values, U.grid.D.all_values, out=self.u1.all_values)
        numpy.divide(V.all_values, V.grid.D.all_values, out=self.v1.all_values)

    def advance_2d_momentum(self, timestep: float, tausx: core.Array, tausy: core.Array, dpdx: core.Array, dpdy: core.Array):
        """Update depth-integrated transports (:attr:`U`, :attr:`V`) and depth-averaged velocities (:attr:`u1`, :attr:`v1`).
        This will also update their halos.
        
        Args:
            timestep: time step (s)
            tausx: surface stress (Pa) in x direction
            tausy: surface stress (Pa) in y direction
            dpdx: surface pressure gradient (dimensionless) in x direction
            dpdx: surface pressure gradient (dimensionless) in y direction
        """
        # Update 2D transports from t-1/2 to t+1/2. This uses advection, diffusion and bottom friction terms
        # (advU, advV, diffu1, diffv1, ru, rv) that were calculated by the call to update_2d_momentum_diagnostics
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

        self.Ui_tmp += self.U.all_values
        self.Vi_tmp += self.V.all_values

    def update_2d_momentum_diagnostics(self, timestep: float, skip_coriolis: bool=False, update_z0b: bool=False):
        """Update 2D momentum diagnostics, including the Coriolis terms that will drive the next 2D update.
        NB the Coriolis update is already done as part of the momentum update itself, so needed only when starting from a restart.
        
        Args:
            timestep: time step (s) to calculate advection of momentum over
            skip_coriolis: flag to indicate that Coriolis terms are already up-to-date and do not need recomputing
        """
        if not skip_coriolis:
            self.coriolis_fu()
            self.coriolis_fv()

        # Calculate sources of transports U and V due to advection (advU, advV) and diffusion (diffu1, diffv1)
        # Transports generally come in at time=-1/2 and are then advanced to time+1/2
        self.transport_2d_momentum(self.U, self.V, timestep, self.advU, self.advV, self.diffu1, self.diffv1, update_z0b)

    def advance_3d_momentum(self, timestep: float, split_factor: int, tausx: core.Array, tausy: core.Array, dpdx: core.Array, dpdy: core.Array, idpdx: core.Array, idpdy: core.Array, viscosity: core.Array):
        """Update depth-explicit transports (:attr:`pk`, :attr:`qk`) and velocities (:attr:`uk`, :attr:`vk`).
        This will also update their halos.
        
        Args:
            timestep: (macro) time step (s)
            split_factor: number of microtimesteps per macrotimestep
            tausx: surface stress (Pa) in x direction
            tausy: surface stress (Pa) in y direction
            dpdx: surface pressure gradient (dimensionless) in x direction
            dpdy: surface pressure gradient (dimensionless) in y direction
            idpdx: internal pressure gradient (m2 s-2) in x direction
            idpdy: internal pressure gradient (m2 s-2) in y direction
            viscosity: turbulent viscosity (m2 s-1)
        """
        # Depth-integrated transports have been summed over all microtimesteps.
        # Average them, then reset depth-integrated transports that will be incremented over the next macrotimestep.
        numpy.multiply(self.Ui_tmp, 1. / split_factor, out=self.Ui.all_values)
        numpy.multiply(self.Vi_tmp, 1. / split_factor, out=self.Vi.all_values)
        self.Ui_tmp.fill(0.)
        self.Vi_tmp.fill(0.)

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

        self.update_3d_momentum_diagnostics(timestep, viscosity, skip_coriolis=True)

    def update_3d_momentum_diagnostics(self, timestep: float, viscosity: core.Array, skip_coriolis: bool=False):
        """Update 3D momentum diagnostics, including the vertical velocity :attr:`ww`,
        the slow terms that will drive the 2D updates over the next macrotimestep,
        and the bottom friction and Coriolis terms that will drive the next 3D update.
        NB the Coriolis update is already done as part of the momentum update itself,
        so needed only when starting from a restart.
        
        Args:
            timestep: time step (s)
            viscosity: turbulent viscosity (T grid, layer interfaces, m2 s-1)
            skip_coriolis: flag to indicate that Coriolis terms are already up-to-date and do not need recomputing
        """
        if not skip_coriolis:
            self.coriolis_fpk()
            self.coriolis_fqk()

        # Infer vertical velocity from horizontal transports and desired layer height change (ho -> hn).
        # This is done at all points surrounding U and V points, so no further halo exchange of w is needed
        # to support interpolation to U and V grids later on. This does require that transports are up to date in halos.
        self.w_3d(timestep)

        itimestep = 1. / timestep

        # Compute 3D velocities (m s-1) from 3D transports (m2 s-1) by dividing by layer heights
        # Both velocities and U/V thicknesses are now at time 1/2
        numpy.divide(self.pk.all_values, self.pk.grid.hn.all_values, where=self.pk.grid.mask.all_values != 0, out=self.uk.all_values)
        numpy.divide(self.qk.all_values, self.qk.grid.hn.all_values, where=self.qk.grid.mask.all_values != 0, out=self.vk.all_values)

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

        # Advect 3D u and v velocity from time=1/2 to 1 1/2 using velocities interpolated to its own advection grids
        # Store the resulting trend, which will be applied as part of the momentum update in the next timestep.
        # They will also be used to calculate the slow advection contribution to depth-integrated momentum equations.
        # JB the alternative would be to interpolate transports and then divide by (colocated) layer heights, like we do for 2D
        self.uadv.apply_3d(self.uua3d, self.uva3d, self.ww.interp(self.uk.grid), timestep, self.uk, new_h=True, skip_initial_halo_exchange=True)
        self.advpk.all_values[...] = (self.uk.all_values * self.uadv.h - self.pk.all_values) * itimestep
        self.vadv.apply_3d(self.vua3d, self.vva3d, self.ww.interp(self.vk.grid), timestep, self.vk, new_h=True, skip_initial_halo_exchange=True)
        self.advqk.all_values[...] = (self.vk.all_values * self.vadv.h - self.qk.all_values) * itimestep

        # Restore velocity at time=1/2 (the final value at the end of the current timestep)
        numpy.divide(self.pk.all_values, self.pk.grid.hn.all_values, where=self.pk.grid.mask.all_values != 0, out=self.uk.all_values)
        numpy.divide(self.qk.all_values, self.qk.grid.hn.all_values, where=self.qk.grid.mask.all_values != 0, out=self.vk.all_values)

        if self.diffuse_momentum:
            # Calculate the momentum trends (diffpk, diffqk) associated with diffusion of 3D u and v velocity
            # between time=1/2 to 1 1/2. Note that thicknesses should be in sync with velocities uk and vk
            # This means they should lag 1/2 a timestep behind the T grid (already the case for X, but for T we use 1/2(ho+hn))
            self.momentum_diffusion_driver(self.domain.h_T_half, self.domain.X.hn, self.uk, self.vk, self.diffpk, self.diffqk)

        # Compute slow (3D) advection and diffusion contribution to to the depth-integrated momentum equations.
        # This is done by comparing the depth-integrated 3D transport calculated above
        # (between centers of the current and next macrotime step) with the newly calculated
        # depth-integrated transport based on accumulated 2D transports (accumulated over the
        # current macrotimestep, and thus representative for its center).
        self.transport_2d_momentum(self.Ui, self.Vi, timestep, self.SxA, self.SyA, self.SxD, self.SyD, False)
        self.SxA.all_values[...] = self.advpk.all_values.sum(axis=0) - self.SxA.all_values
        self.SyA.all_values[...] = self.advqk.all_values.sum(axis=0) - self.SyA.all_values
        self.SxD.all_values[...] = self.diffpk.all_values.sum(axis=0) - self.SxD.all_values
        self.SyD.all_values[...] = self.diffqk.all_values.sum(axis=0) - self.SyD.all_values

        if self.apply_bottom_friction:
            # Note: ru and rv have been updated by transport_2d_momentum, using accumulated transports Ui and Vi (representative for t=1/2, just like uk, vk, rru, rrv)
            # Slow bottom friction (stress/density) is derived by taking the difference between 3D bottom friction and the inferred ru and rv.
            self.SxF.all_values[...] = self.rru.all_values * self.uk.all_values[0, ...] - self.ru.all_values * self.u1.all_values
            self.SyF.all_values[...] = self.rrv.all_values * self.vk.all_values[0, ...] - self.rv.all_values * self.v1.all_values

    @property
    def Ekin(self, rho0: float=RHO0):
        dom = self.domain
        U = self.U.interp(dom.T)
        V = self.V.interp(dom.T)
        vel2_D2 = U**2 + V**2
        return 0.5 * rho0 * dom.T.area * vel2_D2 / dom.T.D

# Expose all Fortran arrays that are a member of Simulation as read-only properties
# The originals are members with and underscore as prefix, ad therefore not visible to the user
# This ensures the user will not accidentally disconnect the Python variable from th underlying Fortran libraries/data
for membername in Simulation._all_fortran_arrays:
    attrs = Simulation._array_args.get(membername[1:], {})
    long_name = attrs.get('long_name')
    units = attrs.get('units')
    doc = ''
    if long_name is not None:
        doc = long_name
        if units:
           doc += ' (%s)' % units 
    prop = property(operator.attrgetter(membername), doc=doc)
    setattr(Simulation, membername[1:], prop)
