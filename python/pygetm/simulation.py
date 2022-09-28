import operator
from typing import Union, Optional, List
import logging
import datetime
import timeit
import functools

import numpy as np
import cftime
import enum

import xarray

from .constants import (
    BAROTROPIC_2D,
    BAROCLINIC,
    INTERFACES,
    FILL_VALUE,
    RHO0,
    CENTERS,
    GRAVITY,
)
from . import _pygetm
from . import core
from . import parallel
from . import output
from . import operators
import pygetm.domain
import pygetm.airsea
import pygetm.ice
import pygetm.density
import pygetm.fabm
import pygetm.input
import pygetm.mixing
import pygetm.momentum
import pygetm.radiation
import pygetm.tracer


class InternalPressure(enum.IntEnum):
    OFF = 0  #: internal pressure is disabled
    BLUMBERG_MELLOR = 1  #: Blumberg and Mellor
    #    BLUMBERG_MELLOR_LIN=2
    #    Z_INTERPOL=3
    #    SONG_WRIGHT=4
    #    CHU_FAN=5
    SHCHEPETKIN_MCWILLIAMS = 6  #: Shchepetkin and McWilliams (2003)
    #    STELLING_VANKESTER=7


def to_cftime(time: Union[datetime.datetime, cftime.datetime]) -> cftime.datetime:
    if isinstance(time, cftime.datetime):
        return time
    elif isinstance(time, datetime.datetime):
        return cftime.datetime(
            time.year,
            time.month,
            time.day,
            time.hour,
            time.minute,
            time.second,
            time.microsecond,
        )
    raise Exception("Unable to convert %r to cftime.datetime" % time)


def log_exceptions(method):
    @functools.wraps(method)
    def wrapper(self, *args, **kwargs):
        try:
            return method(self, *args, **kwargs)
        except Exception as e:
            logger = getattr(self, "logger", None)
            domain = getattr(self, "domain", None)
            if logger is None or domain is None or domain.tiling.n == 1:
                raise
            logger.exception(str(e), stack_info=True, stacklevel=3)
            domain.tiling.comm.Abort(1)

    return wrapper


class Simulation(_pygetm.Simulation):
    _pressure_arrays = "dpdx", "dpdy", "idpdx", "idpdy"
    _sealevel_arrays = ()
    _time_arrays = (
        "timestep",
        "macrotimestep",
        "split_factor",
        "timedelta",
        "time",
        "istep",
        "report",
        "report_totals",
        "default_time_reference",
    )
    _all_fortran_arrays = tuple(
        ["_%s" % name for name in _pressure_arrays + _sealevel_arrays]
    )
    __slots__ = (
        _all_fortran_arrays
        + (
            "output_manager",
            "input_manager",
            "fabm",
            "_yearday",
            "tracers",
            "tracer_totals",
            "logger",
            "momentum",
            "airsea",
            "ice",
            "turbulence",
            "density",
            "buoy",
            "temp",
            "salt",
            "pres",
            "rad",
            "par",
            "par0",
            "rho",
            "sst",
            "temp_sf",
            "salt_sf",
            "ssu_U",
            "ssv_V",
            "ssu",
            "ssv",
            "NN",
            "ustar_s",
            "ustar_b",
            "taub",
            "z0s",
            "z0b",
            "fwf",
            "tausx",
            "tausy",
            "tausxo",
            "tausyo",
            "dpdxo",
            "dpdyo",
            "_cum_river_height_increase",
            "_start_time",
            "_profile",
            "radiation",
            "_initialized_variables",
        )
        + _time_arrays
    )

    _array_args = {
        "dpdx": dict(units="-", long_name="surface pressure gradient in x-direction"),
        "dpdy": dict(units="-", long_name="surface pressure gradient in y-direction"),
        "idpdx": dict(
            units="m2 s-2", long_name="internal pressure gradient in x-direction"
        ),
        "idpdy": dict(
            units="m2 s-2", long_name="internal pressure gradient in y-direction"
        ),
    }

    domain: pygetm.domain.Domain
    runtype: int

    @log_exceptions
    def __init__(
        self,
        domain: pygetm.domain.Domain,
        runtype: int,
        advection_scheme: operators.AdvectionScheme = operators.AdvectionScheme.DEFAULT,
        fabm: Union[pygetm.fabm.FABM, bool, str, None] = None,
        gotm: Union[str, None] = None,
        momentum: Optional[pygetm.momentum.Momentum] = None,
        turbulence: Optional[pygetm.mixing.Turbulence] = None,
        airsea: Optional[pygetm.airsea.Fluxes] = None,
        density: Optional[pygetm.density.Density] = None,
        radiation: Optional[pygetm.radiation.Radiation] = None,
        logger: Optional[logging.Logger] = None,
        log_level: Optional[int] = None,
        internal_pressure_method: InternalPressure = InternalPressure.OFF,
    ):
        """Simulation
        
        Args:
            domain: simulation domain
            """
        self.logger = domain.root_logger
        if log_level is not None:
            self.logger.setLevel(log_level)
        self.output_manager = output.OutputManager(
            domain.fields,
            rank=domain.tiling.rank,
            logger=self.logger.getChild("output_manager"),
        )
        self.input_manager = domain.input_manager

        self.input_manager.set_logger(self.logger.getChild("input_manager"))

        assert not domain._initialized
        super().__init__(
            domain, runtype, internal_pressure_method=internal_pressure_method
        )

        # Flag selected domain/grid fields (elevations, thicknesses, bottom roughness)
        # as part of model state (saved in/loaded from restarts)
        domain.T.z.attrs["_part_of_state"] = True
        domain.T.zo.attrs["_part_of_state"] = True
        if self.runtype > BAROTROPIC_2D:
            domain.T.zio.attrs["_part_of_state"] = True
            domain.T.zin.attrs["_part_of_state"] = True

            # ho cannot be computed from zio,
            # because rivers modify ho-from-zio before it is stored
            domain.T.ho.attrs["_part_of_state"] = True
        if self.runtype == BAROTROPIC_2D:
            domain.U.z0b.attrs["_part_of_state"] = True
            domain.V.z0b.attrs["_part_of_state"] = True

        # Configure momentum provider
        if momentum is None:
            momentum = pygetm.momentum.Momentum()
        self.momentum = momentum
        self.momentum.initialize(
            self.logger.getChild("momentum"), domain, runtype, advection_scheme
        )

        # Make open boundary conditions for elevation and transport/velocity
        # available as part of domain.open_boundaries
        self.wrap(self.domain.open_boundaries.z, b"zbdy", source=3)
        self.momentum.wrap(self.domain.open_boundaries.u, b"bdyu")
        self.momentum.wrap(self.domain.open_boundaries.v, b"bdyv")

        for name in Simulation._pressure_arrays:
            kwargs = dict(fill_value=FILL_VALUE)
            kwargs.update(Simulation._array_args.get(name, {}))
            setattr(
                self,
                "_%s" % name,
                self.wrap(
                    core.Array(name=name, **kwargs), name.encode("ascii"), source=2,
                ),
            )
        for name in Simulation._sealevel_arrays:
            setattr(
                self,
                "_%s" % name,
                self.wrap(
                    core.Array(name=name, **Simulation._array_args.get(name, {})),
                    name.encode("ascii"),
                    source=3,
                ),
            )

        self._cum_river_height_increase = np.zeros((len(self.domain.rivers),))

        #: Provider of air-water fluxes of heat and momentum.
        #: This must inherit from :class:`pygetm.airsea.Fluxes`
        #: and should be provided as argument airsea to :class:`Simulation`.
        self.airsea = airsea or pygetm.airsea.FluxesFromMeteo()
        assert isinstance(self.airsea, pygetm.airsea.Fluxes), (
            "airsea argument should be of type pygetm.airsea.Fluxes, but is %s"
            % type(self.airsea)
        )
        self.airsea.initialize(self.domain.T)

        self.ice = pygetm.ice.Ice()
        self.ice.initialize(self.domain.T)

        # Surface stresses interpolated to U and V grids
        self.tausx = domain.U.array(
            name="tausxu", fill_value=FILL_VALUE, attrs={"_mask_output": True}
        )
        self.tausy = domain.V.array(
            name="tausyv", fill_value=FILL_VALUE, attrs={"_mask_output": True}
        )

        self.fwf = domain.T.array(
            name="fwf",
            units="m s-1",
            long_name="freshwater flux",
            fill_value=FILL_VALUE,
            attrs={"_mask_output": self.airsea.pe.attrs.get("_mask_output", False)},
        )
        self.fwf.fill(0.0)

        #: Collection of tracers that are to be transported.
        #: Optionally they can have sources, open boundary conditions
        #: and riverine concentrations set.
        self.tracers: pygetm.tracer.TracerCollection = pygetm.tracer.TracerCollection(
            self.domain.T, advection_scheme=advection_scheme
        )

        #: List of variables for which the domain-integrated total needs to be reported.
        #: These can be depth-integrated (2D) or depth-explicit (3D).
        self.tracer_totals: List[pygetm.tracer.TracerTotal] = []

        self.fabm = None

        if runtype > BAROTROPIC_2D:
            #: Provider of turbulent viscosity and diffusivity. This must inherit from
            #: :class:`pygetm.mixing.Turbulence` and should be provided as argument
            #: turbulence to :class:`Simulation`.
            self.turbulence = turbulence or pygetm.mixing.GOTM(gotm)
            self.turbulence.initialize(self.domain.T)
            self.NN = domain.T.array(
                z=INTERFACES,
                name="NN",
                units="s-2",
                long_name="buoyancy frequency squared",
                fill_value=FILL_VALUE,
                attrs=dict(
                    standard_name="square_of_brunt_vaisala_frequency_in_sea_water"
                ),
            )
            self.NN.fill(0.0)
            self.ustar_s = domain.T.array(
                fill=0.0,
                name="ustar_s",
                units="m s-1",
                long_name="shear velocity (surface)",
                fill_value=FILL_VALUE,
                attrs=dict(_mask_output=True),
            )
            self.ustar_b = domain.T.array(
                fill=0.0,
                name="ustar_b",
                units="m s-1",
                long_name="shear velocity (bottom)",
                fill_value=FILL_VALUE,
                attrs=dict(_mask_output=True),
            )
            self.z0s = domain.T.array(
                name="z0s",
                units="m",
                long_name="hydrodynamic roughness (surface)",
                fill_value=FILL_VALUE,
            )
            self.z0s.fill(0.1)
            self.taub = domain.T.array(
                fill=0.0,
                name="taub",
                units="Pa",
                long_name="bottom shear stress",
                fill_value=FILL_VALUE,
                fabm_standard_name="bottom_stress",
                attrs=dict(_mask_output=True),
            )

            # Forcing variables for macro/3D momentum update
            # These lag behind the forcing for the micro/2D momentum update
            self.tausxo = domain.U.array()
            self.tausyo = domain.V.array()
            self.dpdxo = domain.U.array()
            self.dpdyo = domain.V.array()

            if fabm:
                if not isinstance(fabm, pygetm.fabm.FABM):
                    fabm = pygetm.fabm.FABM(
                        fabm if isinstance(fabm, str) else "fabm.yaml"
                    )
                self.fabm = fabm
                self.fabm.initialize(
                    self.domain.T,
                    self.tracers,
                    self.tracer_totals,
                    self.logger.getChild("FABM"),
                )

            self.pres = domain.depth
            self.pres.fabm_standard_name = "pressure"

        self.sst = domain.T.array(
            name="sst",
            units="degrees_Celsius",
            long_name="sea surface temperature",
            fill_value=FILL_VALUE,
            attrs=dict(standard_name="sea_surface_temperature", _mask_output=True),
        )

        self.ssu = domain.T.array(fill=0.0)
        self.ssv = domain.T.array(fill=0.0)

        if runtype == BAROCLINIC:
            self.logger.info(
                "Internal pressure method: %s" % internal_pressure_method.name
            )
            self.density = density or pygetm.density.Density()

            self.radiation = radiation or pygetm.radiation.TwoBand()
            self.radiation.initialize(self.domain.T)

            self.temp = self.tracers.add(
                name="temp",
                units="degrees_Celsius",
                long_name="conservative temperature",
                fabm_standard_name="temperature",
                fill_value=FILL_VALUE,
                source=self.radiation.swr_abs,
                surface_flux=self.airsea.shf,
                source_scale=1.0 / (RHO0 * self.density.CP),
                rivers_follow_target_cell=True,
                precipitation_follows_target_cell=True,
                attrs=dict(standard_name="sea_water_conservative_temperature"),
            )
            self.salt = self.tracers.add(
                name="salt",
                units="g kg-1",
                long_name="absolute salinity",
                fabm_standard_name="practical_salinity",
                fill_value=FILL_VALUE,
                attrs=dict(standard_name="sea_water_absolute_salinity"),
            )
            self.pres.saved = True
            self.temp.fill(5.0)
            self.salt.fill(35.0)
            self.rho = domain.T.array(
                z=CENTERS,
                name="rho",
                units="kg m-3",
                long_name="density",
                fabm_standard_name="density",
                fill_value=FILL_VALUE,
                attrs=dict(standard_name="sea_water_density", _mask_output=True),
            )
            self.buoy = domain.T.array(
                z=CENTERS,
                name="buoy",
                units="m s-2",
                long_name="buoyancy",
                attrs=dict(_mask_output=True),
            )
            self.tracer_totals += [
                pygetm.tracer.TracerTotal(
                    self.salt, units="g", per_mass=True, long_name="salt"
                ),
                pygetm.tracer.TracerTotal(
                    self.temp,
                    units="J",
                    per_mass=True,
                    scale_factor=self.density.CP,
                    offset=self.density.CP * 273.15,
                    long_name="heat",
                ),
            ]
            self.temp_sf = self.temp.isel(z=-1)
            self.salt_sf = self.salt.isel(z=-1)
            self.ssu_U = self.momentum.uk.isel(z=-1)
            self.ssv_V = self.momentum.vk.isel(z=-1)
            if internal_pressure_method == InternalPressure.OFF:
                self.idpdx.fill(0.0)
                self.idpdy.fill(0.0)
        else:
            self.temp_sf = None
            self.salt_sf = None

        # Derive old and new elevations, water depths and thicknesses from current
        # surface elevation on T grid. This must be done after self.pres.saved is set
        self.domain.update_depth(_3d=runtype > BAROTROPIC_2D)
        self.domain.update_depth(_3d=runtype > BAROTROPIC_2D)

        self.default_time_reference: Optional[cftime.datetime] = None
        self._initialized_variables = set()

    def __getitem__(self, key: str) -> core.Array:
        return self.output_manager.fields[key]

    @log_exceptions
    def load_restart(
        self, path: str, time: Optional[cftime.datetime] = None, **kwargs
    ) -> cftime.datetime:
        """Load the model state from a restart file. This must be called before :meth:`start`.

        Args:
            path: NetCDF file to load restart state from
            time: time coordinate to take restart information from. This is only
                relevant of the restart file contains the model state at multiple times.
                If not provided, the last time from the file will be used.
            **kwargs: additional keyword arguments passed to :func:`xarray.open_dataset`

        Returns:
            The time from which the restart information was taken.
        """
        kwargs.setdefault("decode_times", True)
        kwargs["use_cftime"] = True
        with xarray.open_dataset(path, **kwargs) as ds:
            timevar = ds["zt"].getm.time
            if timevar.ndim > 1:
                raise Exception(
                    "Time coordinate must be 0D or 1D (%i dimensions found)"
                    % timevar.ndim
                )

            # Use time reference of restart as default time reference for new output
            self.default_time_reference = cftime.num2date(
                0.0,
                units=timevar.encoding["units"],
                calendar=timevar.encoding["calendar"],
            )

            # Determine time index to load restart information from
            time_coord = timevar.values.reshape((-1,))
            itime = -1
            if time is not None:
                # Find time index that matches requested time
                time = to_cftime(time)
                itimes = np.where(time_coord == time)[0]
                if itimes.size == 0:
                    raise Exception(
                        "Requested restart time %s not found in %s, which spans %s - %s"
                        % (
                            time.isoformat(),
                            path,
                            time_coord[0].isoformat(),
                            time_coord[-1].isoformat(),
                        )
                    )
                itime = itimes[0]
            elif time_coord.size > 1:
                self.logger.info(
                    "Restart file %s contains %i time points. Using last: %s"
                    % (path, time_coord.size, time_coord[-1].isoformat())
                )

            # Slice restart file at the required time index
            if time_coord.size > 1:
                ds = ds.isel({timevar.dims[0]: itime})

            # Load all fields that are part of the model state
            for name, field in self.output_manager.fields.items():
                if field.attrs.get("_part_of_state", False):
                    if name not in ds:
                        raise Exception(
                            "Field %s is part of state but not found in %s" % path
                        )
                    field.set(ds[name], on_grid=pygetm.input.OnGrid.ALL, mask=True)
                    self._initialized_variables.add(name)

        if self.runtype > BAROTROPIC_2D:
            # Restore elevation from before open boundary condition was applied
            self.domain.T.z.all_values[...] = self.domain.T.zin.all_values

        return time_coord[itime]

    @log_exceptions
    def start(
        self,
        time: Union[cftime.datetime, datetime.datetime],
        timestep: float,
        split_factor: int = 1,
        report: Union[int, datetime.timedelta] = 10,
        report_totals: Union[int, datetime.timedelta] = datetime.timedelta(days=1),
        profile: Optional[str] = None,
    ):
        """Start a simulation by configuring the time, zeroing velocities, updating
        diagnostics to match the start time, and optionally saving output.

        This should be called after the output configuration is complete
        (because we need to know when variables need to be saved),
        and after the FABM model has been provided with all dependencies.

        Args:
            time (:class:`cftime.datetime`): start time
            timestep: micro time step (s) used for 2D barotropic processes
            split_factor: number of microtimesteps per macrotimestep
            report: time interval or number of microtimesteps between reporting of the
                current time, used as indicator of simulation progress
            report_totals: time interval or number of microtimesteps between reporting
                of integrals over the global domain
            profile: base name for the file to write profiling results to. The proces
                rank and extension ``.prof`` will be appended, so that the final name
                becomes ``<profile>-<rank>.prof``. If the argument is not provided,
                profiling is disabled.
        """
        time = to_cftime(time)
        self.logger.info("Starting simulation at %s" % time)
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

        self.momentum.start()
        self.tracers.start()

        # Ensure U and V points at the land-water interface have non-zero water depth
        # and layer thickness, as (zero) transports at these points will be divided by
        # these quantities
        for grid in (self.domain.U, self.domain.V):
            edges = grid._water_contact & grid._land
            grid.D.all_values[edges] = FILL_VALUE
            grid.hn.all_values[..., edges] = FILL_VALUE

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
            # Since we do not have the preceding (2 time steps before start) zi/h, we
            # explicitly set them (here: T.zio/T.ho) to NaN to make it easier to detect
            # algorithms depending on them.
            # As a result of that, all new metrics on the U, V, X grids will be NaN too!
            self.domain.T.z.all_values[
                ...
            ] = (
                self.domain.T.zio.all_values
            )  # to become T.zin when update_depth is called
            self.domain.T.zio.fill(np.nan)
            self.domain.T.ho.fill(np.nan)
            self.domain.update_depth(_3d=True, timestep=self.macrotimestep)

            # Second 3D depth/thickness update based on zin.
            # Override T.ho with user-provided value if available, since this may
            # incorporate river inflow impacts that our previously calculated ho cannot
            # account for.
            # New metrics for U, V, X grids will be calculated from valid old and new
            # metrics on T grid; therefore they will be valid too. However, old metrics
            # (ho/zio) for U, V, X grids will still be NaN and should not be used.
            self.domain.T.z.all_values[...] = zin_backup

            if "hot" in self._initialized_variables:
                self.domain.T.hn.all_values[...] = ho_T_backup

            # this moves our zin backup into zin, and at the same time moves the
            # current zin (originally zio) to zio
            self.domain.update_depth(_3d=True, timestep=self.macrotimestep)
            self.momentum.update_diagnostics(self.macrotimestep, self.turbulence.num)

        # Update all forcing, which includes the final 2D depth update based on
        # (original) z
        self.domain.T.z.all_values[...] = z_backup
        self.update_forcing(macro_active=self.runtype > BAROTROPIC_2D)

        # Start output manager
        self.output_manager.start(
            self.istep, self.time, default_time_reference=self.default_time_reference,
        )

        # Verify all fields have finite values. Do this after self.output_manager.start
        # so the user can diagnose issues by reviewing the output
        self.check_finite(_3d=self.runtype > BAROTROPIC_2D)

        # Record true start time for performance analysis
        self._start_time = timeit.default_timer()

        # Start profiling if requested
        self._profile = None
        if profile:
            import cProfile

            pr = cProfile.Profile()
            self._profile = (profile, pr)
            pr.enable()

    @log_exceptions
    def advance(self, check_finite: bool = False):
        """Advance the model state by one microtimestep.
        If this completes the current macrotimestep, the part of the state associated
        with that timestep will be advanced too.

        Args:
            check_finite: after the state update, verify that all fields only contain
                finite values
        """
        # Update the time
        macro_updated = self.istep % self.split_factor == 0
        self.time += self.timedelta
        self.istep += 1
        macro_active = self.istep % self.split_factor == 0
        if self.report != 0 and self.istep % self.report == 0:
            self.logger.info(self.time)

        self.output_manager.prepare_save(
            self.timestep * self.istep, self.istep, self.time, macro=macro_updated
        )

        # Update transports U and V from time=-1/2 to +1/2, using surface stresses and
        # pressure gradients defined at time=0
        # Inputs and outputs are on U and V grids. Stresses and pressure gradients have
        # already been updated by the call to update_forcing at the end of the previous
        # time step.
        self.momentum.advance_depth_integrated(
            self.timestep, self.tausx, self.tausy, self.dpdx, self.dpdy
        )

        # Update surface elevation on T grid from time=0 to time=1 using transports
        # U and V at time=1/2 and freshwater fluxes at time=0. This also updates halos
        # so that depths and thicknesses can be computed everywhere without further
        # halo exchange
        self.advance_surface_elevation(
            self.timestep, self.momentum.U, self.momentum.V, self.fwf
        )

        # Track cumulative increase in elevation due to river inflow over the current
        # macrotimestep
        self._cum_river_height_increase += (
            self.domain.rivers.flow * self.domain.rivers.iarea * self.timestep
        )

        if self.runtype > BAROTROPIC_2D and macro_active:
            # Use previous source terms for biogeochemistry (valid for the start of the
            # current macrotimestep) to update tracers. This should be done before the
            # tracer concentrations change due to transport or rivers, as the source
            # terms are only valid for the current tracer concentrations.
            if self.fabm:
                self.fabm.advance(self.macrotimestep)

            # Update layer thicknesses and tracer concentrations to account for
            # precipitation, evaporation and river inflow between start and end of the
            # current macrotimestep.
            self.add_freshwater_inputs(self.macrotimestep)

            # Update 3D elevations and layer thicknesses. New elevation zin on T grid
            # will match elevation at end of the microtimestep, thicknesses on T grid
            # will match. Elevation and thicknesses on U/V grids will be
            # 1/2 MACROtimestep behind, as they are calculated by averaging zio and zin.
            # Old elevations zio and thicknesses ho will store the previous values of
            # zin and hn, and thus are one macrotimestep behind new elevations zin and
            # thicknesses hn.
            self.domain.update_depth(_3d=True, timestep=self.macrotimestep)

            # Update momentum from time=-1/2 to 1/2 of the macrotimestep, using forcing
            # defined at time=0. For this purpose, surface stresses (tausxo, tausyo)
            # and surface pressure gradients (dpdxo, dpdyo) at the end of the previous
            # macrotimestep were saved
            # Internal pressure idpdx and idpdy were calculated at the end of the
            # previous macrotimestep and are therefore ready as-is.
            self.momentum.advance(
                self.macrotimestep,
                self.split_factor,
                self.tausxo,
                self.tausyo,
                self.dpdxo,
                self.dpdyo,
                self.idpdx,
                self.idpdy,
                self.turbulence.num,
            )

            if self.runtype == BAROCLINIC:
                # Update total stresses (x and y combined) and calculate the friction
                # velocities (m s-1) This is for turbulence (GOTM), so all on the T grid
                # Note that surface stress is currently defined at the start of the
                # MICROtimestep (only one microtimestep before the end of the current
                # macrotimestep), whereas bottom stress is at 1/2 of the macrotimestep,
                # since it is computed from velocities that have just been updated.
                self.momentum.update_stresses(self.airsea.taux, self.airsea.tauy)
                np.sqrt(self.momentum.ustar2_s.all_values, out=self.ustar_s.all_values)
                np.sqrt(self.momentum.ustar2_b.all_values, out=self.ustar_b.all_values)
                self.taub.all_values[...] = self.momentum.ustar2_b.all_values * RHO0

                # Update turbulent quantities (T grid - interfaces) from time=0 to
                # time=1 (macrotimestep), using surface/buoyancy-related forcing at
                # time=0, and velocity-related forcing at time=1/2
                # self.domain.T.z0b.all_values[1:, 1:] = 0.5 * (np.maximum(self.domain.U.z0b.all_values[1:, 1:], self.domain.U.z0b.all_values[1:, :-1]) + np.maximum(self.domain.V.z0b.all_values[:-1, 1:], self.domain.V.z0b.all_values[1:, :-1]))
                self.turbulence.advance(
                    self.macrotimestep,
                    self.ustar_s,
                    self.ustar_b,
                    self.z0s,
                    self.domain.T.z0b,
                    self.NN,
                    self.momentum.SS,
                )

                # Advect and diffuse tracers. Source terms are optionally handled too,
                # as part of the diffusion update.
                self.tracers.advance(
                    self.macrotimestep,
                    self.momentum.uk,
                    self.momentum.vk,
                    self.momentum.ww,
                    self.turbulence.nuh,
                )

        if self.report_totals != 0 and self.istep % self.report_totals == 0:
            self.report_domain_integrals()

        # Update all inputs and fluxes that will drive the next state update
        self.update_forcing(
            macro_active,
            skip_2d_coriolis=True,
            update_z0b=self.runtype == BAROTROPIC_2D,
        )

        self.output_manager.save(self.timestep * self.istep, self.istep, self.time)

        if check_finite:
            self.check_finite(_3d=macro_active)

        return macro_active

    def update_forcing(
        self,
        macro_active: bool,
        skip_2d_coriolis: bool = False,
        update_z0b: bool = False,
    ):
        """Update all inputs and fluxes that will drive the next state update.

        Args:
            macro_active: update all quantities associated with the macrotimestep
            skip_2d_coriolis: whether to skip the update of depth-integrated Coriolis
                terms, typically because they have already been recalculated as part of
                the transport update
            update_z0b: update bottom roughness z0b as part of the depth-integrated
                bottom friction calculation
        """
        # Update all inputs.
        self.domain.input_manager.update(self.time, macro=macro_active)

        baroclinic_active = self.runtype == BAROCLINIC and macro_active
        if baroclinic_active:
            # Update tracer values at open boundaries. This must be done after
            # input_manager.update, but before diagnostics/forcing variables derived
            # from the tracers are calculated
            if self.domain.open_boundaries.np:
                for tracer in self.tracers:
                    tracer.open_boundaries.update()

            # Update density, buoyancy and internal pressure to keep them in sync with
            # T and S.
            self.density.get_density(self.salt, self.temp, p=self.pres, out=self.rho)

            # Update density halos: valid rho around all U and V needed for internal
            # pressure; not yet valid because T&S were not valid in halos when rho was
            # calculated. Note BM needs only right/top, SMcW needs left/right/top/bottom
            self.rho.update_halos(parallel.Neighbor.LEFT_AND_RIGHT_AND_TOP_AND_BOTTOM)
            self.buoy.all_values[...] = (-GRAVITY / RHO0) * (self.rho.all_values - RHO0)
            self.update_internal_pressure_gradient(self.buoy)
            self.momentum.SxB.all_values[...] = self.idpdx.all_values.sum(axis=0)
            self.momentum.SyB.all_values[...] = self.idpdy.all_values.sum(axis=0)

            # From conservative temperature to in-situ sea surface temperature,
            # needed to compute heat/momentum fluxes at the surface
            self.density.get_potential_temperature(
                self.salt_sf, self.temp_sf, out=self.sst
            )

            # Calculate squared buoyancy frequency NN
            # (T grid, interfaces between layers)
            self.density.get_buoyancy_frequency(
                self.salt, self.temp, p=self.pres, out=self.NN
            )

            # Interpolate surface velocities to T grid.
            # These are used by airsea as offset for wind speeds
            self.ssu_U.interp(self.ssu)
            self.ssv_V.interp(self.ssv)

        # Update surface elevation z on U, V, X grids and water depth D on all grids
        # This is based on old and new elevation (T grid) for the microtimestep.
        # Thus, for grids lagging 1/2 a timestep behind (U, V, X grids), the elevations
        # and water depths will be representative for 1/2 a MICROtimestep ago.
        # Note that T grid elevations at the open boundary have not yet been updated,
        # so the derived elevations and water depths calculated here will not take
        # those into account. This is intentional: it ensures that water depths on the
        # U and V grids are in sync with the already-updated transports,
        # so that velocities can be calculated correctly.
        # The call to update_surface_elevation_boundaries is made later.
        self.domain.update_depth()

        # Calculate advection and diffusion tendencies of transports, bottom friction
        # and, if needed, Coriolis terms
        self.momentum.update_depth_integrated_diagnostics(
            self.timestep, skip_coriolis=skip_2d_coriolis, update_z0b=update_z0b
        )

        # Update air-sea fluxes of heat and momentum (all on T grid)
        # Note: sst is the in-situ surface temperature, whereas temp_sf is the
        # conservative surface temperature (salt_sf is absolute salinity)
        self.airsea(
            self.time,
            self.sst,
            self.ssu,
            self.ssv,
            calculate_heat_flux=baroclinic_active,
        )
        self.ice(
            baroclinic_active, self.temp_sf, self.salt_sf, self.airsea,
        )

        # Update depth-integrated freshwater fluxes:
        # precipitation, evaporation, condensation, rivers
        self.fwf.all_values[...] = self.airsea.pe.all_values
        self.fwf.all_values[self.domain.rivers.j, self.domain.rivers.i] += (
            self.domain.rivers.flow * self.domain.rivers.iarea
        )

        # Update elevation at the open boundaries. This must be done before
        # update_surface_pressure_gradient
        self.update_surface_elevation_boundaries(self.timestep, self.momentum)

        # Calculate the surface pressure gradient in the U and V points.
        # Note: this requires elevation and surface air pressure (both on T grid) to be
        # valid in the halos, which is guaranteed for elevation (halo exchange happens
        # just after update), and for air pressure if it is managed by the input
        # manager (e.g. read from file)
        self.airsea.sp.update_halos(parallel.Neighbor.TOP_AND_RIGHT)
        self.update_surface_pressure_gradient(self.domain.T.z, self.airsea.sp)

        # Interpolate surface stresses from T to U and V grids
        self.airsea.taux.update_halos(parallel.Neighbor.RIGHT)
        self.airsea.taux.interp(self.tausx)
        self.airsea.tauy.update_halos(parallel.Neighbor.TOP)
        self.airsea.tauy.interp(self.tausy)

        if self.runtype > BAROTROPIC_2D and macro_active:
            # Save surface forcing variables for the next macro momentum update
            self.tausxo.all_values[...] = self.tausx.all_values
            self.tausyo.all_values[...] = self.tausy.all_values
            self.dpdxo.all_values[...] = self.dpdx.all_values
            self.dpdyo.all_values[...] = self.dpdy.all_values

        if baroclinic_active:
            # Update radiation. This must come after the airsea update, which is
            # responsible for calculating swr
            self.radiation(self.airsea.swr)

            # Update source terms of biogeochemistry, using the new tracer
            # concentrations. Do this last because FABM could depend on any of the
            # variables computed before
            if self.fabm:
                self.fabm.update_sources(self.time)

    @log_exceptions
    def finish(self):
        """Clean-up after simulation: save profiling result (if any), write output
        where appropriate (restarts), and close output files
        """
        if self._profile:
            import pstats

            name, pr = self._profile
            pr.disable()
            profile_path = "%s-%03i.prof" % (name, self.domain.tiling.rank)
            self.logger.info("Writing profiling report to %s" % (profile_path,))
            with open(profile_path, "w") as f:
                ps = pstats.Stats(pr, stream=f).sort_stats(pstats.SortKey.TIME)
                ps.print_stats()
        self.logger.info(
            "Time spent in main loop: %.3f s"
            % (timeit.default_timer() - self._start_time,)
        )
        self.output_manager.close(self.timestep * self.istep, self.time)

    def add_freshwater_inputs(self, timestep: float):
        """Update layer thicknesses and tracer concentrations to account for
        precipitation, evaporation and river inflow.
        """
        # Local names for river-related variables
        rivers = self.domain.rivers
        z_increases = self._cum_river_height_increase

        # Depth of layer interfaces for each river cell
        h = self.domain.T.hn.all_values[:, rivers.j, rivers.i].T
        z_if = np.zeros((h.shape[0], h.shape[1] + 1))
        z_if[:, 1:] = -h.cumsum(axis=1)
        z_if -= z_if[:, -1:]

        # Determine the part of every layer over which inflow occurs
        zl = np.minimum(z_if[:, 0], rivers.zl)
        zu = np.minimum(zl - 1e-6, rivers.zu)
        zbot = np.minimum(zl[:, np.newaxis], z_if[:, :-1])
        ztop = np.maximum(zu[:, np.newaxis], z_if[:, 1:])
        h_active = np.maximum(zbot - ztop, 0.0)

        # Change in thickness per layer
        h_increases = h_active * (z_increases / h_active.sum(axis=1))[:, np.newaxis]

        # Calculate the depth-integrated change in tracer, per layer.
        tracer_adds = np.empty((len(self.tracers),) + h.shape)
        for tracer, river_values in zip(self.tracers, tracer_adds):
            follow = tracer.river_follow | (z_increases < 0.0)
            ext = tracer.river_values[:, np.newaxis]
            if follow.any():
                int = tracer.all_values[:, rivers.j, rivers.i].T
                river_values[...] = np.where(follow[:, np.newaxis], int, ext)
            else:
                river_values[...] = ext
        tracer_adds *= h_increases

        # Update thicknesses and tracer values. This must be done iteratively
        # because different rivers can target the same cell (i, j)
        all_tracer_values = [tracer.all_values for tracer in self.tracers]
        for iriver, (i, j) in enumerate(zip(rivers.i, rivers.j)):
            h_old = self.domain.T.hn.all_values[:, j, i]
            h_new = h_old + h_increases[iriver, :]
            h_new_inv = 1.0 / h_new
            for itracer, all_values in enumerate(all_tracer_values):
                tracer_values = all_values[:, j, i]
                add = tracer_adds[itracer, iriver, :]
                tracer_values[:] = (tracer_values * h_old + add) * h_new_inv
            h_old[:] = h_new
            self.domain.T.zin.all_values[j, i] += z_increases[iriver]

        # Precipitation and evaporation (surface layer only)
        # First update halos for the net freshwater flux, as we need to ensure that the
        # layer heights updated as a result remain valid in the halos.
        self.airsea.pe.update_halos()
        unmasked = self.domain.T._water
        h_increase_pe = np.where(unmasked, self.airsea.pe.all_values, 0.0) * timestep
        h = self.domain.T.hn.all_values[-1, :, :]
        h_new = h + h_increase_pe
        for tracer in self.tracers:
            if not tracer.precipitation_follows_target_cell:
                tracer.all_values[-1, :, :] *= h / h_new
        h[:, :] = h_new
        self.domain.T.zin.all_values += h_increase_pe

        # Start tracer halo exchange (to prepare for advection)
        for tracer in self.tracers:
            tracer.update_halos_start(self.tracers._advection.halo1)

        self._cum_river_height_increase.fill(0.0)

    def report_domain_integrals(self):
        """Write totals of selected variables over the global domain
        (those in :attr:`tracer_totals`) to the log.
        """
        unmasked = self.domain.T.mask != 0
        total_volume = (self.domain.T.D * self.domain.T.area).global_sum(where=unmasked)
        if any(tt.per_mass for tt in self.tracer_totals):
            vol = self.domain.T.hn * self.domain.T.area
            vol.all_values *= self.rho.all_values
            total_mass = vol.global_sum(where=unmasked)
        if total_volume is not None:
            self.logger.info("Integrals over global domain:")
            self.logger.info("  volume: %.15e m3" % total_volume)
        if self.fabm:
            self.fabm.update_totals()
        for tt in self.tracer_totals:
            ar = tt.array
            total = ar * ar.grid.area
            if tt.scale_factor != 1.0:
                total.all_values *= tt.scale_factor
            if tt.offset != 0.0:
                total.all_values += tt.offset * ar.grid.area.all_values
            if total.ndim == 3:
                if tt.per_mass:
                    total.all_values *= self.rho.all_values
                total.all_values *= ar.grid.hn.all_values
            total = total.global_sum(where=ar.grid.mask != 0)
            long_name = tt.long_name if tt.long_name is not None else ar.long_name
            units = tt.units if tt.units is not None else f"{ar.units} m3"
            if total is not None:
                ref = total_volume if not tt.per_mass else total_mass
                mean = (total / ref - tt.offset) / tt.scale_factor
                self.logger.info(
                    "  %s: %.15e %s (mean %s: %s %s)"
                    % (long_name, total, units, ar.long_name, mean, ar.units,)
                )

    def advance_surface_elevation(
        self, timestep: float, U: core.Array, V: core.Array, fwf: core.Array
    ):
        """Advance surface elevation (T grid only)

        Args:
            timestep: time step (s)
            U: depth-integrated velocity in x-direction (m2 s-1)
            V: depth-integrated velocity in y-direction (m2 s-1)
            fwf: freshwater flux (m s-1)

        This also updates the surface elevation halos.
        This method does `not` update elevation on the U, V, X grids, nor water depths,
        layer thicknesses or vertical coordinates.
        This is done by :meth:`~pygetm.domain.Domain.update_depth` instead.
        """
        super().advance_surface_elevation(timestep, U, V, fwf)
        self.domain.T.z.update_halos()

    def check_finite(self, _3d: bool = True):
        """Verify that all fields available for output contain finite values.
        Fields with non-finite values are reported in the log as error messages.
        Finally, if any non-finite values were found, an exception is raised.

        Args:
            _3d: also check fields updated on the 3d (macro) timestep, not just
                depth-integrated fields
        """
        nbad = 0
        for field in self.domain.fields.values():
            if field.ndim == 0 or (field.z and not _3d):
                continue
            finite = np.isfinite(field.values)
            unmasked = True if field.on_boundary else field.grid.mask.values > 0
            unmasked = np.broadcast_to(unmasked, field.shape)
            if not finite.all(where=unmasked):
                nbad += 1
                unmasked_count = unmasked.sum()
                bad_count = unmasked_count - finite.sum(where=unmasked)
                self.logger.error(
                    "Field %s has %i non-finite values (out of %i unmasked values)."
                    % (field.name, bad_count, unmasked_count)
                )
        if nbad:
            raise Exception("Non-finite values found in %i fields" % nbad)

    @property
    def Ekin(self, rho0: float = RHO0):
        U = self.momentum.U.interp(self.domain.T)
        V = self.momentum.V.interp(self.domain.T)
        vel2_D2 = U ** 2 + V ** 2
        return 0.5 * rho0 * self.domain.T.area * vel2_D2 / self.domain.T.D


# Expose all Fortran arrays that are a member of Simulation as read-only properties
# The originals are members with and underscore as prefix, ad therefore not visible to
# the user. This ensures the user will not accidentally disconnect the Python variable
# from the underlying Fortran libraries/data
for membername in Simulation._all_fortran_arrays:
    attrs = Simulation._array_args.get(membername[1:], {})
    long_name = attrs.get("long_name")
    units = attrs.get("units")
    doc = ""
    if long_name is not None:
        doc = long_name
        if units:
            doc += " (%s)" % units
    prop = property(operator.attrgetter(membername), doc=doc)
    setattr(Simulation, membername[1:], prop)
