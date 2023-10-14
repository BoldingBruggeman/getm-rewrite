from typing import Union, Optional, List, Tuple, Sequence
import logging
import datetime
import timeit
import functools
import pstats
import enum

import numpy as np
import cftime

import xarray as xr

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
import pygetm.internal_pressure


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
    raise Exception(f"Unable to convert {time!r} to cftime.datetime")


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


class Simulation:
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
    __slots__ = (
        "domain",
        "runtype",
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
        "internal_pressure",
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
        "z0s",
        "z0b",
        "fwf",
        "tausx",
        "tausy",
        "tausxo",
        "tausyo",
        "dpdx",
        "dpdy",
        "dpdxo",
        "dpdyo",
        "_cum_river_height_increase",
        "_start_time",
        "_profile",
        "radiation",
        "_initialized_variables",
        "delay_slow_ip",
        "total_volume_ref",
        "total_area",
        "nuh_ct",
    ) + _time_arrays

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
        internal_pressure: Optional[pygetm.internal_pressure.Base] = None,
        logger: Optional[logging.Logger] = None,
        log_level: Optional[int] = None,
        internal_pressure_method: InternalPressure = InternalPressure.SHCHEPETKIN_MCWILLIAMS,
        delay_slow_ip: bool = False,
    ):
        """Simulation

        Args:
            domain: simulation domain
            runtype: simulation run type (BAROTROPIC_2D, BAROTROPIC_3D, BAROCLINIC)
            delay_slow_ip: let slow internal pressure terms lag one macrotimestep
                behind the 3d internal pressure terms. This can help stabilize
                density-driven flows in deep water
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
        domain.initialize(runtype)
        self.domain = domain
        self.runtype = runtype

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

        unmasked = self.domain.T.mask != 0
        self.total_volume_ref = (self.domain.T.H * self.domain.T.area).global_sum(
            where=unmasked
        )
        self.total_area = self.domain.T.area.global_sum(where=unmasked)

        # Configure momentum provider
        if momentum is None:
            momentum = pygetm.momentum.Momentum()
        self.momentum = momentum
        self.momentum.initialize(
            self.logger.getChild("momentum"), domain, runtype, advection_scheme
        )

        self._cum_river_height_increase = np.zeros((len(self.domain.rivers),))

        #: Provider of air-water fluxes of heat and momentum.
        #: This must inherit from :class:`pygetm.airsea.Fluxes`
        #: and should be provided as argument airsea to :class:`Simulation`.
        self.airsea = airsea or pygetm.airsea.FluxesFromMeteo()
        assert isinstance(self.airsea, pygetm.airsea.Fluxes), (
            "airsea argument should be of type pygetm.airsea.Fluxes,"
            f" but is {type(self.airsea)}"
        )
        self.airsea.initialize(self.domain.T)

        self.ice = pygetm.ice.Ice()
        self.ice.initialize(self.domain.T)

        self.dpdx = domain.U.array(
            name="dpdx",
            units="-",
            long_name="surface pressure gradient in x-direction",
            fill_value=FILL_VALUE,
        )
        self.dpdy = domain.V.array(
            name="dpdy",
            units="-",
            long_name="surface pressure gradient in y-direction",
            fill_value=FILL_VALUE,
        )

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

            self.z0s = domain.T.array(
                name="z0s",
                units="m",
                long_name="hydrodynamic roughness (surface)",
                fill_value=FILL_VALUE,
            )
            self.z0s.fill(0.1)

            # Forcing variables for macro/3D momentum update
            # These lag behind the forcing for the micro/2D momentum update
            self.tausxo = domain.U.array()
            self.tausyo = domain.V.array()
            self.dpdxo = domain.U.array()
            self.dpdyo = domain.V.array()

            self.nuh_ct = None
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

                if self.fabm.has_dependency("vertical_tracer_diffusivity"):
                    self.nuh_ct = domain.T.array(
                        name="nuh_ct",
                        units="m2 s-1",
                        long_name="turbulent diffusivity of heat",
                        z=CENTERS,
                        fill_value=FILL_VALUE,
                        attrs=dict(
                            standard_name="ocean_vertical_heat_diffusivity",
                            _mask_output=True,
                        ),
                    )
                    self.nuh_ct.fabm_standard_name = "vertical_tracer_diffusivity"

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
                molecular_diffusivity=1.4e-7,
                attrs=dict(standard_name="sea_water_conservative_temperature"),
            )
            self.salt = self.tracers.add(
                name="salt",
                units="g kg-1",
                long_name="absolute salinity",
                fabm_standard_name="practical_salinity",
                fill_value=FILL_VALUE,
                molecular_diffusivity=1.1e-9,
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

            if internal_pressure is None:
                if internal_pressure_method == InternalPressure.OFF:
                    internal_pressure = pygetm.internal_pressure.Constant()
                elif internal_pressure_method == InternalPressure.BLUMBERG_MELLOR:
                    internal_pressure = pygetm.internal_pressure.BlumbergMellor()
                else:
                    internal_pressure = pygetm.internal_pressure.ShchepetkinMcwilliams()
            internal_pressure.initialize(self.domain)
            self.logger.info(f"Internal pressure method: {internal_pressure!r}")
            self.internal_pressure = internal_pressure
            self.delay_slow_ip = delay_slow_ip
            if delay_slow_ip:
                self.momentum.SxB.attrs["_part_of_state"] = True
                self.momentum.SyB.attrs["_part_of_state"] = True
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
        """Load the model state from a restart file.
        This must be called before :meth:`start`.

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
        with xr.open_dataset(path, **kwargs) as ds:
            timevar = ds["zt"].getm.time
            if timevar.ndim > 1:
                raise Exception(
                    "Time coordinate must be 0D or 1D"
                    f" ({timevar.ndim} dimensions found)"
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
                        f"Requested restart time {time} not found in {path!r},"
                        f" which spans {time_coord[0]} - {time_coord[-1]}"
                    )
                itime = itimes[0]
            elif time_coord.size > 1:
                self.logger.info(
                    f"Restart file {path!r} contains {time_coord.size} time points."
                    f" Using last: {time_coord[-1]}"
                )

            # Slice restart file at the required time index
            if time_coord.size > 1:
                ds = ds.isel({timevar.dims[0]: itime})

            # Load all fields that are part of the model state
            missing = []
            for name, field in self.output_manager.fields.items():
                if field.attrs.get("_part_of_state", False):
                    if name not in ds:
                        missing.append(name)
                    else:
                        field.set(ds[name], on_grid=pygetm.input.OnGrid.ALL, mask=True)
                        self._initialized_variables.add(name)
            if missing:
                raise Exception(
                    "The following field(s) are part of the model state but not found "
                    f"in {path!r}: {', '.join(missing)}"
                )

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
        self.logger.info(f"Starting simulation at {time}")
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
        self.domain.open_boundaries.start(
            self.momentum.U, self.momentum.V, self.momentum.uk, self.momentum.vk
        )
        # Ensure U and V points at the land-water interface have non-zero water depth
        # and layer thickness, as (zero) transports at these points will be divided by
        # these quantities
        for grid in (self.domain.U, self.domain.V):
            edges = grid._water_contact & grid._land
            grid.D.all_values[edges] = FILL_VALUE
            grid.hn.all_values[..., edges] = FILL_VALUE

        if self.fabm:
            self.fabm.start(self.time)

        # Ensure elevations are valid (not shallower than minimum depth)
        minz = -self.domain.T.H.all_values + self.domain.Dmin
        for zname in ("z", "zo", "zin", "zio"):
            z = getattr(self.domain.T, zname)
            shallow = (z.all_values < minz) & self.domain.T._water
            if shallow.any():
                self.logger.warning(
                    f"Increasing {shallow.sum()} elevations in {zname} to avoid"
                    f" water depths below the minimum depth of {self.domain.Dmin} m."
                )
                np.putmask(z.all_values, shallow, minz)

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
            self.domain.update_depth(_3d=True)

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
            self.istep, self.time, default_time_reference=self.default_time_reference
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
                self.internal_pressure.idpdx,
                self.internal_pressure.idpdy,
                self.turbulence.num,
            )

            if self.runtype == BAROCLINIC:
                # Update turbulent quantities (T grid - interfaces) from time=0 to
                # time=1 (macrotimestep), using surface/buoyancy-related forcing
                # (ustar_s, z0s, NN) at time=0, and bottom/velocity-related forcing
                # (ustar_b, z0b, SS) at time=1/2
                # self.domain.T.z0b.all_values[1:, 1:] = 0.5 * (np.maximum(self.domain.U.z0b.all_values[1:, 1:], self.domain.U.z0b.all_values[1:, :-1]) + np.maximum(self.domain.V.z0b.all_values[:-1, 1:], self.domain.V.z0b.all_values[1:, :-1]))
                self.turbulence.advance(
                    self.macrotimestep,
                    self.ustar_s,
                    self.momentum.ustar_b,
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

                if self.delay_slow_ip:
                    self.internal_pressure.idpdx.all_values.sum(
                        axis=0, out=self.momentum.SxB.all_values
                    )
                    self.internal_pressure.idpdy.all_values.sum(
                        axis=0, out=self.momentum.SyB.all_values
                    )

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
            self.domain.open_boundaries.prepare_depth_explicit()

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
            self.internal_pressure(self.buoy)
            if not self.delay_slow_ip:
                self.internal_pressure.idpdx.all_values.sum(
                    axis=0, out=self.momentum.SxB.all_values
                )
                self.internal_pressure.idpdy.all_values.sum(
                    axis=0, out=self.momentum.SyB.all_values
                )

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
        self.ice(baroclinic_active, self.temp_sf, self.salt_sf, self.airsea)

        # Update depth-integrated freshwater fluxes:
        # precipitation, evaporation, condensation, rivers
        self.fwf.all_values[...] = self.airsea.pe.all_values
        np.add.at(
            self.fwf.all_values,
            (self.domain.rivers.j, self.domain.rivers.i),
            self.domain.rivers.flow * self.domain.rivers.iarea,
        )

        # Update elevation at the open boundaries. This must be done before
        # calculating the surface pressure gradient
        self.domain.T.z.open_boundaries.update()

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
            # Update surface shear velocity (used by GOTM). This requires updated
            # surface stresses and there can only be done after the airesea update.
            _pygetm.surface_shear_velocity(
                self.airsea.taux, self.airsea.tauy, self.ustar_s
            )

            # Update radiation. This must come after the airsea update, which is
            # responsible for calculating swr
            self.radiation(self.airsea.swr, self.fabm.kc if self.fabm else None)

            if self.nuh_ct is not None:
                self.turbulence.nuh.interp(self.nuh_ct)

            # Update source terms of biogeochemistry, using the new tracer
            # concentrations. Do this last because FABM could depend on any of the
            # variables computed before
            if self.fabm:
                self.fabm.update_sources(self.timestep * self.istep, self.time)

        if self.report_totals != 0 and self.istep % self.report_totals == 0:
            self.report_domain_integrals()

    @log_exceptions
    def finish(self):
        """Clean-up after simulation: save profiling result (if any), write output
        where appropriate (restarts), and close output files
        """
        if self._profile:
            name, pr = self._profile
            pr.disable()
            profile_path = f"{name}-{self.domain.tiling.rank:03}.prof"
            self.logger.info(f"Writing profiling report to {profile_path}")
            with open(profile_path, "w") as f:
                ps = pstats.Stats(pr, stream=f).sort_stats(pstats.SortKey.TIME)
                ps.print_stats()
                self.summary_profiling_result(ps)
        nsecs = timeit.default_timer() - self._start_time
        self.logger.info(f"Time spent in main loop: {nsecs:.3f} s")
        self.output_manager.close(self.timestep * self.istep, self.time)

    def summary_profiling_result(self, ps: pstats.Stats):
        if not hasattr(ps, "get_stats_profile"):
            # python < 3.9
            return

        sp = ps.get_stats_profile()
        if "<built-in method Waitall>" not in sp.func_profiles:
            # not a parallel simulation, or advance was never called
            return
        stat = [
            sp.total_tt,
            sp.func_profiles["<built-in method Waitall>"].tottime,
            self.domain.T._water_nohalo.sum(),
        ]
        all_stat = self.domain.tiling.comm.gather(stat)
        if all_stat is not None:
            self.logger.info(
                "Time spent on compute per subdomain (excludes halo exchange):"
            )
            for rank, (tottime, halotime, nwet) in enumerate(all_stat):
                self.logger.info(
                    f"{rank} ({nwet} water points): {tottime - halotime:.3f} s"
                )
            rank = np.argmin([s[1] for s in all_stat])
            self.logger.info(
                f"Most expensive subdomain: {rank}"
                f" (see {self._profile[0]}-{rank:03}.prof)"
            )

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
        zl = np.clip(rivers.zl, 1e-6, z_if[:, 0])
        zu = np.clip(rivers.zu, 0.0, zl - 1e-6)
        zbot = np.minimum(zl[:, np.newaxis], z_if[:, :-1])
        ztop = np.maximum(zu[:, np.newaxis], z_if[:, 1:])
        h_active = np.maximum(zbot - ztop, 0.0)

        # Change in thickness per layer
        h_increase_riv = h_active * (z_increases / h_active.sum(axis=1))[:, np.newaxis]

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
        tracer_adds *= h_increase_riv

        # Precipitation and evaporation (surface layer only)
        # First update halos for the net freshwater flux, as we need to ensure that the
        # layer heights updated as a result remain valid in the halos.
        self.airsea.pe.update_halos()
        unmasked = self.domain.T._water
        z_increase_fwf = np.where(unmasked, self.airsea.pe.all_values, 0.0) * timestep
        h = self.domain.T.hn.all_values[-1, :, :]
        h_new = h + z_increase_fwf
        dilution = h / h_new
        for tracer in self.tracers:
            if not tracer.precipitation_follows_target_cell:
                tracer.all_values[-1, :, :] *= dilution
        h[:, :] = h_new

        # Update thicknesses and tracer values with river inflow. This must be done
        # iteratively because different rivers can target the same cell (i, j)
        all_tracer_values = [tracer.all_values for tracer in self.tracers]
        for iriver, (i, j) in enumerate(zip(rivers.i, rivers.j)):
            h = self.domain.T.hn.all_values[:, j, i]
            h_new = h + h_increase_riv[iriver, :]
            h_new_inv = 1.0 / h_new
            for itracer, all_values in enumerate(all_tracer_values):
                tracer_values = all_values[:, j, i]
                add = tracer_adds[itracer, iriver, :]
                tracer_values[:] = (tracer_values * h + add) * h_new_inv
            h[:] = h_new

        # Update elevation
        np.add.at(z_increase_fwf, (rivers.j, rivers.i), z_increases)
        self.domain.T.zin.all_values += z_increase_fwf

        # Start tracer halo exchange (to prepare for advection)
        for tracer in self.tracers:
            tracer.update_halos_start(self.tracers._advection.halo1)

        self._cum_river_height_increase.fill(0.0)

    @property
    def totals(
        self,
    ) -> Tuple[
        Optional[float],
        Optional[Sequence[Tuple[pygetm.tracer.TracerTotal, float, float]]],
    ]:
        """Global totals of volume and tracers.

        Returns:
            A tuple with total volume and a list with (tracer_total, total, mean)
            tuples on the root subdomains. On non-root subdomains it returns None, None
        """
        unmasked = self.domain.T.mask != 0
        total_volume = (self.domain.T.D * self.domain.T.area).global_sum(where=unmasked)
        if any(tt.per_mass for tt in self.tracer_totals):
            vol = self.domain.T.hn * self.domain.T.area
            vol.all_values *= self.rho.all_values
            total_mass = vol.global_sum(where=unmasked)
        tracer_totals = [] if total_volume is not None else None
        if self.fabm:
            self.fabm.update_totals()
        for tt in self.tracer_totals:
            grid = tt.array.grid
            total = tt.array * grid.area
            if tt.scale_factor != 1.0:
                total.all_values *= tt.scale_factor
            if tt.offset != 0.0:
                total.all_values += tt.offset * grid.area.all_values
            if total.ndim == 3:
                if tt.per_mass:
                    total.all_values *= self.rho.all_values
                total.all_values *= grid.hn.all_values
            total = total.global_sum(where=grid.mask != 0)
            if total is not None:
                ref = total_volume if not tt.per_mass else total_mass
                mean = (total / ref - tt.offset) / tt.scale_factor
                tracer_totals.append((tt, total, mean))
        return total_volume, tracer_totals

    def report_domain_integrals(self):
        """Write totals of selected variables over the global domain
        (those in :attr:`tracer_totals`) to the log.
        """
        total_volume, tracer_totals = self.totals
        if total_volume is not None:
            self.logger.info("Integrals over global domain:")
            mean_z = (total_volume - self.total_volume_ref) / self.total_area
            self.logger.info(
                f"  volume: {total_volume:.15e} m3 (mean elevation: {mean_z} m)"
            )
            for tt, total, mean in tracer_totals:
                ar = tt.array
                long_name = tt.long_name if tt.long_name is not None else ar.long_name
                units = tt.units if tt.units is not None else f"{ar.units} m3"
                self.logger.info(
                    f"  {long_name}: {total:.15e} {units}"
                    f" (mean {ar.long_name}: {mean} {ar.units})"
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
        self.domain.T.zo.all_values[:, :] = self.domain.T.z.all_values
        _pygetm.advance_surface_elevation(timestep, self.domain.T.z, U, V, fwf)
        self.domain.T.z.update_halos()

    def update_surface_pressure_gradient(self, z: core.Array, sp: core.Array):
        _pygetm.surface_pressure_gradient(z, sp, self.dpdx, self.dpdy)

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
            unmasked = True
            if not field.on_boundary:
                unmasked = np.isin(field.grid.mask.values, (1, 2))
            unmasked = np.broadcast_to(unmasked, field.shape)
            if not finite.all(where=unmasked):
                nbad += 1
                unmasked_count = unmasked.sum()
                bad_count = unmasked_count - finite.sum(where=unmasked)
                self.logger.error(
                    f"Field {field.name} has {bad_count} non-finite values"
                    f" (out of {unmasked_count} unmasked values)."
                )
        if nbad:
            raise Exception(f"Non-finite values found in {nbad} fields")

    @property
    def Ekin(self, rho0: float = RHO0):
        U = self.momentum.U.interp(self.domain.T)
        V = self.momentum.V.interp(self.domain.T)
        vel2_D2 = U**2 + V**2
        return 0.5 * rho0 * self.domain.T.area * vel2_D2 / self.domain.T.D
