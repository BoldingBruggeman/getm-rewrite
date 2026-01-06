from typing import Union, Optional, Sequence, Mapping
import logging
import datetime
import timeit
import functools
import pstats

import numpy as np
import cftime

import xarray as xr

from .constants import (
    INTERFACES,
    FILL_VALUE,
    RHO0,
    CENTERS,
    GRAVITY,
    TimeVarying,
    RunType,
    CellType,
)
from . import _pygetm
from . import core
from . import parallel
from . import operators
import pygetm.domain
import pygetm.airsea
import pygetm.ice
import pygetm.density
import pygetm.fabm
import pygetm.input
import pygetm.vertical_mixing
import pygetm.momentum
import pygetm.radiation
import pygetm.tracer
import pygetm.internal_pressure
import pygetm.vertical_coordinates
import pygetm.open_boundaries


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
            tiling = getattr(self, "tiling", None)
            if logger is None or tiling is None or tiling.n == 1:
                raise
            logger.exception(str(e), stack_info=True, stacklevel=3)
            tiling.comm.Abort(1)

    return wrapper


class BaseSimulation:
    __slots__ = (
        "logger",
        "tiling",
        "_fields",
        "output_manager",
        "input_manager",
        "default_time_reference",
        "_initialized_variables",
        "timestep",
        "macrotimestep",
        "split_factor",
        "timedelta",
        "time",
        "istep",
        "report",
        "report_totals",
        "_start_time",
        "_profile",
        "_cached_check_finite_info",
    )

    def __init__(
        self,
        domain: pygetm.domain.Domain,
        *,
        log_level: Optional[int] = None,
        tiling: Optional[parallel.Tiling] = None,
    ):
        self.logger = domain.root_logger
        if log_level is not None:
            self.logger.setLevel(log_level)

        self.tiling = tiling or domain.create_tiling()

        self._fields: Mapping[str, core.Array] = {}

        self.output_manager = pygetm.output.OutputManager(
            self._fields,
            rank=self.tiling.rank,
            logger=self.logger.getChild("output_manager"),
        )

        self.input_manager = pygetm.input.InputManager(
            logger=self.logger.getChild("input_manager")
        )

        self.default_time_reference: Optional[cftime.datetime] = None
        self._initialized_variables = set()
        self._cached_check_finite_info = None

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
        try:
            kwargs.setdefault(
                "decode_times", xr.coders.CFDatetimeCoder(use_cftime=True)
            )
        except AttributeError:
            # xarray < 2025.01.1
            kwargs.setdefault("decode_times", True)
            kwargs["use_cftime"] = True
        kwargs["decode_timedelta"] = False
        with xr.open_dataset(path, **kwargs) as ds:
            timevar = ds["time"]
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
                (itimes,) = (time_coord == time).nonzero()
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
            if timevar.ndim == 1:
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

        self._after_restart()

        return time_coord[itime]

    def _after_restart(self):
        pass

    @log_exceptions
    def start(
        self,
        time: Union[cftime.datetime, datetime.datetime],
        timestep: Union[float, datetime.timedelta],
        split_factor: int = 1,
        report: Union[int, datetime.timedelta] = 10,
        report_totals: Union[int, datetime.timedelta] = datetime.timedelta(days=1),
        profile: Optional[str] = None,
        dump_on_error: bool = False,
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
            profile: base name for the file to write profiling results to. The process
                rank and extension ``.prof`` will be appended, so that the final name
                becomes ``<profile>-<rank>.prof``. If the argument is not provided,
                profiling is disabled.
            dump_on_error: if True, when non-finite values are detected in the model
                state or forcing, the current state will be saved to a NetCDF file
                for debugging purposes.
        """
        self.time = to_cftime(time)
        self.logger.info(f"Starting simulation at {self.time}")

        if isinstance(timestep, datetime.timedelta):
            timestep = timestep.total_seconds()
        self.timestep = timestep
        self.timedelta = datetime.timedelta(seconds=timestep)

        self.split_factor = split_factor
        self.macrotimestep = self.timestep * self.split_factor

        self.istep = 0
        if isinstance(report, datetime.timedelta):
            report = int(round(report.total_seconds() / self.timestep))
        self.report = report
        if isinstance(report_totals, datetime.timedelta):
            report_totals = int(round(report_totals.total_seconds() / self.timestep))
        self.report_totals = report_totals

        self._start()

        # Update all inputs and diagnostics.
        # This includes forcing for the future state update (e.g., tracer source terms)
        self.input_manager.update(self.time, macro=True)
        self._update_forcing_and_diagnostics(macro_active=True)

        # Start output manager
        self.output_manager.start(
            self.istep, self.time, default_time_reference=self.default_time_reference
        )

        # Verify all fields have finite values. Do this after self.output_manager.start
        # so the user can diagnose issues by reviewing the output
        if not self.check_finite(dump=dump_on_error):
            raise Exception("Initial state or forcing is invalid")

        # Record true start time for performance analysis
        self._start_time = timeit.default_timer()

        # Start profiling if requested
        self._profile = None
        if profile:
            import cProfile

            pr = cProfile.Profile()
            self._profile = (profile, pr)
            pr.enable()

    def _start(self):
        pass

    def _update_forcing_and_diagnostics(self, macro_active: bool):
        pass

    @log_exceptions
    def advance(self, check_finite: bool = False):
        """Advance the model state by one microtimestep.
        If this completes the current macrotimestep, the part of the state associated
        with that timestep will be advanced too.

        Args:
            check_finite: after the state update, verify that all fields only contain
                finite values
        """
        # Update time-averaged outputs with the values at the start of the time step
        macro_updated = self.istep % self.split_factor == 0
        self.output_manager.prepare_save(macro=macro_updated)

        # Update the time
        self.time += self.timedelta
        self.istep += 1
        macro_active = self.istep % self.split_factor == 0
        if self.report != 0 and self.istep % self.report == 0:
            self.logger.info(self.time)

        self._advance_state(macro_active)

        # Update all inputs and (state-dependent) diagnostics.
        # This includes forcing for the future state update (e.g., tracer source terms)
        self.input_manager.update(self.time, macro=macro_active)
        self._update_forcing_and_diagnostics(macro_active)

        # Perform output. This is done before the call to check_finite,
        # so that the user can diagnose non-finite values from the output files.
        self.output_manager.save(self.timestep * self.istep, self.istep, self.time)

        if check_finite:
            if not self.check_finite(macro_active):
                raise Exception("Non-finite values found")

    def _advance_state(self, macro_active: bool):
        pass

    @log_exceptions
    def finish(self):
        """Clean-up after simulation: save profiling result (if any), write output
        where appropriate (restarts), and close output files
        """
        if self._profile:
            name, pr = self._profile
            pr.disable()
            prof_fn = f"{name}-{self.tiling.rank:03}.prof"
            self.logger.info(
                f"Writing profiling results to {prof_fn}"
                " (view with snakeviz or similar)"
            )
            pr.dump_stats(prof_fn)
            stats_fn = f"{name}-{self.tiling.rank:03}.prof.stats"
            self.logger.info(f"Writing human-readable profiling summary to {stats_fn}")
            with open(stats_fn, "w") as f:
                ps = pstats.Stats(pr, stream=f).sort_stats(pstats.SortKey.TIME)
                ps.print_stats()
                self._summarize_profiling_result(ps)
        nsecs = timeit.default_timer() - self._start_time
        self.logger.info(f"Time spent in main loop: {nsecs:.3f} s")
        self.output_manager.close(self.timestep * self.istep, self.time)

    def check_finite(self, macro_active: bool = True, dump: bool = True) -> bool:
        """Verify that all fields available for output contain finite values.
        Fields with non-finite values are reported in the log as error messages.
        Finally, if any non-finite values were found, an exception is raised.

        Args:
            macro_active: also check fields updated on the 3d (macro) timestep
        """

        def _collect_info():
            microchecks = []
            macrochecks = []
            for field in self._fields.values():
                if field.ndim == 0:
                    continue
                unmasked = ~field.mask
                default_time_varying = (
                    TimeVarying.MACRO if field.z else TimeVarying.MICRO
                )
                time_varying = field.attrs.get("_time_varying", default_time_varying)
                if time_varying == TimeVarying.MICRO:
                    microchecks.append((field, unmasked))
                macrochecks.append((field, unmasked))
            return (microchecks, macrochecks)

        if self._cached_check_finite_info is None:
            self._cached_check_finite_info = _collect_info()
        checklist = self._cached_check_finite_info[1 if macro_active else 0]
        bad_fields: list[str] = []
        for field, unmasked in checklist:
            finite = np.isfinite(field.values)
            if not finite.all(where=unmasked):
                bad_fields.append(field.name)
                unmasked_count = unmasked.sum()
                bad_count = unmasked_count - finite.sum(where=unmasked)
                self.logger.error(
                    f"Field {field.name} has {bad_count} non-finite values"
                    f" (out of {unmasked_count} unmasked values)."
                )
        nsub = np.empty((self.tiling.n,), dtype=int)
        nsub[self.tiling.rank] = len(bad_fields)
        self.tiling.comm.Allgather(parallel.MPI.IN_PLACE, nsub)
        fail = nsub.any()
        if fail:
            sublist = (f"{i} ({n} fields)" for i, n in enumerate(nsub) if n > 0)
            self.logger.error(
                f"Non-finite values found in {(nsub > 0).sum()} subdomains"
                f" at istep={self.istep}, time={self.time}."
                f" Affected subdomains: {', '.join(sublist)}"
            )
            if dump:
                all_bad_fields = self.tiling.comm.allreduce(bad_fields)
                assert all_bad_fields
                dump_fields = set(all_bad_fields)
                dump_fields = [f for f in self._fields.values() if f.ndim]
                dump_file = f"{parallel.LOGFILE_PREFIX}dump.nc"
                self.logger.info(f"Dumping {len(dump_fields)} fields to {dump_file}")
                self.output_manager.dump(dump_file, *dump_fields)
            self.tiling.comm.Barrier()
        return not fail

    def _summarize_profiling_result(self, ps: pstats.Stats):
        pass


class Simulation(BaseSimulation):
    __slots__ = (
        "runtype",
        "fabm",
        "_yearday",
        "tracers",
        "tracer_totals",
        "momentum",
        "airsea",
        "ice",
        "vertical_mixing",
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
        "_int_river_flow",
        "radiation",
        "delay_slow_ip",
        "total_volume_ref",
        "total_area",
        "nuh_ct",
        "U",
        "V",
        "X",
        "T",
        "Dmin",
        "Dcrit",
        "rivers",
        "open_boundaries",
        "vertical_coordinates",
        "h_T_half",
        "depth",
    )

    @log_exceptions
    def __init__(
        self,
        domain: pygetm.domain.Domain,
        *,
        runtype: RunType = RunType.BAROCLINIC,
        advection_scheme: operators.AdvectionScheme = operators.AdvectionScheme.DEFAULT,
        fabm: Union[pygetm.fabm.FABM, bool, str, None] = None,
        gotm: Union[str, None] = None,
        momentum: Optional[pygetm.momentum.Momentum] = None,
        vertical_mixing: Optional[pygetm.vertical_mixing.VerticalMixing] = None,
        airsea: Optional[pygetm.airsea.Fluxes] = None,
        density: Optional[pygetm.density.Density] = None,
        radiation: Optional[pygetm.radiation.Radiation] = None,
        internal_pressure: Optional[pygetm.internal_pressure.Base] = None,
        vertical_coordinates: Optional[pygetm.vertical_coordinates.Base] = None,
        Dmin: float = 0.02,
        Dcrit: float = 0.1,
        logger: Optional[logging.Logger] = None,
        log_level: Optional[int] = None,
        delay_slow_ip: bool = False,
        tiling: Optional[parallel.Tiling] = None,
    ):
        """Simulation

        Args:
            domain: simulation domain
            runtype: simulation run type
            Dmin: minimum depth (m) for wet points. At this depth, all hydrodynamic
                terms except the pressure gradient and bottom friction are switched off.
            Dcrit: depth (m) at which tapering off of hydrodynamic processes
                (all except pressure gradient and bottom friction) begins.
            delay_slow_ip: let slow internal pressure terms lag one macrotimestep
                behind the 3d internal pressure terms. This can help stabilize
                density-driven flows in deep water
            tiling: subdomain decomposition. If not provided, an optimal subdomain
                division is determined automatically based on the land-sea mask and
                the number of active CPU cores
        """
        super().__init__(domain, log_level=log_level, tiling=tiling)

        if Dmin <= 0.0:
            self.logger.error(f"Dmin ({Dmin} m) must exceed zero")
            raise Exception("Dmin<=0")
        if Dcrit < 2.5 * Dmin:
            self.logger.error(
                f"Dcrit ({Dcrit} m) must equal or exceed 2.5 * Dmin ({Dmin} m)"
                f" = {2.5 * Dmin} m"
            )
            raise Exception("Dcrit < 2.5*Dmin")

        HALO = 2

        self.runtype = runtype

        domain.logger.info(f"Maximum slope factor (rx0): {domain.max_rx0:.3f}")

        maxdt, i, j, depth = domain.cfl_check(return_location=True)
        domain.logger.info(
            f"Maximum timestep for 2D barotropic processes: {maxdt:.3f} s "
            f"(i={i}, j={j}, bathymetric depth={depth:.3f} m)"
        )

        if runtype == RunType.BAROTROPIC_2D:
            vertical_coordinates = None
        elif vertical_coordinates is None:
            self.logger.warn(
                "Argument vertical_coordinates not provided; using a single layer"
            )
            vertical_coordinates = pygetm.vertical_coordinates.Sigma(1)
        self.vertical_coordinates = vertical_coordinates

        self.T = domain.create_grids(
            nz=vertical_coordinates.nz if vertical_coordinates else None,
            halox=HALO,
            haloy=HALO,
            fields=self._fields,
            tiling=self.tiling,
            input_manager=self.input_manager,
            velocity_grids=2,
            t_postfix="t",
        )

        self.open_boundaries = self.T.open_boundaries
        self.rivers = self.T.rivers

        self.U = self.T.ugrid
        self.V = self.T.vgrid
        self.X = self.T.xgrid

        for grid in (self.T, self.U, self.V, self.X):
            self.logger.info(
                f"Number of unmasked {grid.postfix.upper()} points,"
                f" excluding halos: {(grid.mask.values > 0).sum()}"
            )

        self.Dmin = Dmin
        self.Dcrit = Dcrit

        self.T.z = self.T.array(
            name="zt",
            units="m",
            long_name="surface elevation",
            fill_value=FILL_VALUE,
            attrs=dict(
                standard_name="sea_surface_height_above_geopotential_datum",
                _part_of_state=True,
                _minimum=self.Dmin - self.T.H.all_values,
                _valid_at=(CellType.BOUNDARY,),
            ),
        )
        self.T.zo = self.T.array(
            name="zot",
            units="m",
            long_name="surface elevation at previous microtimestep",
            fill_value=FILL_VALUE,
            attrs=dict(_part_of_state=True, _valid_at=self.T.z.attrs["_valid_at"]),
        )
        # Initialize elevation at all water points to 0
        self.T.z.fill(0.0)
        self.T.zo.fill(0.0)

        if runtype > RunType.BAROTROPIC_2D:
            self.T.zin = self.T.array(
                name="zint",
                units="m",
                long_name="surface elevation at macrotimestep",
                fill_value=FILL_VALUE,
                attrs=dict(
                    _part_of_state=True,
                    _time_varying=TimeVarying.MACRO,
                    _valid_at=self.T.z.attrs["_valid_at"],
                ),
            )
            self.T.zio = self.T.array(
                name="ziot",
                units="m",
                long_name="surface elevation at previous macrotimestep",
                fill_value=FILL_VALUE,
                attrs=dict(
                    _part_of_state=True,
                    _time_varying=TimeVarying.MACRO,
                    _valid_at=self.T.zin.attrs["_valid_at"],
                ),
            )

            for grid in (self.T, self.U, self.V):
                grid.ho = grid.array(
                    name="ho" + grid.postfix,
                    z=CENTERS,
                    units="m",
                    long_name="cell thickness at previous time step",
                    fill_value=FILL_VALUE,
                    attrs=dict(_valid_at=grid.hn.attrs["_valid_at"]),
                )

            # On T grid, ho cannot be computed from zio,
            # because rivers modify ho-from-zio before it is stored
            self.T.ho.attrs["_part_of_state"] = True
        else:
            # In 2D barotropic runs, bottom roughness is updated iteratively
            # (new value is a function of the previous value), which requires
            # it being included in restarts as part of the model state
            self.U.z0b.attrs["_part_of_state"] = True
            self.V.z0b.attrs["_part_of_state"] = True

        # Enable open boundaries for surface elevation
        self.T.z.open_boundaries = pygetm.open_boundaries.ArrayOpenBoundaries(self.T.z)
        self.open_boundaries.z = self.T.z.open_boundaries.values

        # Water depths clipped to Dmin (already the default for U,V,X grids)
        self.T.Dclip = self.T.array()

        # Water depth and thicknesses on UU/VV grids will be taken from T grid,
        # which near land has valid values where UU/VV are masked
        self.U.ugrid.D.attrs["_mask_output"] = True
        self.V.vgrid.D.attrs["_mask_output"] = True

        if runtype > RunType.BAROTROPIC_2D:
            self.U.ugrid.hn.attrs["_mask_output"] = True
            self.V.vgrid.hn.attrs["_mask_output"] = True

            # Thicknesses on T grid that lag 1/2 time step behind tracer
            # (i.e., they are in sync with U, V, X grids)
            self.h_T_half = self.T.array(fill=np.nan, z=CENTERS)

            self.depth = self.T.array(
                z=CENTERS,
                name="pres",
                units="dbar",
                long_name="pressure",
                fabm_standard_name="depth",
                fill_value=FILL_VALUE,
                attrs=dict(_valid_at=(CellType.BOUNDARY,)),  # to compute density
            )

        unmasked = self.T.mask != CellType.UNRESOLVED
        self.total_volume_ref = (self.T.H * self.T.area).global_sum(where=unmasked)
        self.total_area = self.T.area.global_sum(where=unmasked)

        # Configure momentum provider
        if momentum is None:
            momentum = pygetm.momentum.Momentum()
        self.momentum = momentum
        self.momentum.initialize(
            self.logger.getChild("momentum"), self.T, runtype, advection_scheme
        )

        self._int_river_flow = np.zeros((len(self.rivers),))

        #: Provider of air-water fluxes of heat and momentum.
        #: This must inherit from :class:`pygetm.airsea.Fluxes`
        #: and should be provided as argument airsea to :class:`Simulation`.
        self.airsea = airsea or pygetm.airsea.FluxesFromMeteo()
        assert isinstance(self.airsea, pygetm.airsea.Fluxes), (
            "airsea argument should be of type pygetm.airsea.Fluxes,"
            f" but is {type(self.airsea)}"
        )
        self.airsea.initialize(self.T, self.logger.getChild("airsea"))

        self.ice = pygetm.ice.Ice()
        self.ice.initialize(self.T, self.logger.getChild("ice"))

        self.dpdx = self.U.array(
            name="dpdx",
            units="-",
            long_name="surface pressure gradient in x-direction",
            fill_value=FILL_VALUE,
        )
        self.dpdy = self.V.array(
            name="dpdy",
            units="-",
            long_name="surface pressure gradient in y-direction",
            fill_value=FILL_VALUE,
        )

        # Surface stresses interpolated to U and V grids
        self.tausx = self.U.array(
            name="tausxu",
            units="Pa",
            long_name="surface stress in x-direction",
            fill_value=FILL_VALUE,
            attrs=dict(_mask_output=True),
        )
        self.tausy = self.V.array(
            name="tausyv",
            units="Pa",
            long_name="surface stress in y-direction",
            fill_value=FILL_VALUE,
            attrs=dict(_mask_output=True),
        )

        self.fwf = self.T.array(
            name="fwf",
            units="m s-1",
            long_name="freshwater flux",
            fill_value=FILL_VALUE,
            attrs=dict(_mask_output=self.airsea.pe.attrs.get("_mask_output", False)),
        )
        self.fwf.fill(0.0)

        #: Collection of tracers that are to be transported.
        #: Optionally they can have sources, open boundary conditions
        #: and riverine concentrations set.
        self.tracers: pygetm.tracer.TracerCollection = pygetm.tracer.TracerCollection(
            self.T, self.logger.getChild("tracers"), advection_scheme=advection_scheme
        )

        #: List of variables for which the domain-integrated total needs to be reported.
        #: These can be depth-integrated (2D) or depth-explicit (3D).
        self.tracer_totals: list[pygetm.tracer.TracerTotal] = []

        self.fabm = None

        # Surface temperature (in-situ) and velocities are needed for all
        # run types as they are used for air-sea exchange
        self.sst = self.T.array(
            name="sst",
            units="degrees_Celsius",
            long_name="sea surface temperature",
            fill_value=FILL_VALUE,
            attrs=dict(
                standard_name="sea_surface_temperature",
                _mask_output=True,
                _valid_at=(CellType.BOUNDARY,),
            ),
        )
        self.ssu = self.T.array(fill=0.0)
        self.ssv = self.T.array(fill=0.0)

        if runtype > RunType.BAROTROPIC_2D:
            #: Provider of turbulent viscosity and diffusivity. This must inherit from
            #: :class:`pygetm.vertical_mixing.VerticalMixing` and should be provided as
            # argument `vertical_mixing` to :class:`Simulation`.
            self.vertical_mixing = vertical_mixing or pygetm.vertical_mixing.GOTM(gotm)
            self.vertical_mixing.initialize(
                self.T, self.logger.getChild("vertical_mixing")
            )
            self.NN = self.T.array(
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
            self.ustar_s = self.T.array(
                fill=0.0,
                name="ustar_s",
                units="m s-1",
                long_name="shear velocity (surface)",
                fill_value=FILL_VALUE,
                attrs=dict(_mask_output=True),
            )

            self.z0s = self.T.array(
                name="z0s",
                units="m",
                long_name="hydrodynamic roughness (surface)",
                fill_value=FILL_VALUE,
                attrs=dict(_time_varying=TimeVarying.MACRO),
            )
            self.z0s.fill(0.1)

            # Forcing variables for macro/3D momentum update
            # These lag behind the forcing for the micro/2D momentum update
            self.tausxo = self.U.array()
            self.tausyo = self.V.array()
            self.dpdxo = self.U.array()
            self.dpdyo = self.V.array()

            self.nuh_ct = None
            if fabm:
                if not isinstance(fabm, pygetm.fabm.FABM):
                    fabm = pygetm.fabm.FABM(
                        fabm if isinstance(fabm, str) else "fabm.yaml"
                    )
                self.fabm = fabm
                self.fabm.initialize(
                    self.T,
                    self.tracers,
                    self.tracer_totals,
                    self.logger.getChild("FABM"),
                )

                if self.fabm.has_dependency("vertical_tracer_diffusivity"):
                    self.nuh_ct = self.T.array(
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

            self.pres = self.depth
            self.pres.fabm_standard_name = "pressure"

            self.ssu_U = self.momentum.uk.isel(z=-1)
            self.ssv_V = self.momentum.vk.isel(z=-1)

        if radiation is None and runtype == RunType.BAROCLINIC:
            radiation = pygetm.radiation.TwoBand()
        self.radiation = radiation
        if self.radiation is not None:
            self.radiation.initialize(self.T, self.logger.getChild("radiation"))

        if runtype == RunType.BAROCLINIC:
            self.density = density or pygetm.density.Density()

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

            # Set initial temperature and salinity to default value throughout
            # the domain
            self.temp.fill(5.0)
            self.salt.fill(35.0)

            # Ensure [approximate] pressure is updated in update_depth as it is
            # needed for equation of state
            self.pres.saved = True

            self.rho = self.T.array(
                z=CENTERS,
                name="rho",
                units="kg m-3",
                long_name="density",
                fabm_standard_name="density",
                fill_value=FILL_VALUE,
                attrs=dict(
                    standard_name="sea_water_density",
                    _mask_output=True,
                    _valid_at=(CellType.BOUNDARY,),
                ),
            )
            self.buoy = self.T.array(
                z=CENTERS,
                name="buoy",
                units="m s-2",
                long_name="buoyancy",
                attrs=dict(_mask_output=True, _valid_at=(CellType.BOUNDARY,)),
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

            # Select surface fields for conservative temperature and absolute
            # salinity, to be used to calculate in-situ surface temperature
            self.temp_sf = self.temp.isel(z=-1)
            self.salt_sf = self.salt.isel(z=-1)

            # Surface temperature will be calculated from 3D temperature and salinity
            # and therefore varies on baroclinic timestep only
            self.sst.attrs.update(_time_varying=TimeVarying.MACRO)
        else:
            self.temp_sf = None
            self.salt_sf = None

        if runtype == RunType.BAROTROPIC_2D:
            internal_pressure = None
        elif runtype == RunType.BAROTROPIC_3D:
            internal_pressure = pygetm.internal_pressure.Prescribed(idpdx=0, idpdy=0)
        elif internal_pressure is None:
            internal_pressure = pygetm.internal_pressure.ShchepetkinMcwilliams()
        self.internal_pressure = internal_pressure
        if self.internal_pressure is not None:
            self.logger.info(f"Internal pressure method: {self.internal_pressure!r}")
            self.internal_pressure.initialize(self.U, self.V)
            self.delay_slow_ip = delay_slow_ip
            if delay_slow_ip:
                self.momentum.SxB.attrs["_part_of_state"] = True
                self.momentum.SyB.attrs["_part_of_state"] = True

        # Initialize vertical coordinates as very last step, so that the underlying
        # logic can access all model fields (e.g., NN and SS for adaptive coordinates)
        if self.vertical_coordinates is not None:
            self.vertical_coordinates.initialize(
                self.T,
                self.U,
                self.V,
                self.X,
                logger=self.logger.getChild("vertical_coordinates"),
            )

        # Derive old and new elevations, water depths and thicknesses from current
        # surface elevation on T grid. This must be done after self.pres.saved is set
        self.update_depth(_3d=runtype > RunType.BAROTROPIC_2D)
        self.update_depth(_3d=runtype > RunType.BAROTROPIC_2D)

    def _after_restart(self):
        if self.runtype > RunType.BAROTROPIC_2D:
            # Restore elevation from before open boundary condition was applied
            self.T.z.all_values = self.T.zin.all_values

    def _start(self):
        self.momentum.start()
        self.tracers.start()
        self.open_boundaries.start(
            self.momentum.U,
            self.momentum.V,
            self.momentum.uk if self.runtype > RunType.BAROTROPIC_2D else None,
            self.momentum.vk if self.runtype > RunType.BAROTROPIC_2D else None,
        )
        # Ensure U and V points at the land-water interface have finite and non-zero
        # water depth and layer thickness, as (zero) transports at these points will
        # be divided by these quantities to compute velocities, which in turn are
        # needed to compute Coriolis and bottom friction terms.
        # NB update_depth will not set these, as they are unresolved (mask=0)
        for grid in (self.U, self.V):
            edges = grid._edge_x | grid._edge_y
            grid.D.all_values[edges] = FILL_VALUE
            if grid.hn is not None:
                grid.hn.all_values[..., edges] = FILL_VALUE

        if self.fabm:
            self.fabm.start(self.macrotimestep, self.time)

        def clip_z(z: core.Array, valid_min: np.ndarray):
            shallow = (z.all_values < valid_min) & self.T._water
            if shallow.any():
                self.logger.warning(
                    f"Increasing {shallow.sum()} elevations in {z.name} to ensure"
                    f" initial water depths equal or exceed the minimum depth of"
                    f" {self.Dmin} m"
                )
                np.putmask(z.all_values, shallow, valid_min)

        # Ensure elevations are valid (not shallower than minimum depth)
        minz = -self.T.H.all_values + self.Dmin
        clip_z(self.T.z, minz)
        clip_z(self.T.zo, minz)
        if self.runtype > RunType.BAROTROPIC_2D:
            clip_z(self.T.zin, minz)
            clip_z(self.T.zio, minz)

        # First (out of two) 2D depth update based on old elevations zo
        z_backup = self.T.z.all_values.copy()
        self.T.z.all_values = self.T.zo.all_values
        self.update_depth(_3d=False)

        if self.runtype > RunType.BAROTROPIC_2D:
            zin_backup = self.T.zin.all_values.copy()
            ho_T_backup = self.T.ho.all_values.copy()

            # First (out of two) 3D depth/thickness update based on zio.
            # This serves to generate T.ho when T.zio is set, but T.ho is not available.
            # Since we do not have the preceding (2 time steps before start) zi/h, we
            # explicitly set them (here: T.zio/T.ho) to NaN to make it easier to detect
            # algorithms depending on them.
            # As a result of that, all new metrics on the U, V, X grids will be NaN too!
            self.T.z.all_values = (
                self.T.zio.all_values
            )  # to become T.zin when update_depth is called
            self.T.zio.fill(np.nan)
            self.T.ho.fill(np.nan)
            self.update_depth(_3d=True)

            # Second 3D depth/thickness update based on zin.
            # Override T.ho with user-provided value if available, since this may
            # incorporate river inflow impacts that our previously calculated ho cannot
            # account for.
            # New metrics for U, V, X grids will be calculated from valid old and new
            # metrics on T grid; therefore they will be valid too. However, old metrics
            # (ho/zio) for U, V, X grids will still be NaN and should not be used.
            self.T.z.all_values = zin_backup

            if "hot" in self._initialized_variables:
                self.T.hn.all_values = ho_T_backup

            # this moves our zin backup into zin, and at the same time moves the
            # current zin (originally zio) to zio
            self.update_depth(_3d=True, timestep=self.macrotimestep)
            self.momentum.update_diagnostics(
                self.macrotimestep, self.vertical_mixing.num
            )

        # Update all forcing, which includes the final 2D depth update based on
        # (original) z
        self.T.z.all_values = z_backup

    def _advance_state(self, macro_active: bool):
        # Update transports U and V from time=-1/2 to +1/2, using surface stresses and
        # pressure gradients defined at time=0
        # Inputs and outputs are on U and V grids. Stresses and pressure gradients have
        # already been updated by the call to _update_forcing_and_diagnostics at the end
        #  of the previous time step.
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

        # Track cumulative river inflow (m3) over the current macrotimestep
        self._int_river_flow += self.rivers.flow * self.timestep

        if self.runtype > RunType.BAROTROPIC_2D and macro_active:
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

            # Update water depth D and layer thicknesses hn on all grids.
            # On the T grid, these will be consistent with surface elevation
            # at the end of the microtimestep, that is, with the result of
            # the call to advance_surface_elevation called above.
            # Water depth and thicknesses on U/V/X grids will be
            # 1/2 MACROtimestep behind.
            # On the T grid, the previous value of surface elevation and
            # thicknesses will be stored in variables zio and ho, respectively.
            # These will thus be a full macrotimestep behind, but do account
            # for freshwater input over the past macrotimestep, as that was
            # added to surface elevation and thicknesses by the call to
            # add_freshwater_inputs above.
            self.update_depth(_3d=True, timestep=self.macrotimestep)

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
                self.vertical_mixing.num,
            )

            # Update turbulent quantities (T grid - interfaces) from time=0 to
            # time=1 (macrotimestep), using surface/buoyancy-related forcing
            # (ustar_s, z0s, NN) at time=0, and bottom/velocity-related forcing
            # (ustar_b, z0b, SS) at time=1/2
            # self.T.z0b.all_values[1:, 1:] = 0.5 * (np.maximum(self.U.z0b.all_values[1:, 1:], self.U.z0b.all_values[1:, :-1]) + np.maximum(self.V.z0b.all_values[:-1, 1:], self.V.z0b.all_values[1:, :-1]))
            self.vertical_mixing.advance(
                self.macrotimestep,
                self.ustar_s,
                self.momentum.ustar_b,
                self.z0s,
                self.T.z0b,
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
                self.vertical_mixing.nuh,
            )

            # If we have to delay slow (2D depth-integrated) terms for internal pressure
            # by one macrotimestep, calculate them now, at the end of state update of the
            # macrotimestep, and just before the new 3D internal pressure is calculated.
            if self.runtype == RunType.BAROCLINIC and self.delay_slow_ip:
                self.internal_pressure.idpdx.all_values.sum(
                    axis=0, out=self.momentum.SxB.all_values
                )
                self.internal_pressure.idpdy.all_values.sum(
                    axis=0, out=self.momentum.SyB.all_values
                )

    def _update_forcing_and_diagnostics(self, macro_active: bool):
        """Update all inputs and fluxes that will drive the next state update.

        Args:
            macro_active: update all quantities associated with the macrotimestep
        """
        starting = self.istep == 0

        if starting:
            self.rivers.flag_prescribed_tracers()

        # Hydrodynamic bottom roughness is updated iteratively in barotropic simulations
        # Do this only at the end of a time step, that is, not at the very start of the
        # simulation.
        update_z0b = self.runtype == RunType.BAROTROPIC_2D and not starting

        update_3d = self.runtype > RunType.BAROTROPIC_2D and macro_active
        update_baroclinic = self.runtype == RunType.BAROCLINIC and macro_active

        # Prepare prescribed open boundaries, colocated model fields, derived metrics
        # For instance, rotate prescribed velocities, extract inward model velocities,
        # and calculate derived quantities specific to some types of boundary condition.
        self.open_boundaries.prepare(update_3d)

        if update_3d:
            # Update tracer values at open boundaries. This must be done after
            # input_manager.update, but before diagnostics/forcing variables derived
            # from the tracers are calculated.
            if self.open_boundaries.np:
                for tracer in self.tracers:
                    tracer.open_boundaries.update()

            # Interpolate surface velocities to T grid.
            # These are used by airsea as offset for wind speeds
            self.ssu_U.interp(self.ssu)
            self.ssv_V.interp(self.ssv)

        if update_baroclinic:
            # Update density, buoyancy and internal pressure to keep them in sync with
            # T and S.
            self.density.get_density(self.salt, self.temp, p=self.pres, out=self.rho)

            # Update density halos: valid rho around all U and V needed for internal
            # pressure; not yet valid because T&S were not valid in halos when rho was
            # calculated. Note BM needs only right/top, SMcW needs left/right/top/bottom
            self.rho.update_halos(parallel.Neighbor.LEFT_AND_RIGHT_AND_TOP_AND_BOTTOM)
            self.buoy.all_values = (-GRAVITY / RHO0) * (self.rho.all_values - RHO0)
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

        # Update water depth D on all grids.
        # For grids lagging 1/2 a timestep behind (U, V, X grids), the
        # water depths will be representative for 1/2 a MICROtimestep ago.
        # They are calculated from old and new elevations on the T grid.
        # Note that T grid elevations at the open boundary have not yet been updated,
        # so water depths calculated here will not take those into account.
        # This is intentional: it ensures that water depths on the
        # U and V grids are in sync with the already-updated transports,
        # so that velocities can be calculated correctly.
        # The call to update_surface_elevation_boundaries is made later.
        self.update_depth()

        # Calculate tendencies of transports (depth-integrated velocities) due
        # to advection and diffusion, bottom friction and Coriolis terms.
        # Only do 2D Coriolis update at the start of the simulation.
        # At subsequent times, this term will already have been updated
        # as part of the momentum update in _advance_state
        self.momentum.update_depth_integrated_diagnostics(
            self.timestep, skip_coriolis=not starting, update_z0b=update_z0b
        )

        # Update air-sea fluxes of heat and momentum (all on T grid)
        # Note: sst is the in-situ surface temperature, whereas temp_sf is the
        # conservative surface temperature (salt_sf is absolute salinity)
        self.airsea(
            self.time,
            self.sst,
            self.ssu,
            self.ssv,
            calculate_heat_flux=update_baroclinic,
        )

        # Update ice coverage. This is done after the airsea update to allow
        # the ice module to manipulate (e.g., suppress) surface fluxes of
        # heat and momentum
        self.ice(update_baroclinic, self.temp_sf, self.salt_sf, self.airsea)

        # Update depth-integrated freshwater fluxes:
        # precipitation/evaporation/condensation from the airsea module, plus rivers
        self.fwf.all_values = self.airsea.pe.all_values
        np.add.at(
            self.fwf.all_values, self.rivers.slice, self.rivers.flow * self.rivers.iarea
        )

        # Update elevation at the open boundaries. This must be done before
        # calculating the surface pressure gradient
        # NB from this moment on, elevations z at the open boundary will be out
        # of sync with water depths D and thicknesses hn (only on T grid).
        # This will last until the next call to update_depth!
        self.T.z.open_boundaries.update()

        # Calculate the surface pressure gradient in the U and V points.
        # Note: this requires elevation and surface air pressure (both on T grid) to be
        # valid in the halos, which is guaranteed for elevation (halo exchange happens
        # just after update), and for air pressure if it is managed by the input
        # manager (e.g. read from file)
        self.airsea.sp.update_halos(parallel.Neighbor.TOP_AND_RIGHT)
        self.update_surface_pressure_gradient(self.T.z, self.airsea.sp)

        # Interpolate surface stresses from T to U and V grids
        self.airsea.taux.update_halos(parallel.Neighbor.RIGHT)
        self.airsea.taux.interp(self.tausx)
        self.airsea.tauy.update_halos(parallel.Neighbor.TOP)
        self.airsea.tauy.interp(self.tausy)

        if update_3d:
            # Save surface forcing variables for the next macro momentum update
            self.tausxo.all_values = self.tausx.all_values
            self.tausyo.all_values = self.tausy.all_values
            self.dpdxo.all_values = self.dpdx.all_values
            self.dpdyo.all_values = self.dpdy.all_values

            # Update surface shear velocity (used by GOTM). This requires updated
            # surface stresses and there can only be done after the airsea update.
            _pygetm.surface_shear_velocity(
                self.airsea.taux, self.airsea.tauy, self.ustar_s
            )

            if self.radiation is not None:
                # Update radiation in the interior.
                # This must come after the airsea update, which is responsible for
                # calculating downwelling shortwave radiation at the water surface (swr)
                self.radiation(self.airsea.swr, self.fabm.kc if self.fabm else None)

            # If we need vertical tracer diffusivity at layer centers (for FABM),
            # calculate it by interpolating diffusivity at the layer interfaces.
            if self.nuh_ct is not None:
                self.turbulence.nuh.interp(self.nuh_ct)

            # Update source terms of biogeochemistry, using the new tracer
            # concentrations. Do this last because FABM could depend on any of the
            # variables computed before (radiation, diffusivity, etc.)
            if self.fabm:
                self.fabm.update_sources(self.timestep * self.istep, self.time)

        if self.report_totals != 0 and self.istep % self.report_totals == 0:
            self.report_domain_integrals()

    def _summarize_profiling_result(self, ps: pstats.Stats):
        sp = ps.get_stats_profile()
        if "<built-in method Waitall>" not in sp.func_profiles:
            # not a parallel simulation, or advance was never called
            return
        stat = [
            sp.total_tt,
            sp.func_profiles["<built-in method Waitall>"].tottime,
            self.T._water_nohalo.sum(),
        ]
        all_stat = self.tiling.comm.gather(stat)
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
        slc = self.rivers.slice

        # Increase in water depth (m)
        z_add = self._int_river_flow * self.rivers.iarea

        # Determine the part of every layer over which inflow or outflow occurs
        h_old = self.T.hn.all_values[slc]
        h_active = self.rivers.get_active_part_of_layers(h_old)

        # Change in thickness per layer
        h_add = h_active * (z_add / h_active.sum(axis=0))
        np.add.at(self.T.hn.all_values, slc, h_add)

        # Update tracers
        # We operate in depth-integrated tracer units to support
        # summing contributions from different rivers to the same cell
        h_new_inv = 1.0 / self.T.hn.all_values[slc]
        for tracer in self.tracers:
            follow = tracer.river_follow | (z_add < 0.0)
            tracer_old = tracer.all_values[slc]
            tracer_add = np.where(follow, tracer_old, tracer.river_values) * h_add
            tracer.all_values[slc] = tracer_old * h_old
            np.add.at(tracer.all_values, slc, tracer_add)
            tracer.all_values[slc] *= h_new_inv

        # Precipitation and evaporation (surface layer only)
        # First update halos for the net freshwater flux, as we need to ensure that the
        # layer heights updated as a result remain valid in the halos.
        self.airsea.pe.update_halos()
        unmasked = self.T._water
        z_add_fwf = np.where(unmasked, self.airsea.pe.all_values, 0.0) * timestep
        h_sf = self.T.hn.all_values[-1, :, :]
        h_sf_new = h_sf + z_add_fwf
        dilution = h_sf / h_sf_new
        for tracer in self.tracers:
            if not tracer.precipitation_follows_target_cell:
                tracer.all_values[-1, :, :] *= dilution
        h_sf[:, :] = h_sf_new

        # Update elevation (first add river contribution to z_add_fwf)
        np.add.at(z_add_fwf, slc, z_add)
        self.T.zin.all_values += z_add_fwf

        # Start tracer halo exchange (to prepare for advection)
        for tracer in self.tracers:
            tracer.update_halos_start(self.tracers._advection.halo1)

        self._int_river_flow.fill(0.0)

    @property
    def totals(
        self,
    ) -> tuple[
        Optional[float],
        Optional[Sequence[tuple[pygetm.tracer.TracerTotal, float, float]]],
    ]:
        """Global totals of volume and tracers.

        Returns:
            A tuple with total volume and a list with (tracer_total, total, mean)
            tuples on the root subdomains. On non-root subdomains it returns None, None
        """
        unmasked = self.T.mask == CellType.ACTIVE
        total_volume = (self.T.D * self.T.area).global_sum(where=unmasked)
        if any(tt.per_mass for tt in self.tracer_totals):
            vol = self.T.hn * self.T.area
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
            total = total.global_sum(where=grid.mask == CellType.ACTIVE)
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
        This is done by :meth:`~update_depth` instead.
        """
        self.T.zo.all_values[:, :] = self.T.z.all_values
        _pygetm.advance_surface_elevation(timestep, self.T.z, U, V, fwf)
        self.T.z.update_halos()

    def update_surface_pressure_gradient(self, z: core.Array, sp: core.Array):
        _pygetm.surface_pressure_gradient(z, sp, self.dpdx, self.dpdy, self.Dmin)

    @property
    def Ekin(self, rho0: float = RHO0):
        U = self.momentum.U.interp(self.T)
        V = self.momentum.V.interp(self.T)
        vel2_D2 = U**2 + V**2
        return 0.5 * rho0 * self.T.area * vel2_D2 / self.T.D

    def update_depth(self, _3d: bool = False, timestep: float = 0.0):
        """Use old and new surface elevation on T grid to calculate elevations on
        U, V, X grids. Subsequently update total water depth ``D`` on all grids.
        If ``_3d`` is True, also update layer thicknesses ``hn``, layer center
        depths ``zc`` and interface depths ``zf`` on all grids.

        Args:
            _3d: update elevations of the macrotimestep (``zin``) rather than
                elevations of the microtimestep (``z``).This first synchronizes the
                elevations of the macrotimestep on the T grid (``self.T.zin``) with
                those of the microtimestep (``self.T.z``). It also updates layer
                thicknesses ``hn``, layer center depths ``zc`` and interface depths
                ``zf`` on all grids.
            timestep: time step (s) for layer thickness change
                if 0, any layer height relaxation is disabled

        This routine will ensure values are up to date in the domain interior and in
        the halos, but that this requires that ``self.T.z`` (and old elevations
        ``self.T.zo`` or ``self.T.zio``) are already up to date in halos.
        """
        if _3d:
            # Store current elevations as previous elevations (on the 3D time step)
            self.T.zio.all_values = self.T.zin.all_values

            # Synchronize new elevations on the 3D time step to those of the 2D time
            # step that has just completed.
            self.T.zin.all_values = self.T.z.all_values

            z_T, zo_T = self.T.zin, self.T.zio
        else:
            z_T, zo_T = self.T.z, self.T.zo

        z_U, z_V, z_X, z_T_half = self.U._work, self.V._work, self.X._work, self.T._work

        # Update total water depth D on T grid
        # This is the only grid where we track raw depth (possibly < Dmin)
        # as well as clipped depth max(D, Dmin). On other grids we only track
        # clipped depth.
        _pygetm.elevation2depth(z_T, self.T.H, -1000.0, self.T.D)
        np.maximum(self.T.D.all_values, self.Dmin, out=self.T.Dclip.all_values)

        # For water depths on U, V, X grids we need elevations that lag 1/2 a
        # timestep behind the T grid. Calculate these by averaging old and
        # new elevations on the T grid.
        np.add(zo_T.all_values, z_T.all_values, out=z_T_half.all_values)
        z_T_half.all_values *= 0.5

        # Total water depth D on U grid
        z_T_half.interp(z_U)
        z_T_half.mirror(z_U)
        _pygetm.elevation2depth(z_U, self.U.H, self.Dmin, self.U.D)

        # Total water depth D on V grid
        z_T_half.interp(z_V)
        z_T_half.mirror(z_V)
        _pygetm.elevation2depth(z_V, self.V.H, self.Dmin, self.V.D)

        # Total water depth D on X grid
        z_T_half.interp(z_X)
        _pygetm.elevation2depth(z_X, self.X.H, self.Dmin, self.X.D)

        # Halo exchange for water depth on U, V grids, needed because the very last
        # points in the halos (x=-1 for U, y=-1 for V) are not valid after
        # interpolating elevation from the T grid above.
        # These depths are needed to later compute velocities from transports
        # These velocities will be advected, and therefore need to be valid througout
        # the halos. We do not need to halo-exchange elevation on the X grid, since
        # that needs to be be valid at the innermost halo point only, which is ensured
        # by z_T exchange.
        self.U.D.update_halos(parallel.Neighbor.RIGHT)
        self.V.D.update_halos(parallel.Neighbor.TOP)

        # Update dampening factor (0-1) for shallow water
        _pygetm.alpha(self.U.D, 2 * self.Dmin, self.Dcrit, self.U.alpha)
        _pygetm.alpha(self.V.D, 2 * self.Dmin, self.Dcrit, self.V.alpha)

        # Update total water depth on advection grids. These must be 1/2 timestep
        # behind the T grid. That's already the case for the X grid, but for the T grid
        # we explicitly compute and use the average of old and new D.
        D_T_half = self.T._work
        _pygetm.elevation2depth(z_T_half, self.T.H, self.Dmin, D_T_half)
        self.U.ugrid.D.all_values[:, :-1] = D_T_half.all_values[:, 1:]
        self.V.vgrid.D.all_values[:-1, :] = D_T_half.all_values[1:, :]
        self.U.vgrid.D.all_values[:, :] = self.X.D.all_values[1:, 1:]
        self.V.ugrid.D.all_values[:, :] = self.U.vgrid.D.all_values

        if _3d:
            # Store previous layer thicknesses
            # NB on U and V grids, ho is needed to estimate thicknesses
            # in between start and stop of the timestep (i.e., in sync with T grid)
            # These are used in the momentum update
            self.T.ho.all_values = self.T.hn.all_values
            self.U.ho.all_values = self.U.hn.all_values
            self.V.ho.all_values = self.V.hn.all_values

            # Update layer thicknesses (hn) on all grids, using bathymetry H and new
            # elevations zin (on the 3D timestep)
            self.vertical_coordinates.update(timestep)

            # Update vertical coordinates, used for e.g., output, internal pressure,
            # vertical interpolation of open boundary forcing of tracers
            for grid in (self.T, self.U, self.V):
                _pygetm.thickness2vertical_coordinates(
                    grid.mask, grid.H, grid.hn, grid.zc, grid.zf
                )

            # Update thicknesses on advection grids. These must be at time=n+1/2
            # That's already the case for the X grid, but for the T grid (now at t=n+1)
            # we explicitly compute thicknesses at time=n+1/2.
            # Note that UU.hn and VV.hn will miss the x=-1 and y=-1 strips,
            # respectively (the last strip of values within their halos);
            # fortunately these values are not needed for advection.
            self.h_T_half.all_values = 0.5 * (
                self.T.ho.all_values + self.T.hn.all_values
            )
            self.U.ugrid.hn.all_values[:, :, :-1] = self.h_T_half.all_values[:, :, 1:]
            self.V.vgrid.hn.all_values[:, :-1, :] = self.h_T_half.all_values[:, 1:, :]
            self.U.vgrid.hn.all_values[:, :, :] = self.X.hn.all_values[:, 1:, 1:]
            self.V.ugrid.hn.all_values[:, :, :] = self.U.vgrid.hn.all_values

            if self.depth.saved:
                # Update depth-below-surface at layer centers.
                # Elsewhere this can be used as approximate pressure in dbar
                _pygetm.thickness2center_depth(self.T.mask, self.T.hn, self.depth)

            # Update vertical coordinate at open boundary, used to interpolate
            # inputs on z grid to dynamic model depths
            if self.open_boundaries.zc.saved:
                self.open_boundaries.zc.all_values = self.T.zc.all_values[
                    :, self.open_boundaries.j, self.open_boundaries.i
                ].T
            if self.open_boundaries.zf.saved:
                self.open_boundaries.zf.all_values = self.T.zf.all_values[
                    :, self.open_boundaries.j, self.open_boundaries.i
                ].T
