from typing import Mapping, Optional, List, Sequence, NamedTuple

import numpy as np

from .constants import ZERO_GRADIENT, CENTERS, INTERFACES, TimeVarying, FILL_VALUE
from . import core
from . import domain
from .open_boundaries import ArrayOpenBoundaries
from . import operators


class TracerTotal(NamedTuple):
    array: core.Array
    scale_factor: float = 1.0
    offset: float = 0.0
    units: Optional[str] = None
    long_name: Optional[str] = None
    per_mass: bool = False


class Tracer(core.Array):
    __slots__ = (
        "source",
        "surface_flux",
        "source_scale",
        "vertical_velocity",
        "open_boundaries",
        "river_values",
        "river_follow",
        "rivers",
        "precipitation_follows_target_cell",
        "molecular_diffusivity",
    )

    def __init__(
        self,
        grid: domain.Grid,
        data: Optional[np.ndarray] = None,
        source: Optional[core.Array] = None,
        surface_flux: Optional[core.Array] = None,
        source_scale: float = 1.0,
        vertical_velocity: Optional[core.Array] = None,
        rivers_follow_target_cell: bool = False,
        precipitation_follows_target_cell: bool = False,
        molecular_diffusivity: float = 0.0,
        **kwargs,
    ):
        """A tracer transported by advection and diffusion, with optional source term,
        surface flux and vertical velocity (e.g., sinking, floating)

        Args:
            grid: the grid on which the tracer is defined
            data: a NumPy array that will hold the values of the tracer.
                If not provided, a new one will be created.
            source: array with source terms that after multiplication with
                ``source_scale`` must have tracer units per time, multiplied by layer
                thickness. Defaults to 0.
            surface_flux: array with surface flux values that after multiplication with
                ``source_scale`` must have tracer units per time, multiplied by layer
                thickness. Defaults to 0.
            source_scale: scale factor for source terms and surface flux
            vertical_velocity: array with vertical velocities (m s-1) describing
                movement through the water, not water movement itself. Positive for
                upward movement (e.g. floating), negative for downward movement
                (e.g. sinking). Defaults to 0 (no movement independent of the water).
            rivers_follow_target_cell: tracer values in river water are assumed equal
                to those in the cell into which the river flows. This can be customized
                further by setting <TRACER>.rivers[<RIVERNAME>].follow_target_cell
                and/or <TRACER>.rivers[<RIVERNAME>].values
            precipitation_follows_target_cell: tracer values in precipitation are
                assumed equal to surface values in the cell into which the
                precipitation falls. If not set, tracer values in precipitation are
                assumed to be 0. This can be adjusted afterwards with
                :attr:`precipitation_follows_target_cell`
            molecular_diffusivity: molecular diffusivity (m2 s-1)
            **kwargs: keyword arguments to be passed to :class:`pygetm.core.Array`
        """
        kwargs.setdefault("attrs", {}).update(
            _part_of_state=True, _time_varying=TimeVarying.MACRO
        )
        super().__init__(grid=grid, shape=grid.hn.all_values.shape, **kwargs)

        if data is None:
            data = np.full_like(grid.hn.all_values, np.nan)
        self.wrap_ndarray(data)

        assert source is None or (source.grid is self.grid and source.z == CENTERS)
        assert surface_flux is None or (
            surface_flux.grid is self.grid and not surface_flux.z
        )
        assert vertical_velocity is None or (
            vertical_velocity.grid is self.grid and vertical_velocity.z == CENTERS
        )

        self.source: Optional[core.Array] = source
        self.surface_flux: Optional[core.Array] = surface_flux
        self.source_scale: float = source_scale
        self.vertical_velocity: Optional[core.Array] = vertical_velocity
        self.molecular_diffusivity: float = molecular_diffusivity
        self.open_boundaries: ArrayOpenBoundaries = ArrayOpenBoundaries(
            self, ZERO_GRADIENT
        )
        self.river_values: np.ndarray = np.zeros((len(grid.domain.rivers),))
        self.river_follow: np.ndarray = np.full(
            (len(grid.domain.rivers),), rivers_follow_target_cell, dtype=bool
        )
        self.precipitation_follows_target_cell: bool = precipitation_follows_target_cell
        self.rivers: Mapping[str, domain.RiverTracer] = {}
        for iriver, river in enumerate(grid.domain.rivers.values()):
            river_tracer = domain.RiverTracer(
                grid,
                river.name,
                self.name,
                self.river_values[..., iriver],
                self.river_follow[..., iriver],
                units=self.units,
                attrs={"_time_varying": TimeVarying.MACRO},
            )
            river._tracers[self.name] = river_tracer
            self.rivers[river.name] = river_tracer


class TracerCollection(Sequence[Tracer]):
    def __init__(
        self,
        grid: domain.Grid,
        advection_scheme: operators.AdvectionScheme = operators.AdvectionScheme.DEFAULT,
        cnpar: float = 1.0,
    ):
        self.logger = grid.domain.root_logger.getChild("tracers")
        self.logger.info(f"Advection scheme: {advection_scheme.name}")
        self.logger.info(f"Crank-Nicolson parameter: {cnpar}")

        self.grid: domain.Grid = grid
        self._tracers: List[Tracer] = []
        self._source = grid.array(z=CENTERS)
        self._advection = operators.Advection(grid, scheme=advection_scheme)
        self._vertical_diffusion = operators.VerticalDiffusion(grid, cnpar=cnpar)
        self._w = grid.array(fill=0.0, z=INTERFACES)

        self.Ah = grid.array(
            name="Ah",
            units="m2 s-1",
            long_name="horizontal diffusivity of tracers",
            fill_value=FILL_VALUE,
            attrs=dict(_require_halos=True, _time_varying=False),
        )
        self.Ah.fill(0.0)
        self.Ah_u = self.Ah_v = None

    def start(self):
        self.Ah.update_halos()
        self.Ah.all_values[self.grid._land] = self.Ah.fill_value
        if (self.Ah.all_values[self.grid._water] == 0.0).all():
            self.logger.info(
                "Disabling horizontal diffusion because Ah is 0 everywhere"
            )
        else:
            self.logger.info(
                f"Horizontal diffusivity Ah ranges between {self.Ah.ma.min()}"
                f" and {self.Ah.ma.max()} m2 s-1"
            )
            self.Ah_u = self.grid.ugrid.array(fill=np.nan)
            self.Ah_v = self.grid.vgrid.array(fill=np.nan)
            self.Ah.interp(self.Ah_u)
            self.Ah.interp(self.Ah_v)
            self.Ah_u.update_halos()
            self.Ah_v.update_halos()

    def __getitem__(self, index: int) -> Tracer:
        return self._tracers[index]

    def __len__(self) -> int:
        return len(self._tracers)

    def add(self, name: str, **kwargs) -> Tracer:
        """Add a tracer that will be subject to advection and diffusion.

        Args:
            name: short name for the tracer (letters, digits, underscores only)
            **kwargs: keyword arguments to be passed to :class:`Tracer`

        Returns:
            tracer instance
        """
        tracer = Tracer(grid=self.grid, name=name, **kwargs)
        self._tracers.append(tracer)
        return tracer

    def advance(
        self,
        timestep: float,
        u: core.Array,
        v: core.Array,
        w: core.Array,
        diffusivity: core.Array,
    ):
        """Advance tracers through advection, diffusion and time integration of
        optional source terms and/or surface fluxes.

        Args:
            timestep: time step (s)
            u: velocity in x-direction (m s-1)
            v: velocity in y-direction (m s-1)
            w: velocity in z-direction (m s-1)
            diffusivity: vertical turbulent diffusivity (m2 s-1)

        This uses operator splitting: advection and horizontal diffusion is done first,
        vertical diffusion (and sources and surface fluxes, if provided) after.

        Velocities ``u``, ``v``, ``w`` are expected to be positioned halfway in time
        between the old and new tracer values (i.e., they should be 1/2 a timestep ahead
        of the tracers upon entry). ``diffusivity`` is expected to be positioned at the
        same time as the new tracer values, and should thus be 1/2 a timestep ahead of
        the advection velocities.
        """

        # Advection of passive tracers (including biogeochemical ones) from time=0 to
        # time=1 (macrotimestep), using velocities defined at time=1/2
        def combined_w(tracer):
            if (
                tracer.vertical_velocity is None
                or not tracer.vertical_velocity.all_values.any()
            ):
                return w
            tracer.vertical_velocity.interp(self._w)
            self._w.all_values += w.all_values
            return self._w

        self._advection.apply_3d_batch(
            u,
            v,
            w,
            timestep,
            self._tracers,
            w_vars=combined_w,
            Ah_u=self.Ah_u,
            Ah_v=self.Ah_v,
        )

        # Vertical diffusion of passive tracers (including biogeochemical ones)
        # This simultaneously time-integrates source terms and surface fluxes, if
        # specified. However, sources of biogeochemical tracers (FABM) are typically
        # dealt with elsewhere; these will not have their source attribute set for
        # use here.
        avmol = -1.0
        for tracer in self._tracers:
            source = None
            if tracer.source is not None or tracer.surface_flux is not None:
                source = self._source  # work array to hold sources
                if tracer.source is not None:
                    source.all_values[...] = tracer.source.all_values
                else:
                    source.all_values.fill(0.0)
                if tracer.surface_flux is not None:
                    source.all_values[-1, ...] += tracer.surface_flux.all_values
                source.all_values *= timestep * tracer.source_scale
            if tracer.molecular_diffusivity != avmol:
                avmol = tracer.molecular_diffusivity
                self._vertical_diffusion.prepare(diffusivity, timestep, avmol)
            self._vertical_diffusion.apply(tracer, ea4=source)
