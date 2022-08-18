from typing import Mapping, Optional, List, Sequence, NamedTuple

import numpy as np

from .constants import ZERO_GRADIENT, CENTERS, INTERFACES, TimeVarying
from . import core
from . import domain
from . import operators


class TracerTotal(NamedTuple):
    array: core.Array
    scale_factor: float = 1.0
    offset: float = 0.0
    units: Optional[str] = None
    long_name: Optional[str] = None
    per_mass: bool = False


class OpenBoundaries:
    __slots__ = "_tracer", "_type", "values"

    def __init__(self, tracer: "Tracer"):
        self._tracer = tracer
        self._type = ZERO_GRADIENT
        self.values: Optional[core.Array] = None

    @property
    def type(self) -> int:
        """Type of open boundary condition"""
        return self._type

    @type.setter
    def type(self, value: int):
        self._type = value
        if self._type == ZERO_GRADIENT:
            self.values = None
        elif self.values is None:
            self.values = self._tracer.grid.array(
                name="%s_bdy" % self._tracer.name,
                z=CENTERS,
                on_boundary=True,
                attrs={"_time_varying": TimeVarying.MACRO},
            )

    def update(self):
        """Update the tracer at the open boundaries"""
        self._tracer.update_boundary(self._type, self.values)


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
        **kwargs
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
        self.open_boundaries: OpenBoundaries = OpenBoundaries(self)
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
        advection_scheme: operators.AdvectionScheme,
        cnpar: float = 1.0,
    ):
        self.grid: domain.Grid = grid
        self._tracers: List[Tracer] = []
        self._source = grid.array(z=CENTERS)
        self._advection = operators.Advection(grid, scheme=advection_scheme)
        self._vertical_diffusion = operators.VerticalDiffusion(grid, cnpar=cnpar)

    def __getitem__(self, index: int) -> Tracer:
        return self._tracers[index]

    def __len__(self) -> int:
        return len(self._tracers)

    def add(self, **kwargs) -> Tracer:
        """Add a tracer that will be subject to advection and diffusion.

        Args:
            **kwargs: keyword arguments to be passed to :class:`Tracer`
        """
        tracer = Tracer(grid=self.grid, **kwargs)
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
            u: velocity in x direction (m s-1)
            v: velocity in y direction (m s-1)
            w: velocity in z direction (m s-1)
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
        w_tracers = []
        for tracer in self._tracers:
            w_tracer = w
            if (
                tracer.vertical_velocity is not None
                and tracer.vertical_velocity.all_values.any()
            ):
                w_tracer = self.grid.array(fill=0.0, z=INTERFACES)
                tracer.vertical_velocity.interp(w_tracer)
                w_tracer.all_values += w.all_values
            w_tracers.append(w_tracer)
        self._advection.apply_3d_batch(
            u, v, w, timestep, self._tracers, w_vars=w_tracers
        )

        # Vertical diffusion of passive tracers (including biogeochemical ones)
        # This simultaneously time-integrates source terms and surface fluxes, if
        # specified. However, sources of biogeochemical tracers (FABM) are typically
        # dealt with elsewhere; these will not have their source attribute set for
        # use here.
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
            self._vertical_diffusion(diffusivity, timestep, tracer, ea4=source)
