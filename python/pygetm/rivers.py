import enum
from typing import Optional, Mapping, Union, Iterable, Any
import operator
import logging

import numpy as np

from . import core
from . import parallel
from .constants import CoordinateType, CellType


class RiverTracer(core.Array):
    """Single tracer in a single river.

    Call :meth:`pygetm.core.Array.set` on this object to prescribe the tracer
    value in the river, or set the :attr:`follow_target_cell` attribute
    to `True` to take the river's tracer value from the model cell it flows
    into.

    If you prescribe the tracer value, the :attr:`follow_target_cell` attribute
    will automatically be set to `False`.

    If you do not prescribe the tracer value and :attr:`follow_target_cell`
    is `False`, the tracer value will default to 0.0.
    """

    __slots__ = ("_follow",)

    def __init__(
        self,
        grid: core.Grid,
        river_name: str,
        tracer_name: str,
        value: np.ndarray,
        follow: np.ndarray,
        **kwargs,
    ):
        super().__init__(
            grid=grid,
            name=f"{tracer_name}_in_river_{river_name}",
            long_name=f"{tracer_name} in river {river_name}",
            **kwargs,
        )
        self.wrap_ndarray(value)
        self._follow = follow

    @property
    def follow_target_cell(self) -> bool:
        """Whether to take the tracer value in the river from the model cell it
        flows into."""
        return bool(self._follow)

    @follow_target_cell.setter
    def follow_target_cell(self, value: bool):
        self._follow[...] = value


class VerticalPosition(enum.Enum):
    DistanceFromSurface = 1
    DistanceFromBottom = 2


class GlobalRiver:
    """Single river in the global domain."""

    def __init__(
        self,
        name: str,
        x: Union[int, float],
        y: Union[int, float],
        coordinate_type: CoordinateType = CoordinateType.IJ,
        **attrs,
    ):
        """
        Args:
            name: unique name for this river
            x: x coordinate of river
            y: y coordinate of river
            coordinate_type: coordinate type of x and y
                (LONLAT spherical, XY for Cartesian coordinates)
            **attrs: additional attributes for this river
        """
        self.name = name
        self.x = x
        self.y = y
        self.coordinate_type = coordinate_type
        self.attrs = attrs
        self.i: Optional[int] = None
        self.j: Optional[int] = None
        if self.coordinate_type == CoordinateType.IJ:
            self.i, self.j = int(round(self.x)), int(round(self.y))

    def locate(self, locator: core.Locator):
        """If this river position is specified by (lon, lat) or (x, y), map it
        to the nearest non-masked grid cell."""
        if self.coordinate_type != CoordinateType.IJ:
            self.i, self.j = locator(
                self.x,
                self.y,
                coordinate_type=self.coordinate_type,
                valid_cell_types=(CellType.ACTIVE,),
            )

    def to_local_grid(self, grid: core.Grid) -> Optional["LocalRiver"]:
        """Map river to local subdomain.

        Args:
            grid: local grid
        Returns:
            local river instance, or None if the river falls outside the local subdomain
        """
        i_loc, j_loc = grid.global_to_local(self.i, self.j, include_halos=True)
        if i_loc is None or j_loc is None:
            return None
        return LocalRiver(grid, self.name, i_loc, j_loc, **self.attrs)

    def __getattr__(self, name: str) -> Any:
        return self.attrs[name]


class LocalRiver(Mapping[str, RiverTracer]):
    """Single river in the local subdomain.

    It acts as a mapping from tracer names to :class:`RiverTracer` instances,
    allowing you to access and control the value of each tracer in this river."""

    def __init__(
        self,
        grid: core.Grid,
        name: str,
        i: int,
        j: int,
        zl: float,
        zu: float,
        vertical_position: VerticalPosition,
        **attrs,
    ):
        self.name = name
        self.i = i
        self.j = j
        self.zl = zl
        self.zu = zu
        self.vertical_position = vertical_position
        self.attrs = attrs
        self._tracers: Mapping[str, RiverTracer] = {}
        self.flow = core.Array(
            grid=grid,
            name=f"river_{name}_flow",
            units="m3 s-1",
            long_name=f"inflow from {name}",
        )

    def __getitem__(self, key) -> RiverTracer:
        return self._tracers[key]

    def __len__(self):
        return len(self._tracers)

    def __iter__(self):
        return iter(self._tracers)

    def __getattr__(self, name: str) -> Any:
        return self.attrs[name]


class LocalRiverCollection(Mapping[str, LocalRiver]):
    """Collection of rivers that fall within the local subdomain.

    It acts as a mapping from river names to :class:`LocalRiver` instances
    """

    def __init__(
        self, grid: core.Grid, rivers: Iterable[LocalRiver], logger: logging.Logger
    ):
        self._rivers = {river.name: river for river in rivers}
        self.logger = logger

        self.flow = np.zeros((len(rivers),))
        self.zl = np.array([river.zl for river in rivers])
        self.zu = np.array([river.zu for river in rivers])
        pos = np.array([river.vertical_position for river in rivers])
        self._relative_to_surface = pos == VerticalPosition.DistanceFromSurface
        for iriver, river in enumerate(rivers):
            river.flow.wrap_ndarray(self.flow[..., iriver])
            river.zl = self.zl[..., iriver]
            river.zu = self.zu[..., iriver]
        self.i = np.array([river.i for river in rivers], dtype=np.intp)
        self.j = np.array([river.j for river in rivers], dtype=np.intp)
        self.slice = (Ellipsis, self.j, self.i)
        self.iarea = grid.iarea.all_values[self.slice]

    def __getitem__(self, key: str) -> LocalRiver:
        return self._rivers[key]

    def __len__(self) -> int:
        return len(self._rivers)

    def __iter__(self):
        return iter(self._rivers)

    def flag_prescribed_tracers(self):
        for river in self._rivers.values():
            for rt in river._tracers.values():
                prescribed = rt.values != rt.fill_value
                if prescribed and rt.follow_target_cell:
                    self.logger.warning(
                        f"Values for {rt.name} are prescribed."
                        " Disabling follow_target_cell."
                    )
                    rt.follow_target_cell = False
                elif not prescribed and not rt.follow_target_cell:
                    self.logger.warning(
                        f"Value for {rt.name} not set. Using default of 0.0"
                    )
                    rt.values[...] = 0.0

    def get_active_part_of_layers(self, h: np.ndarray) -> np.ndarray:
        """For each layer of each river cell, get the part (m) affected by the river

        Args:
            h: layer thicknesses for each river cell (shape: nz x nrivers)
        Returns:
            part of each layer affected by the river (shape: nz x nrivers)
        """
        # Vertical position of layer interfaces (distance from bottom)
        hcum_if = np.zeros((h.shape[0] + 1, h.shape[1]))
        h.cumsum(axis=0, out=hcum_if[1:, :])

        # Lower and upper limit of river penetration for each river cell
        # (distance from bottom, so lower < upper)
        D = hcum_if[-1]
        lower_limit = np.where(self._relative_to_surface, D - self.zl, self.zl)
        upper_limit = np.where(self._relative_to_surface, D - self.zu, self.zu)
        lower_limit.clip(0.0, D - 1e-6, out=lower_limit)
        upper_limit.clip(lower_limit + 1e-6, D, out=upper_limit)

        hcum_bot = np.maximum(lower_limit, hcum_if[:-1, :])
        hcum_top = np.minimum(upper_limit, hcum_if[1:, :])
        return np.maximum(hcum_top - hcum_bot, 0.0)


class GlobalRiverCollection(Mapping[str, GlobalRiver]):
    """Collection of rivers in the global domain.

    It acts as a mapping from river names to :class:`GlobalRiver` instances.
    """

    def __init__(
        self,
        nx: int,
        ny: int,
        default_coordinate_type: CoordinateType,
        logger: logging.Logger,
    ):
        self.nx = nx
        self.ny = ny
        self.default_coordinate_type = default_coordinate_type
        self.logger = logger
        self._rivers: list[GlobalRiver] = []

    def add_by_index(self, name: str, i: int, j: int, **kwargs) -> GlobalRiver:
        """Add a river at a location specified by the indices of a tracer point

        Args:
            name: river name
            i: global domain index in x-direction (0-based)
            j: global domain index in y-direction (0-based)
            **kwargs: additional keyword arguments passed to :meth:`add_by_location`

        Returns:
            river instance
        """
        return self.add_by_location(
            name, i, j, coordinate_type=CoordinateType.IJ, **kwargs
        )

    def add_by_location(
        self,
        name: str,
        x: Union[int, float],
        y: Union[int, float],
        coordinate_type: Optional[CoordinateType] = None,
        zl: Optional[float] = None,
        zu: Optional[float] = None,
        vertical_position: VerticalPosition = VerticalPosition.DistanceFromSurface,
        **attrs,
    ) -> GlobalRiver:
        """Add a river at a location specified by the nearest coordinates.

        By default, the river will spread over the entire water column.
        Distance limits `zl` and `zu` can be specified to control the vertical
        extent of river penetration. They can be identical,
        in which case the river will only affect a single layer. For instance,
        use `zl=zu=0.0` to have the river enter the surface layer, or `zl=zu=0.0`
        with :attr:`vertical_position=VerticalPosition.DistanceFromBottom`
        to have the river enter the bottom layer. If `zu` and/or `zl` exceed
        the current water depth, they will be clipped to that water depth.

        Args:
            name: river name
            x: x coordinate of river
            y: y coordinate of river
            coordinate_type: coordinate type of x and y
                (:attr:`CoordinateType.LONLAT` for spherical, :attr:`CoordinateType.XY`
                for Cartesian coordinates, :attr:`CoordinateType.IJ` for 0-based indices
                into the global tracer grid). If not provided, the default coordinate
                type of the domain will be used.
            zl: lower limit (deepest point) of river penetration (m; >=0).
                This is the distance from the surface or bottom, depending on `vertical_position`.
                Its default value corresponds to the position of the bottom.
            zu: upper limit of river penetration (m; >=0).
                This is the distance from the surface or bottom, depending on `vertical_position`.
                Its default value corresponds to the position of the water surface.
            vertical_position: whether depth limits `zl` and `zu` are distances
                from the surface or from the bottom. If set to
                :attr:`VerticalPosition.DistanceFromSurface`, then `zl` must be >= `zu`.
                If set to :attr:`VerticalPosition.DistanceFromBottom`, then `zl` must be
                <= `zu`.
            **attrs: additional attributes for this river

        Returns:
            river instance
        """
        if coordinate_type is None:
            coordinate_type = self.default_coordinate_type
        if not isinstance(coordinate_type, CoordinateType):
            raise ValueError(
                f"coordinate_type must be a CoordinateType, not {type(coordinate_type)}"
            )
        if coordinate_type == CoordinateType.IJ:
            x = int(round(x))
            y = int(round(y))
            assert x > -self.nx and x < self.nx
            assert y > -self.ny and y < self.ny
            x = x % self.nx
            y = y % self.ny

        if vertical_position == VerticalPosition.DistanceFromBottom:
            zl = 0.0 if zl is None else zl
            zu = np.inf if zu is None else zu
            if zl > zu:
                raise ValueError("For DistanceFromBottom, zl must be <= zu")
        else:
            zl = np.inf if zl is None else zl
            zu = 0.0 if zu is None else zu
            if zl < zu:
                raise ValueError("For DistanceFromSurface, zl must be >= zu")
        if zl < 0.0:
            raise ValueError("zl must be non-negative")
        if zu < 0.0:
            raise ValueError("zu must be non-negative")

        river = GlobalRiver(
            name,
            x,
            y,
            coordinate_type=coordinate_type,
            zl=zl,
            zu=zu,
            vertical_position=vertical_position,
            **attrs,
        )
        self._rivers.append(river)
        return river

    def map_to_grid(self, locator: core.Locator):
        """Map rivers to cell centers.
        This can only be called on MPI nodes that have the full domain
        (typically the root node only).
        """
        for river in self._rivers:
            river.locate(locator)

    def _broadcast_locations(self, comm: parallel.MPI.Comm):
        """Broadcast global river locations (i,j) to all non-root MPI nodes."""
        for river in self._rivers:
            ind = (river.i, river.j) if comm.rank == 0 else None
            river.i, river.j = comm.bcast(ind)

    def initialize(self, grid: core.Grid) -> LocalRiverCollection:
        """Return a collection of only those rivers that fall within the local subdomain."""
        self._broadcast_locations(grid.tiling.comm)

        # Keep only rivers that fall within the local subdomain
        local_rivers = []
        for global_river in self._rivers:
            river = global_river.to_local_grid(grid)
            if river is not None:
                self.logger.info(
                    f"{river.name} at is located at i={river.i}, j={river.j}"
                )
                mask = grid.mask.all_values[river.j, river.i]
                if mask != CellType.ACTIVE:
                    raise Exception(
                        f"{river.name} has been mapped to non-water grid cell"
                        f" (with mask value {mask})."
                    )
                local_rivers.append(river)
            else:
                self.logger.info(f"{global_river.name} falls outside this subdomain")

        return LocalRiverCollection(grid, local_rivers, self.logger)

    def __getitem__(self, key: str) -> GlobalRiver:
        for river in self._rivers:
            if key == river.name:
                return river
        raise KeyError()

    def __len__(self) -> int:
        return len(self._rivers)

    def __iter__(self):
        return map(operator.attrgetter("name"), self._rivers)
