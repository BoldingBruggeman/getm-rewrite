from typing import Optional, Mapping, Union, Iterable
import operator
import logging

import numpy as np

from . import core
from . import parallel
from .constants import CoordinateType, CellType


class RiverTracer(core.Array):
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
        return bool(self._follow)

    @follow_target_cell.setter
    def follow_target_cell(self, value: bool):
        self._follow[...] = value


class GlobalRiver:
    def __init__(
        self,
        name: str,
        x: Union[int, float],
        y: Union[int, float],
        zl: Optional[float] = np.inf,
        zu: Optional[float] = 0.0,
        coordinate_type: CoordinateType = CoordinateType.IJ,
    ):
        """
        Args:
            name: unique name for this river
            x: x coordinate of river
            y: y coordinate of river
            zl: maximum depth to which the river penetrates (non-negative)
            zu: minimum depth from which the river penetrates (non-negative)
            coordinate_type: coordinate type of x and y
                (LONLAT spherical, XY for Cartesian coordinates)
        """
        self.name = name
        self.x = x
        self.y = y
        self.coordinate_type = coordinate_type
        self.zl = zl
        self.zu = zu
        self.i: Optional[int] = None
        self.j: Optional[int] = None

    def locate(self, locator: core.Locator):
        """If this river position is specified by (lon, lat) or (x, y), map it
        to the nearest non-masked grid cell."""
        if self.coordinate_type == CoordinateType.IJ:
            self.i, self.j = int(round(self.x)), int(round(self.y))
        else:
            self.i, self.j = locator(
                self.x,
                self.y,
                coordinate_type=self.coordinate_type,
                valid_cell_types=(CellType.ACTIVE,),
            )

    def to_local_grid(self, grid: core.Grid) -> Optional["LocalRiver"]:
        """Map global river to local subdomain."""
        i_loc, j_loc = grid.global_to_local(self.i, self.j, include_halos=True)
        if i_loc is None or j_loc is None:
            return None
        river = LocalRiver(grid, self.name, i_loc, j_loc, zl=self.zl, zu=self.zu)
        for att in ("original_name", "split"):
            if hasattr(self, att):
                setattr(river, att, getattr(self, att))
        return river


class LocalRiver(Mapping[str, RiverTracer]):
    def __init__(
        self, grid: core.Grid, name: str, i: int, j: int, zl: float, zu: float
    ):
        self.name = name
        self.i = i
        self.j = j
        self.zl = zl
        self.zu = zu
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


class LocalRiverCollection(Mapping[str, LocalRiver]):
    def __init__(
        self, grid: core.Grid, rivers: Iterable[LocalRiver], logger: logging.Logger
    ):
        self._rivers = {river.name: river for river in rivers}
        self.logger = logger

        self.flow = np.zeros((len(rivers),))
        self.zl = np.array([river.zl for river in rivers])
        self.zu = np.array([river.zu for river in rivers])
        for iriver, river in enumerate(rivers):
            river.flow.wrap_ndarray(self.flow[..., iriver])
            river.zl = self.zl[..., iriver]
            river.zu = self.zu[..., iriver]
        self.i = np.array([river.i for river in rivers], dtype=np.intp)
        self.j = np.array([river.j for river in rivers], dtype=np.intp)
        self.iarea = grid.iarea.all_values[self.j, self.i]

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


class GlobalRiverCollection(Mapping[str, GlobalRiver]):
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
            **kwargs: additional keyword arguments passed to :class:`GlobalRiver`

        Returns:
            river instance
        """
        assert i >= 0 and i < self.nx
        assert j >= 0 and j < self.ny
        return self.add_by_location(
            name, i, j, coordinate_type=CoordinateType.IJ, **kwargs
        )

    def add_by_location(
        self,
        name: str,
        x: Union[int, float],
        y: Union[int, float],
        coordinate_type: Optional[CoordinateType] = None,
        **kwargs,
    ) -> GlobalRiver:
        """Add a river at a location specified by the nearest coordinates

        Args:
            name: river name
            x: x coordinate of river
            y: y coordinate of river
            coordinate_type: coordinate type of x and y
                (LONLAT for spherical, XY for Cartesian coordinates,
                IJ for indices into the global tracer grid)
            **kwargs: additional keyword arguments passed to :class:`GlobalRiver`

        Returns:
            river instance
        """
        if coordinate_type is None:
            coordinate_type = self.default_coordinate_type
        river = GlobalRiver(name, x, y, coordinate_type=coordinate_type, **kwargs)
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
        """Freeze the river collection. Drop those outside the current subdomain
        and verify the remaining ones are on unmasked T points.
        """
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
