from typing import Optional, Mapping, List
import operator
import logging

import numpy as np

from . import core
from . import parallel
from .constants import CoordinateType


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


class River:
    flow: core.Array  #: flow rate in m3 s-1

    def __init__(
        self,
        name: str,
        i: int,
        j: int,
        zl: Optional[float] = np.inf,
        zu: Optional[float] = 0.0,
        x: Optional[float] = None,
        y: Optional[float] = None,
        coordinate_type: CoordinateType = CoordinateType.IJ,
    ):
        """
        Args:
            name: unique name for this river
            i: global index in x-direction (0-based)
            j: global index in y-direction (0-based)
            zl: maximum depth to which the river penetrates (non-negative)
            zu: minimum depth from which the river penetrates (non-negative)
            x: x coordinate of river
            y: y coordinate of river
            coordinate_type: coordinate type of x and y
                (LONLAT spherical, XY for Cartesian coordinates)
        """
        self.name = name
        self.i_glob = i
        self.j_glob = j
        self.x = x
        self.y = y
        self.coordinate_type = coordinate_type
        self.zl = zl
        self.zu = zu
        self.i_loc: Optional[int] = None
        self.j_loc: Optional[int] = None
        self._tracers: Mapping[str, RiverTracer] = {}

    def locate(self, locator: core.Locator):
        """If this river position is specified by (lon, lat) or (x, y), map it
        to the nearest non-masked grid cell."""
        if self.x is not None and self.y is not None:
            self.i_glob, self.j_glob = locator(
                self.x, self.y, coordinate_type=self.coordinate_type
            )

    def to_local_grid(self, grid: core.Grid, logger: logging.Logger) -> bool:
        """Map global river position (i,j) to local subdomain."""
        self.i_loc, self.j_loc = grid.global_to_local(
            self.i_glob, self.j_glob, include_halos=True
        )
        if self.i_loc is None or self.j_loc is None:
            logger.info(f"{self.name} falls outside this subdomain")
            return False
        logger.info(f"{self.name} at is located at i={self.i_loc}, j={self.j_loc}")

        mask = grid.mask.all_values[self.j_loc, self.i_loc]
        if mask != 1:
            raise Exception(
                f"{self.name} has been mapped to non-water grid cell"
                f" (with mask value {mask})."
            )

        return True

    def initialize(self, grid: core.Grid, flow: np.ndarray):
        self.flow = core.Array(
            grid=grid,
            name="river_" + self.name + "_flow",
            units="m3 s-1",
            long_name=f"inflow from {self.name}",
        )
        self.flow.wrap_ndarray(flow)

    def __getitem__(self, key) -> RiverTracer:
        return self._tracers[key]

    def __len__(self):
        return len(self._tracers)

    def __iter__(self):
        return iter(self._tracers)


class Rivers(Mapping[str, River]):
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
        self._rivers: List[River] = []
        self.global_rivers = self._rivers
        self._frozen = False

    def add_by_index(self, name: str, i: int, j: int, **kwargs):
        """Add a river at a location specified by the indices of a tracer point

        Args:
            name: river name
            i: global domain index in x-direction (0-based)
            j: global domain index in y-direction (0-based)
            **kwargs: additional keyword arguments passed to :class:`River`

        Returns:
            river instance
        """
        assert not self._frozen, (
            "The river collection has already been initialized"
            " and can no longer be modified."
        )

        assert i >= 0 and i < self.nx
        assert j >= 0 and j < self.ny
        river = River(name, i, j, **kwargs)
        self._rivers.append(river)
        return river

    def add_by_location(
        self,
        name: str,
        x: float,
        y: float,
        coordinate_type: Optional[CoordinateType] = None,
        **kwargs,
    ):
        """Add a river at a location specified by the nearest coordinates

        Args:
            name: river name
            x: x coordinate of river
            y: y coordinate of river
            **kwargs: additional keyword arguments passed to :class:`River`
        """
        if coordinate_type is None:
            coordinate_type = self.default_coordinate_type
        river = River(
            name, None, None, x=x, y=y, coordinate_type=coordinate_type, **kwargs
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
            ind = (river.i_glob, river.j_glob) if comm.rank == 0 else None
            river.i_glob, river.j_glob = comm.bcast(ind)

    def initialize(self, grid: core.Grid):
        """Freeze the river collection. Drop those outside the current subdomain
        and verify the remaining ones are on unmasked T points.
        """
        assert not self._frozen, "The river collection has already been initialized"
        self._frozen = True

        self._broadcast_locations(grid.tiling.comm)

        # Keep only rivers that fall within the local subdomain
        self._rivers = [
            river for river in self._rivers if river.to_local_grid(grid, self.logger)
        ]

        self.flow = np.zeros((len(self._rivers),))
        for iriver, river in enumerate(self._rivers):
            river.initialize(grid, self.flow[..., iriver])
        self.i = np.array([river.i_loc for river in self._rivers], dtype=int)
        self.j = np.array([river.j_loc for river in self._rivers], dtype=int)
        self.iarea = grid.iarea.all_values[self.j, self.i]
        self.zl = np.array([river.zl for river in self._rivers])
        self.zu = np.array([river.zu for river in self._rivers])

    def flag_prescribed_tracers(self):
        for river in self._rivers:
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

    def __getitem__(self, key: str) -> River:
        for river in self._rivers:
            if key == river.name:
                return river
        raise KeyError()

    def __len__(self) -> int:
        return len(self._rivers)

    def __iter__(self):
        return map(operator.attrgetter("name"), self._rivers)
