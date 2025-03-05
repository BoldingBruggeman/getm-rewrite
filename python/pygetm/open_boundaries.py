import enum
from typing import (
    Optional,
    Sequence,
    List,
    Mapping,
    Callable,
    Iterable,
    Tuple,
    Union,
    Set,
    TYPE_CHECKING,
)
import functools
import logging

import numpy as np

from pygetm import core
from .constants import (
    CENTERS,
    GRAVITY,
    INTERFACES,
    ZERO_GRADIENT,
    SPONGE,
    FLATHER_ELEV,
    FLATHER_TRANSPORT,
    CLAMPED,
    TimeVarying,
    FILL_VALUE,
)

if TYPE_CHECKING:
    from .domain import Domain


class Side(enum.IntEnum):
    WEST = 1
    NORTH = 2
    EAST = 3
    SOUTH = 4


class OpenBoundary:
    def __init__(
        self,
        name: str,
        side: Side,
        l: int,
        mstart: int,
        mstop: int,
        type_2d: int,
        type_3d: int,
    ):
        side = Side(side)

        self.name = name
        self.side = side
        self.l_glob = l
        self.mstart_glob = mstart
        self.mstop_glob = mstop
        self.mstep = 1 if mstop > mstart else -1
        self.type_2d = type_2d
        self.type_3d = type_3d
        self.inflow_sign = 1.0 if side in (Side.WEST, Side.SOUTH) else -1.0

    def __repr__(self) -> str:
        return (
            f"{__class__.__name__}({self.name!r}, {self.side.name},"
            f" {self.l_glob}, {self.mstart_glob}, {self.mstop_glob})"
        )

    def to_local_grid(self, grid: core.Grid):
        # NB below we convert to indices in the T grid of the current subdomain
        # INCLUDING halos
        # We also limit the indices to the range valid for the current subdomain.
        side = self.side
        xoffset = grid.tiling.xoffset - grid.halox
        yoffset = grid.tiling.yoffset - grid.haloy
        if side in (Side.WEST, Side.EAST):
            l_offset, m_offset, l_max, m_max = (xoffset, yoffset, grid.nx_, grid.ny_)
        else:
            l_offset, m_offset, l_max, m_max = (yoffset, xoffset, grid.ny_, grid.nx_)
        self.l = self.l_glob - l_offset
        self.mstart_ = self.mstart_glob - m_offset
        self.mstop_ = self.mstop_glob - m_offset
        self.mstart = min(max(0, self.mstart_), m_max)
        self.mstop = min(max(0, self.mstop_), m_max)
        if self.l < 0 or self.l >= l_max or self.mstart == self.mstop:
            # Boundary lies completely outside current subdomain. Record it anyway,
            # so we can later set up a global -> local map of open boundary points
            self.l, self.mstart, self.mstop = None, None, None
            return

        mslice = slice(self.mstart, self.mstop, self.mstep)
        ms = np.arange(self.mstart, self.mstop, self.mstep)
        self.np = ms.size
        l = self.l
        ls = np.full_like(ms, l)
        if side in (Side.WEST, Side.EAST):
            self.i, self.j = ls, ms
            self.slice_t = (Ellipsis, mslice, l)
            self.slice_uv_in = (Ellipsis, mslice, l if side == Side.WEST else l - 1)
        else:
            self.i, self.j = ms, ls
            self.slice_t = (Ellipsis, l, mslice)
            self.slice_uv_in = (Ellipsis, l if side == Side.SOUTH else l - 1, mslice)

    def extract_inward(
        self, values: np.ndarray, start: int, stop: Optional[int] = None
    ) -> np.ndarray:
        """Extract one or more rows/columns on the T grid, inward from the open boundary.

        Args:
            values: 2D or 3D array of values on the T grid
            start: index of the first row/column to extract
            stop: index of the last row/column to extract (exclusive)

        Returns:
            Extracted of values on the T grid.
            If a single row/column is extracted, the shape of the extracted
            values will be ``(..., np)``, where ``np`` is the number of points of the
            open boundary, and ``...`` represents any non-xy dimensions of the
            input array. If ``n  = stop - start`` rows/columns are extracted,
            the shape will be ``(..., np, n)``
        """
        l_inward = {Side.WEST: 1, Side.EAST: -1, Side.SOUTH: 1, Side.NORTH: -1}[
            self.side
        ]
        mslice = slice(self.mstart, self.mstop, self.mstep)
        ldim = -1 if self.side in (Side.WEST, Side.EAST) else -2
        lstart = self.l + l_inward * start
        assert lstart >= 0 and lstart < values.shape[ldim]
        if stop is None:
            # Extract a single row/column
            lslice = lstart
        else:
            # Extract a range of rows/columns
            llast = self.l + l_inward * (stop - 1)
            lstop = None if llast == 0 else llast + l_inward
            lslice = slice(lstart, lstop, l_inward)
        if self.side in (Side.WEST, Side.EAST):
            return values[Ellipsis, mslice, lslice]
        else:
            if stop is not None:
                return np.swapaxes(values[Ellipsis, lslice, mslice], -1, -2)
            return values[Ellipsis, lslice, mslice]

    def extract_uv_in(self, u: np.ndarray, v: np.ndarray) -> np.ndarray:
        uv = u if self.side in (Side.WEST, Side.EAST) else v
        return uv[self.slice_uv_in]


class BoundaryCondition:
    def initialize(self, grid: core.Grid):
        pass

    def get_updater(
        self, boundary: OpenBoundary, array: core.Array, bdy: np.ndarray
    ) -> Callable[[], None]:
        """Return a function that updates the model values at the open boundary.

        Args:
            boundary: the open boundary
            array: array with all model values
            bdy: prescribed values at the open boundary
        """
        raise NotImplementedError

    def prepare_depth_explicit(self):
        pass

    def get_clipper(
        self, array: core.Array, slc: Tuple[Union[slice, int], ...]
    ) -> Callable[[np.ndarray, np.ndarray], None]:
        valid_min = array.attrs.get("_minimum")
        if isinstance(valid_min, np.ndarray):
            valid_min = valid_min[slc]
        valid_max = array.attrs.get("_maximum")
        if isinstance(valid_max, np.ndarray):
            valid_max = valid_max[slc]
        if valid_min is not None and valid_max is not None:
            return functools.partial(np.clip, a_min=valid_min, a_max=valid_max)
        elif valid_min is not None:
            return functools.partial(np.maximum, valid_min)
        elif valid_max is not None:
            return functools.partial(np.minimum, valid_min)

        def no_op(x: np.ndarray, out: np.ndarray):
            out[...] = x

        return no_op


class Sponge(BoundaryCondition):
    def __init__(
        self,
        n: int = 3,
        tmrlx: bool = False,
        tmrlx_max: float = 0.25,
        tmrlx_min: float = 0.0,
        tmrlx_ucut: float = 0.02,
        tmrlx_umin: Optional[float] = None,
    ):
        """Create a sponge boundary condition, with ``n`` points inward from
        the open boundary being relaxed to values at the boundary.

        If a flow-dependent relaxation coefficient is used (``tmrlx=True``),
        values at the boundary itself will blend the prescribed values
        and a weighted mean over the sponge zone. The fraction based on
        prescribed values (the relaxation coefficient) will increase with the
        inward flow velocity.

        Args:
            n: number of points in the sponge zone; the boundary itself is not
                 included, so ``n=0`` is equivalent to a clamped boundary.
            tmrlx: use a flow-dependent relaxation coefficient
            tmrlx_max: maximum relaxation coefficient (fraction per timestep).
                Only used if ``tmrlx`` is ``True``.
            tmrlx_min: minimum relaxation coefficient (fraction per timestep).
                Only used if ``tmrlx`` is ``True``.
            tmrlx_ucut: inward flow velocity (m s-1) at which relaxation
                coefficient reaches its maximum value. Only used if ``tmrlx``
                is ``True``.
            tmrlx_umin: inward flow velocity (m s-1) at which relaxation
                coefficient reaches its minimum value. If ``None``,
                ``tmrlx_umin`` is set to ``-0.25 * tmrlx_ucut``.
                Only used if ``tmrlx`` is ``True``.
        """
        self.n = n
        self.tmrlx = tmrlx
        self.tmrlx_max = tmrlx_max
        self.tmrlx_min = tmrlx_min
        self.tmrlx_ucut = tmrlx_ucut
        if tmrlx_umin is None:
            tmrlx_umin = -0.25 * tmrlx_ucut
        self.tmrlx_umin = tmrlx_umin
        self.rlxcoef: Optional[core.Array] = None
        self.inflow: Optional[np.ndarray] = None

        # relaxation coefficient as per Martinsen & Engedahl (1987), Eq 3
        # https://doi.org/10.1016/0378-3839(87)90028-7
        # These represent the fraction of the cell's value that is replaced
        # by the value at the boundary, per timestep, as we move inward from
        # the boundary.
        self.sponge_coef = ((n - np.arange(n)) / (n + 1.0)) ** 2

    def initialize(self, grid: core.Grid):
        """Create the relaxation coefficient array if using flow-dependent relaxation.
        This will be called for every open boundary of every array that uses this
        boundary condition, so we take care to only create the array once.
        """
        if self.tmrlx and self.rlxcoef is None:
            self.rlxcoef = grid.array(z=CENTERS, on_boundary=True)
            grid.open_boundaries.velocity_3d_in.saved = True
            self.inflow = grid.open_boundaries.velocity_3d_in.all_values

    def prepare_depth_explicit(self):
        """Update the relaxation coefficient based on the current flow velocity.
        This is done once per timestep and processes all open boundaries at once.
        """
        if self.tmrlx:
            self.rlxcoef.all_values[...] = (self.tmrlx_max - self.tmrlx_min) * np.clip(
                (self.inflow - self.tmrlx_umin) / (self.tmrlx_ucut - self.tmrlx_umin),
                0.0,
                1.0,
            ) + self.tmrlx_min

    def get_updater(
        self, boundary: OpenBoundary, array: core.Array, bdy: np.ndarray
    ) -> Callable[[], None]:
        mask = boundary.extract_inward(
            array.grid.mask.all_values, start=1, stop=self.n + 1
        )
        n = mask.shape[-1]
        assert self.n == n

        sp = np.empty((boundary.np, n), dtype=float)
        sp[...] = self.sponge_coef

        if self.tmrlx:
            # The flow-dependent relaxation rlxcoef varies in time
            # and will be updated by prepare_depth_explicit.
            # It covers all boundaries; here we obtain its slice
            # corresponding to the current boundary.
            sp[mask == 0] = 0.0  # only include water points
            w = sp / sp.sum(axis=1, keepdims=True)
            rlxcoef = self.rlxcoef.all_values[boundary.slice_bdy].T
        else:
            # Values at the boundary will be set to prescribed values
            w = None
            rlxcoef = None
        sp[mask != 1] = 0.0  # only include updatable water points (exclude bdy)

        return functools.partial(
            self.update_boundary,
            array.all_values[boundary.slice_t],
            boundary.extract_inward(array.all_values, start=1, stop=n + 1),
            bdy.T,
            sp,
            rlxcoef,
            w,
        )

    @staticmethod
    def update_boundary(
        values: np.ndarray,
        sponge_values: np.ndarray,
        bdy_values: np.ndarray,
        sp: np.ndarray,
        rlxcoef: Optional[np.ndarray],
        w: Optional[np.ndarray],
    ):
        """Update model values at open boundary and in sponge zone,
        for a single model variable and single open boundary.

        Args:
            values: current model values at the open boundary (nz x nbdy)
            sponge_values: current model values in the sponge zone
                (inward from open boundary). (nz x nbdy x nsponge)
            bdy_values: prescribed values at open boundary (nz x nbdy)
            sp: fraction of model values in sponge zone to be replaced by
                prescribed values (decreasing inward from boundary) (nbdy x nsponge)
            rlxcoef: fraction of model values at the open boundary to be taken
                from prescribed values. The remainder will be based on a
                weighted mean over the sponge. (nz x nbdy)
            w: weight for sponge mean used as part of new model values at the
                open boundary. (nbdy x nsponge)
        """
        if w is not None:
            # note: where=w != 0.0 is used to avoid mixing in NaNs from areas where w=0
            sponge_mean = (w * sponge_values).sum(axis=-1, where=w != 0.0)
            bdy_values = rlxcoef * bdy_values + (1.0 - rlxcoef) * sponge_mean
        values[...] = bdy_values
        sponge_values += sp * (bdy_values[..., np.newaxis] - sponge_values)


class ZeroGradient(BoundaryCondition):
    """Zero-gradient boundary condition

    The model values at the open boundary are set to the values at the
    first point inward from the boundary.
    """

    def get_updater(
        self,
        boundary: OpenBoundary,
        array: core.Array,
        bdy: np.ndarray,
    ) -> Callable[[], None]:
        return functools.partial(
            self.get_clipper(array, boundary.slice_t),
            boundary.extract_inward(array.all_values, start=1),
            out=array.all_values[boundary.slice_t],
        )


class Clamped(BoundaryCondition):
    """Clamped boundary condition

    The model values at the open boundary are set to the prescribed values.
    """

    def get_updater(
        self,
        boundary: OpenBoundary,
        array: core.Array,
        bdy: np.ndarray,
    ):
        return functools.partial(
            self.get_clipper(array, boundary.slice_t),
            bdy.T,
            out=array.all_values[boundary.slice_t],
        )


class Flather(BoundaryCondition):
    def __init__(self, transport: bool = False):
        self.transport = transport

    def get_updater(self, boundary: OpenBoundary, array: core.Array, bdy: np.ndarray):
        return functools.partial(
            self.update_transport if self.transport else self.update_velocity,
            array.all_values[boundary.slice_t],
            bdy,
            boundary.tp,
            boundary.flow_ext,
            array.grid.Dclip.all_values[boundary.slice_t],
            boundary.inflow_sign,
            self.get_clipper(array, boundary.slice_t),
        )

    @staticmethod
    def update_velocity(
        z: np.ndarray,
        z_ext: np.ndarray,
        tp: np.ndarray,
        vel_ext: np.ndarray,
        D: np.ndarray,
        inflow_sign: float,
        clipper: Callable[[np.ndarray, np.ndarray], None],
    ):
        clipper(z_ext - inflow_sign * (tp - vel_ext * D) / np.sqrt(D * GRAVITY), out=z)

    @staticmethod
    def update_transport(
        z: np.ndarray,
        z_ext: np.ndarray,
        tp: np.ndarray,
        tp_ext: np.ndarray,
        D: np.ndarray,
        inflow_sign: float,
        clipper: Callable[[np.ndarray, np.ndarray], None],
    ):
        clipper(z_ext - inflow_sign * (tp - tp_ext) / np.sqrt(D * GRAVITY), out=z)


class ArrayOpenBoundary:
    def __init__(
        self,
        open_boundaries: "OpenBoundaries",
        boundary: OpenBoundary,
        prescribed_values: np.ndarray,
        model_values: np.ndarray,
        type: BoundaryCondition,
    ):
        self._make_bc = open_boundaries._make_bc
        self._boundary = boundary
        self._prescribed_values = prescribed_values
        self._model_values = model_values
        self._type = type

    @property
    def type(self) -> BoundaryCondition:
        return self._type

    @type.setter
    def type(self, value: int):
        self._type = self._make_bc(value)

    @property
    def values(self) -> int:
        return self._prescribed_values


class ArrayOpenBoundaries:
    __slots__ = "_array", "values", "_bdy", "updaters"

    def __init__(self, array: core.Array, type=None):
        self._array = array
        self.values = array.grid.array(
            name=f"{array.name}_bdy",
            z=array.z,
            on_boundary=True,
            fill_value=array.fill_value,
            attrs={
                "_time_varying": array.attrs.get("_time_varying", TimeVarying.MICRO)
            },
        )
        if type is not None:
            type = array.grid.open_boundaries._make_bc(type)
        self._bdy: List[ArrayOpenBoundary] = []
        for bdy in array.grid.open_boundaries.active:
            self._bdy.append(
                ArrayOpenBoundary(
                    array.grid.open_boundaries,
                    bdy,
                    self.values.all_values[bdy.slice_bdy],
                    array.all_values[bdy.slice_t],
                    type or (bdy.type_2d if array.ndim == 2 else bdy.type_3d),
                )
            )
        self.updaters: List[Callable[[], None]] = []

    def _set_type(self, value: int):
        value = self._array.grid.open_boundaries._make_bc(value)
        for bdy in self._bdy:
            bdy.type = value

    type = property(fset=_set_type)

    def initialize(self) -> Iterable[BoundaryCondition]:
        bcs: Set[BoundaryCondition] = set()
        for bdy in self._bdy:
            bdy._type.initialize(self._array.grid)
            self.updaters.append(
                bdy._type.get_updater(
                    bdy._boundary, self._array, bdy._prescribed_values
                )
            )
            bcs.add(bdy._type)
        return bcs

    def update(self):
        """Update the tracer at the open boundaries"""
        for updater in self.updaters:
            updater()


class OpenBoundaries(Sequence[OpenBoundary]):
    __slots__ = (
        "nx",
        "ny",
        "logger",
        "np",
        "np_glob",
        "i",
        "j",
        "i_glob",
        "j_glob",
        "z",
        "u",
        "v",
        "u_rot",
        "v_rot",
        "lon",
        "lat",
        "zc",
        "zf",
        "local_to_global",
        "_boundaries",
        "_frozen",
        "_rotator",
        "sponge",
        "zero_gradient",
        "active",
        "mirror_U",
        "mirror_V",
        "mirror_TU",
        "mirror_TV",
        "velocity_3d_in",
        "bcs",
    )

    def __init__(self, nx: int, ny: int, logger: logging.Logger):
        self.nx = nx
        self.ny = ny
        self.logger = logger
        self._boundaries: List[OpenBoundary] = []
        self.sponge = Sponge()
        self.zero_gradient = ZeroGradient()
        self.active: List[OpenBoundary] = []
        self._frozen = False
        self.bcs: Set[BoundaryCondition] = set()

    def _make_bc(self, value) -> BoundaryCondition:
        if isinstance(value, BoundaryCondition):
            return value
        if value == ZERO_GRADIENT:
            return self.zero_gradient
        elif value == SPONGE:
            return self.sponge
        elif value == CLAMPED:
            return Clamped()
        elif value == FLATHER_ELEV:
            return Flather()
        elif value == FLATHER_TRANSPORT:
            return Flather(transport=True)
        else:
            raise Exception(f"Unknown boundary type {value} specified")

    def add_by_index(
        self,
        side: Side,
        l: int,
        mstart: int,
        mstop: int,
        type_2d: int,
        type_3d: int,
        name: Optional[str] = None,
    ):
        """Note that l, mstart, mstop are 0-based indices of a T point in the global domain.
        mstop indicates the upper limit of the boundary - it is the first index that is
        EXcluded.
        """

        assert (
            not self._frozen
        ), "The open boundary collection has already been initialized"
        along_y = side in (Side.WEST, Side.EAST)
        if along_y:
            l_max, m_max = (self.nx, self.ny)
        else:
            l_max, m_max = (self.ny, self.nx)

        mstep = 1 if mstop > mstart else -1
        assert mstart >= 0 and mstart < m_max
        assert mstop - mstep > 0 and mstop - mstep < m_max
        assert mstart != mstop
        assert l >= 0 and l < l_max

        if name is None:
            name = str(len(self._boundaries))

        mmin = min(mstart, mstop - mstep)
        mmax = max(mstart, mstop - mstep)
        for b in self._boundaries:
            current_along_y = b.side in (Side.WEST, Side.EAST)
            if along_y == current_along_y:
                if l == b.l_glob:
                    ostart = max(mmin, min(b.mstart_glob, b.mstop_glob - b.mstep))
                    ostop = min(mmax, max(b.mstart_glob, b.mstop_glob - b.mstep)) + 1
                    overlap = ostop - ostart
                    if overlap > 0:
                        raise Exception(
                            f"New boundary {name} overlaps in {overlap} points"
                            f" with existing boundary {b.name}"
                        )
            else:
                cross = (
                    l >= min(b.mstart_glob, b.mstop_glob - b.mstep)
                    and l <= max(b.mstart_glob, b.mstop_glob - b.mstep)
                    and b.l_glob >= mmin
                    and b.l_glob <= mmax
                )
                if cross:
                    i, j = (l, b.l_glob) if along_y else (b.l_glob, l)
                    raise Exception(
                        f"New boundary {name} crosses existing boundary {b.name}"
                        f" at i={i}, j={j}"
                    )

        self._boundaries.append(
            OpenBoundary(name, side, l, mstart, mstop, type_2d, type_3d)
        )

    def add_top_boundary(
        self, name: str, j: int, istart: int, istop: int, type_2d: int, type_3d: int
    ):
        """Add an open boundary with the model exterior above and the model
        interior below. This is a northern open boundary in a spherical domain.
        """
        self.add_by_index(Side.NORTH, j, istart, istop, type_2d, type_3d, name=name)

    def add_bottom_boundary(
        self, name: str, j: int, istart: int, istop: int, type_2d: int, type_3d: int
    ):
        """Add an open boundary with the model exterior below and the model
        interior above. This is a southern open boundary in a spherical domain.
        """
        self.add_by_index(Side.SOUTH, j, istart, istop, type_2d, type_3d, name=name)

    def add_left_boundary(
        self, name: str, i: int, jstart: int, jstop: int, type_2d: int, type_3d: int
    ):
        """Add an open boundary with the model exterior to the left and the
        model interior to the right. This is a western open boundary in a
        spherical domain.
        """
        self.add_by_index(Side.WEST, i, jstart, jstop, type_2d, type_3d, name=name)

    def add_right_boundary(
        self, name: str, i: int, jstart: int, jstop: int, type_2d: int, type_3d: int
    ):
        """Add an open boundary with the model exterior to the right and the
        model interior to the left. This is an eastern open boundary in a
        spherical domain.
        """
        self.add_by_index(Side.EAST, i, jstart, jstop, type_2d, type_3d, name=name)

    def clear(self):
        """Delete all open boundaries."""
        self._boundaries.clear()

    def adjust_mask(self, mask: np.ndarray):
        """Adjust the global mask for open boundaries

        This assigns mask values 2 (masked T point), 3 (velocity point in
        between open boundary T points) and 4 (velocity point just
        outside the open boundary)

        Args:
            mask: writeable mask defined on the supergrid, with shape (1+2\*ny, 1+2\*nx)
        """
        assert mask.shape == (1 + 2 * self.ny, 1 + 2 * self.nx)
        umask = mask[1::2, 0::2]
        vmask = mask[0::2, 1::2]
        tmask = mask[1::2, 1::2]
        for boundary in self._boundaries:
            l = boundary.l_glob
            mslice = slice(boundary.mstart_glob, boundary.mstop_glob, boundary.mstep)
            np = boundary.mstop_glob - boundary.mstart_glob
            if boundary.side in (Side.WEST, Side.EAST):
                bdy_mask = tmask[mslice, l]
            else:
                bdy_mask = tmask[l, mslice]

            bdy_on_land = bdy_mask == 0
            if bdy_on_land.any():
                self.logger.error(
                    f"Open boundary {boundary.name}: {bdy_on_land.sum()} of"
                    f" {np} points of this {boundary.side.name.capitalize()}ern"
                    " boundary are on land"
                )
                raise Exception()

            bdy_mask[:] = 2

            # Velocity points that lie just outside the T points of the open boundary
            if boundary.side in (Side.WEST, Side.EAST):
                umask[mslice, l if boundary.side == Side.WEST else l + 1] = 4
            else:
                vmask[l if boundary.side == Side.SOUTH else l + 1, mslice] = 4

        # Velocity points that lie in between open boundary T points
        umask[:, 1:-1][(tmask[:, :-1] == 2) & (tmask[:, 1:] == 2)] = 3
        vmask[1:-1, :][(tmask[:-1, :] == 2) & (tmask[1:, :] == 2)] = 3

    def initialize(self, grid: core.Grid):
        """Freeze the open boundary collection. Drop those outside the current
        subdomain.
        """
        assert (
            not self._frozen
        ), "The open boundary collection has already been initialized"

        for boundary in self._boundaries:
            boundary.to_local_grid(grid)

        nbdyp = 0
        nbdyp_glob = 0

        # Indices of T point in open boundary into arrays that INCLUDE halos.
        # We start with an empty array to support later concatenation in
        # subdomains without any open boundary
        bdy_i = [np.empty((0,), dtype=np.intc)]
        bdy_j = [np.empty((0,), dtype=np.intc)]
        side2count = {Side.WEST: 0, Side.NORTH: 0, Side.EAST: 0, Side.SOUTH: 0}
        self.local_to_global: List[slice] = []
        bdy_types_2d = {}
        for boundary in self._boundaries:
            if boundary.l is not None:
                mskip = (boundary.mstart - boundary.mstart_) // boundary.mstep
                assert mskip >= 0
                boundary.start = nbdyp
                boundary.stop = nbdyp + boundary.np
                boundary.slice_bdy = (
                    slice(boundary.start, boundary.stop),
                    Ellipsis,
                )
                if boundary.type_2d not in bdy_types_2d:
                    bdy_types_2d[boundary.type_2d] = self._make_bc(boundary.type_2d)
                boundary.type_2d = bdy_types_2d[boundary.type_2d]
                nbdyp += boundary.np

                bdy_i.append(boundary.i)
                bdy_j.append(boundary.j)

                start_glob = nbdyp_glob + mskip
                stop_glob = start_glob + boundary.np
                if self.local_to_global and self.local_to_global[-1].stop == start_glob:
                    # Attach to previous boundary by dropping previous slice
                    # and adopting its start index
                    start_glob = self.local_to_global.pop().start
                self.local_to_global.append(slice(start_glob, stop_glob))
                side2count[boundary.side] += 1
                self.active.append(boundary)
            nbdyp_glob += (boundary.mstop_ - boundary.mstart_) // boundary.mstep

        # Number of open boundary points (local and global)
        self.np = nbdyp
        self.np_glob = nbdyp_glob
        grid.open_boundaries = self

        # Local indices of open boundary points within local subdomain
        # (T grid, for arrays that INclude halos)
        self.i = np.concatenate(bdy_i, dtype=np.intc)
        self.j = np.concatenate(bdy_j, dtype=np.intc)
        assert (grid.mask.all_values[self.j, self.i] == 2).all()

        # Global indices of open boundary points within local subdomain
        # (T grid, for arrays that EXclude halos)
        self.i_glob = self.i - grid.halox + grid.tiling.xoffset
        self.j_glob = self.j - grid.haloy + grid.tiling.yoffset

        self.logger.info(
            f"{sum(side2count.values())} open boundaries"
            f" ({side2count[Side.WEST]} West, {side2count[Side.NORTH]} North,"
            f" {side2count[Side.EAST]} East, {side2count[Side.SOUTH]} South)"
        )

        if self.np > 0:
            if self.np == self.np_glob:
                # This subdomain includes all open boundaries in full
                assert (
                    len(self.local_to_global) == 1
                    and self.local_to_global[0].start == 0
                    and self.local_to_global[0].stop == self.np_glob
                )
                self.local_to_global = None
            else:
                # This subdomain includes part of the open boundaries
                slices = ", ".join(
                    [f"[{s.start}:{s.stop}]" for s in self.local_to_global]
                )
                self.logger.info(f"global-to-local open boundary map: {slices}")

        if grid.ugrid:
            self._configure_mirroring(grid)

        # Horizontal coordinates of open boundary points
        if grid.lon is not None:
            self.lon = grid.array(
                on_boundary=True, fill=grid.lon.all_values[self.j, self.i]
            )
        if grid.lat is not None:
            self.lat = grid.array(
                on_boundary=True, fill=grid.lat.all_values[self.j, self.i]
            )

        if grid.nz:
            # Vertical coordinates of open boundary points.
            # These will be updated by indexing into the full zc and zf
            # after every update of surface elevation/water depth/cell thickness
            kwargs = dict(fill_value=FILL_VALUE, units="m", on_boundary=True)
            self.zc = grid.array(
                name="zc_bdy",
                z=CENTERS,
                long_name="height",
                attrs=grid.zc.attrs,
                **kwargs,
            )
            self.zf = grid.array(
                name="zf_bdy",
                z=INTERFACES,
                long_name="interface height",
                attrs=grid.zf.attrs,
                **kwargs,
            )

        # Prescribed depth-averaged or depth-integrated velocity at the open boundaries
        self.u = grid.array(name="u_bdy", on_boundary=True)
        self.v = grid.array(name="v_bdy", on_boundary=True)

        if grid.rotation is not None and self.np > 0:
            self._rotator = core.Rotator(grid.rotation.all_values[self.j, self.i])
            self.u_rot = grid.array(name="u_rot_bdy", on_boundary=True)
            self.v_rot = grid.array(name="v_rot_bdy", on_boundary=True)
        else:
            self._rotator = None
            self.u_rot = self.u
            self.v_rot = self.v

        self.velocity_3d_in = grid.array(z=CENTERS, on_boundary=True)

        self._frozen = True

    def _configure_mirroring(self, grid: core.Grid):
        """Set up indices for mirroring across boundaries on U and V grids"""
        mirror_U = []
        mirror_V = []
        mirror_TU = []
        mirror_TV = []
        umask = grid.ugrid.mask.all_values
        vmask = grid.vgrid.mask.all_values
        for boundary in self._boundaries:
            if boundary.l is None:
                continue

            # Identify velocity points that lie within the open boundary in between
            # tracer points with mask=2. Their indices will be (i_velout, j_velout)
            # on the corresponding velocity grid (U or V). Values at these points
            # will be mirrored from the interior velocity point (i_velin, j_velin)
            # We look 1 point beyond the start and stop of the boundary within the
            # current subdomain to catch cases where (a) the boundary extends into the
            # next subdomain and (b) the current boundary neighbors another boundary
            # We then only include points with mask value 3 (see "select" below).
            if boundary.side in (Side.WEST, Side.EAST):
                j_velout = np.arange(max(boundary.j.min() - 1, 0), boundary.j.max() + 1)
                i_velout = np.full_like(j_velout, boundary.i[0])
                mirror_mask = vmask[j_velout, i_velout]
            else:
                i_velout = np.arange(max(boundary.i.min() - 1, 0), boundary.i.max() + 1)
                j_velout = np.full_like(i_velout, boundary.j[0])
                mirror_mask = umask[j_velout, i_velout]
            select = mirror_mask == 3
            i_velout = i_velout[select]
            j_velout = j_velout[select]
            i_velin = i_velout + {Side.EAST: -1, Side.WEST: 1}.get(boundary.side, 0)
            j_velin = j_velout + {Side.NORTH: -1, Side.SOUTH: 1}.get(boundary.side, 0)

            # Identify velocity points along an open boundary
            # The indices of inner points will be (i_in, j_in);
            # those indices of outer points will be (i_out, j_out).
            # Values at outer points will be mirrored from either the neighboring
            # inner T point (boundary.i, boundary.j) [e.g., elevations at mask=2]
            # or the neighboring inner velocity point.
            i_in = boundary.i + {Side.EAST: -1}.get(boundary.side, 0)
            j_in = boundary.j + {Side.NORTH: -1}.get(boundary.side, 0)
            i_out = boundary.i + {Side.WEST: -1}.get(boundary.side, 0)
            j_out = boundary.j + {Side.SOUTH: -1}.get(boundary.side, 0)
            select = (i_in >= 0) & (j_in >= 0) & (i_out >= 0) & (j_out >= 0)
            i_in, i_out, j_in, j_out = (a[select] for a in (i_in, i_out, j_in, j_out))

            if boundary.side in (Side.WEST, Side.EAST):
                assert (grid.ugrid.mask.all_values[j_out, i_out] == 4).all()
                mirror_U.append([i_in, j_in, i_out, j_out])
                mirror_V.append([i_velin, j_velin, i_velout, j_velout])
                mirror_TU.append([boundary.i, boundary.j, i_out, j_out])
            else:
                assert (grid.vgrid.mask.all_values[j_out, i_out] == 4).all()
                mirror_V.append([i_in, j_in, i_out, j_out])
                mirror_U.append([i_velin, j_velin, i_velout, j_velout])
                mirror_TV.append([boundary.i, boundary.j, i_out, j_out])

        def pack_mirror(indices: Sequence[Sequence[np.ndarray]]):
            if indices:
                i_in = np.concatenate([i[0] for i in indices], dtype=np.intp)
                j_in = np.concatenate([i[1] for i in indices], dtype=np.intp)
                i_out = np.concatenate([i[2] for i in indices], dtype=np.intp)
                j_out = np.concatenate([i[3] for i in indices], dtype=np.intp)
                return ((Ellipsis, j_in, i_in), (Ellipsis, j_out, i_out))

        grid.ugrid._mirrors[grid.ugrid] = pack_mirror(mirror_U)
        grid.vgrid._mirrors[grid.vgrid] = pack_mirror(mirror_V)
        grid._mirrors[grid.ugrid] = pack_mirror(mirror_TU)
        grid._mirrors[grid.vgrid] = pack_mirror(mirror_TV)

    def start(
        self,
        U: core.Array,
        V: core.Array,
        uk: Optional[core.Array],
        vk: Optional[core.Array],
        fields: Mapping[str, core.Array],
    ):
        for boundary in self.active:
            boundary.tp = boundary.extract_uv_in(U.all_values, V.all_values)
            if uk is not None:
                boundary.vel = boundary.extract_uv_in(uk.all_values, vk.all_values)
            if boundary.side in (Side.EAST, Side.WEST):
                boundary.flow_ext = self.u_rot[boundary.slice_bdy]
            else:
                boundary.flow_ext = self.v_rot[boundary.slice_bdy]
            boundary.velocity_3d_in = self.velocity_3d_in.all_values[boundary.slice_bdy]
        for field in fields.values():
            if hasattr(field, "open_boundaries"):
                self.bcs.update(field.open_boundaries.initialize())

    def __getitem__(self, key: int) -> OpenBoundary:
        return self._boundaries[key]

    def __len__(self) -> int:
        return len(self._boundaries)

    def __iter__(self):
        return iter(self._boundaries)

    @property
    def local_to_global_indices(self) -> np.ndarray:
        indices = np.arange(self.np, dtype=int)
        if self.local_to_global is not None:
            i = 0
            for s in self.local_to_global:
                indices[i : i + s.stop - s.start] = np.arange(s.start, s.stop)
                i += s.stop - s.start
            assert i == self.np
        return indices

    def prepare(self, macro_active: bool):
        """Prepare open boundary forcing. This includes rotating prescribed velocity
        fields, as well as extracting 3D velocities from which each class of
        open boundary condition can calculate derived metrics."""
        if self._rotator:
            u_rot, v_rot = self._rotator(self.u.all_values, self.v.all_values)
            self.u_rot.all_values[:] = u_rot
            self.v_rot.all_values[:] = v_rot
        if macro_active:
            if self.velocity_3d_in.saved:
                for boundary in self.active:
                    boundary.velocity_3d_in[:] = boundary.inflow_sign * boundary.vel.T
            for bc in self.bcs:
                bc.prepare_depth_explicit()
