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
    NamedTuple,
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


class Side(enum.IntEnum):
    WEST = 1
    NORTH = 2
    EAST = 3
    SOUTH = 4


class OpenBoundary:
    """A single open boundary.
    This defines the location of the open boundary, but not the type of
    boundary condition to apply to the different variables. For that, use
    :attr:`core.Array.open_boundaries`.
    """

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
        """Create an open boundary.
        The supplied indices are global indices in the T grid.

        Args:
            name: unique name for this boundary
            side: side of the domain
            l: fixed index (i for left/right boundaries, j for top/bottom)
            mstart: start index (j for left/right boundaries, i for top/bottom)
            mstop: stop index (exclusive)
            type_2d: default type of boundary condition for 2D variables
            type_3d: default type of boundary condition for 3D variables
        """
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
        """Convert global indices to local indices of the current subdomain."""
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
        mlast = self.mstop_ - self.mstep
        if (
            self.l < 0
            or self.l >= l_max
            or max(self.mstart_, mlast) < 0
            or min(self.mstart_, mlast) >= m_max
        ):
            # Boundary lies completely outside current subdomain. Record it anyway,
            # so we can later set up a global -> local map of open boundary points
            self.l, self.mstart, self.mstop = None, None, None
            return
        self.mstart = min(max(0, self.mstart_), m_max - 1)
        self.mstop = min(max(0, mlast), m_max - 1) + self.mstep

        mslice = slice(self.mstart, self.mstop, self.mstep)
        ms = np.arange(self.mstart, self.mstop, self.mstep)
        self.np = ms.size
        l = self.l
        ls = np.broadcast_to(l, ms.shape)
        if side in (Side.WEST, Side.EAST):
            self.i, self.j = ls, ms
            self.slice_t = (Ellipsis, mslice, l)
            self.slice_uv_in = (Ellipsis, mslice, l if side == Side.WEST else l - 1)
        else:
            self.i, self.j = ms, ls
            self.slice_t = (Ellipsis, l, mslice)
            self.slice_uv_in = (Ellipsis, l if side == Side.SOUTH else l - 1, mslice)
        assert (self.i >= 0).all() and (self.i < grid.nx_).all()
        assert (self.j >= 0).all() and (self.j < grid.ny_).all()

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
        """Extract velocity or transport across the boundary from components at
        velocity points (U or V) just inward from the boundary.

        Args:
            u: 2D or 3D component of velocity or transport in x-direction.
                It must be defined on the U grid.
            v: 2D or 3D component of velocity or transport in y-direction.
                It must be defined on the V grid.

        Returns:
            2D or 3D array of velocity or transport across the boundary.
        """
        uv = u if self.side in (Side.WEST, Side.EAST) else v
        return uv[self.slice_uv_in]


class BoundaryCondition:
    """Base class for an open boundary condition that can be applied to any
    open boundary of any array. Specific types of open boundary condition
    (clamped, zero-gradient, etc.) must subclass this base class and
    implement :meth:`get_updater`
    """

    def initialize(self, grid: core.Grid):
        """Initialize the boundary condition.

        Subclasses for specific types of open boundary condition can perform
        their own initialization here, e.g., by creating auxiliary arrays.
        This method is called once for each open boundary of each array that
        uses this boundary condition.

        Args:
            grid: the T grid on which the boundaries are defined
        """
        pass

    def get_updater(
        self, boundary: OpenBoundary, array: core.Array, bdy: np.ndarray
    ) -> Callable[[], None]:
        """Return a function that updates the model values at the open
        boundary.

        Boundary conditions that use nearby interior values (e.g., zero
        gradient, flow-dependent relaxation) can only use points with mask=1
        (wet and active). They must skip mask=2 to avoid dependence on the
        order in which boundaries are updated.

        Args:
            boundary: the open boundary
            array: array with all model values
            bdy: prescribed values at the open boundary.
                If depth-explicit, its shape is (np, nz), otherwise (np,)
        """
        raise NotImplementedError

    def prepare_depth_explicit(self):
        pass

    @staticmethod
    def get_setter(
        array: core.Array,
        slc: Tuple[Union[slice, int], ...],
        where: Union[bool, np.ndarray] = True,
    ) -> Callable[[np.ndarray], None]:
        """Takes an array object and a slice, and returns a function that sets
        values in the array at the slice, respecting the array's mask and valid
        range.

        Args:
            array: array that will be updated
            slc: slice of the values in the array that will be updated
            where: boolean array indicating which points to update
        """
        out = array.all_values[slc]
        kwargs = dict()
        where = (array.grid.mask.all_values[slc] != 0) & where
        if not where.all():
            kwargs["where"] = np.broadcast_to(where, out.shape)
        valid_min = array.attrs.get("_minimum")
        if isinstance(valid_min, np.ndarray):
            valid_min = valid_min[slc]
        valid_max = array.attrs.get("_maximum")
        if isinstance(valid_max, np.ndarray):
            valid_max = valid_max[slc]
        if valid_min is not None and valid_max is not None:
            return functools.partial(np.clip, valid_min, valid_max, out, **kwargs)
        elif valid_min is not None:
            return functools.partial(np.maximum, valid_min, out=out, **kwargs)
        elif valid_max is not None:
            return functools.partial(np.minimum, valid_min, out=out, **kwargs)
        elif "where" in kwargs:
            return functools.partial(np.copyto, out, **kwargs)
        else:
            return functools.partial(out.__setitem__, (Ellipsis,))


class Clamped(BoundaryCondition):
    """Clamped boundary condition

    The model values at the open boundary are set to the prescribed values.
    """

    def get_updater(
        self, boundary: OpenBoundary, array: core.Array, bdy: np.ndarray
    ) -> Callable[[], None]:
        return functools.partial(self.get_setter(array, boundary.slice_t), bdy.T)


class ZeroGradient(BoundaryCondition):
    """Zero-gradient boundary condition

    The model values at the open boundary are set to the values at the
    first point inward from the boundary.
    """

    def get_updater(
        self, boundary: OpenBoundary, array: core.Array, bdy: np.ndarray
    ) -> Callable[[], None]:
        where = boundary.extract_inward(array.grid.mask.all_values, start=1) == 1
        return functools.partial(
            self.get_setter(array, boundary.slice_t, where),
            boundary.extract_inward(array.all_values, start=1),
        )


class Sponge(Clamped):
    """Sponge boundary condition, with ``n`` points inward from
    the open boundary being relaxed to values at the boundary.

    If flow-dependent relaxation is used (``tmrlx=True``),
    values at the boundary itself will blend the prescribed values
    and a weighted mean over the sponge zone. The fraction taken from the
    prescribed values will increase with increasing inward flow velocity.

    If flow-dependent relaxation is not used (``tmrlx=False``), values at
    the boundary will be set to the prescribed values.

    Relaxation rates for the ``n`` points are available as :attr:`sp`.
    They default to the formulation of `Martinsen & Engedahl (1987)
    <https://doi.org/10.1016/0378-3839(87)90028-7>`_ (Eq 3)."""

    def __init__(
        self,
        n: int = 3,
        tmrlx: bool = False,
        tmrlx_max: float = 0.25,
        tmrlx_min: float = 0.0,
        tmrlx_ucut: float = 0.02,
        tmrlx_umin: Optional[float] = None,
    ):
        """
        Args:
            n: number of points in the sponge zone; the boundary itself is not
                 included, so ``n=0`` is equivalent to a clamped boundary.
            tmrlx: make use of prescribed values at the open boundary
                dependent on flow direction. This is typically used to
                gradually disable use of prescribed values as the current
                switches from inflow to outflow over the boundary.
            tmrlx_max: maximum fraction to take from prescribed values at the
                open boundary. The remaining fraction will be taken from a
                weighted mean over the sponge zone. This maximum is reached
                at strong inflow (>= ``tmrlx_ucut``).
                Only used if ``tmrlx`` is ``True``.
            tmrlx_min: minimum fraction to take from prescribed values at the
                open boundary. The remaining fraction will be taken from a
                weighted mean over the sponge zone. This minimum is reached
                when water flows outward over the boundary (inflow <=
                ``tmrlx_umin``). Only used if ``tmrlx`` is ``True``.
            tmrlx_ucut: inward flow velocity (m s-1) at which the use of
                prescribed boundary values becomes highest. At this inflow, or
                higher, boundary values are for fraction ``tmrlx_max`` based on
                prescribed values; the remaining ``1 - tmrlx_max`` is taken
                from a weighted mean over the sponge zone.
                Only used if ``tmrlx`` is ``True``.
            tmrlx_umin: inward flow velocity (m s-1) at which the use of
                prescribed boundary values becomes minimal. Use negative values
                to indicate outflow. At the specified inflow value, or lower,
                boundary values are for fraction ``tmrlx_min`` based on
                prescribed values; the remaining ``1 - tmrlx_min`` is taken
                from a weighted mean over the sponge zone. If ``None``,
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
        self.sp = ((n - np.arange(n)) / (n + 1.0)) ** 2

    def initialize(self, grid: core.Grid):
        """Create the relaxation coefficient array if using flow-dependent relaxation.
        This will be called for every open boundary of every array that uses this
        boundary condition, so we take care to only create the array once.
        """
        if self.tmrlx and self.rlxcoef is None:
            self.rlxcoef = grid.array(z=CENTERS, on_boundary=True)
            grid.open_boundaries.uv_in.saved = True
            self.inflow = grid.open_boundaries.uv_in.all_values

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
        slicer = functools.partial(boundary.extract_inward, start=1, stop=self.n + 1)

        # Select sponge points that are wet and active (mask=1) and connected
        # to the open boundary (mask=2) via similar such points.
        sponge_active = slicer(array.grid.mask.all_values) == 1
        assert sponge_active.shape[-1] == self.n
        sponge_active[:, 0] &= array.grid.mask.all_values[boundary.slice_t] == 2
        np.logical_and.accumulate(sponge_active, axis=-1, out=sponge_active)

        # Relax sponge zone to boundary values using coefficients sp
        sponge_target = array.all_values[boundary.slice_t][..., np.newaxis]
        collection = array.open_boundaries
        collection.add_relaxation(slicer, self.sp, sponge_target, sponge_active)

        if self.tmrlx:
            # Boundary values will blend prescribed values and a weighted mean
            # over the sponge zone. The fraction taken from the prescribed
            # values will increase with increasing inward flow velocity.
            w = np.where(sponge_active, self.sp, 0.0)
            wsum = w.sum(axis=1, keepdims=True)
            np.divide(w, wsum, out=w, where=wsum != 0.0)

            # The flow-dependent relaxation rlxcoef varies in time and will be
            # updated by prepare_depth_explicit. It covers all boundaries; here
            # we obtain its slice corresponding to the current boundary.
            rlxcoef = self.rlxcoef.all_values[boundary.slice_bdy].T

            return functools.partial(
                self.update_boundary,
                slicer(array.all_values),
                bdy.T,
                rlxcoef,
                w,
                w != 0.0,
                (w == 0.0).all(axis=1),
                self.get_setter(array, boundary.slice_t),
            )
        else:
            # Clamp boundary to prescribed values
            return super().get_updater(boundary, array, bdy)

    @staticmethod
    def update_boundary(
        sponge_values: np.ndarray,
        bdy_values: np.ndarray,
        rlxcoef: np.ndarray,
        w: np.ndarray,
        nonzerow: np.ndarray,
        clamp: np.ndarray,
        bdy_setter: Callable[[np.ndarray], None],
    ):
        """Update model values at open boundary.
        The sponge zone update is done from :meth:`ArrayOpenBoundary.update`,
        via the relaxation term added in :meth:`get_updater`

        Args:
            sponge_values: current model values in the sponge zone
                (inward from open boundary). (nz x np x nsponge)
            bdy_values: prescribed values at open boundary (nz x np)
            rlxcoef: fraction of model values at the open boundary to be taken
                from prescribed values. The remainder will be based on a
                weighted mean over the sponge. (nz x np)
            w: weight for sponge mean used as part of new model values at the
                open boundary. (np x nsponge)
            nonzerow: boolean array indicating where w is non-zero.
                Cells with zero weight will be excluded from the sponge mean
                in order to avoid blending in NaNs.
            clamp: boolean array indicating where no sponge mean is available.
                In these cells, the model values will be clamped to the
                prescribed values.
            bdy_setter: function that updates values at the boundary,
                respecting their valid range and mask
        """
        sponge_mean = (w * sponge_values).sum(axis=-1, where=nonzerow)
        blend = rlxcoef * bdy_values + (1.0 - rlxcoef) * sponge_mean
        bdy_setter(np.where(clamp, bdy_values, blend))


class Flather(BoundaryCondition):
    """Flather boundary condition for elevation.

    This infers the elevation at the open boundary from prescribed
    elevation and from the difference between prescribed and modelled
    velocity across the boundary.

    See also:

    Flather, R.A. (1976) A tidal model of the north-west European
    continental shelf. Mem. Soc. R. Sci. Liege 6 (10), 141-164.

    Blayo, E., & Debreu, L. (2005). Revisiting open boundary conditions
    from the point of view of characteristic variables. Ocean Modelling,
    9(3), 231–252. `10.1016/j.ocemod.2004.07.001 <https://doi.org/10.1016/j.ocemod.2004.07.001>`_
    """

    def __init__(self, transport: bool = False):
        """
        Args:
            transport: use prescribed transport at the boundary instead of
                prescribed velocity
        """
        self.transport = transport

    def get_updater(
        self, boundary: OpenBoundary, array: core.Array, bdy: np.ndarray
    ) -> Callable[[], None]:
        return functools.partial(
            self.update_transport if self.transport else self.update_velocity,
            bdy,
            boundary.UV,
            boundary.flow_ext,
            array.grid.Dclip.all_values[boundary.slice_t],
            boundary.inflow_sign,
            self.get_setter(array, boundary.slice_t),
        )

    @staticmethod
    def update_velocity(
        z_ext: np.ndarray,
        tp: np.ndarray,
        vel_ext: np.ndarray,
        D: np.ndarray,
        inflow_sign: float,
        setter: Callable[[np.ndarray], None],
    ):
        setter(z_ext - inflow_sign * (tp - vel_ext * D) / np.sqrt(D * GRAVITY))

    @staticmethod
    def update_transport(
        z_ext: np.ndarray,
        tp: np.ndarray,
        tp_ext: np.ndarray,
        D: np.ndarray,
        inflow_sign: float,
        setter: Callable[[np.ndarray, np.ndarray], None],
    ):
        setter(z_ext - inflow_sign * (tp - tp_ext) / np.sqrt(D * GRAVITY))


class ArrayOpenBoundary:
    """Single open boundary for a single array.
    The type of open boundary condition is accessible as :attr:`type`.
    The values prescribed at the boundary are accessible as :attr:`values`.
    """

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
    def type(self, value: Union[int, BoundaryCondition]):
        self._type = self._make_bc(value)

    @property
    def values(self) -> np.ndarray:
        return self._prescribed_values


class Relaxation(NamedTuple):
    slicer: Callable[[np.ndarray], np.ndarray]
    weights: np.ndarray
    target: np.ndarray
    local: np.ndarray


class ArrayOpenBoundaries:
    """All open boundaries for a single array.
    The values prescribed at all boundaries combined are accessible
    as :attr:`values`. The type of open boundary condition can be
    set for all boundaries at once via the :attr:`type` property.
    """

    __slots__ = "_array", "values", "_bdy", "updaters", "relaxation"

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
        self.relaxation: List[Relaxation] = []

    def _set_type(self, value: Union[int, BoundaryCondition]):
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

        if self.relaxation:
            # Enumerate relaxation terms beyond the open boundary itself
            # (e.g., sponge zones) and normalize the weights so that the
            # combined relaxation from all open boundaries can at most exactly
            # replace the model values at the boundary (max total weight = 1)
            weight_sum = np.zeros(self.values.grid.mask.all_values.shape)
            for relax in self.relaxation:
                local_weight_sum = relax.slicer(weight_sum)
                local_weight_sum += relax.weights
            for relax in self.relaxation:
                relax.weights[...] /= np.maximum(relax.slicer(weight_sum), 1.0)

        return bcs

    def add_relaxation(
        self,
        slicer: Callable[[np.ndarray], np.ndarray],
        weights: np.ndarray,
        target: np.ndarray,
        where: np.ndarray,
    ):
        """Add a boundary-influenced relaxation term, e.g., a sponge zone.

        Args:
            slicer: function that extracts the relaxation area from the model
                values. It must take a 2D or 3D array of model values and
                return the slice corrresponding to the relaxation area.
            weights: fraction of the model values replaced by relaxation.
            target: target values to relax towards.
            where: boolean array indicating where the relaxation term should
                be applied.
        """
        local = slicer(self._array.all_values)
        where = (slicer(self._array.grid.mask.all_values) == 1) & where
        weights = np.where(where, weights, 0.0)
        self.relaxation.append(Relaxation(slicer, weights, target, local))

    def update(self):
        """Update the tracer at the open boundaries"""
        for updater in self.updaters:
            updater()
        olds = [relax.local.copy() for relax in self.relaxation]
        for relax, old in zip(self.relaxation, olds):
            relax.local[...] += relax.weights * (relax.target - old)


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
        "uv_in",
        "bcs",
        "allow_on_land",
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
        self.allow_on_land = False

    def _make_bc(self, value: Union[int, BoundaryCondition]) -> BoundaryCondition:
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
        bdy = OpenBoundary(name, side, l, mstart, mstop, type_2d, type_3d)
        self._boundaries.append(bdy)

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
            mask: writeable mask defined on the supergrid, with shape (1+2*ny, 1+2*nx)
        """
        assert mask.shape == (1 + 2 * self.ny, 1 + 2 * self.nx)
        umask = mask[1::2, 0::2]
        vmask = mask[0::2, 1::2]
        tmask = mask[1::2, 1::2]
        for boundary in self._boundaries:
            l = boundary.l_glob
            mslice = slice(boundary.mstart_glob, boundary.mstop_glob, boundary.mstep)
            if boundary.side in (Side.WEST, Side.EAST):
                bdy_tmask = tmask[mslice, l]
                bdy_uvmask = umask[mslice, l if boundary.side == Side.WEST else l + 1]
            else:
                bdy_tmask = tmask[l, mslice]
                bdy_uvmask = vmask[l if boundary.side == Side.SOUTH else l + 1, mslice]

            wet = bdy_tmask != 0
            if not wet.all():
                self.logger.log(
                    logging.WARNING if self.allow_on_land else logging.ERROR,
                    f"Open boundary {boundary.name}: {wet.size - wet.sum()} of"
                    f" {wet.size} points of this {boundary.side.name.capitalize()}ern"
                    " boundary are on land",
                )
                if not self.allow_on_land:
                    raise Exception()

            # The open boundary points themselves (T grid)
            bdy_tmask[wet] = 2

            # Velocity points that lie just outside the T points of the open boundary
            bdy_uvmask[wet] = 4

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

        self.np = 0
        self.np_glob = 0

        # Indices of T point in open boundary into arrays that INCLUDE halos.
        # We start with an empty array to support later concatenation in
        # subdomains without any open boundary
        all_i = [np.empty((0,), dtype=np.intp)]
        all_j = [np.empty((0,), dtype=np.intp)]
        side2count = {Side.WEST: 0, Side.NORTH: 0, Side.EAST: 0, Side.SOUTH: 0}
        self.local_to_global: List[slice] = []
        bdy_types_2d = {}
        for boundary in self._boundaries:
            boundary.to_local_grid(grid)

            if boundary.l is not None:
                # This boundary falls at least partially in the local subdomain
                self.active.append(boundary)
                side2count[boundary.side] += 1

                # Determine start/stop indices in arrays with all local boundary points
                mskip = (boundary.mstart - boundary.mstart_) // boundary.mstep
                assert mskip >= 0
                boundary.start = self.np
                boundary.stop = self.np + boundary.np
                boundary.slice_bdy = (slice(boundary.start, boundary.stop), Ellipsis)

                # Determine start/stop indices in arrays with all global boundary points
                start_glob = self.np_glob + mskip
                stop_glob = start_glob + boundary.np
                if self.local_to_global and self.local_to_global[-1].stop == start_glob:
                    # Attach to previous boundary by dropping previous slice
                    # and adopting its start index
                    start_glob = self.local_to_global.pop().start
                self.local_to_global.append(slice(start_glob, stop_glob))

                all_i.append(boundary.i)
                all_j.append(boundary.j)

                if boundary.type_2d not in bdy_types_2d:
                    bdy_types_2d[boundary.type_2d] = self._make_bc(boundary.type_2d)
                boundary.type_2d = bdy_types_2d[boundary.type_2d]

                self.np += boundary.np

            self.np_glob += (boundary.mstop_ - boundary.mstart_) // boundary.mstep

        # Number of open boundary points (local and global)
        grid.open_boundaries = self

        # Local indices of open boundary points within local subdomain
        # (T grid, for arrays that INclude halos)
        self.i = np.concatenate(all_i, dtype=np.intp)
        self.j = np.concatenate(all_j, dtype=np.intp)
        assert self.allow_on_land or (grid.mask.all_values[self.j, self.i] == 2).all()

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

        # Longitude/latitude of open boundary points
        if grid.lon is not None:
            self.lon = grid.array(
                name="lon_bdy",
                on_boundary=True,
                fill=grid.lon.all_values[self.j, self.i],
                attrs=dict(_time_varying=False),
            )
        if grid.lat is not None:
            self.lat = grid.array(
                name="lat_bdy",
                on_boundary=True,
                fill=grid.lat.all_values[self.j, self.i],
                attrs=dict(_time_varying=False),
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

            # Depth-explicit inward velocity at the open boundaries
            self.uv_in = grid.array(z=CENTERS, on_boundary=True)

        # Prescribed depth-averaged or depth-integrated velocity at the open boundaries
        self.u = grid.array(name="u_bdy", on_boundary=True)
        self.v = grid.array(name="v_bdy", on_boundary=True)

        if grid.rotation is not None and self.np > 0:
            # Prescribed geocentric velocity (eastwards/northwards) needs to be rotated
            # to model grid. Create rotator and arrays to hold the rotated velocity.
            self._rotator = core.Rotator(grid.rotation.all_values[self.j, self.i])
            self.u_rot = grid.array(name="u_rot_bdy", on_boundary=True)
            self.v_rot = grid.array(name="v_rot_bdy", on_boundary=True)
        else:
            # Model grid is geocentric - no rotation of prescribed velocities necessary
            self._rotator = None
            self.u_rot = self.u
            self.v_rot = self.v

        self._frozen = True

    def _configure_mirroring(self, grid: core.Grid):
        """Set up indices for mirroring across boundaries on U and V grids"""

        assert grid.ugrid is not None and grid.vgrid is not None
        nx, ny = grid.nx_, grid.ny_

        class Mirror:
            def __init__(self):
                self.i_in = []
                self.j_in = []
                self.i_out = []
                self.j_out = []

            def append(self, i_in, j_in, i_out, j_out, where, check=None, value=None):
                i_in, j_in, i_out, j_out = np.broadcast_arrays(i_in, j_in, i_out, j_out)
                keep = (i_in >= 0) & (j_in >= 0) & (i_in < nx) & (j_in < ny)
                keep &= (i_out >= 0) & (j_out >= 0) & (i_out < nx) & (j_out < ny)
                keep &= where
                if keep.any():
                    self.i_in.append(i_in[keep])
                    self.j_in.append(j_in[keep])
                    self.i_out.append(i_out[keep])
                    self.j_out.append(j_out[keep])
                    if check is not None:
                        assert (check[self.j_out[-1], self.i_out[-1]] == value).all()

            def get_slices(self) -> Optional[Tuple[Tuple, Tuple]]:
                if self.i_in:
                    i_in = np.concatenate(self.i_in, dtype=np.intp)
                    j_in = np.concatenate(self.j_in, dtype=np.intp)
                    i_out = np.concatenate(self.i_out, dtype=np.intp)
                    j_out = np.concatenate(self.j_out, dtype=np.intp)
                    return (Ellipsis, j_in, i_in), (Ellipsis, j_out, i_out)

        mirror_U = Mirror()
        mirror_V = Mirror()
        mirror_TU = Mirror()
        mirror_TV = Mirror()
        tmask = grid.mask.all_values
        umask = grid.ugrid.mask.all_values
        vmask = grid.vgrid.mask.all_values
        for boundary in self.active:
            # Select the T points from the open boundary that are actually wet
            # (some may be on land with mask = 0)
            tsel = tmask[boundary.j, boundary.i] == 2

            # Velocity points that lie WITHIN the open boundary (mask=3),
            # that is, they lie in between tracer points with mask=2.
            # Values at these points will be mirrored from the interior velocity point.
            # We look 1 point beyond the start and stop of the boundary within the
            # current subdomain to catch cases where (a) the boundary extends into the
            # next subdomain and (b) the current boundary neighbors another boundary
            # Note that the values in boundary.{i,j} may have stride -1, which is why
            # we use min()/max() to determine the range.
            if boundary.side in (Side.WEST, Side.EAST):
                j = np.arange(max(boundary.j.min() - 1, 0), boundary.j.max() + 1)
                i_out = boundary.i[0]
                i_in = i_out + {Side.EAST: -1, Side.WEST: 1}[boundary.side]
                mirror_V.append(i_in, j, i_out, j, vmask[j, i_out] == 3)
            else:
                i = np.arange(max(boundary.i.min() - 1, 0), boundary.i.max() + 1)
                j_out = boundary.j[0]
                j_in = j_out + {Side.NORTH: -1, Side.SOUTH: 1}[boundary.side]
                mirror_U.append(i, j_in, i, j_out, umask[j_out, i] == 3)

            # Velocity points OUTSIDE AND ALONG the open boundary (mask=4)
            # Values at these points will be mirrored from either the neighboring
            # inner velocity point, or from inner T point (boundary.i, boundary.j)
            # [e.g., elevations at mask=2]
            if boundary.side in (Side.WEST, Side.EAST):
                j = boundary.j
                i_in = boundary.i[0] + {Side.WEST: 0, Side.EAST: -1}[boundary.side]
                i_out = boundary.i[0] + {Side.WEST: -1, Side.EAST: 0}[boundary.side]
                mirror_U.append(i_in, j, i_out, j, tsel, umask, 4)
                mirror_TU.append(boundary.i, j, i_out, j, tsel, umask, 4)
            else:
                i = boundary.i
                j_in = boundary.j[0] + {Side.SOUTH: 0, Side.NORTH: -1}[boundary.side]
                j_out = boundary.j[0] + {Side.SOUTH: -1, Side.NORTH: 0}[boundary.side]
                mirror_V.append(i, j_in, i, j_out, tsel, vmask, 4)
                mirror_TV.append(i, boundary.j, i, j_out, tsel, vmask, 4)

        grid.ugrid._mirrors[grid.ugrid] = mirror_U.get_slices()
        grid.vgrid._mirrors[grid.vgrid] = mirror_V.get_slices()
        grid._mirrors[grid.ugrid] = mirror_TU.get_slices()
        grid._mirrors[grid.vgrid] = mirror_TV.get_slices()

    def start(
        self,
        U: core.Array,
        V: core.Array,
        uk: Optional[core.Array],
        vk: Optional[core.Array],
        fields: Mapping[str, core.Array],
    ):
        for boundary in self.active:
            # Set up depth-integrated/depth-averaged velocity perpendicular
            # to the open boundary for Flather boundary conditions

            # Modelled depth-integrated velocity (slice from full U/V)
            boundary.UV = boundary.extract_uv_in(U.all_values, V.all_values)

            # Prescribed depth-integrated or depth-averaged velocity
            if boundary.side in (Side.EAST, Side.WEST):
                boundary.flow_ext = self.u_rot[boundary.slice_bdy]
            else:
                boundary.flow_ext = self.v_rot[boundary.slice_bdy]

            if uk is not None:
                # Depth-explicit horizontal model velocities for sponge BCs

                # Velocity perpendicular to the open boundary (slice from full u/v)
                boundary.uv = boundary.extract_uv_in(uk.all_values, vk.all_values).T

                # Inward velocity (to be calculated from uv and inward sign)
                boundary.uv_in = self.uv_in.all_values[boundary.slice_bdy]

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
        """Indices of local open boundary points in a global open boundary array."""
        indices = np.arange(self.np, dtype=np.intp)
        if self.local_to_global is not None:
            i = 0
            for s in self.local_to_global:
                indices[i : i + s.stop - s.start] = np.arange(s.start, s.stop)
                i += s.stop - s.start
            assert i == self.np
        return indices

    def prepare(self, macro_active: bool):
        """Prepare open boundary forcing. This includes rotating prescribed velocity
        fields, as well as extracting depth-explicit horizontal velocities from which
        each class of open boundary condition can calculate derived metrics."""

        if self._rotator:
            # Rotate prescribed geocentric velocity to model grid
            u_rot, v_rot = self._rotator(self.u.all_values, self.v.all_values)
            self.u_rot.all_values[:] = u_rot
            self.v_rot.all_values[:] = v_rot

        if macro_active:
            if self.uv_in.saved:
                # Set modelled depth-explicit inward velocity for all boundaries
                # combined by looping over individual boundaries and assigning
                # to their velocity slices
                for boundary in self.active:
                    np.multiply(boundary.uv, boundary.inflow_sign, out=boundary.uv_in)

            for bc in self.bcs:
                bc.prepare_depth_explicit()
