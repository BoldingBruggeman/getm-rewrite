import enum
from typing import Optional, Sequence, List, TYPE_CHECKING
import functools

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
        mstart_: int,
        mstop_: int,
        type_2d: int,
        type_3d: int,
    ):
        side = Side(side)

        self.name = name
        self.side = side
        self.l = l
        self.mstart = mstart
        self.mstop = mstop
        self.mstart_ = mstart_
        self.mstop_ = mstop_
        self.type_2d = type_2d
        self.type_3d = type_3d
        self.inflow_sign = 1.0 if side in (Side.WEST, Side.SOUTH) else -1.0

        if l is None:
            return

        self.np = mstop - mstart
        mslice = slice(self.mstart, self.mstop)
        if side in (Side.WEST, Side.EAST):
            self.i = np.repeat(l, mstop - mstart)
            self.j = np.arange(mstart, mstop)
            self.slice_t = (Ellipsis, mslice, l)
            self.slice_uv_in = (Ellipsis, mslice, l if side == Side.WEST else l - 1)
            self.slice_uv_out = (Ellipsis, mslice, l - 1 if side == Side.WEST else l)
        else:
            self.i = np.arange(mstart, mstop)
            self.j = np.repeat(l, mstop - mstart)
            self.slice_t = (Ellipsis, l, mslice)
            self.slice_uv_in = (Ellipsis, l if side == Side.SOUTH else l - 1, mslice)
            self.slice_uv_out = (Ellipsis, l - 1 if side == Side.SOUTH else l, mslice)

    def extract_inward(
        self, values: np.ndarray, start: int, stop: Optional[int] = None
    ):
        l_inward = {Side.WEST: 1, Side.EAST: -1, Side.SOUTH: 1, Side.NORTH: -1}[
            self.side
        ]
        mslice = slice(self.mstart, self.mstop)
        if stop is None:
            lslice = self.l + l_inward * start
        else:
            lslice = slice(
                self.l + l_inward * start, self.l + l_inward * stop, l_inward
            )
        if self.side in (Side.WEST, Side.EAST):
            return values[Ellipsis, mslice, lslice]
        else:
            if stop is not None:
                return np.swapaxes(values[Ellipsis, lslice, mslice], -1, -2)
            return values[Ellipsis, lslice, mslice]

    def extract_uv_in(self, u: np.ndarray, v: np.ndarray):
        uv = u if self.side in (Side.WEST, Side.EAST) else v
        return uv[self.slice_uv_in]


class BoundaryCondition:
    def initialize(self, open_boundaries: "OpenBoundaries"):
        pass

    def get_updater(
        self,
        domain: "Domain",
        boundary: "OpenBoundary",
        array: core.Array,
        bdy: np.ndarray,
    ):
        raise NotImplementedError

    def prepare_depth_explicit(self):
        pass


class Sponge(BoundaryCondition):
    def __init__(self, n: int = 3):
        self.n = n
        self.tmrlx = False
        self.tmrlx_max = 0.25
        self.tmrlx_min = 0.0
        self.tmrlx_ucut = 0.02
        self.tmrlx_umin = -0.25 * self.tmrlx_ucut
        self.rlxcoef = None

    def initialize(self, open_boundaries: "OpenBoundaries"):
        if self.tmrlx and self.rlxcoef is None:
            self.rlxcoef = open_boundaries.domain.T.array(z=CENTERS, on_boundary=True)
            open_boundaries.velocity_3d_in.saved = True
            self.inflow = open_boundaries.velocity_3d_in.all_values

    def prepare_depth_explicit(self):
        if self.tmrlx:
            self.rlxcoef.all_values[...] = (self.tmrlx_max - self.tmrlx_min) * np.clip(
                (self.inflow - self.tmrlx_umin) / (self.tmrlx_ucut - self.tmrlx_umin),
                0.0,
                1.0,
            ) + self.tmrlx_min

    def get_updater(
        self, boundary: "OpenBoundary", array: core.Array, bdy: np.ndarray,
    ):
        mask = boundary.extract_inward(
            array.grid.mask.all_values, start=1, stop=self.n + 1
        )
        n = mask.shape[-1]

        # relaxation coeffiicent as per Martinsen & Engedahl (1987), Eq 3
        # https://doi.org/10.1016/0378-3839(87)90028-7
        sp = np.empty((boundary.np, n), dtype=float)
        sp[...] = ((self.n - np.arange(n)) / (self.n + 1.0)) ** 2

        if self.tmrlx:
            sp[mask == 0] = 0.0  # only include water points
            w = sp / sp.sum(axis=1, keepdims=True)
            rlxcoef = self.rlxcoef.all_values[boundary.slice_bdy]
        else:
            w = None
            rlxcoef = None
        sp[mask != 1] = 0.0  # only include updatable water points (exclude bdy)

        return (
            self.update_boundary,
            (
                array.all_values[boundary.slice_t],
                boundary.extract_inward(array.all_values, start=1, stop=n + 1),
                bdy,
                sp,
                rlxcoef,
                w,
            ),
        )

    @staticmethod
    def update_boundary(
        values: np.ndarray,
        sponge_values: np.ndarray,
        bdy_values: np.ndarray,
        sp: np.ndarray,
        rlxcoef: np.ndarray,
        w: Optional[np.ndarray],
    ):
        if w is not None:
            # note: where=w != 0.0 is used to avoid mixing in NaNs from areas where w=0
            sponge_mean = (w * sponge_values).sum(axis=-1, where=w != 0.0).T
            bdy_values = rlxcoef * bdy_values + (1.0 - rlxcoef) * sponge_mean
        bdy_values = bdy_values.T
        blend = sp * bdy_values[..., np.newaxis] + (1.0 - sp) * sponge_values
        sponge_values[...] = blend
        values[...] = bdy_values


class ZeroGradient(BoundaryCondition):
    def get_updater(
        self, boundary: "OpenBoundary", array: core.Array, bdy: np.ndarray,
    ):
        return (
            self.update,
            (
                array.all_values[boundary.slice_t],
                boundary.extract_inward(array.all_values, start=1),
            ),
        )

    @staticmethod
    def update(values: np.ndarray, inward_values: np.ndarray):
        values[:] = inward_values


class Clamped(BoundaryCondition):
    def get_updater(
        self, boundary: "OpenBoundary", array: core.Array, bdy: np.ndarray,
    ):
        return self.update, (array.all_values[boundary.slice_t], bdy.T)

    @staticmethod
    def update(values: np.ndarray, prescribed_values: np.ndarray):
        values[:] = prescribed_values


class Flather(BoundaryCondition):
    def __init__(self, transport: bool = False):
        self.update = self.update_transport if transport else self.update_velocity

    def get_updater(
        self, boundary: "OpenBoundary", array: core.Array, bdy: np.ndarray,
    ):
        return (
            self.update,
            (
                array.all_values[boundary.slice_t],
                bdy,
                boundary.tp,
                boundary.flow_ext,
                array.grid.D.all_values[boundary.slice_t],
                boundary.inflow_sign,
            ),
        )

    @staticmethod
    def update_velocity(
        z: np.ndarray,
        z_ext: np.ndarray,
        tp: np.ndarray,
        vel_ext: np.ndarray,
        D: np.ndarray,
        inflow_sign: float,
    ):
        z[:] = z_ext - inflow_sign * (tp - vel_ext * D) / np.sqrt(D * GRAVITY)

    @staticmethod
    def update_transport(
        z: np.ndarray,
        z_ext: np.ndarray,
        tp: np.ndarray,
        tp_ext: np.ndarray,
        D: np.ndarray,
        inflow_sign: float,
    ):
        z[:] = z_ext - inflow_sign * (tp - tp_ext) / np.sqrt(D * GRAVITY)


def make_bc(open_boundaries: "OpenBoundaries", value) -> BoundaryCondition:
    if isinstance(value, BoundaryCondition):
        return value
    if value == ZERO_GRADIENT:
        return open_boundaries.zero_gradient
    elif value == SPONGE:
        return open_boundaries.sponge
    elif value == CLAMPED:
        return Clamped()
    elif value == FLATHER_ELEV:
        return Flather()
    elif value == FLATHER_TRANSPORT:
        return Flather(transport=True)
    else:
        raise Exception(f"Unknown boundary type {value} specified")


class ArrayOpenBoundary:
    def __init__(
        self,
        domain: "Domain",
        boundary: OpenBoundary,
        prescribed_values: np.ndarray,
        model_values: np.ndarray,
        type: BoundaryCondition,
    ):
        self.domain = domain
        self._boundary = boundary
        self._prescribed_values = prescribed_values
        self._model_values = model_values
        self._type = type

    @property
    def type(self) -> BoundaryCondition:
        return self._type

    @type.setter
    def type(self, value: int):
        self._type = make_bc(self.domain.open_boundaries, value)

    @property
    def values(self) -> int:
        return self._prescribed_values


class ArrayOpenBoundaries:
    __slots__ = "_array", "values", "_bdy", "updaters"

    def __init__(self, array: core.Array, type=None):
        self._array = array
        self.values = array.grid.array(
            name="%s_bdy" % array.name,
            z=array.z,
            on_boundary=True,
            fill_value=array.fill_value,
            attrs={
                "_time_varying": array.attrs.get("_time_varying", TimeVarying.MICRO)
            },
        )
        if type is not None:
            type = make_bc(array.grid.domain.open_boundaries, type)
        self._bdy = []
        for bdy in array.grid.domain.open_boundaries.active:
            self._bdy.append(
                ArrayOpenBoundary(
                    array.grid.domain,
                    bdy,
                    self.values.all_values[bdy.slice_bdy],
                    array.all_values[bdy.slice_t],
                    type or (bdy.type_2d if array.ndim == 2 else bdy.type_3d),
                )
            )
        self.updaters = []

    def _set_type(self, value: int):
        value = make_bc(self._array.grid.domain.open_boundaries, value)
        for bdy in self._bdy:
            bdy.type = value

    type = property(fset=_set_type)

    def initialize(self):
        bcs = set()
        for bdy in self._bdy:
            bdy._type.initialize(self._array.grid.domain.open_boundaries)
            fn, args = bdy._type.get_updater(
                bdy._boundary, self._array, bdy._prescribed_values
            )
            self.updaters.append(functools.partial(fn, *args))
            bcs.add(bdy._type)
        return bcs

    def update(self):
        """Update the tracer at the open boundaries"""
        for updater in self.updaters:
            updater()


class OpenBoundaries(Sequence[OpenBoundary]):
    __slots__ = (
        "domain",
        "np",
        "np_glob",
        "i",
        "j",
        "i_glob",
        "j_glob",
        "z",
        "u",
        "v",
        "lon",
        "lat",
        "zc",
        "zf",
        "local_to_global",
        "_boundaries",
        "_frozen",
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

    def __init__(self, domain: "Domain"):
        self.domain = domain
        self._boundaries: List[OpenBoundary] = []
        self.sponge = Sponge()
        self.zero_gradient = ZeroGradient()
        self.active: List[OpenBoundary] = []
        self._frozen = False
        self.bcs = set()

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
        # NB below we convert to indices in the T grid of the current subdomain
        # INCLUDING halos
        # We also limit the indices to the range valid for the current subdomain.
        xoffset, yoffset = (
            self.domain.tiling.xoffset - self.domain.halox,
            self.domain.tiling.yoffset - self.domain.haloy,
        )
        if side in (Side.WEST, Side.EAST):
            l_offset, m_offset, l_max, m_max = (
                xoffset,
                yoffset,
                self.domain.T.nx_,
                self.domain.T.ny_,
            )
        else:
            l_offset, m_offset, l_max, m_max = (
                yoffset,
                xoffset,
                self.domain.T.ny_,
                self.domain.T.nx_,
            )
        l_loc = l - l_offset
        mstart_loc_ = mstart - m_offset
        mstop_loc_ = mstop - m_offset
        mstart_loc = min(max(0, mstart_loc_), m_max)
        mstop_loc = min(max(0, mstop_loc_), m_max)
        if l_loc < 0 or l_loc >= l_max or mstop_loc <= mstart_loc:
            # Boundary lies completely outside current subdomain. Record it anyway,
            # so we can later set up a global -> local map of open boundary points
            l_loc, mstart_loc, mstop_loc = None, None, None
        if name is None:
            name = str(len(self._boundaries))
        self._boundaries.append(
            OpenBoundary(
                name,
                side,
                l_loc,
                mstart_loc,
                mstop_loc,
                mstart_loc_,
                mstop_loc_,
                type_2d,
                type_3d,
            )
        )

        if self.domain.glob is not None and self.domain.glob is not self.domain:
            self.domain.glob.open_boundaries.add_by_index(
                side, l, mstart, mstop, type_2d, type_3d
            )

    def initialize(self):
        """Freeze the open boundary collection. Drop those outside the current
        subdomain.
        """
        assert (
            not self._frozen
        ), "The open boundary collection has already been initialized"

        nbdyp = 0
        nbdyp_glob = 0
        bdy_i, bdy_j = [], []
        side2count = {}
        self.local_to_global = []
        umask = self.domain.mask_[1::2, 2::2]
        vmask = self.domain.mask_[2::2, 1::2]
        tmask = self.domain.mask_[1::2, 1::2]
        bdy_types_2d = {}
        for side in (Side.WEST, Side.NORTH, Side.EAST, Side.SOUTH):
            n = 0
            for boundary in [b for b in self._boundaries if b.side == side]:
                if boundary.l is not None:
                    mskip = boundary.mstart - boundary.mstart_
                    assert mskip >= 0
                    boundary.start = nbdyp
                    boundary.stop = nbdyp + boundary.np
                    boundary.slice_bdy = (
                        slice(boundary.start, boundary.stop),
                        Ellipsis,
                    )
                    if boundary.type_2d not in bdy_types_2d:
                        bdy_types_2d[boundary.type_2d] = make_bc(self, boundary.type_2d)
                    boundary.type_2d = bdy_types_2d[boundary.type_2d]
                    nbdyp += boundary.np

                    bdy_mask = tmask[boundary.slice_t]
                    if (bdy_mask == 0).any():
                        self.domain.logger.error(
                            "Open boundary %s: %i of %i points of this %sern boundary"
                            " are on land"
                            % (
                                boundary.name,
                                (bdy_mask == 0).sum(),
                                boundary.np,
                                side.name.capitalize(),
                            )
                        )
                        raise Exception()

                    bdy_mask[:] = 2
                    if side in (Side.WEST, Side.EAST):
                        umask[boundary.slice_uv_out] = 4
                    else:
                        vmask[boundary.slice_uv_out] = 4
                    bdy_i.append(boundary.i)
                    bdy_j.append(boundary.j)

                    start_glob = nbdyp_glob + mskip
                    if (
                        self.local_to_global
                        and self.local_to_global[-1][1] == start_glob
                    ):
                        # attach to previous boundary
                        self.local_to_global[-1][1] += boundary.np
                    else:
                        # gap; add new slice
                        self.local_to_global.append(
                            [start_glob, start_glob + boundary.np]
                        )
                    n += 1
                    self.active.append(boundary)
                nbdyp_glob += boundary.mstop_ - boundary.mstart_
            side2count[side] = n

        umask[:, :-1][(tmask[:, :-1] == 2) & (tmask[:, 1:] == 2)] = 3
        vmask[:-1, :][(tmask[:-1, :] == 2) & (tmask[1:, :] == 2)] = 3

        self.mirror_U = []
        self.mirror_V = []
        self.mirror_TU = []
        self.mirror_TV = []
        for boundary in self._boundaries:
            if boundary.l is None:
                continue

            # Identify velocity points that lie within the open boundary in between
            # tracer points with mask=2. Their indices will be (i_velout, j_velout)
            # on the corresponding velocity grid (U or V). Values at these points
            # will be mirrored from the interior velocity point (i_velin, j_velin)
            if boundary.side in (Side.WEST, Side.EAST):
                i_velout = np.repeat(boundary.i[0], boundary.i.size + 1)
                j_velout = np.arange(boundary.j[0] - 1, boundary.j[-1] + 1)
                mirror_mask = self.domain.mask_[2 + j_velout * 2, 1 + i_velout * 2]
            else:
                i_velout = np.arange(boundary.i[0] - 1, boundary.i[-1] + 1)
                j_velout = np.repeat(boundary.j[0], boundary.j.size + 1)
                mirror_mask = self.domain.mask_[1 + j_velout * 2, 2 + i_velout * 2]
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
                self.mirror_U.append([i_in, j_in, i_out, j_out])
                self.mirror_V.append([i_velin, j_velin, i_velout, j_velout])
                self.mirror_TU.append([boundary.i, boundary.j, i_out, j_out])
            else:
                self.mirror_V.append([i_in, j_in, i_out, j_out])
                self.mirror_U.append([i_velin, j_velin, i_velout, j_velout])
                self.mirror_TV.append([boundary.i, boundary.j, i_out, j_out])

        def pack_mirror(indices):
            if indices:
                i_in = np.concatenate([i[0] for i in indices], dtype=np.intp)
                j_in = np.concatenate([i[1] for i in indices], dtype=np.intp)
                i_out = np.concatenate([i[2] for i in indices], dtype=np.intp)
                j_out = np.concatenate([i[3] for i in indices], dtype=np.intp)
                return ((Ellipsis, j_in, i_in), (Ellipsis, j_out, i_out))

        self.mirror_U = pack_mirror(self.mirror_U)
        self.mirror_V = pack_mirror(self.mirror_V)
        self.mirror_TU = pack_mirror(self.mirror_TU)
        self.mirror_TV = pack_mirror(self.mirror_TV)

        self.np = nbdyp
        self.np_glob = nbdyp_glob
        self.i = (
            np.empty((0,), dtype=np.intc)
            if self.np == 0
            else np.concatenate(bdy_i, dtype=np.intc)
        )
        self.j = (
            np.empty((0,), dtype=np.intc)
            if self.np == 0
            else np.concatenate(bdy_j, dtype=np.intc)
        )

        self.i_glob = self.i - self.domain.halox + self.domain.tiling.xoffset
        self.j_glob = self.j - self.domain.haloy + self.domain.tiling.yoffset
        self.domain.logger.info(
            "%i open boundaries (%i West, %i North, %i East, %i South)"
            % (
                sum(side2count.values()),
                side2count[Side.WEST],
                side2count[Side.NORTH],
                side2count[Side.EAST],
                side2count[Side.SOUTH],
            )
        )
        if self.np > 0:
            if self.np == self.np_glob:
                assert (
                    len(self.local_to_global) == 1
                    and self.local_to_global[0][0] == 0
                    and self.local_to_global[0][1] == self.np_glob
                )
                self.local_to_global = None
            else:
                self.domain.logger.info(
                    "global-to-local open boundary map: %s" % (self.local_to_global,)
                )

        # Coordinates of open boundary points
        self.zc = self.domain.T.array(z=CENTERS, on_boundary=True)
        self.zf = self.domain.T.array(z=INTERFACES, on_boundary=True)
        if self.domain.lon is not None:
            self.lon = self.domain.T.array(
                on_boundary=True, fill=self.domain.T.lon.all_values[self.j, self.i]
            )
        if self.domain.lat is not None:
            self.lat = self.domain.T.array(
                on_boundary=True, fill=self.domain.T.lat.all_values[self.j, self.i]
            )

        # The arrays below are placeholders that will be assigned data
        # (from momentum/sealevel Fortran modules) when linked to the Simulation
        self.u = self.domain.T.array(name="u_bdy", on_boundary=True)
        self.v = self.domain.T.array(name="v_bdy", on_boundary=True)

        self.velocity_3d_in = self.domain.T.array(z=CENTERS, on_boundary=True)

        self._frozen = True

    def start(self, U: core.Array, V: core.Array, uk: core.Array, vk: core.Array):
        for boundary in self.active:
            boundary.tp = boundary.extract_uv_in(U.all_values, V.all_values)
            if uk is not None:
                boundary.vel = boundary.extract_uv_in(uk.all_values, vk.all_values)
            if boundary.side in (Side.EAST, Side.WEST):
                boundary.flow_ext = self.u[boundary.slice_bdy]
            else:
                boundary.flow_ext = self.v[boundary.slice_bdy]
            boundary.velocity_3d_in = self.velocity_3d_in.all_values[boundary.slice_bdy]
        for field in self.domain.fields.values():
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
            for start, stop in self.local_to_global:
                indices[i : i + stop - start] = np.arange(start, stop)
                i += stop - start
            assert i == self.np
        return indices

    def prepare_depth_explicit(self):
        if self.velocity_3d_in.saved:
            for boundary in self.active:
                boundary.velocity_3d_in[:] = boundary.inflow_sign * boundary.vel.T
        for bc in self.bcs:
            bc.prepare_depth_explicit()
