from typing import Callable, MutableMapping, Optional, Tuple, List, Mapping
import operator
import logging
import os.path
import enum
import functools

import numpy as np
from numpy.typing import ArrayLike
import xarray
import netCDF4

from . import _pygetm
from . import core
from . import parallel
from . import input
from .constants import FILL_VALUE, CENTERS, GRAVITY, INTERFACES


class VerticalCoordinates(enum.IntEnum):
    SIGMA = 1
    #    Z = 2
    GVC = 3
    #    HYBRID = 4
    #    ADAPTIVE = 5


class Grid(_pygetm.Grid):
    _coordinate_arrays = "x", "y", "lon", "lat"
    _readonly_arrays = _coordinate_arrays + (
        "dx",
        "dy",
        "idx",
        "idy",
        "dlon",
        "dlat",
        "area",
        "iarea",
        "cor",
    )
    _fortran_arrays = _readonly_arrays + (
        "H",
        "D",
        "mask",
        "z",
        "zo",
        "ho",
        "hn",
        "zc",
        "zf",
        "z0b",
        "z0b_min",
        "zio",
        "zin",
        "alpha",
    )
    _all_arrays = tuple(
        ["_%s" % n for n in _fortran_arrays]
        + ["_%si" % n for n in _coordinate_arrays]
        + ["_%si_" % n for n in _coordinate_arrays]
    )
    __slots__ = _all_arrays + (
        "halo",
        "type",
        "ioffset",
        "joffset",
        "postfix",
        "ugrid",
        "vgrid",
        "_sin_rot",
        "_cos_rot",
        "rotation",
        "nbdyp",
        "overlap",
        "_interpolators",
        "rotated",
        "_water_contact",
        "_land",
        "_water",
    )

    _array_args = {
        "x": dict(units="m", attrs=dict(_time_varying=False)),
        "y": dict(units="m", attrs=dict(_time_varying=False)),
        "lon": dict(
            units="degrees_east",
            long_name="longitude",
            attrs=dict(standard_name="longitude", axis="X", _time_varying=False),
        ),
        "lat": dict(
            units="degrees_north",
            long_name="latitude",
            attrs=dict(standard_name="latitude", axis="Y", _time_varying=False),
        ),
        "dx": dict(units="m", attrs=dict(_time_varying=False)),
        "dy": dict(units="m", attrs=dict(_time_varying=False)),
        "idx": dict(units="m-1", attrs=dict(_time_varying=False)),
        "idy": dict(units="m-1", attrs=dict(_time_varying=False)),
        "dlon": dict(units="degrees_east", attrs=dict(_time_varying=False)),
        "dlat": dict(units="degrees_north", attrs=dict(_time_varying=False)),
        "H": dict(
            units="m", long_name="water depth at rest", attrs=dict(_time_varying=False)
        ),
        "D": dict(
            units="m",
            long_name="water depth",
            attrs=dict(standard_name="sea_floor_depth_below_sea_surface"),
        ),
        "mask": dict(attrs=dict(_time_varying=False), fill_value=0),
        "z": dict(units="m", long_name="elevation"),
        "zo": dict(units="m", long_name="elevation at previous microtimestep"),
        "zin": dict(units="m", long_name="elevation at macrotimestep"),
        "zio": dict(units="m", long_name="elevation at previous macrotimestep"),
        "area": dict(
            units="m2",
            long_name="cell area",
            attrs=dict(standard_name="cell_area", _time_varying=False),
        ),
        "iarea": dict(
            units="m-2",
            long_name="inverse of cell area",
            attrs=dict(_time_varying=False),
        ),
        "cor": dict(
            units="1", long_name="Coriolis parameter", attrs=dict(_time_varying=False)
        ),
        "ho": dict(units="m", long_name="cell thickness at previous time step"),
        "hn": dict(
            units="m",
            long_name="cell thickness",
            attrs=dict(standard_name="cell_thickness"),
        ),
        "zc": dict(
            units="m",
            long_name="height",
            attrs=dict(
                axis="Z", positive="up", standard_name="height_above_mean_sea_level"
            ),
        ),
        "zf": dict(
            units="m",
            long_name="interface height",
            attrs=dict(
                axis="Z", positive="up", standard_name="height_above_mean_sea_level"
            ),
        ),
        "z0b": dict(units="m", long_name="hydrodynamic bottom roughness"),
        "z0b_min": dict(
            units="m",
            long_name="physical bottom roughness",
            attrs=dict(_time_varying=False),
        ),
        "alpha": dict(units="1", long_name="dampening"),
    }
    domain: "Domain"

    def __init__(
        self,
        domain: "Domain",
        grid_type: int,
        ioffset: int,
        joffset: int,
        overlap: int = 0,
        ugrid: Optional["Grid"] = None,
        vgrid: Optional["Grid"] = None,
    ):
        _pygetm.Grid.__init__(self, domain, grid_type)
        self.halo = domain.halo
        self.type = grid_type
        self.ioffset = ioffset
        self.joffset = joffset
        self.overlap = overlap
        self.postfix = {
            _pygetm.TGRID: "t",
            _pygetm.UGRID: "u",
            _pygetm.VGRID: "v",
            _pygetm.XGRID: "x",
            _pygetm.UUGRID: "_uu_adv",
            _pygetm.VVGRID: "_vv_adv",
            _pygetm.UVGRID: "_uv_adv",
            _pygetm.VUGRID: "_vu_adv",
        }[grid_type]
        self.ugrid: Optional[Grid] = ugrid
        self.vgrid: Optional[Grid] = vgrid
        self._sin_rot: Optional[np.ndarray] = None
        self._cos_rot: Optional[np.ndarray] = None
        self._interpolators = {}

        for name in self._readonly_arrays:
            self._setup_array(name)
        with np.errstate(divide="ignore"):
            self._iarea.all_values[...] = 1.0 / self._area.all_values
            self._idx.all_values[...] = 1.0 / self._dx.all_values
            self._idy.all_values[...] = 1.0 / self._dy.all_values
        for name in self._readonly_arrays:
            getattr(self, name).all_values.flags.writeable = False

    def initialize(self, nbdyp: int):
        for name in self._fortran_arrays:
            if not hasattr(self, "_%s" % name):
                self._setup_array(name)
        self.rotation = core.Array.create(
            grid=self,
            dtype=self.x.dtype,
            name="rotation" + self.postfix,
            units="rad",
            long_name="grid rotation with respect to true North",
            fill_value=np.nan,
        )
        self._setup_array("rotation", self.rotation)

        self._land = self.mask.all_values == 0
        self._water = ~self._land
        self.zc.all_values[...] = -self.H.all_values
        self.zf.all_values[...] = -self.H.all_values

        self.H.all_values[self._land] = FILL_VALUE
        self.z0b_min.all_values[self._land] = FILL_VALUE
        self.H.all_values.flags.writeable = False
        self.z0b_min.all_values.flags.writeable = False
        self.z0b.all_values[...] = self.z0b_min.all_values

        self.z.all_values[self._water] = 0.0
        self.zo.all_values[...] = self.z.all_values
        self.zio.all_values[...] = self.z.all_values
        self.zin.all_values[...] = self.z.all_values

        self.nbdyp = nbdyp
        self.rotated = self.rotation.all_values[self._water].any()

    def _setup_array(
        self, name: str, array: Optional[core.Array] = None, from_supergrid: bool = True
    ) -> core.Array:
        if array is None:
            # No array provided, so it must live in Fortran; retrieve it
            kwargs = dict(fill_value=FILL_VALUE)
            kwargs.update(self._array_args[name])
            array = core.Array(name=name + self.postfix, **kwargs)
            setattr(self, "_%s" % name, self.wrap(array, name.encode("ascii")))

        # Obtain corresponding array on the supergrid.
        # If this does not exist, we are done
        source = getattr(self.domain, name + "_", None)
        if source is None or not from_supergrid:
            return array

        imax, jmax = self.ioffset + 2 * self.nx_, self.joffset + 2 * self.ny_
        values = source[(slice(self.joffset, jmax, 2), slice(self.ioffset, imax, 2))]
        slc = (Ellipsis,)
        if name in ("z0b_min",):
            slc = self.mask.all_values[: values.shape[0], : values.shape[1]] > 0
        array.all_values[: values.shape[0], : values.shape[1]][slc] = values[slc]
        if values.shape != array.all_values.shape:
            # supergrid does not span entire grid; fill remainder by exchanging halos
            array.update_halos()

        has_bounds = (
            self.ioffset > 0
            and self.joffset > 0
            and source.shape[-1] >= imax
            and source.shape[-2] >= jmax
        )
        if has_bounds and name in self._coordinate_arrays:
            # Generate interface coordinates. These are not represented in Fortran as
            # they are only needed for plotting. The interface coordinates are slices
            # that point to the supergrid data; they thus do not consume additional
            # memory.
            values_i = source[self.joffset - 1 : jmax : 2, self.ioffset - 1 : imax : 2]
            setattr(self, "_%si_" % name, values_i)
            setattr(
                self,
                "_%si" % name,
                values_i[self.halo : -self.halo, self.halo : -self.halo],
            )

        return array

    def interpolator(self, target: "Grid") -> Callable[[np.ndarray, np.ndarray], None]:
        ip = self._interpolators.get(target)
        if ip:
            return ip
        assert self.domain is target.domain
        if self.ioffset == target.ioffset + 1 and self.joffset == target.joffset:
            ip = functools.partial(_pygetm.interp_x, offset=1)
        elif self.ioffset == target.ioffset - 1 and self.joffset == target.joffset:
            ip = functools.partial(_pygetm.interp_x, offset=0)
        elif self.joffset == target.joffset + 1 and self.ioffset == target.ioffset:
            ip = functools.partial(_pygetm.interp_y, offset=1)
        elif self.joffset == target.joffset - 1 and self.ioffset == target.ioffset:
            ip = functools.partial(_pygetm.interp_y, offset=0)
        elif self.ioffset == target.ioffset - 1 and self.joffset == target.joffset - 1:
            ip = functools.partial(_pygetm.interp_xy, ioffset=0, joffset=0)
        elif self.ioffset == target.ioffset + 1 and self.joffset == target.joffset + 1:
            ip = functools.partial(_pygetm.interp_xy, ioffset=1, joffset=1)
        else:
            raise NotImplementedError(
                "Cannot interpolate from grid type %s to grid type %s"
                % (self.postfix, target.postfix)
            )
        self._interpolators[target] = ip
        return ip

    def rotate(
        self, u: ArrayLike, v: ArrayLike, to_grid: bool = True
    ) -> Tuple[ArrayLike, ArrayLike]:
        """Rotate a geocentric velocity field to the model coordinate system,
        or a model velocity field to the geocentric coordinate system.

        Args:
            u: velocity in x-direction in source coordinate system
                (Eastward velocity if the source is a geocentric velocity field)
            v: velocity in y-direction in source coordinate system
                (Northward velocity if the source is a geocentric velocity field)
            to_grid: rotate from geocentric to model coordinate system, not vice versa
        """
        if not self.rotated:
            return u, v
        elif self._sin_rot is None:
            self._sin_rot = np.sin(self.rotation.all_values)
            self._cos_rot = np.cos(self.rotation.all_values)

            # hardcode cos(0.5*pi)=0 to increase precision in 90 degree rotaton tests
            self._cos_rot[self.rotation.all_values == 0.5 * np.pi] = 0
        sin_rot = -self._sin_rot if to_grid else self._sin_rot
        u_new = u * self._cos_rot - v * sin_rot
        v_new = u * sin_rot + v * self._cos_rot
        return u_new, v_new

    def array(self, *args, **kwargs) -> core.Array:
        return core.Array.create(self, *args, **kwargs)

    def add_to_netcdf(self, nc: netCDF4.Dataset, postfix: str = ""):
        xdim, ydim = "x" + postfix, "y" + postfix

        def save(name, units="", long_name=None):
            data = getattr(self, name)
            ncvar = nc.createVariable(name + postfix, data.dtype, (ydim, xdim))
            ncvar[...] = data
            ncvar.units = units
            ncvar.long_name = long_name or name

        ny, nx = self.x.shape
        nc.createDimension(xdim, nx)
        nc.createDimension(ydim, ny)
        save("dx", "m")
        save("dy", "m")
        save("H", "m", "undisturbed water depth")
        save("mask")
        save("area", "m2")
        save("cor", "s-1", "Coriolis parameter")

    def nearest_point(
        self,
        x: float,
        y: float,
        mask: Optional[Tuple[int]] = None,
        include_halos: bool = False,
    ) -> Optional[Tuple[int, int]]:
        """Return index (i,j) of point nearest to specified coordinate."""
        if not self.domain.contains(x, y, include_halos=include_halos):
            return None
        local_slice, _, _, _ = self.domain.tiling.subdomain2slices(
            halo_sub=self.halo,
            halo_glob=self.halo,
            share=self.overlap,
            exclude_halos=not include_halos,
            exclude_global_halos=True,
        )
        allx, ally = (self.lon, self.lat) if self.domain.spherical else (self.x, self.y)
        actx, acty = allx.all_values[local_slice], ally.all_values[local_slice]
        dist = (actx - x) ** 2 + (acty - y) ** 2
        if mask is not None:
            if isinstance(mask, int):
                mask = (mask,)
            invalid = np.ones(dist.shape, dtype=bool)
            for mask_value in mask:
                invalid &= self.mask.all_values[local_slice] != mask_value
            dist[invalid] = np.inf
        idx = np.nanargmin(dist)
        j, i = np.unravel_index(idx, dist.shape)
        return j + local_slice[-2].start, i + local_slice[-1].start


for membername in Grid._all_arrays:
    info = Grid._array_args.get(membername[1:], {})
    long_name = info.get("long_name")
    units = info.get("units")
    doc = ""
    if long_name:
        doc = long_name
        if units:
            doc += " (%s)" % units
    setattr(Grid, membername[1:], property(operator.attrgetter(membername), doc=doc))


class RiverTracer(core.Array):
    __slots__ = ("_follow",)

    def __init__(
        self,
        grid: Grid,
        river_name: str,
        tracer_name: str,
        value: np.ndarray,
        follow: np.ndarray,
        **kwargs
    ):
        super().__init__(
            grid=grid,
            name=tracer_name + "_in_river_" + river_name,
            long_name="%s in river %s" % (tracer_name, river_name),
            **kwargs
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
    def __init__(
        self,
        name: str,
        i: int,
        j: int,
        zl: Optional[float] = np.inf,
        zu: Optional[float] = 0.0,
        x: Optional[float] = None,
        y: Optional[float] = None,
    ):
        self.name = name
        self.i_glob = i
        self.j_glob = j
        self.x = x
        self.y = y
        self.zl = zl
        self.zu = zu
        self.i = None
        self.j = None
        self._tracers: Mapping[str, RiverTracer] = {}

    def locate(self, grid: Grid) -> bool:
        if self.x is not None:
            # Location is specified by x, y coordinate.
            # Look up nearest unmasked grid cell.
            ind = grid.nearest_point(self.x, self.y, mask=1, include_halos=True)
            if ind is None:
                grid.domain.logger.info(
                    "River %s at x=%s, y=%s not present in this subdomain"
                    % (self.name, self.x, self.y)
                )
                if grid.domain is grid.domain.glob:
                    raise Exception(
                        "River %s is located at x=%s, y=%s, "
                        "which does not fall within the global model domain"
                        % (self.name, self.x, self.y)
                    )
                return False
            self.j, self.i = ind
            grid.domain.logger.info(
                "River %s at x=%s, y=%s is located at i=%i, j=%i in this subdomain"
                % (self.name, self.x, self.y, self.i, self.j)
            )
        else:
            # Location is specified by global i, j. Map to subdomain.
            i_loc = self.i_glob - grid.domain.tiling.xoffset + grid.domain.halox
            j_loc = self.j_glob - grid.domain.tiling.yoffset + grid.domain.haloy
            if i_loc < 0 or j_loc < 0 or i_loc >= grid.nx_ or j_loc >= grid.ny_:
                return False
            self.i, self.j = i_loc, j_loc
        return True

    def initialize(self, grid: Grid, flow: np.ndarray):
        assert self.i is not None and self.j is not None
        self.flow = core.Array(
            grid=grid,
            name="river_" + self.name + "_flow",
            units="m3 s-1",
            long_name="inflow from %s" % self.name,
        )
        self.flow.wrap_ndarray(flow)

    def __getitem__(self, key) -> RiverTracer:
        return self._tracers[key]

    def __len__(self):
        return len(self._tracers)

    def __iter__(self):
        return iter(self._tracers)


class Rivers(Mapping[str, River]):
    def __init__(self, grid: Grid):
        self.grid = grid
        self._rivers: List[River] = []
        self._frozen = False

    def add_by_index(self, name: str, i: int, j: int, **kwargs):
        """Add a river at a location specified by the indices of a tracer point

        Args:
            name: river name
            i: index in x-direction (0-based)
            j: index in y-direction (0-based)
            **kwargs: additional keyword arguments passed to :class:`River`
        """
        assert not self._frozen, (
            "The river collection has already been initialized"
            " and can no longer be modified."
        )

        domain = self.grid.domain
        if domain.glob is not None and domain.glob is not domain:
            domain.glob.rivers.add_by_index(name, i, j, **kwargs)

        river = River(name, i, j, **kwargs)
        self._rivers.append(river)
        return river

    def add_by_location(self, name: str, x: float, y: float, **kwargs):
        """Add a river at a location specified by the nearest coordinates
        (longitude and latitude on a spherical grid)
        """
        if (
            self.grid.domain.glob is not None
            and self.grid.domain.glob is not self.grid.domain
        ):
            self.grid.domain.glob.rivers.add_by_location(name, x, y, **kwargs)
        river = River(name, None, None, x=x, y=y, **kwargs)
        self._rivers.append(river)
        return river

    def initialize(self):
        """Freeze the river collection. Drop those outside the current subdomain
        and verify the remaining ones are on unmasked T points.
        """
        assert not self._frozen, "The river collection has already been initialized"
        self._frozen = True
        self._rivers = [river for river in self._rivers if river.locate(self.grid)]
        self.flow = np.zeros((len(self._rivers),))
        for iriver, river in enumerate(self._rivers):
            mask = self.grid.mask.all_values[river.j, river.i]
            if mask != 1:
                raise Exception(
                    "River %s is located at i=%i, j=%i, which is not water"
                    " (it has mask value %i)."
                    % (river.name, river.i_glob, river.j_glob, mask)
                )
            river.initialize(self.grid, self.flow[..., iriver])
        self.i = np.array([river.i for river in self._rivers], dtype=int)
        self.j = np.array([river.j for river in self._rivers], dtype=int)
        self.iarea = self.grid.iarea.all_values[self.j, self.i]
        self.zl = np.array([river.zl for river in self._rivers])
        self.zu = np.array([river.zu for river in self._rivers])

    def __getitem__(self, key) -> River:
        for river in self._rivers:
            if key == river.name:
                return river
        raise KeyError()

    def __len__(self) -> int:
        return len(self._rivers)

    def __iter__(self):
        return map(operator.attrgetter("name"), self._rivers)


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
        self.name = name
        self.side = Side(side)
        self.l = l
        self.mstart = mstart
        self.mstop = mstop
        self.mstart_ = mstart_
        self.mstop_ = mstop_
        self.type_2d = type_2d
        self.type_3d = type_3d


class BoundaryCondition:
    def initialize(self, open_boundaries: "OpenBoundaries"):
        pass

    def __call__(self, array: core.Array, bdy: Optional[core.Array] = None):
        raise NotImplementedError


class Sponge(BoundaryCondition):
    def __init__(self, n: int = 3):
        self.n = n
        self.tmrlx = False
        self.tmrlx_max = 0.25
        self.tmrlx_min = 0.0
        self.tmrlx_ucut = 0.02
        self.tmrlx_umin = -0.25 * self.tmrlx_ucut

    def initialize(self, open_boundaries: "OpenBoundaries"):
        tmask = open_boundaries.domain.mask_[1::2, 1::2]
        self.i = np.empty((open_boundaries.np, self.n), dtype=np.intc)
        self.j = np.empty((open_boundaries.np, self.n), dtype=np.intc)
        for boundary in open_boundaries.active:
            i_inward = {Side.WEST: 1, Side.EAST: -1}.get(boundary.side, 0)
            j_inward = {Side.SOUTH: 1, Side.NORTH: -1}.get(boundary.side, 0)
            for n in range(self.n):
                self.i[boundary.start : boundary.stop, n] = (
                    boundary.i + (n + 1) * i_inward
                )
                self.j[boundary.start : boundary.stop, n] = (
                    boundary.j + (n + 1) * j_inward
                )
        i_valid = (self.i >= 0) & (self.i < tmask.shape[1])
        j_valid = (self.j >= 0) & (self.j < tmask.shape[0])
        water = tmask[self.j, self.i] != 0  # can include other bdy points
        sponge_valid = i_valid & j_valid & water
        self.sp = np.empty(self.i.shape, dtype=float)
        self.sp[...] = ((self.n - np.arange(self.n)) / (self.n + 1.0)) ** 2
        self.sp[~sponge_valid] = 0.0
        self.w = self.sp / self.sp.sum(axis=1, keepdims=True)
        self.sp[tmask[self.j, self.i] != 1] = 0.0  # only relax water points
        self.open_boundaries = open_boundaries

        self.inflow = open_boundaries.domain.T.array(z=CENTERS, on_boundary=True)
        self.rlxcoef = open_boundaries.domain.T.array(z=CENTERS, on_boundary=True)

    def update(self, u: core.Array, v: core.Array):
        if not self.tmrlx:
            return

        for boundary in self.open_boundaries.active:
            if boundary.side == Side.EAST:
                inflow = -u.all_values[:, boundary.j, boundary.i - 1]
            elif boundary.side == Side.WEST:
                inflow = u.all_values[:, boundary.j, boundary.i]
            elif boundary.side == Side.NORTH:
                inflow = -v.all_values[:, boundary.j - 1, boundary.i]
            elif boundary.side == Side.SOUTH:
                inflow = v.all_values[:, boundary.j, boundary.i]
            self.inflow.all_values[boundary.start : boundary.stop, :] = inflow.T

        self.rlxcoef.all_values[...] = (self.tmrlx_max - self.tmrlx_min) * np.clip(
            (self.inflow.all_values - self.tmrlx_umin)
            / (self.tmrlx_ucut - self.tmrlx_umin),
            0.0,
            1.0,
        ) + self.tmrlx_min

    def __call__(self, array: core.Array, bdy: core.Array):
        bdy_values = bdy.all_values
        values = array.all_values
        sponge_values = values[..., self.j, self.i]
        if self.tmrlx:
            r = self.rlxcoef.all_values
            sponge_mean = (self.w * sponge_values).sum(axis=-1).T
            bdy_values = r * bdy_values + (1.0 - r) * sponge_mean
        bdy_values = bdy_values.T
        blend = self.sp * bdy_values[..., np.newaxis] + (1.0 - self.sp) * sponge_values
        values[..., self.j, self.i] = blend
        values[..., self.open_boundaries.j, self.open_boundaries.i] = bdy_values


class ZeroGradient(BoundaryCondition):
    def initialize(self, open_boundaries: "OpenBoundaries"):
        tmask = open_boundaries.domain.mask_[1::2, 1::2]
        i = np.empty_like(open_boundaries.i)
        j = np.empty_like(open_boundaries.j)
        for boundary in open_boundaries.active:
            i_inward = {Side.WEST: 1, Side.EAST: -1}.get(boundary.side, 0)
            j_inward = {Side.SOUTH: 1, Side.NORTH: -1}.get(boundary.side, 0)
            i[boundary.start : boundary.stop] = boundary.i + i_inward
            j[boundary.start : boundary.stop] = boundary.j + j_inward
        assert (tmask[j, i] != 0).all(), "Land at boundary interior"
        self.source_slice = (Ellipsis, j, i)
        self.target_slice = (Ellipsis, open_boundaries.j, open_boundaries.i)

    def __call__(self, array: core.Array, bdy: Optional[core.Array] = None):
        array.all_values[self.target_slice] = array.all_values[self.source_slice]


class OpenBoundaries(Mapping):
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
    )

    def __init__(self, domain: "Domain"):
        self.domain = domain
        self._boundaries: List[OpenBoundary] = []
        self.sponge = Sponge()
        self.zero_gradient = ZeroGradient()
        self.active = []
        self._frozen = False

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

        HALO = 2
        nbdyp = 0
        nbdyp_glob = 0
        bdyinfo, bdy_i, bdy_j = [], [], []
        side2count = {}
        self.local_to_global = []
        umask = self.domain.mask_[1::2, 2::2]
        vmask = self.domain.mask_[2::2, 1::2]
        tmask = self.domain.mask_[1::2, 1::2]
        for side in (Side.WEST, Side.NORTH, Side.EAST, Side.SOUTH):
            n = 0
            for boundary in [b for b in self._boundaries if b.side == side]:
                if boundary.l is not None:
                    mskip = boundary.mstart - boundary.mstart_
                    assert mskip >= 0
                    # Note that bdyinfo needs indices into the T grid EXCLUDING halos
                    bdyinfo.append(
                        np.array(
                            (
                                boundary.l - HALO,
                                boundary.mstart - HALO,
                                boundary.mstop - HALO,
                                boundary.type_2d,
                                boundary.type_3d,
                                nbdyp,
                            ),
                            dtype=np.intc,
                        )
                    )
                    len_bdy = boundary.mstop - boundary.mstart
                    boundary.start = nbdyp
                    boundary.stop = nbdyp + len_bdy
                    nbdyp += len_bdy

                    if side in (Side.WEST, Side.EAST):
                        bdy_mask = tmask[boundary.mstart : boundary.mstop, boundary.l]
                        boundary.i = np.repeat(boundary.l, len_bdy)
                        boundary.j = np.arange(boundary.mstart, boundary.mstop)
                    else:
                        bdy_mask = tmask[boundary.l, boundary.mstart : boundary.mstop]
                        boundary.i = np.arange(boundary.mstart, boundary.mstop)
                        boundary.j = np.repeat(boundary.l, len_bdy)

                    if (bdy_mask == 0).any():
                        self.domain.logger.error(
                            "Open boundary %s: %i of %i points of this %sern boundary"
                            " are on land"
                            % (
                                boundary.name,
                                (bdy_mask == 0).sum(),
                                len_bdy,
                                side.name.capitalize(),
                            )
                        )
                        raise Exception()

                    bdy_mask[:] = 2
                    bdy_i.append(boundary.i)
                    bdy_j.append(boundary.j)

                    start_glob = nbdyp_glob + mskip
                    if (
                        self.local_to_global
                        and self.local_to_global[-1][1] == start_glob
                    ):
                        # attach to previous boundary
                        self.local_to_global[-1][1] += len_bdy
                    else:
                        # gap; add new slice
                        self.local_to_global.append([start_glob, start_glob + len_bdy])
                    n += 1
                    self.active.append(boundary)
                nbdyp_glob += boundary.mstop_ - boundary.mstart_
            side2count[side] = n

        umask[:, :-1][(tmask[:, :-1] == 2) & (tmask[:, 1:] == 2)] = 3
        vmask[:-1, :][(tmask[:-1, :] == 2) & (tmask[1:, :] == 2)] = 3

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
                len(bdyinfo),
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
            bdyinfo = np.stack(bdyinfo, axis=-1)
            self.domain.initialize_open_boundaries(
                nwb=side2count[Side.WEST],
                nnb=side2count[Side.NORTH],
                neb=side2count[Side.EAST],
                nsb=side2count[Side.SOUTH],
                nbdyp=self.np,
                bdy_i=self.i - HALO,
                bdy_j=self.j - HALO,
                bdy_info=bdyinfo,
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
        self.z = self.domain.T.array(name="z_bdy", on_boundary=True, register=False)
        self.u = self.domain.T.array(name="u_bdy", on_boundary=True, register=False)
        self.v = self.domain.T.array(name="v_bdy", on_boundary=True, register=False)

        self.sponge.initialize(self)
        self.zero_gradient.initialize(self)

        self._frozen = True

    def __getitem__(self, key) -> OpenBoundary:
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


def find_interfaces(c: ArrayLike) -> np.ndarray:
    c_if = np.empty((c.size + 1),)
    d = np.diff(c)
    c_if[1:-1] = c[:-1] + 0.5 * d
    c_if[0] = c[0] - 0.5 * d[0]
    c_if[-1] = c[-1] + 0.5 * d[-1]
    return c_if


DEG2RAD = np.pi / 180  # degree to radian conversion
R_EARTH = 6378815.0  # radius of the earth (m)
OMEGA = (
    2.0 * np.pi / 86164.0
)  # rotation rate of the earth (rad/s), 86164 is number of seconds in a sidereal day


def coriolis(lat: ArrayLike) -> ArrayLike:
    """Calculate Coriolis parameter f for the given latitude.

    Args:
        lat: latitude in degrees North

    Returns:
        Coriolis parameter f
    """
    return 2.0 * OMEGA * np.sin(DEG2RAD * lat)


def centers_to_supergrid_1d(data: ArrayLike) -> np.ndarray:
    assert data.ndim == 1, "data must be one-dimensional"
    assert data.size > 1, "data must have at least 2 elements"
    data_sup = np.empty((data.size * 2 + 1,))
    data_sup[1::2] = data
    data_sup[2:-2:2] = 0.5 * (data[1:] + data[:-1])
    data_sup[0] = 2 * data_sup[1] - data_sup[2]
    data_sup[-1] = 2 * data_sup[-2] - data_sup[-3]
    return data_sup


def centers_to_supergrid_2d(
    source: ArrayLike,
    ioffset: int,
    joffset: int,
    nx: int,
    ny: int,
    dtype=None,
    missing_values=(),
):
    if dtype is None:
        dtype = source.dtype
    data = np.ma.masked_all(source.shape[:-2] + (ny * 2 + 1, nx * 2 + 1), dtype=dtype)

    # Create an array to hold data at centers (T points),
    # with strips of size 1 on all sides to support interpolation to interfaces
    data_centers = np.ma.masked_all(source.shape[:-2] + (ny + 2, nx + 2), dtype=dtype)

    # Extend the read domain (T grid) by 1 each side, where possible
    # That will allow us to interpolate (rater than extrapolate) to values at the
    # interfaces
    ex_imin = 0 if ioffset == 0 else 1
    ex_imax = 0 if ioffset + nx == source.shape[-1] else 1
    ex_jmin = 0 if joffset == 0 else 1
    ex_jmax = 0 if joffset + ny == source.shape[-2] else 1
    data_centers[
        ..., 1 - ex_jmin : 1 + ny + ex_jmax, 1 - ex_imin : 1 + nx + ex_imax
    ] = source[
        ...,
        joffset - ex_jmin : joffset + ny + ex_jmax,
        ioffset - ex_imin : ioffset + nx + ex_imax,
    ]
    for missing_value in missing_values:
        data_centers = np.ma.masked_equal(data_centers, missing_value, copy=False)

    data_if_ip = np.ma.masked_all(
        (4,) + source.shape[:-2] + (ny + 1, nx + 1), dtype=dtype
    )
    data_if_ip[0, ...] = data_centers[..., :-1, :-1]
    data_if_ip[1, ...] = data_centers[..., 1:, :-1]
    data_if_ip[2, ...] = data_centers[..., :-1, 1:]
    data_if_ip[3, ...] = data_centers[..., 1:, 1:]
    data[..., 1::2, 1::2] = data_centers[1:-1, 1:-1]
    data[..., ::2, ::2] = data_if_ip.mean(axis=0)

    data_if_ip[0, ..., :-1] = data_centers[..., :-1, 1:-1]
    data_if_ip[1, ..., :-1] = data_centers[..., 1:, 1:-1]
    data[..., ::2, 1::2] = data_if_ip[:2, ..., :-1].mean(axis=0)

    data_if_ip[0, ..., :-1, :] = data_centers[..., 1:-1, :-1]
    data_if_ip[1, ..., :-1, :] = data_centers[
        ..., 1:-1, 1:,
    ]
    data[..., 1::2, ::2] = data_if_ip[:2, ..., :-1, :].mean(axis=0)

    return data


def interfaces_to_supergrid_1d(
    data: ArrayLike, out: Optional[np.ndarray] = None
) -> np.ndarray:
    assert data.ndim == 1, "data must be one-dimensional"
    assert data.size > 1, "data must have at least 2 elements"
    if out is None:
        out = np.empty((data.size * 2 - 1,))
    out[0::2] = data
    out[1::2] = 0.5 * (data[1:] + data[:-1])
    return out


def interfaces_to_supergrid_2d(
    data: ArrayLike, out: Optional[np.ndarray] = None
) -> np.ndarray:
    assert data.ndim == 2, "data must be two-dimensional"
    assert data.shape[0] > 1 and data.shape[1] > 1, "dimensions must have length >= 2"
    if out is None:
        out = np.empty((data.shape[0] * 2 - 1, data.shape[1] * 2 - 1))
    out[0::2, 0::2] = data
    out[1::2, 0::2] = 0.5 * (data[:-1, :] + data[1:, :])
    out[0::2, 1::2] = 0.5 * (data[:, :-1] + data[:, 1:])
    out[1::2, 1::2] = 0.25 * (
        data[:-1, :-1] + data[:-1, 1:] + data[1:, :-1] + data[1:, 1:]
    )
    return out


def create_cartesian(
    x: ArrayLike, y: ArrayLike, nz: int, interfaces=False, **kwargs
) -> "Domain":
    """Create Cartesian domain from 1D arrays with x coordinates and  y coordinates.

    Args:
        x: 1d array with x coordinates
            (at cell interfaces if `interfaces=True`, else at cell centers)
        y: 1d array with y coordinates
            (at cell interfaces if `interfaces=True`, else at cell centers)
        nz: number of vertical layers
        interfaces: coordinates are given at cell interfaces, rather than cell centers.
        **kwargs: additional arguments passed to :func:`create`
    """
    assert x.ndim == 1, "x coordinate must be one-dimensional"
    assert y.ndim == 1, "y coordinate must be one-dimensional"

    if interfaces:
        nx, ny = x.size - 1, y.size - 1
        x, y = interfaces_to_supergrid_1d(x), interfaces_to_supergrid_1d(y)
    else:
        nx, ny = x.size, y.size
        x, y = centers_to_supergrid_1d(x), centers_to_supergrid_1d(y)
    return create(nx, ny, nz, x=x, y=y[:, np.newaxis], **kwargs)


def create_spherical(
    lon: ArrayLike, lat: ArrayLike, nz: int, interfaces=False, **kwargs
) -> "Domain":
    """Create spherical domain from 1D arrays with longitudes and latitudes.

    Args:
        lon: 1d array with longitude coordinates
            (at cell interfaces if `interfaces=True`, else at cell centers)
        lat: 1d array with latitude coordinates
            (at cell interfaces if `interfaces=True`, else at cell centers)
        nz: number of vertical layers
        interfaces: coordinates are given at cell interfaces, rather than cell centers.
        **kwargs: additional arguments passed to :func:`create`
    """
    assert lon.ndim == 1, "longitude coordinate must be one-dimensional"
    assert lat.ndim == 1, "latitude coordinate must be one-dimensional"

    if interfaces:
        nx, ny = lon.size - 1, lat.size - 1
        lon, lat = interfaces_to_supergrid_1d(lon), interfaces_to_supergrid_1d(lat)
    else:
        nx, ny = lon.size, lat.size
        lon, lat = centers_to_supergrid_1d(lon), centers_to_supergrid_1d(lat)
    return create(nx, ny, nz, lon=lon, lat=lat[:, np.newaxis], spherical=True, **kwargs)


def create_spherical_at_resolution(
    minlon: float,
    maxlon: float,
    minlat: float,
    maxlat: float,
    resolution: float,
    nz: int,
    **kwargs
) -> "Domain":
    """Create spherical domain encompassing the specified longitude range and latitude
    range and desired resolution in m.

    Args:
        minlon: minimum longitude
        maxlon: maximum longitude
        minlat: minimum latitude
        maxlat: maximum latitude
        resolution: maximum grid cell length and width (m)
        nz: number of vertical layers
        **kwargs: additional arguments passed to :func:`create`
    """
    assert maxlon > minlon, (
        "Maximum longitude %s must be greater than minimum longitude %s"
        % (maxlon, minlon)
    )
    assert (
        maxlat > minlat
    ), "Maximum latitude %s must be greater than minimum latitude %s" % (maxlat, minlat)
    assert resolution > 0, "Desired resolution must be greater than 0, but is %s m" % (
        resolution,
    )
    dlat = resolution / (DEG2RAD * R_EARTH)
    minabslat = min(abs(minlat), abs(maxlat))
    dlon = resolution / (DEG2RAD * R_EARTH) / np.cos(DEG2RAD * minabslat)
    nx = int(np.ceil((maxlon - minlon) / dlon)) + 1
    ny = int(np.ceil((maxlat - minlat) / dlat)) + 1
    return create_spherical(
        np.linspace(minlon, maxlon, nx),
        np.linspace(minlat, maxlat, ny),
        nz=nz,
        interfaces=True,
        **kwargs
    )


def create(
    nx: int,
    ny: int,
    nz: int,
    lon: Optional[np.ndarray] = None,
    lat: Optional[np.ndarray] = None,
    x: Optional[np.ndarray] = None,
    y: Optional[np.ndarray] = None,
    spherical: bool = False,
    mask: Optional[np.ndarray] = 1,
    H: Optional[np.ndarray] = None,
    z0: Optional[np.ndarray] = 0.0,
    f: Optional[np.ndarray] = None,
    tiling: Optional[parallel.Tiling] = None,
    periodic_x: bool = False,
    periodic_y: bool = False,
    logger: Optional[logging.Logger] = None,
    glob: bool = False,
    **kwargs
) -> "Domain":
    """Create new domain. By default, the local subdomain is returned if
    multiple processing cores are active (and `glob` is `False`)

    Args:
        nx: number of tracer points in x-direction
        ny: number of tracer points in y-direction
        nz: number of vertical layers
        lon: longitude (degrees East)
        lat: latitude (degrees North)
        x: x coordinate (m)
        y: y coordinate (m)
        spherical: grid is spherical (as opposed to Cartesian). If True, at least
            ``lon`` and ``lat`` must be provided. Otherwise at least ``x`` and ``y``
            must be provided.
        mask: initial mask (0: land, 1: water)
        H: initial bathymetric depth. This is the distance between the bottom and
            some arbitrary depth reference (m, positive if bottom lies below the
            depth reference). Typically the depth reference is mean sea level.
        z0: initial bottom roughness (m)
        f: Coriolis parameter. By default this is calculated from latitude ``lat``
            if provided.
        tiling: subdomain decomposition
        periodic_x: use periodic boundary in x-direction (left == right)
            This is only used if `tiling` is not provided.
        periodic_y: use periodic boundary in y-direction (top == bottom)
            This is only used if `tiling` is not provided.
        logger: target for log messages
        glob: return the global domain rather than the local subdomain
        **kwargs: additional keyword arguments passed to :class:`Domain`
    """
    global_domain = None
    logger = logger or parallel.get_logger()
    parlogger = logger.getChild("parallel")

    # Determine subdomain division

    if glob:
        # We need to return the global domain. Create a dummy subdomain layout.
        tiling = parallel.Tiling(
            nrow=1, ncol=1, ncpus=1, periodic_x=periodic_x, periodic_y=periodic_y
        )
        tiling.rank = 0

    if tiling is None:
        # No tiling provided - autodetect
        mask = np.broadcast_to(mask, (1 + 2 * ny, 1 + 2 * nx))
        tiling = parallel.Tiling.autodetect(
            mask=mask[1::2, 1::2],
            logger=parlogger,
            periodic_x=periodic_x,
            periodic_y=periodic_y,
        )
    elif isinstance(tiling, str):
        # Path to dumped Tiling object provided
        if not os.path.isfile(tiling):
            logger.critical(
                "Cannot find file %s. If tiling is a string, it must be the path to"
                " an existing file with a pickled tiling object." % tiling
            )
            raise Exception()
        tiling = parallel.Tiling.load(tiling)
    else:
        # Existing tiling object provided
        # Transfer extent of global domain to determine subdomain sizes
        if isinstance(tiling, tuple):
            tiling = parallel.Tiling(
                nrow=tiling[0],
                ncol=tiling[1],
                periodic_x=periodic_x,
                periodic_y=periodic_y,
            )
        tiling.set_extent(nx, ny)
    tiling.report(parlogger)

    global_tiling = tiling
    if tiling.n > 1:
        # The global tiling object is a simple 1x1 partition
        global_tiling = parallel.Tiling(
            nrow=1, ncol=1, ncpus=1, periodic_x=periodic_x, periodic_y=periodic_y
        )
        global_tiling.set_extent(nx, ny)

    # If on master node (possibly only node), create global domain object
    if tiling.rank == 0:
        global_domain = Domain(
            nx,
            ny,
            nz,
            lon,
            lat,
            x,
            y,
            spherical,
            tiling=global_tiling,
            mask=mask,
            H=H,
            z0=z0,
            f=f,
            logger=logger,
            **kwargs
        )

    # If there is only one node, return the global domain immediately
    if tiling.n == 1:
        return global_domain

    # Create the subdomain, and (if on root) attach a pointer to the global domain
    subdomain = Domain.partition(
        tiling,
        nx,
        ny,
        nz,
        global_domain,
        has_xy=x is not None,
        has_lonlat=lon is not None,
        spherical=spherical,
        logger=logger,
        **kwargs
    )
    subdomain.glob = global_domain

    return subdomain


def load(path: str, nz: int, **kwargs) -> "Domain":
    """Load domain from file. Typically this is a file created by :meth:`Domain.save`.

    Args:
        path: NetCDF file to load from
        nz: number of vertical layers
        **kwargs: additional keyword arguments to pass to :class:`Domain`
    """
    with netCDF4.Dataset(path) as nc:
        for name in ("lon", "lat", "x", "y", "H", "mask", "z0b_min", "cor"):
            if name in nc.variables and name not in kwargs:
                kwargs[name] = nc.variables[name][...]
    spherical = kwargs.get("spherical", "lon" in kwargs)
    ny, nx = kwargs["lon" if spherical else "x"].shape
    nx = (nx - 1) // 2
    ny = (ny - 1) // 2
    return create(nx, ny, nz, spherical=spherical, **kwargs)


class Domain(_pygetm.Domain):
    @staticmethod
    def partition(
        tiling: parallel.Tiling,
        nx: int,
        ny: int,
        nz: int,
        global_domain: Optional["Domain"],
        halo: int = 2,
        has_xy: bool = True,
        has_lonlat: bool = True,
        logger: Optional[logging.Logger] = None,
        **kwargs
    ):
        assert nx == tiling.nx_glob and ny == tiling.ny_glob, (
            "Extent of global domain (%i, %i) does not match that of tiling (%i, %i)."
            % (ny, nx, tiling.ny_glob, tiling.nx_glob)
        )
        assert tiling.n == tiling.comm.Get_size(), (
            "Number of active cores in subdomain decompositon (%i) does not match "
            "available number of cores (%i)." % (tiling.n, tiling.comm.Get_size())
        )

        # supergrid metrics: double extent/halos, one point overlap between subdomains
        SCALE = 2
        HALO = SCALE * halo
        SHARE = 1

        local_slice, _, _, _ = tiling.subdomain2slices(
            halo_sub=HALO,
            halo_glob=HALO,
            scale=SCALE,
            share=SHARE,
            exclude_halos=False,
            exclude_global_halos=True,
        )

        coordinates = {"f": "cor"}
        if has_xy:
            coordinates["x"] = "x"
            coordinates["y"] = "y"
        if has_lonlat:
            coordinates["lon"] = "lon"
            coordinates["lat"] = "lat"
        for name, att in coordinates.items():
            c = np.empty(
                (
                    SCALE * tiling.ny_sub + 2 * HALO + SHARE,
                    SCALE * tiling.nx_sub + 2 * HALO + SHARE,
                )
            )
            scatterer = parallel.Scatter(
                tiling, c, halo=HALO, share=SHARE, scale=SCALE, fill_value=np.nan
            )
            scatterer(
                None if global_domain is None else getattr(global_domain, att + "_")
            )
            assert not np.isnan(
                c[local_slice]
            ).any(), "Subdomain %s contains NaN after initial scatter"
            kwargs[name] = c

        domain = Domain(
            tiling.nx_sub, tiling.ny_sub, nz, tiling=tiling, logger=logger, **kwargs
        )

        parallel.Scatter(tiling, domain.mask_, halo=HALO, share=SHARE, scale=SCALE)(
            None if global_domain is None else global_domain.mask_
        )
        parallel.Scatter(tiling, domain.H_, halo=HALO, share=SHARE, scale=SCALE)(
            None if global_domain is None else global_domain.H_
        )
        parallel.Scatter(tiling, domain.z0b_min_, halo=HALO, share=SHARE, scale=SCALE)(
            None if global_domain is None else global_domain.z0b_min_
        )

        return domain

    def _exchange_metric(
        self, data, relative_in_x: bool = False, relative_in_y: bool = False
    ):
        if not self.tiling:
            return

        expected_shape = (
            1 + 2 * (self.ny + 2 * self.haloy),
            1 + 2 * (self.nx + 2 * self.halox),
        )
        assert data.shape == expected_shape, "Wrong shape: got %s, expected %s." % (
            data.shape,
            expected_shape,
        )

        HALO = 4  # supergrid
        fill_value = 0 if np.issubdtype(data.dtype, np.integer) else FILL_VALUE
        valid_before = data != fill_value
        assert np.isfinite(data).all(), str(data)

        # Expand the data array one each side
        shape_ext = (data.shape[0] + 2, data.shape[1] + 2)
        data_ext = np.full(shape_ext, fill_value, dtype=data.dtype)
        data_ext[1 + HALO : -1 - HALO, 1 + HALO : -1 - HALO] = data[
            HALO:-HALO, HALO:-HALO
        ]
        if relative_in_x or relative_in_y:
            # Pre-fill the halo zones with existing values
            # This is needed in cases where some neighbors are missing
            data_ext[: HALO + 1, : HALO + 1] = data[: HALO + 1, : HALO + 1]
            data_ext[: HALO + 1, HALO + 1 : -HALO - 1] = data[: HALO + 1, HALO:-HALO]
            data_ext[: HALO + 1, -HALO - 1 :] = data[: HALO + 1, -HALO - 1 :]
            data_ext[HALO + 1 : -HALO - 1, : HALO + 1] = data[HALO:-HALO, : HALO + 1]
            data_ext[HALO + 1 : -HALO - 1, -HALO - 1 :] = data[HALO:-HALO, -HALO - 1 :]
            data_ext[-HALO - 1 :, : HALO + 1] = data[-HALO - 1 :, : HALO + 1]
            data_ext[-HALO - 1 :, HALO + 1 : -HALO - 1] = data[-HALO - 1 :, HALO:-HALO]
            data_ext[-HALO - 1 :, -HALO - 1 :] = data[-HALO - 1 :, -HALO - 1 :]
        self.tiling.wrap(data_ext, HALO + 1).update_halos()

        # For values in the halo, compute their difference with the outer boundary of
        # the subdomain we exchanged with (now the innermost halo point). Then use that
        # difference plus the value on our own boundary as values inside the halo.
        # This ensures coordinate variables are monotonically increasing in interior
        # AND halos, even if periodic boundary conditions are used.
        if relative_in_x:
            data_ext[:, : HALO + 1] += (
                data_ext[:, HALO + 1 : HALO + 2] - data_ext[:, HALO : HALO + 1]
            )
            data_ext[:, -HALO - 1 :] += (
                data_ext[:, -HALO - 2 : -HALO - 1] - data_ext[:, -HALO - 1 : -HALO]
            )
        if relative_in_y:
            data_ext[: HALO + 1, :] += (
                data_ext[HALO + 1 : HALO + 2, :] - data_ext[HALO : HALO + 1, :]
            )
            data_ext[-HALO - 1 :, :] += (
                data_ext[-HALO - 2 : -HALO - 1, :] - data_ext[-HALO - 1 : -HALO, :]
            )

        # Since subdomains share the outer boundary, that boundary will be replicated
        # in the outermost interior point and in the innermost halo point
        # We move the outer part of the halos (all but their innermost point) one point
        # inwards to eliminate that overlapping point
        # Where we do not have a subdomain neighbor, we keep the original values.
        if self.tiling.bottomleft != -1:
            data[:HALO, :HALO] = data_ext[:HALO, :HALO]
        if self.tiling.bottom != -1:
            data[:HALO, HALO:-HALO] = data_ext[:HALO, HALO + 1 : -HALO - 1]
        if self.tiling.bottomright != -1:
            data[:HALO, -HALO:] = data_ext[:HALO, -HALO:]
        if self.tiling.left != -1:
            data[HALO:-HALO, :HALO] = data_ext[HALO + 1 : -HALO - 1, :HALO]
        if self.tiling.right != -1:
            data[HALO:-HALO, -HALO:] = data_ext[HALO + 1 : -HALO - 1, -HALO:]
        if self.tiling.topleft != -1:
            data[-HALO:, :HALO] = data_ext[-HALO:, :HALO]
        if self.tiling.top != -1:
            data[-HALO:, HALO:-HALO] = data_ext[-HALO:, HALO + 1 : -HALO - 1]
        if self.tiling.topright != -1:
            data[-HALO:, -HALO:] = data_ext[-HALO:, -HALO:]

        # Values in halos where there is no matching neighbor should have been preserved
        # Therefore we cannot have gained invalid values anywhere - verify this.
        assert np.isfinite(data).all()
        valid_after = True if isinstance(fill_value, int) else data != fill_value
        still_ok = np.where(valid_before, valid_after, True)
        assert still_ok.all(), "Rank %i: _exchange_metric corrupted %i values: %s." % (
            self.tiling.rank,
            still_ok.size - still_ok.sum(),
            still_ok,
        )

    def _map_array(self, source: ArrayLike, target: np.ndarray):
        supergrid_shape = (
            (self.ny + 2 * self.haloy) * 2 + 1,
            (self.nx + 2 * self.halox) * 2 + 1,
        )
        assert target.shape[-2:] == supergrid_shape
        nx_glob, ny_glob = self.tiling.nx_glob, self.tiling.ny_glob
        source_shape = np.shape(source)

        target_slice, source_slice, _, _ = self.tiling.subdomain2slices(
            halo_sub=4,
            halo_glob=0,
            scale=2,
            share=1,
            exclude_halos=False,
            exclude_global_halos=True,
        )

        def cast(target_shape: Tuple[int]) -> bool:
            nonlocal source
            try:
                np.broadcast_shapes(source_shape, target_shape)
            except ValueError:
                return False
            source = np.broadcast_to(source, target_shape)
            return True

        # Note: the first of these clauses will catch assignment to a scalar
        # - a very common case!
        if cast((self.ny * 2 + 1, self.nx * 2 + 1)):
            # local domain, supergrid EXcluding halos
            target_slice = (Ellipsis, slice(4, -4), slice(4, -4))
            source_slice = (Ellipsis,)
        elif cast(target.shape):
            # local domain, supergrid INcluding halos
            target_slice = (Ellipsis,)
            source_slice = (Ellipsis,)
        elif cast((self.ny, self.nx)):
            # local domain, T grid, no halos
            source = centers_to_supergrid_2d(source, 0, 0, self.nx, self.ny)
            target_slice = (Ellipsis, slice(4, -4), slice(4, -4))
            source_slice = (Ellipsis,)
        elif cast((self.ny + 1, self.nx + 1)):
            # local domain, X grid, no halos
            source = interfaces_to_supergrid_2d(source)
            target_slice = (Ellipsis, slice(4, -4), slice(4, -4))
            source_slice = (Ellipsis,)
        elif cast((ny_glob, nx_glob)):
            # global domain, T grid
            source = centers_to_supergrid_2d(source, 0, 0, nx_glob, ny_glob)
        elif cast((ny_glob + 1, nx_glob + 1)):
            # global domain, X grid
            source = interfaces_to_supergrid_2d(source)
        else:
            raise Exception(
                "Cannot map array with shape %s to local supergrid with shape %s."
                % (source_shape, target.shape)
            )

        target[target_slice] = np.nan_to_num(source[source_slice], nan=FILL_VALUE)

    def __init__(
        self,
        nx: int,
        ny: int,
        nz: int,
        lon: Optional[np.ndarray] = None,
        lat: Optional[np.ndarray] = None,
        x: Optional[np.ndarray] = None,
        y: Optional[np.ndarray] = None,
        spherical: bool = False,
        mask: Optional[np.ndarray] = 1,
        H: Optional[np.ndarray] = None,
        z0: Optional[np.ndarray] = 0.0,
        f: Optional[np.ndarray] = None,
        tiling: Optional[parallel.Tiling] = None,
        logger: Optional[logging.Logger] = None,
        Dmin: float = 1.0,
        Dcrit: float = 2.0,
        vertical_coordinate_method: VerticalCoordinates = VerticalCoordinates.SIGMA,
        ddl: float = 0.0,
        ddu: float = 0.0,
        Dgamma: float = 0.0,
        gamma_surf: bool = True,
        **kwargs
    ):
        """Create domain with coordinates, bathymetry, mask defined on the supergrid.

        Args:
            nx: number of tracer points in x-direction
            ny: number of tracer points in y-direction
            nz: number of vertical layers
            lon: longitude (degrees East)
            lat: latitude (degrees North)
            x: x coordinate (m)
            y: y coordinate (m)
            spherical: grid is spherical (as opposed to Cartesian). If True, at least
                ``lon`` and ``lat`` must be provided. Otherwise at least ``x`` and ``y``
                must be provided.
            mask: initial mask (0: land, 1: water)
            H: initial bathymetric depth. This is the distance between the bottom and
                some arbitrary depth reference (m, positive if bottom lies below the
                depth reference). Typically the depth reference is mean sea level.
            z0: initial bottom roughness (m)
            f: Coriolis parameter. By default this is calculated from latitude ``lat``
                if provided.
            tiling: subdomain decomposition
            Dmin: minimum depth (m) for wet points. At this depth, all hydrodynamic
                terms except the pressure gradient and bottom friction are switched off.
            Dcrit: depth (m) at which tapering of processes (all except pressure
                gradient and bottom friction) begins.
            vertical_coordinate_method: type of vertical coordinate to use
            ddl: dimensionless factor for zooming towards the bottom (0: no zooming,
                > 2: strong zooming)
            ddl: dimensionless factor for zooming towards the surface (0: no zooming,
                > 2: strong zooming)
            Dgamma: depth (m) range over which z-like coordinates should be used
            gamma_surf: use z-like coordinates in surface layer (as opposed to bottom
                layer)
            **kwargs: additional keyword arguments passed to
                :class:`pygetm.parallel.Tiling`
        """
        if nx <= 0:
            raise Exception("Number of x points is %i but must be > 0" % nx)
        if ny <= 0:
            raise Exception("Number of y points is %i but must be > 0" % ny)
        if nz <= 0:
            raise Exception("Number of z points is %i but must be > 0" % nz)
        if lat is None and f is None:
            raise Exception(
                "Either lat of f must be provided to determine the Coriolis parameter."
            )

        # Loggers
        if logger is None:
            logger = parallel.get_logger()
        self.root_logger: logging.Logger = logger
        self.logger: logging.Logger = self.root_logger.getChild("domain")

        #: collection of all model fields
        self.fields: MutableMapping[str, core.Array] = {}

        #: input manager responsible for reading from NetCDF and for spatial and
        # temporal interpolation
        self.input_manager: input.InputManager = input.InputManager()

        #: global model domain, set on root (rank=0) only
        self.glob: Optional["Domain"] = self

        self.logger.info("Domain size (T grid): %i x %i (%i cells)" % (nx, ny, nx * ny))

        self.imin, self.imax = 1, nx
        self.jmin, self.jmax = 1, ny
        self.kmin, self.kmax = 1, nz

        super().__init__(
            self.imin, self.imax, self.jmin, self.jmax, self.kmin, self.kmax
        )

        halo = 2

        shape = (2 * ny + 1, 2 * nx + 1)
        superhalo = 2 * halo
        shape_ = (shape[0] + 2 * superhalo, shape[1] + 2 * superhalo)

        # Set up subdomain partition information to enable halo exchanges
        if tiling is None:
            tiling = parallel.Tiling(**kwargs)
        elif kwargs:
            raise Exception("Encountered unexpected arguments: %s" % ", ".join(kwargs))
        self.tiling = tiling

        def setup_metric(
            source: Optional[np.ndarray] = None,
            optional: bool = False,
            fill_value=FILL_VALUE,
            relative_in_x: bool = False,
            relative_in_y: bool = False,
            dtype: np.typing.DTypeLike = float,
            writeable: bool = True,
        ) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
            if optional and source is None:
                return None, None
            data = np.full(shape_, fill_value, dtype)
            data_int = data[superhalo:-superhalo, superhalo:-superhalo]
            if source is not None:
                self._map_array(source, data)
                self._exchange_metric(data, relative_in_x, relative_in_y)
            data.flags.writeable = data_int.flags.writeable = writeable
            return data_int, data

        # Supergrid metrics:
        # without underscore=interior only, with underscore=including halos
        self.x, self.x_ = setup_metric(
            x, optional=True, relative_in_x=True, writeable=False
        )
        self.y, self.y_ = setup_metric(
            y, optional=True, relative_in_y=True, writeable=False
        )
        self.lon, self.lon_ = setup_metric(
            lon, optional=True, relative_in_x=True, writeable=False
        )
        self.lat, self.lat_ = setup_metric(
            lat, optional=True, relative_in_y=True, writeable=False
        )
        self.H, self.H_ = setup_metric(H)
        self.z0b_min, self.z0b_min_ = setup_metric(z0)
        self.mask, self.mask_ = setup_metric(mask, dtype=np.intc, fill_value=0)

        cor = f if f is not None else coriolis(lat)
        self.cor, self.cor_ = setup_metric(cor, writeable=False)

        # Compute dx, dy from Cartesian or spherical coordinates
        # These have had their halo exchanges as part of setup_metric and are therefore
        # defined inside the halos too.
        self.dx, self.dx_ = setup_metric()
        self.dy, self.dy_ = setup_metric()
        if spherical:
            dlon_x = self.lon_[:, 2:] - self.lon_[:, :-2]
            dlat_x = self.lat_[:, 2:] - self.lat_[:, :-2]
            dx_x = DEG2RAD * dlon_x * R_EARTH * np.cos(DEG2RAD * self.lat_[:, 1:-1])
            dy_x = DEG2RAD * dlat_x * R_EARTH
            dlon_y = self.lon_[2:, :] - self.lon_[:-2, :]
            dlat_y = self.lat_[2:, :] - self.lat_[:-2, :]
            dx_y = DEG2RAD * dlon_y * R_EARTH * np.cos(DEG2RAD * self.lat_[1:-1, :])
            dy_y = DEG2RAD * dlat_y * R_EARTH
        else:
            dx_x = self.x_[:, 2:] - self.x_[:, :-2]
            dy_x = self.y_[:, 2:] - self.y_[:, :-2]
            dx_y = self.x_[2:, :] - self.x_[:-2, :]
            dy_y = self.y_[2:, :] - self.y_[:-2, :]
        self.dx_[:, 1:-1] = np.hypot(dx_x, dy_x)
        self.dy_[1:-1, :] = np.hypot(dx_y, dy_y)

        # Halo exchange for dx, dy, needed to ensure the outer strips of the halos are
        # valid. Those outermost strips could not be computed by central-differencing
        # the coordinates as that would require points outside the domain.
        self._exchange_metric(self.dx_)
        self._exchange_metric(self.dy_)

        self.dx_.flags.writeable = self.dx.flags.writeable = False
        self.dy_.flags.writeable = self.dy.flags.writeable = False

        self.rotation, self.rotation_ = setup_metric()

        def supergrid_rotation(x, y):
            # For each point, draw lines to the nearest neighbor (1/2 a grid cell) on
            # the left, right, top and bottom.
            rotation_left = np.arctan2(y[:, 1:-1] - y[:, :-2], x[:, 1:-1] - x[:, :-2])
            rotation_right = np.arctan2(y[:, 2:] - y[:, 1:-1], x[:, 2:] - x[:, 1:-1])
            rotation_bot = (
                np.arctan2(y[1:-1, :] - y[:-2, :], x[1:-1, :] - x[:-2, :]) - 0.5 * np.pi
            )
            rotation_top = (
                np.arctan2(y[2:, :] - y[1:-1, :], x[2:, :] - x[1:-1, :]) - 0.5 * np.pi
            )
            x_dum = (
                np.cos(rotation_left[1:-1, :])
                + np.cos(rotation_right[1:-1, :])
                + np.cos(rotation_bot[:, 1:-1])
                + np.cos(rotation_top[:, 1:-1])
            )
            y_dum = (
                np.sin(rotation_left[1:-1, :])
                + np.sin(rotation_right[1:-1, :])
                + np.sin(rotation_bot[:, 1:-1])
                + np.sin(rotation_top[:, 1:-1])
            )
            return np.arctan2(y_dum, x_dum)

        if self.lon_ is not None and self.lat_ is not None:
            # Proper rotation with respect to true North
            self.rotation_[1:-1, 1:-1] = supergrid_rotation(self.lon_, self.lat_)
        else:
            # Rotation with respect to y axis - assumes y axis always point to true
            # North (can be valid only for infinitesimally small domain)
            self.rotation_[1:-1, 1:-1] = supergrid_rotation(self.x_, self.y_)
        self._exchange_metric(self.rotation_)
        self.rotation.flags.writeable = False

        self.area, self.area_ = setup_metric(self.dx * self.dy, writeable=False)

        self.spherical = spherical

        # Determine if we have simple coordinates (useful for xarray and plotting in
        # general)
        self.lon_is_1d = (
            self.lon is not None and (self.lon[:1, :] == self.lon[:, :]).all()
        )
        self.lat_is_1d = (
            self.lat is not None and (self.lat[:, :1] == self.lat[:, :]).all()
        )
        self.x_is_1d = self.x is not None and (self.x[:1, :] == self.x[:, :]).all()
        self.y_is_1d = self.y is not None and (self.y[:, :1] == self.y[:, :]).all()

        self.halo = self.halox
        self.shape = (nz, ny + 2 * self.halo, nx + 2 * self.halo)

        # Advection grids (two letters: first for advected quantity, second for
        # advection direction)
        self.UU = Grid(self, _pygetm.UUGRID, ioffset=3, joffset=1)
        self.VV = Grid(self, _pygetm.VVGRID, ioffset=1, joffset=3)
        self.UV = Grid(self, _pygetm.UVGRID, ioffset=2, joffset=2)
        self.VU = Grid(self, _pygetm.VUGRID, ioffset=2, joffset=2)

        # Create grids
        self.U = Grid(
            self, _pygetm.UGRID, ioffset=2, joffset=1, ugrid=self.UU, vgrid=self.UV
        )
        self.V = Grid(
            self, _pygetm.VGRID, ioffset=1, joffset=2, ugrid=self.VU, vgrid=self.VV
        )
        self.T = Grid(
            self, _pygetm.TGRID, ioffset=1, joffset=1, ugrid=self.U, vgrid=self.V
        )
        self.X = Grid(self, _pygetm.XGRID, ioffset=0, joffset=0, overlap=1)

        self.Dmin = Dmin
        self.Dcrit = Dcrit
        self.ddl = ddl
        self.ddu = ddu
        self.Dgamma = Dgamma
        self.gamma_surf = gamma_surf
        self.vertical_coordinate_method = vertical_coordinate_method

        self._initialized = False
        self.open_boundaries = OpenBoundaries(self)
        self.rivers = Rivers(self.T)

    def cfl_check(self, z: float = 0.0, log: bool = True) -> float:
        """Determine maximum time step for depth-integrated equations

        Args:
            z: surface elevation (m) at rest
            log: whether to write the maximum time step and its location to the log
        """
        mask = self.mask[1::2, 1::2] > 0
        dx = self.dx[1::2, 1::2]
        dy = self.dy[1::2, 1::2]
        H = self.H[1::2, 1::2]
        maxdts = (
            dx
            * dy
            / np.sqrt(
                2.0 * GRAVITY * (H + z) * (dx ** 2 + dy ** 2),
                where=mask,
                out=np.ones_like(H),
            )
        )
        maxdts[~mask] = np.inf
        maxdt = maxdts.min()
        if log:
            j, i = np.unravel_index(np.argmin(maxdts), maxdts.shape)
            self.logger.info(
                "Maximum dt = %.3f s (i=%i, j=%i, bathymetric depth=%.3f m)"
                % (maxdt, i, j, H[j, i])
            )
        return maxdt

    def initialize(self, runtype: int):
        """Initialize the domain. This updates the mask in order for it to be
        consistent across T, U, V, X grids. Values for the mask, bathymetry, and
        bottom roughness are subsequently read-only.
        """
        assert not self._initialized, "Domain has already been initialized"
        if self.glob is not None and self.glob is not self:
            self.glob.initialize(runtype)

        # The user could modify mask, bathymetry, bottom roughness freely until now.
        # Ensure they are marked invalid inside halos and outside the global domain.
        local_slice, _, _, _ = self.tiling.subdomain2slices(
            halo_sub=4, scale=2, share=1, exclude_halos=True, exclude_global_halos=True,
        )
        outside = np.full(self.H_.shape, True)
        outside[local_slice] = False
        self.mask_[outside] = 0

        # Now update halos.
        # This ensures the mask in halos is valid only if there is an actual neighbor.
        self._exchange_metric(self.mask_)
        self._exchange_metric(self.H_)
        self._exchange_metric(self.z0b_min_)

        # Mask U,V,X points without any valid T neighbor - this mask will be maintained
        # by the domain to be used for e.g. plotting
        tmask = self.mask_[1::2, 1::2]
        umask = self.mask_[1::2, 2::2]
        vmask = self.mask_[2::2, 1::2]
        xmask = self.mask_[::2, ::2]

        backup = self.mask_.copy()
        umask[:, :-1][(tmask[:, 1:] != 0) | (tmask[:, :-1] != 0)] = 1
        vmask[:-1, :][(tmask[1:, :] != 0) | (tmask[:-1, :] != 0)] = 1
        xmask[1:-1, 1:-1][
            (tmask[1:, 1:] != 0)
            | (tmask[:-1, 1:] != 0)
            | (tmask[1:, :-1] != 0)
            | (tmask[:-1, :-1] != 0)
        ] = 0
        self._exchange_metric(self.mask_)

        self.T._water_contact = tmask != 0
        self.U._water_contact = umask != 0
        self.V._water_contact = vmask != 0
        self.X._water_contact = xmask != 0
        self.mask_[...] = backup

        self.open_boundaries.initialize()

        # At this point we could make a copy of self.mask_ that would be useful for
        # plotting, since it does not hide U, V, X points at the edges of valid T cells.
        # Now mask U,V,X points unless all their T neighbors are valid - this mask will
        # be sent to Fortran and determine which points are computed
        umask[:, :-1][(tmask[:, 1:] == 0) | (tmask[:, :-1] == 0)] = 0
        vmask[:-1, :][(tmask[1:, :] == 0) | (tmask[:-1, :] == 0)] = 0
        xmask[1:-1, 1:-1][
            (tmask[1:, 1:] == 0)
            | (tmask[:-1, 1:] == 0)
            | (tmask[1:, :-1] == 0)
            | (tmask[:-1, :-1] == 0)
        ] = 0
        self._exchange_metric(self.mask_)

        uumask = self.UU._setup_array("mask", from_supergrid=False).all_values
        uvmask = self.UV._setup_array("mask", from_supergrid=False).all_values
        vumask = self.VU._setup_array("mask", from_supergrid=False).all_values
        vvmask = self.VV._setup_array("mask", from_supergrid=False).all_values
        uumask[:, :-1][(umask[:, :-1] > 0) & (umask[:, 1:] > 0)] = 1
        uvmask[:-1, :][(umask[:-1, :] > 0) & (umask[1:, :] > 0)] = 1
        vumask[:, :-1][(vmask[:, :-1] > 0) & (vmask[:, 1:] > 0)] = 1
        vvmask[:-1, :][(vmask[:-1, :] > 0) & (vmask[1:, :] > 0)] = 1

        for grid in self.grids.values():
            grid.initialize(self.open_boundaries.np)

        self.logger.info(
            "Number of unmasked points excluding halos: %i on T grid, %i on U grid,"
            " %i on V grid, %i on X grid"
            % (
                (self.T.mask.values > 0).sum(),
                (self.U.mask.values > 0).sum(),
                (self.V.mask.values > 0).sum(),
                (self.X.mask.values > 0).sum(),
            )
        )

        # Water depth and thicknesses on UU/VV grids will be taken from T grid,
        # which near land has valid values where UU/VV are masked
        self.UU.D.attrs["_mask_output"] = True
        self.VV.D.attrs["_mask_output"] = True
        self.UU.hn.attrs["_mask_output"] = True
        self.VV.hn.attrs["_mask_output"] = True

        # Elevation on U/V/X will be interpolated from the T grid and therefore
        # contain values other than fill_value in masked points
        for grid in (self.U, self.V, self.X):
            grid.z.attrs["_mask_output"] = True
            grid.zio.attrs["_mask_output"] = True
            grid.zin.attrs["_mask_output"] = True

        self.H_.flags.writeable = self.H.flags.writeable = False
        self.z0b_min_.flags.writeable = self.z0b_min.flags.writeable = False
        self.mask_.flags.writeable = self.mask.flags.writeable = False

        super().initialize(
            runtype,
            Dmin=self.Dmin,
            method_vertical_coordinates=self.vertical_coordinate_method,
            ddl=self.ddl,
            ddu=self.ddu,
            Dgamma=self.Dgamma,
            gamma_surf=self.gamma_surf,
        )

        self.rivers.initialize()

        # Water depth and thicknesses on T grid that lag 1/2 time step behind tracer
        # (i.e., they are in sync with U, V, X grids)
        self.D_T_half = self.T.array(fill=np.nan)
        self.h_T_half = self.T.array(fill=np.nan, z=CENTERS)
        self.depth = self.T.array(
            z=CENTERS,
            name="pres",
            units="dbar",
            long_name="pressure",
            fabm_standard_name="depth",
            fill_value=FILL_VALUE,
        )

        self.cfl_check()

        self._initialized = True

    def set_bathymetry(
        self, depth: ArrayLike, scale_factor=None, periodic_lon: bool = False
    ):
        """Set bathymetric depth on supergrid. The bathymetric depth is the distance
        between some arbitrary depth reference (often mean sea level) and the bottom,
        positive for greater depth.

        Args:
            depth: depth values for the global domain or the local subdomain.
                They can be provided on the supergrid, or on the T or X grid.
            scale_factor: apply scale factor to provided depths
            periodic_lon: depth source spans entire globe in longitude
        """
        assert (
            not self._initialized
        ), "set_bathymetry cannot be called after the domain has been initialized."

        if not isinstance(depth, xarray.DataArray):
            self._map_array(depth, self.H_)
        else:
            # Depth is provided as xarray object that includes coordinates (we require
            # CF compliant longitude, latitude). Interpolate to target grid.
            depth = input.limit_region(
                depth,
                self.lon.min(),
                self.lon.max(),
                self.lat.min(),
                self.lat.max(),
                periodic_lon=periodic_lon,
            )
            depth = input.horizontal_interpolation(depth, self.lon, self.lat)
            self.H[...] = depth.values
        if scale_factor is not None:
            self.H *= scale_factor

    def mask_shallow(self, minimum_depth: float):
        """Mask all points shallower less the specified value.

        Args:
            minimum_depth: minimum bathmetric depth :attr:`H`; points that are
                shallower will be masked
        """
        self.mask[self.H < minimum_depth] = 0

    def limit_velocity_depth(self, critical_depth: Optional[float] = None):
        """Decrease bathymetric depth of velocity (U, V) points to the minimum of the
        bathymetric depth of both neighboring T points, wherever one of these two
        points is shallower than the specified critical depth.

        Args:
            critical_depth: neighbor depth at which the limiting starts. If either
                neighbor (T grid) is shallower than this value, the depth of velocity
                point (U or V grid) is restricted. If not provided, ``self.Dcrit`` is
                used.
        """
        assert (
            not self._initialized
        ), "limit_velocity_depth cannot be called after the domain has been initialized"
        if critical_depth is None:
            critical_depth = self.Dcrit
        tdepth = self.H_[1::2, 1::2]
        Vsel = (tdepth[1:, :] <= critical_depth) | (tdepth[:-1, :] <= critical_depth)
        self.H_[2:-2:2, 1::2][Vsel] = np.minimum(tdepth[1:, :], tdepth[:-1, :])[Vsel]
        Usel = (tdepth[:, 1:] <= critical_depth) | (tdepth[:, :-1] <= critical_depth)
        self.H_[1::2, 2:-2:2][Usel] = np.minimum(tdepth[:, 1:], tdepth[:, :-1])[Usel]
        self.logger.info(
            "limit_velocity_depth has decreased depth in %i U points (%i currently"
            " unmasked), %i V points (%i currently unmasked)."
            % (
                Usel.sum(),
                Usel.sum(where=self.mask_[1::2, 2:-2:2] > 0),
                Vsel.sum(),
                Vsel.sum(where=self.mask_[2:-2:2, 1::2] > 0),
            )
        )

    def mask_rectangle(
        self,
        xmin: Optional[float] = None,
        xmax: Optional[float] = None,
        ymin: Optional[float] = None,
        ymax: Optional[float] = None,
        value: int = 0,
    ):
        """Mask all points that fall within the specified rectangle.

        Args:
            xmin: lower native x coordinate of the rectangle to mask
                (default: left boundary of the domain)
            xmax: upper native x coordinate of the rectangle to mask
                (default: right boundary of the domain)
            ymin: lower native y coordinate of the rectangle to mask
                (default: bottom boundary of the domain)
            ymax: upper native y coordinate of the rectangle to mask
                (default: top boundary of the domain)

        Coordinates will be interpreted as longitude, latitude if the domain is
        configured as spherical; otherwise they will be interpreted as Cartesian
        x and y (m).
        """
        assert (
            not self._initialized
        ), "adjust_mask cannot be called after the domain has been initialized."
        selected = np.ones(self.mask.shape, dtype=bool)
        x, y = (self.lon, self.lat) if self.spherical else (self.x, self.y)
        if xmin is not None:
            selected &= x >= xmin
        if xmax is not None:
            selected &= x <= xmax
        if ymin is not None:
            selected &= y >= ymin
        if ymax is not None:
            selected &= y <= ymax
        self.mask[selected] = value

    def mask_indices(self, istart, istop, jstart, jstop, value: int = 0):
        """Mask all points that fall within the specified rectangle.

        Args:
            istart: lower x index (first that is included)
            istop: upper x index (first that is EXcluded)
            jstart: lower y index (first that is included)
            jstop: upper y index (first that is EXcluded)
        """
        istart = istart - self.tiling.xoffset + self.halox
        istop = istop - self.tiling.xoffset + self.halox
        jstart = jstart - self.tiling.yoffset + self.haloy
        jstop = jstop - self.tiling.yoffset + self.haloy
        istart = max(0, istart)
        istop = max(istart, min(self.nx + 2 * self.halox, istop))
        jstart = max(0, jstart)
        jstop = max(jstart, min(self.ny + 2 * self.haloy, jstop))
        self.mask_[
            1 + 2 * jstart : 1 + 2 * jstop, 1 + 2 * istart : 1 + 2 * istop
        ] = value

    def rotate(self) -> "Domain":
        def tp(array):
            return None if array is None else np.transpose(array)[::-1, :]

        return create(
            self.ny,
            self.nx,
            self.nz,
            tp(self.lon),
            tp(self.lat),
            tp(self.x),
            tp(self.y),
            self.spherical,
            tp(self.mask),
            tp(self.H),
            tp(self.z0b_min),
            tp(self.cor),
            Dmin=self.Dmin,
            vertical_coordinate_method=self.vertical_coordinate_method,
            ddl=self.ddl,
            ddu=self.ddu,
            Dgamma=self.Dgamma,
            gamma_surf=self.gamma_surf,
        )

    def plot(
        self,
        fig=None,
        show_bathymetry: bool = True,
        show_mask: bool = False,
        show_mesh: bool = True,
        show_rivers: bool = True,
        editable: bool = False,
        sub: bool = False,
        spherical: Optional[bool] = None,
    ):
        """Plot the domain, optionally including bathymetric depth, mesh and river positions.

        Args:
            fig: :class:`matplotlib.figure.Figure` instance to plot to. If not provided,
                a new figure is created.
            show_bathymetry: show bathymetry as color map
            show_mask: show mask as color map (this disables ``show_bathymetry``)
            show_mesh: show model grid
            show_rivers: show rivers with position and name
            editable: allow interactive selection of rectangular regions in the domain
                plot that are subsequently masked out
            sub: plot the subdomain, not the global domain

        Returns:
            :class:`matplotlib.figure.Figure` instance for processes with rank 0 or if
                ``sub`` is ``True``, otherwise ``None``
        """
        import matplotlib.pyplot
        import matplotlib.collections
        import matplotlib.widgets

        if self.glob is not self and not sub:
            # We need to plot the global domain; not the current subdomain.
            # If we are the root, divert the plot command to the global domain.
            # Otherwise just ignore this and return.
            if self.glob:
                return self.glob.plot(
                    fig, show_bathymetry, show_mask, show_mesh, show_rivers, editable
                )
            return

        if fig is None:
            fig, ax = matplotlib.pyplot.subplots(
                figsize=(0.15 * self.nx, 0.15 * self.ny)
            )
        else:
            ax = fig.gca()

        if spherical is None:
            spherical = self.spherical
        x, y = (self.lon, self.lat) if spherical else (self.x, self.y)

        local_slice, _, _, _ = self.tiling.subdomain2slices(
            halo_sub=0, halo_glob=4, scale=2, share=1, exclude_global_halos=True
        )
        if show_mask:
            c = ax.pcolormesh(
                x[local_slice],
                y[local_slice],
                self.mask[local_slice],
                alpha=0.5 if show_mesh else 1,
                shading="auto",
            )
        elif show_bathymetry:
            import cmocean

            cm = cmocean.cm.deep
            cm.set_bad("gray")
            c = ax.pcolormesh(
                x[local_slice],
                y[local_slice],
                np.ma.array(self.H[local_slice], mask=self.mask[local_slice] == 0),
                alpha=0.5 if show_mesh else 1,
                shading="auto",
                cmap=cm,
            )
            # c = ax.contourf(x, y, np.ma.array(self.H, mask=self.mask==0), 20, alpha=0.5 if show_mesh else 1)
            cb = fig.colorbar(c)
            cb.set_label("undisturbed water depth (m)")

        if show_rivers:
            for river in self.rivers.values():
                iloc, jloc = 1 + river.i * 2, 1 + river.j * 2
                lon = self.lon_[jloc, iloc]
                lat = self.lat_[jloc, iloc]
                ax.plot([lon], [lat], ".r")
                ax.text(lon, lat, river.name, color="r")

        def plot_mesh(ax, x, y, **kwargs):
            segs1 = np.stack((x, y), axis=2)
            segs2 = segs1.transpose(1, 0, 2)
            ax.add_collection(matplotlib.collections.LineCollection(segs1, **kwargs))
            ax.add_collection(matplotlib.collections.LineCollection(segs2, **kwargs))

        if show_mesh:
            plot_mesh(
                ax, x[::2, ::2], y[::2, ::2], colors="k", linestyle="-", linewidth=0.3
            )
            # ax.pcolor(x[1::2, 1::2], y[1::2, 1::2], np.ma.array(x[1::2, 1::2], mask=True), edgecolors='k', linestyles='--', linewidth=.2)
            # pc = ax.pcolormesh(x[1::2, 1::2], y[1::2, 1::2],  np.ma.array(x[1::2, 1::2], mask=True), edgecolor='gray', linestyles='--', linewidth=.2)
            ax.plot(x[::2, ::2], y[::2, ::2], "xk", markersize=3.0)
            ax.plot(x[1::2, 1::2], y[1::2, 1::2], ".k", markersize=2.5)
        ax.set_xlabel("longitude (degrees East)" if spherical else "x (m)")
        ax.set_ylabel("latitude (degrees North)" if spherical else "y (m)")
        if not spherical:
            ax.axis("equal")
        xmin, xmax = np.nanmin(x), np.nanmax(x)
        ymin, ymax = np.nanmin(y), np.nanmax(y)
        xmargin = 0.05 * (xmax - xmin)
        ymargin = 0.05 * (ymax - ymin)
        ax.set_xlim(xmin - xmargin, xmax + xmargin)
        ax.set_ylim(ymin - ymargin, ymax + ymargin)

        def on_select(eclick, erelease):
            xmin, xmax = (
                min(eclick.xdata, erelease.xdata),
                max(eclick.xdata, erelease.xdata),
            )
            ymin, ymax = (
                min(eclick.ydata, erelease.ydata),
                max(eclick.ydata, erelease.ydata),
            )
            self.mask_rectangle(xmin, xmax, ymin, ymax)
            c.set_array(np.ma.array(self.H, mask=self.mask == 0).ravel())
            fig.canvas.draw()
            # self.sel.set_active(False)
            # self.sel = None
            # ax.draw()
            # fig.clf()
            # self.plot(fig=fig, show_mesh=show_mesh)

        if editable:
            self.sel = matplotlib.widgets.RectangleSelector(
                ax, on_select, useblit=True, button=[1], interactive=False
            )
        return fig

    def save(self, path: str, full: bool = False, sub: bool = False):
        """Save grid to a NetCDF file that can be interpreted by :func:`load`.

        Args:
            path: NetCDF file to save to
            full: also save every field including its halos separately
                This can be useful for debugging.
            sub: save the local subdomain, not the global domain
        """
        if self.glob is not self and not sub:
            # We need to save the global domain; not the current subdomain.
            # If we are the root, divert the plot command to the global domain.
            # Otherwise just ignore this and return.
            if self.glob:
                self.glob.save(path, full)
            return

        with netCDF4.Dataset(path, "w") as nc:

            def create(
                name, units, long_name, values, coordinates: str, dimensions=("y", "x")
            ):
                fill_value = None
                if np.ma.getmask(values) is not np.ma.nomask:
                    fill_value = values.fill_value
                ncvar = nc.createVariable(
                    name, values.dtype, dimensions, fill_value=fill_value
                )
                ncvar.units = units
                ncvar.long_name = long_name
                # ncvar.coordinates = coordinates
                ncvar[...] = values

            def create_var(name, units, long_name, values, values_):
                if values is None:
                    return
                create(
                    name,
                    units,
                    long_name,
                    values,
                    coordinates="lon lat" if self.spherical else "x y",
                )
                if full:
                    create(
                        name + "_",
                        units,
                        long_name,
                        values_,
                        dimensions=("y_", "x_"),
                        coordinates="lon_ lat_" if self.spherical else "x_ y_",
                    )

            nc.createDimension("x", self.H.shape[1])
            nc.createDimension("y", self.H.shape[0])
            if full:
                nc.createDimension("x_", self.H_.shape[1])
                nc.createDimension("y_", self.H_.shape[0])
                create_var("dx", "m", "dx", self.dx, self.dx_)
                create_var("dy", "m", "dy", self.dy, self.dy_)
                create_var("area", "m2", "area", self.dx * self.dy, self.dx_ * self.dy_)
            create_var("lat", "degrees_north", "latitude", self.lat, self.lat_)
            create_var("lon", "degrees_east", "longitude", self.lon, self.lon_)
            create_var("x", "m", "x", self.x, self.x_)
            create_var("y", "m", "y", self.y, self.y_)
            create_var("H", "m", "undisturbed water depth", self.H, self.H_)
            create_var("mask", "", "mask", self.mask, self.mask_)
            create_var("z0b_min", "m", "bottom roughness", self.z0b_min, self.z0b_min_)
            create_var("cor", "", "Coriolis parameter", self.cor, self.cor_)

    def save_grids(self, path: str):
        with netCDF4.Dataset(path, "w") as nc:
            self.T.add_to_netcdf(nc)
            self.U.add_to_netcdf(nc, postfix="u")
            self.V.add_to_netcdf(nc, postfix="v")
            self.X.add_to_netcdf(nc, postfix="x")

    def contains(self, x: float, y: float, include_halos: bool = False) -> bool:
        """Determine whether the domain contains the specified point.

        Args:
            x: native x coordinate (longitude for spherical grids, Cartesian coordinate
                in m otherwise)
            y: native y coordinate (latitude for spherical grids, Cartesian coordinate
                in m otherwise)
            include_halos: whether to also search the halos

        Returns:
            True if the point falls within the domain, False otherwise
        """
        local_slice, _, _, _ = self.tiling.subdomain2slices(
            halo_sub=4,
            halo_glob=4,
            scale=2,
            share=1,
            exclude_halos=not include_halos,
            exclude_global_halos=True,
        )
        allx, ally = (self.lon_, self.lat_) if self.spherical else (self.x_, self.y_)
        allx, ally = allx[local_slice], ally[local_slice]
        ny, nx = allx.shape

        # Determine whether point falls within current subdomain
        # based on https://wrf.ecse.rpi.edu/Research/Short_Notes/pnpoly.html
        x_bnd = np.concatenate(
            (
                allx[0, :-1],
                allx[:-1, -1],
                allx[-1, nx - 1 : 0 : -1],
                allx[ny - 1 : 0 : -1, 0],
            )
        )
        y_bnd = np.concatenate(
            (
                ally[0, :-1],
                ally[:-1, -1],
                ally[-1, nx - 1 : 0 : -1],
                ally[ny - 1 : 0 : -1, 0],
            )
        )
        assert not np.isnan(x_bnd).any(), "Invalid x boundary: %s." % (x_bnd,)
        assert not np.isnan(y_bnd).any(), "Invalid y boundary: %s." % (y_bnd,)
        assert x_bnd.size == 2 * ny + 2 * nx - 4
        inside = False
        for i, (vertxi, vertyi) in enumerate(zip(x_bnd, y_bnd)):
            vertxj, vertyj = x_bnd[i - 1], y_bnd[i - 1]
            if (vertyi > y) != (vertyj > y) and x < (vertxj - vertxi) * (y - vertyi) / (
                vertyj - vertyi
            ) + vertxi:
                inside = not inside
        return inside

    def update_depth(self, _3d: bool = False, timestep: float = 1.0):
        """Use old and new surface elevation on T grid to update elevations on U, V, X grids
        and subsequently update total water depth ``D`` on all grids.

        Args:
            _3d: update elevations of the macrotimestep (``zin``) rather than
                elevations of the microtimestep (``z``).This first synchronizes the
                elevations of the macrotimestep on the T grid (``self.T.zin``) with
                those of the microtimestep (``self.T.z``). It also updates layer
                thicknesses ``hn``, layer center depths ``zc`` and interface depths
                ``zf`` on all grids.
            timestep: time step (s) for layer thickness change

        This routine will ensure values are up to date in the domain interior and in
        the halos, but that this requires that ``self.T.z`` (and old elevations
        ``self.T.zo`` or ``self.T.zio``) are already up to date in halos.
        """
        if _3d:
            # Store current elevations as previous elevations (on the 3D time step)
            self.T.zio.all_values[...] = self.T.zin.all_values
            self.U.zio.all_values[...] = self.U.zin.all_values
            self.V.zio.all_values[...] = self.V.zin.all_values
            self.X.zio.all_values[...] = self.X.zin.all_values

            # Synchronize new elevations on the 3D time step to those of the 2D time
            # step that has just completed.
            self.T.zin.all_values[...] = self.T.z.all_values

            z_T, z_U, z_V, z_X, zo_T = (
                self.T.zin,
                self.U.zin,
                self.V.zin,
                self.X.zin,
                self.T.zio,
            )
        else:
            z_T, z_U, z_V, z_X, zo_T = self.T.z, self.U.z, self.V.z, self.X.z, self.T.zo

        # Compute surface elevation on U, V, X grids.
        # These must lag 1/2 a timestep behind the T grid.
        # They are therefore calculated from the average of old and new elevations on
        # the T grid.
        z_T_half = 0.5 * (zo_T + z_T)
        z_T_half.interp(z_U)
        z_T_half.interp(z_V)
        z_T_half.interp(z_X)
        _pygetm.clip_z(z_U, self.Dmin)
        _pygetm.clip_z(z_V, self.Dmin)
        _pygetm.clip_z(z_X, self.Dmin)

        # Halo exchange for elevation on U, V grids, needed because the very last
        # points in the halos (x=-1 for U, y=-1 for V) are not valid after
        # interpolating from the T grid above.
        # These elevations are needed to later compute velocities from transports
        # (by dividing by layer thicknesses, which are computed from elevation)
        # These velocities will be advected, and therefore need to be valid througout
        # the halos. We do not need to halo-exchange elevation on the X grid, since
        # that needs to be be valid at the innermost halo point only, which is ensured
        # by z_T exchange.
        z_U.update_halos(parallel.Neighbor.RIGHT)
        z_V.update_halos(parallel.Neighbor.TOP)

        # Update total water depth D on T, U, V, X grids
        # This also processes the halos; no further halo exchange needed.
        np.add(
            self.T.H.all_values,
            z_T.all_values,
            where=self.T._water,
            out=self.T.D.all_values,
        )
        np.add(
            self.U.H.all_values,
            z_U.all_values,
            where=self.U._water,
            out=self.U.D.all_values,
        )
        np.add(
            self.V.H.all_values,
            z_V.all_values,
            where=self.V._water,
            out=self.V.D.all_values,
        )
        np.add(
            self.X.H.all_values,
            z_X.all_values,
            where=self.X._water,
            out=self.X.D.all_values,
        )

        # Update dampening factor (0-1) for shallow water
        _pygetm.alpha(self.U.D, self.Dmin, self.Dcrit, self.U.alpha)
        _pygetm.alpha(self.V.D, self.Dmin, self.Dcrit, self.V.alpha)

        # Update total water depth on advection grids. These must be 1/2 timestep
        # behind the T grid. That's already the case for the X grid, but for the T grid
        # we explicitly compute and use the average of old and new D.
        np.add(self.T.H.all_values, z_T_half.all_values, out=self.D_T_half.all_values)
        self.UU.D.all_values[:, :-1] = self.D_T_half.all_values[:, 1:]
        self.VV.D.all_values[:-1, :] = self.D_T_half.all_values[1:, :]
        self.UV.D.all_values[:, :] = self.VU.D.all_values[:, :] = self.X.D.all_values[
            1:, 1:
        ]

        if _3d:
            # Store previous layer thicknesses
            self.T.ho.all_values[...] = self.T.hn.all_values
            self.U.ho.all_values[...] = self.U.hn.all_values
            self.V.ho.all_values[...] = self.V.hn.all_values
            self.X.ho.all_values[...] = self.X.hn.all_values

            # Update layer thicknesses (hn) on all grids, using bathymetry H and new
            # elevations zin (on the 3D timestep)
            self.do_vertical(timestep)

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
            self.h_T_half.all_values[...] = 0.5 * (
                self.T.ho.all_values + self.T.hn.all_values
            )
            self.UU.hn.all_values[:, :, :-1] = self.h_T_half.all_values[:, :, 1:]
            self.VV.hn.all_values[:, :-1, :] = self.h_T_half.all_values[:, 1:, :]
            self.UV.hn.all_values[:, :, :] = self.VU.hn.all_values[
                :, :, :
            ] = self.X.hn.all_values[:, 1:, 1:]

            if self.depth.saved:
                # Update depth-below-surface at layer centers.
                # Elsewhere this can be used as approximate pressure in dbar
                _pygetm.thickness2center_depth(self.T.mask, self.T.hn, self.depth)

            if self.open_boundaries.zc.saved:
                # Update vertical coordinate at open boundary, used to interpolate
                # inputs on z grid to dynamic model depths
                self.open_boundaries.zc.all_values[...] = self.T.zc.all_values[
                    :, self.open_boundaries.j, self.open_boundaries.i
                ].T
