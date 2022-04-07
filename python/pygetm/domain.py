from typing import Optional, Tuple, List, Mapping
import operator
import logging
import os.path
import collections

import numpy
import numpy.typing
import xarray
import netCDF4

from . import _pygetm
from . import core
from . import parallel
from . import output
from . import input
from .constants import FILL_VALUE, CENTERS, INTERFACES

WEST  = 1
NORTH = 2
EAST  = 3
SOUTH = 4

def find_interfaces(c: numpy.ndarray):
    c_if = numpy.empty((c.size + 1),)
    d = numpy.diff(c)
    c_if[1:-1] = c[:-1] + 0.5 * d
    c_if[0] = c[0] - 0.5 * d[0]
    c_if[-1] = c[-1] + 0.5 * d[-1]
    return c_if

class Grid(_pygetm.Grid):
    _coordinate_arrays = 'x', 'y', 'lon', 'lat'
    _readonly_arrays = _coordinate_arrays + ('dx', 'dy', 'idx', 'idy', 'dlon', 'dlat', 'area', 'iarea', 'cor')
    _fortran_arrays = _readonly_arrays + ('H', 'D', 'mask', 'z', 'zo', 'ho', 'hn', 'zc', 'zf', 'z0b', 'z0b_min', 'zio', 'zin')
    _all_arrays = tuple(['_%s' % n for n in _fortran_arrays] + ['_%si' % n for n in _coordinate_arrays] + ['_%si_' % n for n in _coordinate_arrays])
    __slots__ = _all_arrays + ('halo', 'type', 'ioffset', 'joffset', 'postfix', 'ugrid', 'vgrid', '_sin_rot', '_cos_rot', 'rotation', 'nbdyp', 'overlap')

    array_args = {
        'x': dict(units='m', constant=True, fill_value=FILL_VALUE),
        'y': dict(units='m', constant=True, fill_value=FILL_VALUE),
        'lon': dict(units='degrees_north', long_name='longitude', constant=True, fill_value=FILL_VALUE),
        'lat': dict(units='degrees_east', long_name='latitude', constant=True, fill_value=FILL_VALUE),
        'dx': dict(units='m', constant=True, fill_value=FILL_VALUE),
        'dy': dict(units='m', constant=True, fill_value=FILL_VALUE),
        'idx': dict(units='m-1', constant=True, fill_value=FILL_VALUE),
        'idy': dict(units='m-1', constant=True, fill_value=FILL_VALUE),
        'dlon': dict(units='degrees_north', constant=True, fill_value=FILL_VALUE),
        'dlat': dict(units='degrees_east', constant=True, fill_value=FILL_VALUE),
        'H': dict(units='m', long_name='water depth at rest', constant=True, fill_value=FILL_VALUE),
        'D': dict(units='m', long_name='water depth', fill_value=FILL_VALUE),
        'mask': dict(constant=True, fill_value=0),
        'z': dict(units='m', long_name='elevation', fill_value=FILL_VALUE),
        'zo': dict(units='m', long_name='elevation at previous time step', fill_value=FILL_VALUE),
        'zin': dict(units='m', long_name='elevation at 3d time step', fill_value=FILL_VALUE),
        'zio': dict(units='m', long_name='elevation at previous 3d time step', fill_value=FILL_VALUE),
        'area': dict(units='m2', long_name='grid cell area', constant=True, fill_value=FILL_VALUE),
        'iarea': dict(units='m-2', long_name='inverse of grid cell area', constant=True, fill_value=FILL_VALUE),
        'cor': dict(units='1', long_name='Coriolis parameter', constant=True, fill_value=FILL_VALUE),
        'ho': dict(units='m', long_name='layer heights at previous time step', fill_value=FILL_VALUE),
        'hn': dict(units='m', long_name='layer heights', fill_value=FILL_VALUE),
        'zc': dict(units='m', long_name='depth', fill_value=FILL_VALUE),
        'zf': dict(units='m', long_name='interface depth', fill_value=FILL_VALUE),
        'z0b': dict(units='m', long_name='hydrodynamic bottom roughness', fill_value=FILL_VALUE),
        'z0b_min': dict(units='m', long_name='physical bottom roughness', constant=True, fill_value=FILL_VALUE),
    }

    def __init__(self, domain: 'Domain', grid_type: int, ioffset: int, joffset: int, overlap: int=0, ugrid: Optional['Grid']=None, vgrid: Optional['Grid']=None):
        _pygetm.Grid.__init__(self, domain, grid_type)
        self.halo = domain.halo
        self.type = grid_type
        self.ioffset = ioffset
        self.joffset = joffset
        self.overlap = overlap
        self.postfix = {_pygetm.TGRID: 't', _pygetm.UGRID: 'u', _pygetm.VGRID: 'v', _pygetm.XGRID: 'x', _pygetm.UUGRID: '_uu_adv', _pygetm.VVGRID: '_vv_adv', _pygetm.UVGRID: '_uv_adv', _pygetm.VUGRID: '_vu_adv'}[grid_type]
        self.ugrid: Optional[Grid] = ugrid
        self.vgrid: Optional[Grid] = vgrid
        self._sin_rot: Optional[numpy.ndarray] = None
        self._cos_rot: Optional[numpy.ndarray] = None

        for name in self._readonly_arrays:
            self.setup_array(name, register=False)
        self._iarea.all_values[...] = 1. / self._area.all_values
        self._idx.all_values[...] = 1. / self._dx.all_values
        self._idy.all_values[...] = 1. / self._dy.all_values
        for name in self._readonly_arrays:
            getattr(self, name).all_values.flags.writeable = True

    def initialize(self, nbdyp: int):
        for name in self._fortran_arrays:
            if name in self._readonly_arrays:
                getattr(self, name).register()
            else:
                self.setup_array(name)
        self.rotation = core.Array.create(grid=self, dtype=self.x.dtype, name='rotation' + self.postfix, units='rad', long_name='grid rotation with respect to true North')
        self.setup_array(name, self.rotation)
        self.zc.all_values.fill(0.)
        self.zf.all_values.fill(0.)
        self.z0b.all_values[...] = self.z0b_min.all_values
        self.nbdyp = nbdyp

    def setup_array(self, name: str, array: Optional[core.Array]=None, register: bool=True):
        if array is None:
            # No array provided, so it must live in Fortran; retrieve it
            array = core.Array(name=name + self.postfix, **self.array_args[name])
            setattr(self, '_%s' % name, self.wrap(array, name.encode('ascii'), register=register))

        # Obtain corresponding array on the supergrid. If this does not exist, we are done
        source = getattr(self.domain, name + '_', None)
        if source is None:
            return

        nj, ni = self.ny_, self.nx_
        has_bounds = self.ioffset > 0 and self.joffset > 0 and self.domain.H_.shape[-1] >= self.ioffset + 2 * ni and self.domain.H_.shape[-2] >= self.joffset + 2 * nj
        valid = self.domain.mask_[self.joffset:self.joffset + 2 * nj:2, self.ioffset:self.ioffset + 2 * ni:2] > 0
        values = source[self.joffset:self.joffset + 2 * nj:2, self.ioffset:self.ioffset + 2 * ni:2]
        array.all_values.fill(numpy.nan)
        slc = valid if name in ('z', 'zo',  'z0b_min') else (Ellipsis,)
        array.all_values[:values.shape[0], :values.shape[1]][slc] = values[slc]
        if has_bounds and name in self._coordinate_arrays:
            # Generate interface coordinates. These are not represented in Fortran as they are only needed for plotting.
            # The interface coordinates are slices that point to the supergrid data; they thus do not consume additional memory.
            values_i = source[self.joffset - 1:self.joffset + 2 * nj + 1:2, self.ioffset - 1:self.ioffset + 2 * ni + 1:2]
            setattr(self, '_%si_' % name, values_i)
            setattr(self, '_%si' % name, values_i[self.halo:-self.halo, self.halo:-self.halo])

    def rotate(self, u: numpy.typing.ArrayLike, v: numpy.typing.ArrayLike, to_grid: bool= True) -> Tuple[numpy.typing.ArrayLike, numpy.typing.ArrayLike]:
        if self._sin_rot is None:
            self._sin_rot = numpy.sin(self.rotation.all_values)
            self._cos_rot = numpy.cos(self.rotation.all_values)
        sin_rot = -self._sin_rot if to_grid else self._sin_rot
        u_new = u * self._cos_rot - v * sin_rot
        v_new = u * sin_rot + v * self._cos_rot
        return u_new, v_new

    def array(self, *args, **kwargs) -> core.Array:
        return core.Array.create(self,  *args, **kwargs)

    def add_to_netcdf(self, nc: netCDF4.Dataset, postfix: str=''):
        xdim, ydim = 'x' + postfix, 'y' + postfix
        def save(name, units='', long_name=None):
            data = getattr(self, name)
            ncvar = nc.createVariable(name + postfix, data.dtype, (ydim, xdim))
            ncvar[...] = data
            ncvar.units = units
            ncvar.long_name = long_name or name
        ny, nx = self.x.shape
        nc.createDimension(xdim, nx)
        nc.createDimension(ydim, ny)
        save('dx', 'm')
        save('dy', 'm')
        save('H', 'm', 'undisturbed water depth')
        save('mask')
        save('area', 'm2')
        save('cor', 's-1', 'Coriolis parameter')

    def nearest_point(self, x: float, y: float, mask: Optional[Tuple[int]]=None) -> Optional[Tuple[int, int]]:
        """Return index (i,j) of point nearest to specified coordinate."""
        if not self.domain.contains(x, y):
            return None
        allx, ally = (self.lon, self.lat) if self.domain.spherical else (self.x, self.y)
        dist = (allx.values - x)**2 + (ally.values - y)**2
        if mask is not None:
            if isinstance(mask, int):
                mask = (mask,)
            invalid = numpy.ones(dist.shape, dtype=bool)
            for mask_value in mask:
                invalid = numpy.logical_and(invalid, self.mask.values != mask_value)
            dist[invalid] = numpy.inf
        idx = numpy.nanargmin(dist)
        return numpy.unravel_index(idx, dist.shape)

for membername in Grid._all_arrays:
    setattr(Grid, membername[1:], property(operator.attrgetter(membername)))

def read_centers_to_supergrid(ncvar, ioffset: int, joffset: int, nx: int, ny: int, dtype=None):
    if dtype is None:
        dtype = ncvar.dtype
    data = numpy.ma.masked_all(ncvar.shape[:-2] + (ny * 2 + 1, nx * 2 + 1), dtype=dtype)

    # Create an array to data at centers (T points),
    # with strips of size 1 on all sides to support interpolation to interfaces
    data_centers = numpy.ma.masked_all(ncvar.shape[:-2] + (ny + 2, nx + 2), dtype=dtype)

    masked_values = []
    if hasattr(ncvar, 'missing_value'):
        masked_values.append(numpy.array(ncvar.missing_value, ncvar.dtype))

    # Extend the read domain (T grid) by 1 each side, where possible
    # That will allow us to interpolate (rater than extrapolate) to values at the interfaces
    ex_imin = 0 if ioffset == 0 else 1
    ex_imax = 0 if ioffset + nx == ncvar.shape[-1] else 1
    ex_jmin = 0 if joffset == 0 else 1
    ex_jmax = 0 if joffset + ny == ncvar.shape[-2] else 1
    data_centers[..., 1 - ex_jmin:1 + ny + ex_jmax, 1 - ex_imin:1 + nx + ex_imax] = ncvar[..., joffset - ex_jmin:joffset + ny + ex_jmax, ioffset - ex_imin:ioffset + nx + ex_imax]
    for value in masked_values:
        data_centers = numpy.ma.masked_equal(data_centers, value, copy=False)

    data_if_ip = numpy.ma.masked_all((4,) + ncvar.shape[:-2] + (ny + 1, nx + 1), dtype=dtype)
    data_if_ip[0, ...] = data_centers[...,  :-1,  :-1]
    data_if_ip[1, ...] = data_centers[..., 1:,    :-1]
    data_if_ip[2, ...] = data_centers[...,  :-1, 1:  ]
    data_if_ip[3, ...] = data_centers[..., 1:,   1:  ]
    data[..., 1::2, 1::2] = data_centers[1:-1, 1:-1]
    data[..., ::2, ::2] = data_if_ip.mean(axis=0)

    data_if_ip[0, ..., :-1] = data_centers[...,  :-1,  1:-1]
    data_if_ip[1, ..., :-1] = data_centers[..., 1:,    1:-1]
    data[..., ::2, 1::2] = data_if_ip[:2, ..., :-1].mean(axis=0)

    data_if_ip[0, ..., :-1, :] = data_centers[..., 1:-1,  :-1]
    data_if_ip[1, ..., :-1, :] = data_centers[..., 1:-1, 1:, ]
    data[..., 1::2, ::2] = data_if_ip[:2, ..., :-1, :].mean(axis=0)

    if len(masked_values) > 0:
        data.set_fill_value(masked_values[0])

    return data

deg2rad = numpy.pi / 180        # degree to radian conversion
rearth = 6378815.               # radius of the earth (m)
omega = 2. * numpy.pi / 86164.  # rotation rate of the earth (rad/s), 86164 is number of seconds in a sidereal day

def center_to_supergrid_1d(data) -> numpy.ndarray:
    assert data.ndim == 1, 'data must be one-dimensional'
    assert data.size > 1, 'data must have at least 2 elements'
    data_sup = numpy.empty((data.size * 2 + 1,))
    data_sup[1::2] = data
    data_sup[2:-2:2] = 0.5 * (data[1:] + data[:-1])
    data_sup[0] = 2 * data_sup[1] - data_sup[2]
    data_sup[-1] = 2 * data_sup[-2] - data_sup[-3]
    return data_sup

def interfaces_to_supergrid_1d(data, out: Optional[numpy.ndarray]=None) -> numpy.ndarray:
    assert data.ndim == 1, 'data must be one-dimensional'
    assert data.size > 1, 'data must have at least 2 elements'
    if out is None:
        out = numpy.empty((data.size * 2 - 1,))
    out[0::2] = data
    out[1::2] = 0.5 * (data[1:] + data[:-1])
    return out

def interfaces_to_supergrid_2d(data, out: Optional[numpy.ndarray]=None) -> numpy.ndarray:
    assert data.ndim == 2, 'data must be two-dimensional'
    assert data.shape[0] > 1 and data.shape[1] > 1, 'data must have at least 2 elements'
    if out is None:
        out = numpy.empty((data.shape[1] * 2 - 1, data.shape[0] * 2 - 1))
    out[0::2, 0::2] = data
    out[1::2, 0::2] = 0.5 * (data[:-1, :] + data[1:, :])
    out[0::2, 1::2] = 0.5 * (data[:,:-1] + data[:,1:])
    out[1::2, 1::2] = 0.25 * (data[:-1,:-1] + data[:-1,1:] + data[1:,:-1] + data[1:,1:])
    return out

def create_cartesian(x, y, nz: int, interfaces=False, **kwargs) -> 'Domain':
    """Create Cartesian domain from 1D arrays with longitudes and latitudes."""
    assert x.ndim == 1, 'x coordinate must be one-dimensional'
    assert y.ndim == 1, 'y coordinate must be one-dimensional'

    if interfaces:
        nx, ny = x.size - 1, y.size - 1
        x, y = interfaces_to_supergrid_1d(x), interfaces_to_supergrid_1d(y)
    else:
        nx, ny = x.size, y.size
        x, y = center_to_supergrid_1d(x), center_to_supergrid_1d(y)
    return Domain.create(nx, ny, nz, x=x, y=y[:, numpy.newaxis], **kwargs)

def create_spherical(lon, lat, nz: int, interfaces=False, **kwargs) -> 'Domain':
    """Create spherical domain from 1D arrays with longitudes and latitudes."""
    assert lon.ndim == 1, 'longitude coordinate must be one-dimensional'
    assert lat.ndim == 1, 'latitude coordinate must be one-dimensional'

    if interfaces:
        nx, ny = lon.size - 1, lat.size - 1
        x, y = interfaces_to_supergrid_1d(lon), interfaces_to_supergrid_1d(lat)
    else:
        nx, ny = lon.size, lat.size
        x, y = center_to_supergrid_1d(lon), center_to_supergrid_1d(lat)
    return Domain.create(nx, ny, nz, lon=lon, lat=lat[:, numpy.newaxis], spherical=True, **kwargs)

def create_spherical_at_resolution(minlon: float, maxlon: float, minlat: float, maxlat: float, resolution: float, nz: int, **kwargs) -> 'Domain':
    """Create spherical domain longitude range, latitude range and desired resolution in m."""
    assert maxlon > minlon, 'Maximum longitude %s must be greater than minimum longitude %s' % (maxlon, minlon)
    assert maxlat > minlat, 'Maximum latitude %s must be greater than minimum latitude %s' % (maxlat, minlat)
    assert resolution > 0, 'Desired resolution must be greater than 0, but is %s m' % (resolution,)
    dlat = resolution / (deg2rad * rearth)
    minabslat = min(abs(minlat), abs(maxlat))
    dlon = resolution / (deg2rad * rearth) / numpy.cos(deg2rad * minabslat)
    nx = int(numpy.ceil((maxlon - minlon) / dlon)) + 1
    ny = int(numpy.ceil((maxlat - minlat) / dlat)) + 1
    return create_spherical(numpy.linspace(minlon, maxlon, nx), numpy.linspace(minlat, maxlat, ny), nz=nz, interfaces=True, **kwargs)

class RiverTracer(core.Array):
    __slots__ = ('_follow',)
    def __init__(self, grid, river_name: str, tracer_name: str, value: numpy.ndarray, follow: numpy.ndarray, **kwargs):
        super().__init__(grid=grid, name=tracer_name + '_in_river_' + river_name, long_name='%s in river %s' % (tracer_name, river_name), **kwargs)
        self.wrap_ndarray(value)
        self._follow = follow

    @property
    def follow_target_cell(self) -> bool:
        return self._follow

    @follow_target_cell.setter
    def follow_target_cell(self, value: bool):
        self._follow[...] = value

class River:
    def __init__(self, name: str, i: int, j: int, zl: Optional[float]=None, zu: Optional[float]=None, x: Optional[float]=None, y: Optional[float]=None):
        self.name = name
        self.i = i
        self.j = j
        self.x = x
        self.y = y
        self.zl = zl
        self.zu = zu
        self.active = True
        self._tracers: Mapping[str, RiverTracer] = {}

    def locate(self,  grid: Grid):
        if self.x is not None:
            ind = grid.nearest_point(self.x, self.y, mask=1)
            if ind is None:
                self.active = False
            else:
                self.i = ind[1] + grid.domain.halox
                self.j = ind[0] + grid.domain.haloy
        return self.active

    def initialize(self, grid: Grid, flow: numpy.ndarray):
        assert self.active
        self.flow = core.Array(grid=grid, name='river_' + self.name + '_flow', units='m3 s-1', long_name='inflow from %s' % self.name)
        self.flow.wrap_ndarray(flow)

    def __getitem__(self, key) -> RiverTracer:
        return self._tracers[key]

    def __len__(self):
        return len(self._tracers)

    def __iter__(self):
        return iter(self._tracers)

class Rivers(collections.Mapping):
    def __init__(self, grid: Grid):
        self.grid = grid
        self._rivers: List[River] = []
        self._tracers = []
        self._frozen = False

    def add_by_index(self, name: str, i: int, j: int, **kwargs):
        """Add a river at a location specified by the indices of a tracer point"""
        assert not self._frozen, 'The river collection has already been initialized and can no longer be modified.'

        i_loc, j_loc =  i - self.grid.domain.tiling.xoffset, j - self.grid.domain.tiling.yoffset

        #if self.grid.domain.glob is not None and self.grid.domain is not self.grid.domain.glob:
        #    self.grid.domain.glob.rivers.add_by_index(name, i, j, **kwargs)

        if i_loc >= 0 and j_loc >= 0 and i_loc < self.grid.domain.T.nx and j_loc < self.grid.domain.T.ny:
            river = River(name, i_loc + self.grid.domain.halox, j_loc + self.grid.domain.haloy, **kwargs)
            self._rivers.append(river)
            return river

    def add_by_location(self, name: str, x: float, y: float, **kwargs):
        """Add a river at a location specified by the nearest coordinates (longitude and latitide on a spherical grid)"""
        #if self.grid.domain.glob is not None and self.grid.domain is not self.grid.domain.glob:
        #    self.grid.domain.glob.rivers.add_by_location(name, x, y, **kwargs)
        river = River(name, None, None, x=x, y=y, **kwargs)
        self._rivers.append(river)
        return river

    def add_tracer(self, name: str, units: str, follow_target_cell: bool=False) -> Tuple[numpy.ndarray, numpy.ndarray, List[RiverTracer]]:
        assert self._frozen, 'Tracers can be added only after the river collection has been initialized.'
        assert name not in self._tracers, 'A tracer with name %s has already been added.' % name
        values = numpy.zeros((len(self._rivers),))
        follow = numpy.full((len(self._rivers),), follow_target_cell, dtype=bool)
        self._tracers.append((values, follow))
        river_tracers: List[RiverTracer] = []
        for iriver, river in enumerate(self._rivers):
            river_tracer = RiverTracer(self.grid, river.name, name, values[..., iriver], follow[..., iriver], units=units, attrs={'_3d_only': True})
            river_tracers.append(river_tracer)
            river._tracers[name] = river_tracer
        return values, follow, river_tracers

    def initialize(self):
        assert not self._frozen
        self._frozen = True
        self._rivers = [river for river in self._rivers if river.locate(self.grid)]
        self.flow = numpy.zeros((len(self._rivers),))
        self.i = numpy.empty((len(self._rivers),), dtype=int)
        self.j = numpy.empty((len(self._rivers),), dtype=int)
        self.iarea = numpy.empty((len(self._rivers),))
        for iriver, river in enumerate(self._rivers):
            assert self.grid.mask.all_values[river.j, river.i] == 1, 'River %s is located at i=%i, j=%i, which is not water (it has mask value %i).' % (river.name, river.i, river.j, self.grid.mask.all_values[river.j, river.i])
            river.initialize(self.grid, self.flow[..., iriver])
            self.i[iriver] = river.i
            self.j[iriver] = river.j
            self.iarea[iriver] = self.grid.iarea.all_values[river.j, river.i]

    def __getitem__(self, key) -> River:
        for river in self._rivers:
            if key == river.name:
                return river
        raise KeyError()

    def __len__(self) -> int:
        return len(self._rivers)

    def __iter__(self):
        return map(operator.attrgetter('name'), self._rivers)

class Domain(_pygetm.Domain):
    @staticmethod
    def partition(tiling: parallel.Tiling, nx: int, ny: int, nz: int, global_domain: Optional['Domain'], halo: int=2, has_xy: bool=True, has_lonlat: bool=True, logger: Optional[logging.Logger]=None, **kwargs):
        assert nx == tiling.nx_glob and ny == tiling.ny_glob, 'Extent of global domain (%i, %i) does not match that of tiling (%i, %i).' % (ny, nx, tiling.ny_glob, tiling.nx_glob)
        assert tiling.n == tiling.comm.Get_size(), 'Number of active cores in subdomain decompositon (%i) does not match available number of cores (%i).' % (tiling.n, tiling.comm.Get_size())
        assert global_domain is None or global_domain.initialized

        halo = 4   # coordinates are scattered without their halo - Domain object will update halos upon creation
        share = 1  # one X point overlap in both directions between subdomains for variables on the supergrid

        coordinates = {'f': 'cor'}
        if has_xy:
            coordinates['x'] = 'x'
            coordinates['y'] = 'y'
        if has_lonlat:
            coordinates['lon'] = 'lon'
            coordinates['lat'] = 'lat'
        for name, att in coordinates.items():
            c = numpy.empty((2 * tiling.ny_sub + share, 2 * tiling.nx_sub + share))
            scatterer = parallel.Scatter(tiling, c, halo=0, share=share, scale=2, fill_value=numpy.nan)
            scatterer(None if global_domain is None else getattr(global_domain, att))
            kwargs[name] = c

        domain = Domain(tiling.nx_sub, tiling.ny_sub, nz, tiling=tiling, logger=logger, **kwargs)

        halo = 4
        parallel.Scatter(tiling, domain.mask_, halo=halo, share=share, scale=2)(None if global_domain is None else global_domain.mask_)
        parallel.Scatter(tiling, domain.H_, halo=halo, share=share, scale=2)(None if global_domain is None else global_domain.H_)
        parallel.Scatter(tiling, domain.z0b_min_, halo=halo, share=share, scale=2)(None if global_domain is None else global_domain.z0b_min_)
        parallel.Scatter(tiling, domain.z_, halo=halo, share=share, scale=2)(None if global_domain is None else global_domain.z_)
        parallel.Scatter(tiling, domain.zo_, halo=halo, share=share, scale=2)(None if global_domain is None else global_domain.zo_)

        return domain

    def exchange_metric(self, data, relative_in_x: bool=False, relative_in_y: bool=False, fill_value=numpy.nan):
        if not self.tiling:
            return

        halo = 2
        superhalo = 2 * halo

        # Expand the data array one each side
        data_ext = numpy.full((data.shape[0] + 2, data.shape[1] + 2), fill_value, dtype=data.dtype)
        data_ext[1:-1, 1:-1] = data
        self.tiling.wrap(data_ext, superhalo + 1).update_halos()

        # For values in the halo, compute their difference with the outer boundary of the subdomain we exchanged with (now the innermost halo point).
        # Then use that difference plus the value on our own boundary as values inside the halo.
        # This ensures coordinate variables are monotonically increasing in interior AND halos, even if periodic boundary conditions are used.
        if relative_in_x:
            data_ext[:, :superhalo + 1] += data_ext[:, superhalo + 1:superhalo + 2] - data_ext[:, superhalo:superhalo + 1]
            data_ext[:, -superhalo - 1:] += data_ext[:, -superhalo - 2:-superhalo - 1] - data_ext[:, -superhalo - 1:-superhalo]
        if relative_in_y:
            data_ext[:superhalo + 1, :] += data_ext[superhalo + 1:superhalo + 2, :] - data_ext[superhalo:superhalo + 1, :]
            data_ext[-superhalo - 1:, :] += data_ext[-superhalo - 2:-superhalo - 1, :] - data_ext[-superhalo - 1:-superhalo, :]

        # Since subdomains share the outer boundary, that boundary will be replicated in the outermost interior point and in the innermost halo point
        # We move the outer part of the halos (all but their innermost point) one point inwards to eliminate that overlapping point
        data[:superhalo,              : superhalo] = data_ext[             :superhalo,                    : superhalo    ]
        data[:superhalo,   superhalo  :-superhalo] = data_ext[             :superhalo,       superhalo + 1:-superhalo - 1]
        data[:superhalo,  -superhalo  :          ] = data_ext[             :superhalo,      -superhalo    :              ]
        data[ superhalo:  -superhalo, :superhalo]  = data_ext[superhalo + 1:-superhalo - 1,               : superhalo    ]
        data[ superhalo:  -superhalo, -superhalo:] = data_ext[superhalo + 1:-superhalo - 1, -superhalo    :              ]
        data[-superhalo:,             : superhalo] = data_ext[-superhalo:,                                : superhalo    ]
        data[-superhalo:,  superhalo  :-superhalo] = data_ext[-superhalo:,                   superhalo + 1:-superhalo - 1]
        data[-superhalo:, -superhalo  :          ] = data_ext[-superhalo:,                      -superhalo:              ]

    @staticmethod
    def create(nx: int, ny: int, nz: int, runtype: int=1, lon: Optional[numpy.ndarray]=None, lat: Optional[numpy.ndarray]=None, x: Optional[numpy.ndarray]=None, y: Optional[numpy.ndarray]=None, spherical: bool=False, mask: Optional[numpy.ndarray]=1, H: Optional[numpy.ndarray]=None, z0: Optional[numpy.ndarray]=0., f: Optional[numpy.ndarray]=None, tiling: Optional[parallel.Tiling]=None, z: Optional[numpy.ndarray]=0., zo: Optional[numpy.ndarray]=0., logger: Optional[logging.Logger]=None, **kwargs):
        global_domain = None
        logger = logger or parallel.getLogger()
        parlogger = logger.getChild('parallel')

        # Determine subdomain division
        if tiling is None:
            # No tiling provided - autodetect
            mask = numpy.broadcast_to(mask, (1 + 2 * ny, 1 + 2 * nx))
            tiling = parallel.Tiling.autodetect(mask=mask[1::2, 1::2], logger=parlogger, **kwargs)
        elif isinstance(tiling, str):
            # Path to dumped Tiling object provided
            if not os.path.isfile(tiling):
                logger.critical('Cannot find file %s. If tiling is a string, it must be the path to an existing file with a pickled tiling object.' % tiling)
                raise Exception()
            tiling = parallel.Tiling.load(tiling)
        else:
            # Existing tiling object provided - transfer extent of global domain to determine subdomain sizes
            if isinstance(tiling, tuple):
                tiling = parallel.Tiling(nrow=tiling[0], ncol=tiling[1], **kwargs)
            tiling.set_extent(nx, ny)
        tiling.report(parlogger)

        global_tiling = tiling
        if tiling.n > 1:
            # The global tiling object is a simple 1x1 partition
            global_tiling = parallel.Tiling(nrow=1, ncol=1, ncpus=1, **kwargs)
            global_tiling.set_extent(nx, ny)

        # If on master node (possibly only node), create global domain object
        if tiling.rank == 0:
            global_domain = Domain(nx, ny, nz, lon, lat, x, y, spherical, tiling=global_tiling, mask=mask, H=H, z0=z0, f=f, z=z, zo=zo, logger=logger)

        # If there is only one node, return the global domain immediately
        if tiling.n == 1:
            return global_domain

        # If on root, initialize the global domain [JB 2021-11-23 Why? Is this really needed?]
        if global_domain is not None:
            global_domain.initialize(runtype=runtype)

        # Create the subdomain, and (if on root) attach a pointer to the global domain
        subdomain = Domain.partition(tiling, nx, ny, nz, global_domain, runtype=runtype, has_xy=x is not None, has_lonlat=lon is not None, spherical=spherical, logger=logger)
        subdomain.glob = global_domain

        return subdomain

    def __init__(self, nx: int, ny: int, nz: int, lon: Optional[numpy.ndarray]=None, lat: Optional[numpy.ndarray]=None, x: Optional[numpy.ndarray]=None, y: Optional[numpy.ndarray]=None, spherical: bool=False, mask: Optional[numpy.ndarray]=1, H: Optional[numpy.ndarray]=None, z0: Optional[numpy.ndarray]=0., f: Optional[numpy.ndarray]=None, tiling: Optional[parallel.Tiling]=None, z: Optional[numpy.ndarray]=0., zo: Optional[numpy.ndarray]=0., logger: Optional[logging.Logger]=None, **kwargs):
        assert nx > 0, 'Number of x points is %i but must be > 0' % nx
        assert ny > 0, 'Number of y points is %i but must be > 0' % ny
        assert nz > 0, 'Number of z points is %i but must be > 0' % nz
        assert lat is not None or f is not None, 'Either lat of f must be provided to determine the Coriolis parameter.'

        self.root_logger = logger if logger is not None else parallel.getLogger()
        self.logger = self.root_logger.getChild('domain')
        self.field_manager: Optional[output.FieldManager] = None
        self.input_manager = input.InputManager()
        self.glob: Optional['Domain'] = self

        self.logger.info('Domain size (T grid): %i x %i (%i cells)' % (nx, ny, nx * ny))

        halo = 2

        shape = (2 * ny + 1, 2 * nx + 1)
        superhalo = 2 * halo
        shape_ = (shape[0] + 2 * superhalo, shape[1] + 2 * superhalo)

        # Set up subdomain partition information to enable halo exchanges
        if tiling is None:
            tiling = parallel.Tiling(**kwargs)
        self.tiling = tiling

        def setup_metric(source: Optional[numpy.ndarray]=None, optional: bool=False, fill_value=numpy.nan, relative_in_x: bool=False, relative_in_y: bool=False, dtype: numpy.typing.DTypeLike=float, writeable: bool=True) -> Tuple[Optional[numpy.ndarray], Optional[numpy.ndarray]]:
            if optional and source is None:
                return None, None
            data = numpy.full(shape_, fill_value, dtype)
            data_int = data[superhalo:-superhalo, superhalo:-superhalo]
            if source is not None:
                try:
                    # First try if data has been provided on the supergrid
                    data_int[...] = source
                except ValueError:
                    try:
                        # Now try if data has been provided on the X (corners) grid
                        source_on_X = numpy.broadcast_to(source, (ny + 1, nx + 1))
                        interfaces_to_supergrid_2d(source_on_X, out=data_int)
                    except ValueError:
                        raise Exception('Cannot array broadcast to supergrid (%i x %i) or X grid (%i x %i)' % ((ny * 2 + 1, nx * 2 + 1, ny + 1, nx + 1)))
                self.exchange_metric(data, relative_in_x, relative_in_y, fill_value=fill_value)
            data.flags.writeable = data_int.flags.writeable = writeable
            return data_int, data

        # Supergrid metrics (without underscore=interior only, with underscore=including halos)
        self.x, self.x_ = setup_metric(x, optional=True, relative_in_x=True, writeable=False)
        self.y, self.y_ = setup_metric(y, optional=True, relative_in_y=True, writeable=False)
        self.lon, self.lon_ = setup_metric(lon, optional=True, relative_in_x=True, writeable=False)
        self.lat, self.lat_ = setup_metric(lat, optional=True, relative_in_y=True, writeable=False)
        self.H, self.H_ = setup_metric(H)
        self.z, self.z_ = setup_metric(z)       # elevation
        self.zo, self.zo_ = setup_metric(zo)    # elevaton on previous time step
        self.z0b_min, self.z0b_min_ = setup_metric(z0)
        self.mask, self.mask_ = setup_metric(mask, dtype=numpy.intc, fill_value=0)

        cor = f if f is not None else 2. * omega * numpy.sin(deg2rad * lat)
        self.cor, self.cor_ = setup_metric(cor, writeable=False)

        # Compute dx, dy from Cartesian or spherical coordinates
        # These have had their halo exchanges as part of setup_metric and are therefore defined inside the halos too.
        self.dx, self.dx_ = setup_metric()
        self.dy, self.dy_ = setup_metric()
        if spherical:
            dlon_ = self.lon_[:, 2:] - self.lon_[:, :-2]
            dlat_ = self.lat_[2:, :] - self.lat_[:-2, :]
            self.dx_[:, 1:-1] = deg2rad * dlon_ * rearth * numpy.cos(deg2rad * self.lat_[:, 1:-1])
            self.dy_[1:-1, :] = deg2rad * dlat_ * rearth
        else:
            self.dx_[:, 1:-1] = self.x_[:, 2:] - self.x_[:, :-2]
            self.dy_[1:-1, :] = self.y_[2:, :] - self.y_[:-2, :]

        # Halo exchange for dx, dy, needed to ensure the outer strips of the halos are valid
        # Those outermost strips could not be computed by central-differencing the coordinates as that would require points outside the domain.
        self.exchange_metric(self.dx_)
        self.exchange_metric(self.dy_)

        self.dx_.flags.writeable = self.dx.flags.writeable = False
        self.dy_.flags.writeable = self.dy.flags.writeable = False

        self.rotation, self.rotation_ = setup_metric()
        def supergrid_rotation(x, y):
            # For each point, draw lines to the nearest neighbor (1/2 a grid cell) on the left, right, top and bottom.
            # Determinethe angle bey
            rotation_left = numpy.arctan2(y[:, 1:-1] - y[:, :-2], x[:, 1:-1] - x[:, :-2])
            rotation_right = numpy.arctan2(y[:, 2:] - y[:, 1:-1], x[:, 2:] - x[:, 1:-1])
            rotation_bot = numpy.arctan2(y[1:-1, :] - y[:-2, :], x[1:-1, :] - x[:-2, :]) - 0.5 * numpy.pi
            rotation_top = numpy.arctan2(y[2:, :] - y[1:-1, :], x[2:, :] - x[1:-1, :]) - 0.5 * numpy.pi
            x_dum = numpy.cos(rotation_left[1:-1,:]) + numpy.cos(rotation_right[1:-1,:]) + numpy.cos(rotation_bot[:,1:-1]) + numpy.cos(rotation_top[:,1:-1])
            y_dum = numpy.sin(rotation_left[1:-1,:]) + numpy.sin(rotation_right[1:-1,:]) + numpy.sin(rotation_bot[:,1:-1]) + numpy.sin(rotation_top[:,1:-1])
            return numpy.arctan2(y_dum, x_dum)
        if self.lon_ is not None and self.lat_ is not None:
            # Proper rotation with respect to true North
            self.rotation_[1:-1,1:-1] = supergrid_rotation(self.lon_, self.lat_)
        else:
            # Rotation with respect to y axis - assumes y axis always point to true North (can be valid on for infinitesimally small domain)
            self.rotation_[1:-1,1:-1] = supergrid_rotation(self.x_, self.y_)
        self.exchange_metric(self.rotation)
        self.rotation.flags.writeable = False

        self.area, self.area_ = setup_metric(self.dx * self.dy, writeable=False)

        self.spherical = spherical

        # Determine if we have simple coordinates (useful for xarray and plotting in general)
        self.lon_is_1d = self.lon is not None and (self.lon[:1, :] == self.lon[:, :]).all()
        self.lat_is_1d = self.lat is not None and (self.lat[:, :1] == self.lat[:, :]).all()
        self.x_is_1d = self.x is not None and (self.x[:1, :] == self.x[:, :]).all()
        self.y_is_1d = self.y is not None and (self.y[:, :1] == self.y[:, :]).all()

        self.imin, self.imax = 1, nx
        self.jmin, self.jmax = 1, ny
        self.kmin, self.kmax = 1, nz

        _pygetm.Domain.__init__(self, self.imin, self.imax, self.jmin, self.jmax, self.kmin, self.kmax)
        self.halo = self.halox
        self.shape = (nz, ny + 2 * self.halo, nx + 2 * self.halo)

        # Advection grids (two letters: first for advected quantity, second for advection direction)
        self.UU = Grid(self, _pygetm.UUGRID, ioffset=3, joffset=1)
        self.VV = Grid(self, _pygetm.VVGRID, ioffset=1, joffset=3)
        self.UV = Grid(self, _pygetm.UVGRID, ioffset=2, joffset=2)
        self.VU = Grid(self, _pygetm.VUGRID, ioffset=2, joffset=2)

        # Create grids
        self.U = Grid(self, _pygetm.UGRID, ioffset=2, joffset=1, ugrid=self.UU, vgrid=self.UV)
        self.V = Grid(self, _pygetm.VGRID, ioffset=1, joffset=2, ugrid=self.VU, vgrid=self.VV)
        self.T = Grid(self, _pygetm.TGRID, ioffset=1, joffset=1, ugrid=self.U, vgrid=self.V)
        self.X = Grid(self, _pygetm.XGRID, ioffset=0, joffset=0, overlap=1)

        self.Dmin = 1.

        self.initialized = False
        self.open_boundaries = {}
        self.rivers = Rivers(self.T)

    def add_open_boundary(self, side: int, l: int, mstart: int, mstop: int, type_2d: int, type_3d: int):
        """Note that l, mstart, mstop are 0-based indices of a T point in the global domain.
        mstop indicates the upper limit of the boundary - it is the first index that is EXcluded."""
        # NB below we convert to indices in the T grid of the current subdomain INCLUDING halos
        # We also limit the indices to the range valid for the current subdomain.
        HALO = 2
        if side in (WEST, EAST):
            l_offset, m_offset, l_max, m_max = self.tiling.xoffset, self.tiling.yoffset, self.T.nx_, self.T.ny_
        else:
            l_offset, m_offset, l_max, m_max = self.tiling.yoffset, self.tiling.xoffset, self.T.ny_, self.T.nx_
        l_loc = HALO + l - l_offset
        mstart_loc = min(max(0, HALO + mstart - m_offset), m_max)
        mstop_loc = min(max(0, HALO + mstop - m_offset), m_max)
        if l_loc >= 0 and l_loc < l_max and mstop_loc > mstart_loc:
            # Boundary lies at least partially within current subdomain
            if side in (WEST, EAST):
                mask = self.mask_[1 + 2 * mstart_loc:2 * mstop_loc:2, 1 + l_loc * 2]
            else:
                mask = self.mask_[1 + l_loc * 2, 1 + 2 * mstart_loc:2 * mstop_loc:2]
            if (mask == 0).any():
                self.logger.error('%i of %i points of this open boundary are on land' % ((mask == 0).sum(), mstop_loc - mstart_loc))
                raise Exception()
            m_skip = mstart_loc - (HALO + mstart - m_offset)
            assert m_skip >= 0
        else:
            # Boundary lies completely outside current subdomain. Record it anyway, so we can later set up a
            # global -> local map of open bounfary points
            l_loc, mstart_loc, mstop_loc, m_skip = None, None, None, mstop - mstart
        self.open_boundaries.setdefault(side, []).append((l_loc, mstart_loc, mstop_loc, m_skip, mstop - mstart, type_2d, type_3d))

    def initialize(self, runtype: int, field_manager: Optional[output.FieldManager]=None):
        assert not self.initialized, 'Domain has already been initialized'
        # Mask U,V,X points without any valid T neighbor - this mask will be maintained by the domain to be used for e.g. plotting
        tmask = self.mask_[1::2, 1::2]
        self.mask_[2:-2:2, 1::2][numpy.logical_and(tmask[1:, :] == 0, tmask[:-1, :] == 0)] = 0
        self.mask_[1::2, 2:-2:2][numpy.logical_and(tmask[:, 1:] == 0, tmask[:, :-1] == 0)] = 0
        self.mask_[2:-2:2, 2:-2:2][numpy.logical_and(numpy.logical_and(tmask[1:, 1:] == 0, tmask[:-1, 1:] == 0), numpy.logical_and(tmask[1:, :-1] == 0, tmask[:-1, :-1] == 0))] = 0
        self.exchange_metric(self.mask_, fill_value=0)

        if field_manager is not None:
            self.field_manager = field_manager
        self.field_manager = self.field_manager or output.FieldManager()

        HALO = 2
        nbdyp = 0
        nbdyp_glob = 0
        bdyinfo, bdy_i, bdy_j = [], [], []
        side2count = {}
        self.local_to_global_ob = []
        for side in (WEST, NORTH, EAST, SOUTH):
            bounds = self.open_boundaries.get(side, [])
            kept_bounds = []
            for l, mstart, mstop, mskip, len_glob, type_2d, type_3d in bounds:
                if l is not None:
                    # Note that bdyinfo needs indices into the T grid EXCLUDING halos
                    bdyinfo.append(numpy.array((l - HALO, mstart - HALO, mstop - HALO, type_2d, type_3d, nbdyp), dtype=numpy.intc))
                    nbdyp += mstop - mstart
                    if side == WEST:
                        self.mask_[1 + 2 * mstart:2 * mstop:2, 1 + l * 2] = 2
                        self.mask_[2 + 2 * mstart:2 * mstop:2, 1 + l * 2] = 3
                    elif side == EAST:
                        self.mask_[1 + 2 * mstart:2 * mstop:2, 1 + l * 2] = 2
                        self.mask_[2 + 2 * mstart:2 * mstop:2, 1 + l * 2] = 3
                    elif side == SOUTH:
                        self.mask_[1 + l * 2, 1 + 2 * mstart:2 * mstop:2] = 2
                        self.mask_[1 + l * 2, 2 + 2 * mstart:2 * mstop:2] = 3
                    elif side == NORTH:
                        self.mask_[1 + l * 2, 1 + 2 * mstart:2 * mstop:2] = 2
                        self.mask_[1 + l * 2, 2 + 2 * mstart:2 * mstop:2] = 3

                    if side in (WEST, EAST):
                        bdy_i.append(numpy.repeat(l, mstop - mstart))
                        bdy_j.append(numpy.arange(mstart, mstop))
                    else:
                        bdy_i.append(numpy.arange(mstart, mstop))
                        bdy_j.append(numpy.repeat(l, mstop - mstart))
                    if self.local_to_global_ob and self.local_to_global_ob[-1][1] == nbdyp_glob + mskip:
                        # attach to previous boundary
                        self.local_to_global_ob[-1][1] += mstop - mstart
                    else:
                        # gap; add new slice
                        self.local_to_global_ob.append([nbdyp_glob + mskip, nbdyp_glob + mskip + mstop - mstart])
                    kept_bounds.append(bounds)
                nbdyp_glob += len_glob
            side2count[side] = len(kept_bounds)
            if not kept_bounds:
                # No open boundaries on this side fall within the current subdomain
                # Delete any reference to that side from our open_boundaries list.
                self.open_boundaries.pop(side, None)
            else:
                # Replace the list of open boundaries on this side with only those that fall within the current subdomain
                self.open_boundaries[side] = kept_bounds
        self.bdy_i = numpy.empty((0,), dtype=numpy.intc) if nbdyp == 0 else numpy.concatenate(bdy_i, dtype=numpy.intc)
        self.bdy_j = numpy.empty((0,), dtype=numpy.intc) if nbdyp == 0 else numpy.concatenate(bdy_j, dtype=numpy.intc)
        self.nbdyp_glob = nbdyp_glob
        self.logger.info('%i open boundaries (%i West, %i North, %i East, %i South)' % (len(bdyinfo), side2count[WEST], side2count[NORTH], side2count[EAST], side2count[SOUTH]))
        if nbdyp > 0:
            if nbdyp == nbdyp_glob:
                assert len(self.local_to_global_ob) == 1 and self.local_to_global_ob[0][0] == 0 and self.local_to_global_ob[0][1] == nbdyp_glob
                self.local_to_global_ob = None
            else:
                self.logger.info('global-to-local open boundary map: %s' % (self.local_to_global_ob,))
            bdyinfo = numpy.stack(bdyinfo, axis=-1)
            self.initialize_open_boundaries(nwb=side2count[WEST], nnb=side2count[NORTH], neb=side2count[EAST], nsb=side2count[SOUTH], nbdyp=nbdyp, bdy_i=self.bdy_i - HALO, bdy_j=self.bdy_j - HALO, bdy_info=bdyinfo)

        # Mask U,V,X points unless all their T neighbors are valid - this mask will be sent to Fortran and determine which points are computed
        mask_ = numpy.array(self.mask_, copy=True)
        tmask = mask_[1::2, 1::2]
        mask_[2:-2:2, 1::2][numpy.logical_or(tmask[1:, :] == 0, tmask[:-1, :] == 0)] = 0
        mask_[1::2, 2:-2:2][numpy.logical_or(tmask[:, 1:] == 0, tmask[:, :-1] == 0)] = 0
        mask_[2:-2:2, 2:-2:2][numpy.logical_or(numpy.logical_or(tmask[1:, 1:] == 0, tmask[:-1, 1:] == 0), numpy.logical_or(tmask[1:, :-1] == 0, tmask[:-1, :-1] == 0))] = 0
        self.exchange_metric(mask_, fill_value=0)
        self.mask_[...] = mask_

        for grid in self.grids.values():
            grid.initialize(nbdyp)
        self.UU.mask.all_values.fill(0)
        self.UV.mask.all_values.fill(0)
        self.VU.mask.all_values.fill(0)
        self.VV.mask.all_values.fill(0)
        self.UU.mask.all_values[:, :-1][numpy.logical_and(self.U.mask.all_values[:, :-1], self.U.mask.all_values[:,1:])] = 1
        self.UV.mask.all_values[:-1, :][numpy.logical_and(self.U.mask.all_values[:-1, :], self.U.mask.all_values[1:,:])] = 1
        self.VU.mask.all_values[:, :-1][numpy.logical_and(self.V.mask.all_values[:, :-1], self.V.mask.all_values[:,1:])] = 1
        self.VV.mask.all_values[:-1, :][numpy.logical_and(self.V.mask.all_values[:-1, :], self.V.mask.all_values[1:,:])] = 1

        self.zc_bdy = self.T.array(z = CENTERS, on_boundary=True)
        self.zf_bdy = self.T.array(z = INTERFACES, on_boundary=True)

        self.logger.info('Number of unmasked points excluding halos: %i on T grid, %i on U grid, %i on V grid, %i on X grid' % ((self.T.mask.values > 0).sum(), (self.U.mask.values > 0).sum(), (self.V.mask.values > 0).sum(), (self.X.mask.values > 0).sum()))

        self.H_.flags.writeable = self.H.flags.writeable = False
        self.z0b_min_.flags.writeable = self.z0b_min.flags.writeable = False
        self.mask_.flags.writeable = self.mask.flags.writeable = False
        self.z_.flags.writeable = self.z.flags.writeable = False
        self.zo_.flags.writeable = self.zo.flags.writeable = False

        _pygetm.Domain.initialize(self, runtype, Dmin=self.Dmin)

        self.rivers.initialize()

        self.initialized = True

    def set_bathymetry(self, depth, scale_factor=None):
        assert not self.initialized, 'set_bathymetry cannot be called after the domain has been initialized.'
        if not isinstance(depth, xarray.DataArray):
            # Depth is provided as raw data and therefore must be already on the supergrid
            self.H[...] = depth
        else:
            # Depth is provided as xarray object that includes coordinates (we require CF compliant longitude, latitude)
            # Interpolate to target grid.
            depth = input.limit_region(depth, self.lon.min(), self.lon.max(), self.lat.min(), self.lat.max())
            depth = input.horizontal_interpolation(depth, self.lon, self.lat)            
            self.H[...] = depth.values
        if scale_factor is not None:
            self.H *= scale_factor

    def mask_shallow(self, minimum_depth: float):
        self.mask[self.H < minimum_depth] = 0

    def mask_rectangle(self, xmin: Optional[float]=None, xmax: Optional[float]=None, ymin: Optional[float]=None, ymax: Optional[float]=None, value: int=0):
        assert not self.initialized, 'adjust_mask cannot be called after the domain has been initialized.'
        selected = numpy.ones(self.mask.shape, dtype=bool)
        x, y = (self.lon, self.lat) if self.spherical else (self.x, self.y)
        if xmin is not None: selected = numpy.logical_and(selected, x >= xmin)
        if xmax is not None: selected = numpy.logical_and(selected, x <= xmax)
        if ymin is not None: selected = numpy.logical_and(selected, y >= ymin)
        if ymax is not None: selected = numpy.logical_and(selected, y <= ymax)
        self.mask[selected] = value

    def plot(self, fig=None, show_H: bool=True, show_mesh: bool=True, show_rivers: bool=True, editable: bool=False):
        import matplotlib.pyplot
        import matplotlib.collections
        import matplotlib.widgets
        if fig is None:
            fig, ax = matplotlib.pyplot.subplots(figsize=(0.15 * self.nx, 0.15 * self.ny))
        else:
            ax = fig.gca()

        x, y = (self.lon, self.lat) if self.spherical else (self.x, self.y)

        local_slice, _, _, _ = self.tiling.subdomain2slices(halo_sub=0, halo_glob=4, scale=2, share=1, exclude_global_halos=True)
        if show_H:
            import cmocean
            cm = cmocean.cm.deep
            cm.set_bad('gray')
            c = ax.pcolormesh(x[local_slice], y[local_slice], numpy.ma.array(self.H[local_slice], mask=self.mask[local_slice]==0), alpha=0.5 if show_mesh else 1, shading='auto', cmap=cm)
            #c = ax.contourf(x, y, numpy.ma.array(self.H, mask=self.mask==0), 20, alpha=0.5 if show_mesh else 1)
            cb = fig.colorbar(c)
            cb.set_label('undisturbed water depth (m)')

        if show_rivers:
            for river in self.rivers.values():
                iloc, jloc = 1 + river.i * 2, 1 + river.j * 2
                lon = self.lon_[jloc, iloc]
                lat = self.lat_[jloc, iloc]
                ax.plot([lon], [lat], '.r')
                ax.text(lon, lat, river.name, color='r')

        def plot_mesh(ax, x, y, **kwargs):
            segs1 = numpy.stack((x, y), axis=2)
            segs2 = segs1.transpose(1, 0, 2)
            ax.add_collection(matplotlib.collections.LineCollection(segs1, **kwargs))
            ax.add_collection(matplotlib.collections.LineCollection(segs2, **kwargs))

        if show_mesh:
            plot_mesh(ax, x[::2, ::2], y[::2, ::2], colors='k', linestyle='-', linewidth=.3)
            #ax.pcolor(x[1::2, 1::2], y[1::2, 1::2], numpy.ma.array(x[1::2, 1::2], mask=True), edgecolors='k', linestyles='--', linewidth=.2)
            #pc = ax.pcolormesh(x[1::2, 1::2], y[1::2, 1::2],  numpy.ma.array(x[1::2, 1::2], mask=True), edgecolor='gray', linestyles='--', linewidth=.2)
            ax.plot(x[::2, ::2], y[::2, ::2], '.k', markersize=3.)
            ax.plot(x[1::2, 1::2], y[1::2, 1::2], 'xk', markersize=2.5)
        ax.set_xlabel('longitude (degrees East)' if self.spherical else 'x (m)')
        ax.set_ylabel('latitude (degrees North)' if self.spherical else 'y (m)')
        if not self.spherical:
            ax.axis('equal')
        xmin, xmax = numpy.nanmin(x), numpy.nanmax(x)
        ymin, ymax = numpy.nanmin(y), numpy.nanmax(y)
        xmargin = 0.05 * (xmax - xmin)
        ymargin = 0.05 * (ymax - ymin)
        ax.set_xlim(xmin - xmargin, xmax + xmargin)
        ax.set_ylim(ymin - ymargin, ymax + ymargin)

        def on_select(eclick, erelease):
            xmin, xmax = min(eclick.xdata, erelease.xdata), max(eclick.xdata, erelease.xdata)
            ymin, ymax = min(eclick.ydata, erelease.ydata), max(eclick.ydata, erelease.ydata)
            self.mask_rectangle(xmin, xmax, ymin, ymax)
            c.set_array(numpy.ma.array(self.H, mask=self.mask==0).ravel())
            fig.canvas.draw()
            #self.sel.set_active(False)
            #self.sel = None
            #ax.draw()
            #fig.clf()
            #self.plot(fig=fig, show_mesh=show_mesh)
        if editable:
            self.sel = matplotlib.widgets.RectangleSelector(
                    ax, on_select,
                    useblit=True,
                    button=[1],
                    interactive=False)
        return fig

    def save(self, path: str, full: bool=False):
        with netCDF4.Dataset(path, 'w') as nc:
            def create(name, units, long_name, values, coordinates: str, dimensions=('y', 'x')):
                fill_value = None
                if numpy.ma.getmask(values) is not numpy.ma.nomask:
                    fill_value = values.fill_value
                ncvar = nc.createVariable(name, values.dtype, dimensions, fill_value=fill_value)
                ncvar.units = units
                ncvar.long_name = long_name
                #ncvar.coordinates = coordinates
                ncvar[...] = values

            def create_var(name, units, long_name, values, values_):
                if values is None:
                    return
                create(name, units, long_name, values, coordinates='lon lat' if self.spherical else 'x y')
                if full: create(name + '_', units, long_name, values_, dimensions=('y_', 'x_'), coordinates='lon_ lat_' if self.spherical else 'x_ y_')

            nc.createDimension('x', self.H.shape[1])
            nc.createDimension('y', self.H.shape[0])
            if full:
                nc.createDimension('x_', self.H_.shape[1])
                nc.createDimension('y_', self.H_.shape[0])
                create_var('dx', 'm', 'dx', self.dx, self.dx_)
                create_var('dy', 'm', 'dy', self.dy, self.dy_)
                create_var('area', 'm2', 'area', self.dx * self.dy, self.dx_ * self.dy_)
            create_var('lat', 'degrees_north', 'latitude', self.lat, self.lat_)
            create_var('lon', 'degrees_east', 'longitude', self.lon, self.lon_)
            create_var('x', 'm', 'x', self.x, self.x_)
            create_var('y', 'm', 'y', self.y, self.y_)
            create_var('H', 'm', 'undisturbed water depth', self.H, self.H_)
            create_var('mask', '', 'mask', self.mask, self.mask_)

    def save_grids(self, path: str):
        with netCDF4.Dataset(path, 'w') as nc:
            self.T.add_to_netcdf(nc)
            self.U.add_to_netcdf(nc, postfix='u')
            self.V.add_to_netcdf(nc, postfix='v')
            self.X.add_to_netcdf(nc, postfix='x')

    def contains(self, x: float, y: float):
        allx, ally = (self.lon, self.lat) if self.spherical else (self.x, self.y)
        ny, nx = allx.shape

        # Determine whether point falls within current subdomain
        # based on https://wrf.ecse.rpi.edu/Research/Short_Notes/pnpoly.html
        x_bnd = numpy.concatenate((allx[0, :-1], allx[:-1, -1], allx[-1, nx-1:0:-1], allx[ny-1:0:-1, 0]))
        y_bnd = numpy.concatenate((ally[0, :-1], ally[:-1, -1], ally[-1, nx-1:0:-1], ally[ny-1:0:-1, 0]))
        assert x_bnd.size == 2 * ny + 2 * nx - 4
        inside = False
        for i, (vertxi, vertyi) in enumerate(zip(x_bnd, y_bnd)):
            vertxj, vertyj = x_bnd[i - 1], y_bnd[i - 1]
            if (vertyi > y) != (vertyj > y) and x < (vertxj - vertxi) * (y - vertyi) / (vertyj - vertyi) + vertxi:
                inside = not inside
        return inside