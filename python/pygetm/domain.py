from typing import Optional, Tuple
import operator
import logging
import os.path

import numpy
import numpy.typing
import xarray
import netCDF4

from . import _pygetm
from . import core
from . import parallel
from . import output
from . import input

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
    _fortran_arrays = _coordinate_arrays + ('dx', 'dy', 'dlon', 'dlat', 'H', 'D', 'mask', 'z', 'zo', 'area', 'iarea', 'cor', 'ho', 'hn', 'zc', 'z0b', 'z0b_min')
    _all_fortran_arrays = tuple(['_%s' % n for n in _fortran_arrays] + ['_%si' % n for n in _coordinate_arrays] + ['_%si_' % n for n in _coordinate_arrays])
    __slots__ = _all_fortran_arrays + ('halo', 'type', 'ioffset', 'joffset', 'postfix', 'xypostfix', 'zpostfix', 'ugrid', 'vgrid', 'wgrid', '_sin_rot', '_cos_rot', 'rotation', 'nbdyp')

    def __init__(self, domain: 'Domain', grid_type: int, ioffset: int, joffset: int, ugrid: Optional['Grid']=None, vgrid: Optional['Grid']=None, wgrid: Optional['Grid']=None):
        _pygetm.Grid.__init__(self, domain, grid_type)
        self.halo = domain.halo
        self.type = grid_type
        self.ioffset = ioffset
        self.joffset = joffset
        self.postfix = {_pygetm.TGRID: 't', _pygetm.UGRID: 'u', _pygetm.VGRID: 'v', _pygetm.XGRID: 'x', _pygetm.WGRID: 'w', _pygetm.UUGRID: '_uu_adv', _pygetm.VVGRID: '_vv_adv', _pygetm.UVGRID: '_uv_adv', _pygetm.VUGRID: '_vu_adv'}[grid_type]
        self.xypostfix = {_pygetm.WGRID: 't'}.get(grid_type, self.postfix)
        self.zpostfix = {_pygetm.WGRID: 'w'}.get(grid_type, '')
        self.ugrid: Optional[Grid] = ugrid
        self.vgrid: Optional[Grid] = vgrid
        self.wgrid: Optional[Grid] = wgrid
        self._sin_rot: Optional[numpy.ndarray] = None
        self._cos_rot: Optional[numpy.ndarray] = None

    def initialize(self, nbdyp: int):
        for name in self._fortran_arrays:
            setattr(self, '_%s' % name, self.wrap(core.Array(name=name + self.postfix), name.encode('ascii')))
        self.rotation = core.Array.create(grid=self, dtype=self.x.dtype, name='rotation' + self.postfix, units='rad', long_name='grid rotation with respect to true North')
        self.fill()
        self.z0b.all_values[...] = self.z0b_min.all_values
        self.nbdyp = nbdyp

    def fill(self):
        read_only = ('dx', 'dy', 'lon', 'lat', 'x', 'y', 'cor', 'area', 'rotation', 'z0b_min')
        nj, ni = self.H.all_values.shape
        has_bounds = self.ioffset > 0 and self.joffset > 0 and self.domain.H_.shape[-1] >= self.ioffset + 2 * ni and self.domain.H_.shape[-2] >= self.joffset + 2 * nj
        for name in ('H', 'mask', 'dx', 'dy', 'lon', 'lat', 'x', 'y', 'cor', 'area', 'z', 'zo', 'rotation', 'z0b_min'):
            source = getattr(self.domain, name + '_')
            if source is not None:
                target = getattr(self, name).all_values
                values = source[self.joffset:self.joffset + 2 * nj:2, self.ioffset:self.ioffset + 2 * ni:2]
                target.fill(numpy.nan)
                target[:values.shape[0], :values.shape[1]] = values
                target.flags.writeable = name not in read_only
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

    def array(self, fill=None, is_3d: bool=False, dtype: numpy.typing.DTypeLike=float, **kwargs) -> core.Array:
        return core.Array.create(self, fill, is_3d, dtype, **kwargs)

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

for membername in Grid._all_fortran_arrays:
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

class Domain(_pygetm.Domain):
    @staticmethod
    def create_cartesian(x, y, nz: int, interfaces=False, **kwargs) -> 'Domain':
        assert x.ndim == 1, 'x coordinate must be one-dimensional'
        assert y.ndim == 1, 'y coordinate must be one-dimensional'

        if interfaces:
            nx, ny = x.size - 1, y.size - 1
            x, y = interfaces_to_supergrid_1d(x), interfaces_to_supergrid_1d(y)
        else:
            nx, ny = x.size, y.size
            x, y = center_to_supergrid_1d(x), center_to_supergrid_1d(y)
        return Domain.create(nx, ny, nz, x=x, y=y[:, numpy.newaxis], **kwargs)

    @staticmethod
    def partition(tiling, nx: int, ny: int, nz: int, global_domain: Optional['Domain'], halo: int=2, has_xy: bool=True, has_lonlat: bool=True, logger: Optional[logging.Logger]=None, **kwargs):
        assert nx == tiling.nx_glob and ny == tiling.ny_glob, 'Extent of global domain (%i, %i) does not match that of tiling (%i, %i).' % (ny, nx, tiling.ny_glob, tiling.nx_glob)
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
            scatterer = parallel.Scatter(tiling, c, halo=0, share=share, scale=2)
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
        self.W = Grid(self, _pygetm.WGRID, ioffset=1, joffset=1)
        self.T = Grid(self, _pygetm.TGRID, ioffset=1, joffset=1, ugrid=self.U, vgrid=self.V, wgrid=self.W)
        self.X = Grid(self, _pygetm.XGRID, ioffset=0, joffset=0)

        self.Dmin = 1.

        self.initialized = False
        self.open_boundaries = {}

    def add_open_boundary(self, side: int, l: int, mstart: int, mstop: int, type_2d: int, type_3d: int):
        """Note that l, mstart, mstop are 0-based indices in the global domain.
        mstop indicates the upper limit of the boundary - it is the first index that is EXcluded."""
        # NB below we convert to indices in the T grid of the current subdomain INCLUDING halos
        # We also limit the indices to the range valid for the current subdomain.
        HALO = 2
        if side in (WEST, EAST):
            l = HALO + l - self.tiling.xoffset
            mstart = min(max(0, HALO + mstart - self.tiling.yoffset), self.T.ny_)
            mstop = min(max(0, HALO + mstop - self.tiling.yoffset), self.T.ny_)
            lmax = self.T.nx_
        else:
            l = HALO + l - self.tiling.yoffset
            mstart = min(max(0, HALO + mstart - self.tiling.xoffset), self.T.nx_)
            mstop = min(max(0, HALO + mstop - self.tiling.xoffset), self.T.nx_)
            lmax = self.T.ny_
        if l >= 0 and l < lmax and mstop > mstart:
            # Boundary lies at least partially within current subdomain
            if side in (WEST, EAST):
                mask = self.mask_[1 + 2 * mstart:2 * mstop:2, 1 + l * 2]
            else:
                mask = self.mask_[1 + l * 2, 1 + 2 * mstart:2 * mstop:2]
            if (mask == 0).any():
                self.logger.error('%i of %i points of this open boundary are on land' % ((mask == 0).sum(), mstop - mstart))
                raise Exception()
            self.open_boundaries.setdefault(side, []).append((l, mstart, mstop, type_2d, type_3d))

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
        bdyinfo, bdy_i, bdy_j = [], [], []
        side2count = {}
        for side in (WEST, NORTH, EAST, SOUTH):
            bounds = self.open_boundaries.get(side, [])
            side2count[side] = len(bounds)
            for l, mstart, mstop, type_2d, type_3d in bounds:
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
        self.bdy_i = numpy.empty((0,), dtype=numpy.intc) if nbdyp == 0 else numpy.concatenate(bdy_i, dtype=numpy.intc)
        self.bdy_j = numpy.empty((0,), dtype=numpy.intc) if nbdyp == 0 else numpy.concatenate(bdy_j, dtype=numpy.intc)
        self.logger.info('%i open boundaries (%i West, %i North, %i East, %i South)' % (len(bdyinfo), side2count[WEST], side2count[NORTH], side2count[EAST], side2count[SOUTH]))
        if nbdyp > 0:
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
        self.UU.mask.all_values[...] = 0
        self.UV.mask.all_values[...] = 0
        self.VU.mask.all_values[...] = 0
        self.VV.mask.all_values[...] = 0
        self.UU.mask.all_values[:, :-1][numpy.logical_and(self.U.mask.all_values[:, :-1], self.U.mask.all_values[:,1:])] = 1
        self.UV.mask.all_values[:-1, :][numpy.logical_and(self.U.mask.all_values[:-1, :], self.U.mask.all_values[1:,:])] = 1
        self.VU.mask.all_values[:, :-1][numpy.logical_and(self.V.mask.all_values[:, :-1], self.V.mask.all_values[:,1:])] = 1
        self.VV.mask.all_values[:-1, :][numpy.logical_and(self.V.mask.all_values[:-1, :], self.V.mask.all_values[1:,:])] = 1

        self.logger.info('Number of unmasked points excluding halos: %i on T grid, %i on U grid, %i on V grid, %i on X grid' % ((self.T.mask.values > 0).sum(), (self.U.mask.values > 0).sum(), (self.V.mask.values > 0).sum(), (self.X.mask.values > 0).sum()))

        self.H_.flags.writeable = self.H.flags.writeable = False
        self.z0b_min_.flags.writeable = self.z0b_min.flags.writeable = False
        self.mask_.flags.writeable = self.mask.flags.writeable = False
        self.z_.flags.writeable = self.z.flags.writeable = False
        self.zo_.flags.writeable = self.zo.flags.writeable = False

        _pygetm.Domain.initialize(self, runtype, Dmin=self.Dmin)

        self.initialized = True

    def set_bathymetry(self, depth, scale_factor=None, minimum_depth=None):
        if not isinstance(depth, xarray.DataArray):
            # Depth is provided as raw data and therefore must be already on the supergrid
            self.H[...] = depth
        else:
            # Depth is provided as xarray object that includes coordinates (we require CF compliant longitude, latitude)
            # Interpolate to target grid.
            self.H[...] = input.SpatialInterpolation(input.Variable(depth), self.lon, self.lat).x.values
        if scale_factor is not None:
            self.H *= scale_factor
        if minimum_depth is not None:
            self.H = numpy.ma.masked_less(self.H, minimum_depth)

    def distribute(self, field):
        return self.tiling.wrap(field, halo=self.halo)

    def plot(self, fig=None, show_H: bool=True):
        import matplotlib.pyplot
        import matplotlib.collections
        if fig is None:
            fig, ax = matplotlib.pyplot.subplots(figsize=(0.15 * self.nx, 0.15 * self.ny))
        else:
            ax = fig.gca()
        x, y = (self.lon, self.lat) if self.spherical else (self.x, self.y)
        if show_H:
            c = ax.contourf(x, y, self.H, 20, alpha=0.5)
            cb = fig.colorbar(c)
            cb.set_label('undisturbed water depth (m)')

        def plot_mesh(ax, x, y, **kwargs):
            segs1 = numpy.stack((x, y), axis=2)
            segs2 = segs1.transpose(1, 0, 2)
            ax.add_collection(matplotlib.collections.LineCollection(segs1, **kwargs))
            ax.add_collection(matplotlib.collections.LineCollection(segs2, **kwargs))

        plot_mesh(ax, x[::2, ::2], y[::2, ::2], colors='k', linestyle='-', linewidth=.3)
        #ax.pcolor(x[1::2, 1::2], y[1::2, 1::2], numpy.ma.array(x[1::2, 1::2], mask=True), edgecolors='k', linestyles='--', linewidth=.2)
        #pc = ax.pcolormesh(x[1::2, 1::2], y[1::2, 1::2],  numpy.ma.array(x[1::2, 1::2], mask=True), edgecolor='gray', linestyles='--', linewidth=.2)
        ax.plot(x[::2, ::2], y[::2, ::2], '.k', markersize=3.)
        ax.plot(x[1::2, 1::2], y[1::2, 1::2], 'xk', markersize=2.5)
        ax.set_xlabel('longitude (degrees East)' if self.spherical else 'x (m)')
        ax.set_ylabel('latitude (degrees North)' if self.spherical else 'y (m)')
        if not self.spherical:
            ax.axis('equal')

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
