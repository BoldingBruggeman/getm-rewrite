from typing import Optional, Tuple

import numpy
import numpy.typing
import xarray
import netCDF4

from . import _pygetm, core, parallel

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

class Grid(_pygetm.Grid, core.FortranObject):
    def __init__(self, domain: 'Domain', grid_type: bytes, ioffset: int, joffset: int):
        _pygetm.Grid.__init__(self, domain, ord(grid_type))
        self.halo = domain.halo
        self.domain = domain
        self.type = grid_type
        self.ioffset = ioffset
        self.joffset = joffset

    def initialize(self):
        self.get_arrays(('x', 'y', 'dx', 'dy', 'lon', 'lat', 'dlon', 'dlat', 'H', 'D', 'mask', 'z', 'zo', 'area', 'iarea', 'cor'), dtypes={'mask': int}, halo=self.halo)
        self.fill()

    def fill(self):
        read_only = ('dx', 'dy', 'lon', 'lat', 'x', 'y', 'cor', 'area')
        nj, ni = self.H_.shape
        has_bounds = self.ioffset > 0 and self.joffset > 0 and self.domain.H_.shape[-1] >= self.ioffset + 2 * ni and self.domain.H_.shape[-2] >= self.joffset + 2 * nj
        for name in ('H', 'mask', 'dx', 'dy', 'lon', 'lat', 'x', 'y', 'cor', 'area', 'z', 'zo'):
            source = getattr(self.domain, name + '_')
            if source is not None:
                target = getattr(self, name + '_')
                target[:, :] = source[self.joffset:self.joffset + 2 * nj:2, self.ioffset:self.ioffset + 2 * ni:2]
                target.flags.writeable = name not in read_only
                if has_bounds:
                    # Generate interface coordinates. These are not represented in Fortran as they are only needed for plotting.
                    # The interface coordinates are slices that point to the supergrid data; they thus do not consume additional memory.
                    values_i = source[self.joffset - 1:self.joffset + 2 * nj + 1:2, self.ioffset - 1:self.ioffset + 2 * ni + 1:2]
                    setattr(self, name + 'i_', values_i)
                    setattr(self, name + 'i', values_i[self.halo:-self.halo, self.halo:-self.halo])
        self.iarea_[:, :] = 1. / self.area_[:, :]

    def array(self, fill=None, dtype=float):
        data = numpy.empty(self.H_.shape, dtype=dtype)
        if fill is not None:
            data[...] = fill
        return data[self.halo:-self.halo, self.halo:-self.halo], data

    def as_xarray(self, data):
        assert data.shape == self.x.shape
        x = xarray.DataArray(self.c1, dims=['x',])
        y = xarray.DataArray(self.c2, dims=['y',])
        x.attrs['units'] = 'm'
        y.attrs['units'] = 'm'
        return xarray.DataArray(data, coords={'x': x, 'y': y}, dims=['y', 'x'])

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

def interfaces_to_supergrid_1d(data) -> numpy.ndarray:
    assert data.ndim == 1, 'data must be one-dimensional'
    assert data.size > 1, 'data must have at least 2 elements'
    data_sup = numpy.empty((data.size * 2 - 1,))
    data_sup[0::2] = data
    data_sup[1::2] = 0.5 * (data[1:] + data[:-1])
    return data_sup

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
        domain = Domain(nx, ny, nz, x=x, y=y[:, numpy.newaxis], **kwargs)
        return domain

    @staticmethod
    def partition(tiling, nx: int, ny: int, nlev: int, global_domain, runtype):
        assert nx % tiling.ncol == 0
        assert ny % tiling.nrow == 0
        assert global_domain is None or global_domain.initialized

        domain = Domain(1, nlev, 1, ny // tiling.nrow, 1, nx // tiling.ncol, tiling=tiling)

        # Scatter coordinates, bathymetry and mask to subdomains
        domain.distribute(domain.T.x_).scatter(None if global_domain is None else global_domain.T.x)
        domain.distribute(domain.T.y_).scatter(None if global_domain is None else global_domain.T.y)
        domain.distribute(domain.T.H_).scatter(None if global_domain is None else global_domain.T.H)
        domain.distribute(domain.T.mask_).scatter(None if global_domain is None else global_domain.T.mask)

        # Extract 1D coordinate vectors per subdomain
        domain.T.c1_[:] = domain.T.x_[domain.halo, :]
        domain.T.c2_[:] = domain.T.y_[:, domain.halo]

        domain.initialize(runtype)

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
        data[:superhalo, :] = data_ext[:superhalo, 1:-1]
        data[-superhalo:, :] = data_ext[-superhalo:, 1:-1]
        data[:, :superhalo] = data_ext[1:-1, :superhalo]
        data[:, -superhalo:] = data_ext[1:-1, -superhalo:]

    def __init__(self, nx: int, ny: int, nz: int, lon: Optional[numpy.ndarray]=None, lat: Optional[numpy.ndarray]=None, x: Optional[numpy.ndarray]=None, y: Optional[numpy.ndarray]=None, spherical: bool=False, mask: Optional[numpy.ndarray]=1, H: Optional[numpy.ndarray]=None, z0: Optional[numpy.ndarray]=None, f: Optional[numpy.ndarray]=None, tiling: Optional[parallel.Tiling]=None, z: Optional[numpy.ndarray]=0., zo: Optional[numpy.ndarray]=0., **kwargs):
        assert nx > 0, 'Number of x points is %i but must be > 0' % nx
        assert ny > 0, 'Number of y points is %i but must be > 0' % ny
        assert nz > 0, 'Number of z points is %i but must be > 0' % nz
        assert lat is not None or f is not None, 'Either lat of f must be provided to determine the Coriolis parameter.'

        halo = 2

        shape = (2 * ny + 1, 2 * nx + 1)
        superhalo = 2 * halo
        shape_ = (shape[0] + 2 * superhalo, shape[1] + 2 * superhalo)

        # Set up subdomain partition information to enable halo exchanges
        if tiling is None and kwargs:
            tiling = parallel.Tiling(1, 1, **kwargs)
        self.tiling = tiling

        def setup_metric(source: Optional[numpy.ndarray]=None, optional: bool=False, fill_value=numpy.nan, relative_in_x: bool=False, relative_in_y: bool=False, dtype: numpy.typing.DTypeLike=float, writeable: bool=True) -> Tuple[Optional[numpy.ndarray], Optional[numpy.ndarray]]:
            if optional and source is None:
                return None, None
            data = numpy.full(shape_, fill_value, dtype)
            data_int = data[superhalo:-superhalo, superhalo:-superhalo]
            if source is not None:
                data_int[...] = source
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
        self.z0, self.z0_ = setup_metric(z0)
        self.mask, self.mask_ = setup_metric(mask, dtype=int, fill_value=0)

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

        # Create grids
        self.T = Grid(self, b'T', ioffset=1, joffset=1)
        self.U = Grid(self, b'U', ioffset=2, joffset=1)
        self.V = Grid(self, b'V', ioffset=1, joffset=2)
        self.X = Grid(self, b'X', ioffset=0, joffset=0)

        self.initialized = False
        self.open_boundaries = {}

    def add_open_boundary(self, side: int, l: int, mstart: int, mstop: int, type_2d: int, type_3d: int):
        self.open_boundaries.setdefault(side, []).append((l, mstart, mstop, type_2d, type_3d))

    def initialize(self, runtype):
        assert not self.initialized, 'Domain has already been initialized'
        # Mask U,V,X points without any valid T neighbor - this mask will be maintained by the domain to be used for e.g. plotting
        tmask = self.mask_[1::2, 1::2]
        self.mask_[2:-2:2, 1::2][numpy.logical_and(tmask[1:, :] == 0, tmask[:-1, :] == 0)] = 0
        self.mask_[1::2, 2:-2:2][numpy.logical_and(tmask[:, 1:] == 0, tmask[:, :-1] == 0)] = 0
        self.mask_[2:-2:2, 2:-2:2][numpy.logical_and(numpy.logical_and(tmask[1:, 1:] == 0, tmask[:-1, 1:] == 0), numpy.logical_and(tmask[1:, :-1] == 0, tmask[:-1, :-1] == 0))] = 0
        self.exchange_metric(self.mask_)

        nbdyp = 0
        bdy_i, bdy_j = [], []
        for side in (WEST, NORTH, EAST, SOUTH):
            bounds = self.open_boundaries.get(side, [])
            for l, mstart, mstop, type_2d, type_3d in bounds:
                nbdyp += mstop - mstart + 1
                if side == WEST:
                    self.mask[-1 + 2 * mstart:2 * mstop:2, l * 2 - 1] = 2
                    self.mask[2 * mstart:2 * mstop:2, l * 2 - 1] = 3
                elif side == EAST:
                    self.mask[-1 + 2 * mstart:2 * mstop:2, l * 2 - 1] = 2
                    self.mask[2 * mstart:2 * mstop:2, l * 2 - 1] = 3
                elif side == SOUTH:
                    self.mask[l * 2 - 1, -1 + 2 * mstart:2 * mstop:2] = 2
                    self.mask[l * 2 - 1, 2 * mstart:2 * mstop:2] = 3
                elif side == NORTH:
                    self.mask[l * 2 - 1, -1 + 2 * mstart:2 * mstop:2] = 2
                    self.mask[l * 2 - 1, 2 * mstart:2 * mstop:2] = 3

                if side in (WEST, EAST):
                    bdy_i.append(numpy.repeat(l - 1, mstop - mstart + 1))
                    bdy_j.append(numpy.arange(mstart - 1, mstop))
                else:
                    bdy_i.append(numpy.arange(mstart - 1, mstop))
                    bdy_j.append(numpy.repeat(l - 1, mstop - mstart + 1))
        self.bdy_i = [] if nbdyp == 0 else numpy.concatenate(bdy_i)
        self.bdy_j = [] if nbdyp == 0 else numpy.concatenate(bdy_j)

        # Mask U,V,X points unless all their T neighbors are valid - this mask will be sent to Fortran and determine which points are computed
        mask_ = numpy.array(self.mask_, copy=True)
        tmask = mask_[1::2, 1::2]
        mask_[2:-2:2, 1::2][numpy.logical_or(tmask[1:, :] == 0, tmask[:-1, :] == 0)] = 0
        mask_[1::2, 2:-2:2][numpy.logical_or(tmask[:, 1:] == 0, tmask[:, :-1] == 0)] = 0
        mask_[2:-2:2, 2:-2:2][numpy.logical_or(numpy.logical_or(tmask[1:, 1:] == 0, tmask[:-1, 1:] == 0), numpy.logical_or(tmask[1:, :-1] == 0, tmask[:-1, :-1] == 0))] = 0
        self.exchange_metric(mask_)
        self.mask_[...] = mask_

        for grid in (self.T, self.U, self.V, self.X):
            grid.initialize()

        self.H_.flags.writeable = self.H.flags.writeable = False
        self.z0_.flags.writeable = self.z0.flags.writeable = False
        self.mask_.flags.writeable = self.mask.flags.writeable = False
        self.z_.flags.writeable = self.z.flags.writeable = False
        self.zo_.flags.writeable = self.zo.flags.writeable = False

        _pygetm.Domain.initialize(self, runtype)

        self.initialized = True

    def set_bathymetry(self, depth, scale_factor=None, minimum_depth=None):
        if not isinstance(depth, xarray.DataArray):
            # Depth is provided as raw data and therefore must be already on the supergrid
            self.H[...] = depth
        else:
            # Depth is provided as xarray object that includes coordinates (we require CF compliant longitude, latitude)
            # Interpolate to target grid.
            self.H[...] = depth.getm.interp(self.lon, self.lat).values
        if scale_factor is not None:
            self.H *= scale_factor
        if minimum_depth is not None:
            self.H = numpy.ma.masked_less(self.H, minimum_depth)

    def distribute(self, field):
        return self.tiling.wrap(field, halo=self.halo)

    def plot(self):
        from matplotlib import pyplot
        fig, ax = pyplot.subplots(figsize=(12,12))
        c = ax.contourf(self.lon, self.lat, self.H, 20)
        ax.pcolormesh(self.lon[::2, ::2], self.lat[::2, ::2], numpy.ma.array(self.lon[::2, ::2], mask=True), edgecolors='k', linestyles='-', linewidth=.2)
        pc = ax.pcolormesh(self.lon[1::2, 1::2], self.lat[1::2, 1::2],  numpy.ma.array(self.lon[1::2, 1::2], mask=True), edgecolor='gray', linestyles='--', linewidth=.2)
        cb = fig.colorbar(c)
        cb.set_label('undisturbed water depth (m)')
        fig.savefig('H.png', dpi=300)

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
