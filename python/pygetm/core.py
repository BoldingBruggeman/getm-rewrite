import numbers
from typing import Optional

import numpy, numpy.lib.mixins, numpy.typing
import xarray

from . import _pygetm
from . import parallel

class Array(_pygetm.Array, numpy.lib.mixins.NDArrayOperatorsMixin):
    def __init__(self, name: Optional[str]=None, units: Optional[str]=None, long_name: Optional[str]=None):
        self._xarray: Optional[xarray.DataArray] = None
        self._scatter: Optional[parallel.Scatter] = None
        self._gather: Optional[parallel.Gather] = None
        self._dist: Optional[parallel.DistributedArray] = None
        self._name = name
        self._units = units
        self._long_name = long_name

    def __repr__(self) -> str:
        return super().__repr__() + self.grid.postfix

    def update_halos(self, *args, **kwargs):
        if not self.grid.domain.tiling:
            return
        if self._dist is None:
            self._dist = parallel.DistributedArray(self.grid.domain.tiling, self.all_values, self.grid.halo)
        return self._dist.update_halos(*args, **kwargs)

    def compare_halos(self, *args, **kwargs):
        if not self.grid.domain.tiling:
            return True
        if self._dist is None:
            self._dist = parallel.DistributedArray(self.grid.domain.tiling, self.all_values, self.grid.halo)
        return self._dist.compare_halos(*args, **kwargs)

    def scatter(self, global_data: Optional['Array']):
        if self.grid.domain.tiling.n == 1:
            self.values[...] = global_data
            return
        if self._scatter is None:
            self._scatter = parallel.Scatter(self.grid.domain.tiling, self.all_values, halo=self.grid.halo)
        self._scatter(None if global_data is None else global_data.all_values)

    def gather(self, out: Optional['Array']=None, slice_spec=()):
        if self.grid.domain.tiling.n == 1:
            if out is not None:
                out[slice_spec + (Ellipsis,)] = self.values
            return self
        if self._gather is None:
            self._gather = parallel.Gather(self.grid.domain.tiling, self.values)
        result = self._gather(out.values if isinstance(out, Array) else out, slice_spec=slice_spec)
        if result is not None and out is None:
            out = self.grid.domain.glob.grids[self.grid.type].array(dtype=self.dtype)
            out[...] = result
        return out

    def global_sum(self, reproducible: bool=False, where: Optional['Array']=None) -> Optional[numpy.ndarray]:
        if reproducible:
            all = self.gather()
            if where is not None:
                where = where.gather()
            if all is not None:
                return all.values.sum(where=numpy._NoValue if where is None else where.values)
        else:
            local_sum = self.values.sum(where=numpy._NoValue if where is None else where.values)
            return parallel.Sum(self.grid.domain.tiling, local_sum)()

    def global_mean(self, reproducible: bool=False, where: Optional['Array']=None) -> Optional[numpy.ndarray]:
        sum = self.global_sum(reproducible=reproducible, where=where)
        if where is not None:
            count = where.global_sum()
        else:
            count = parallel.Sum(self.grid.domain.tiling, self.values.size)()
        if sum is not None:
            return sum / count

    def finish_initialization(self):
        assert self.grid is not None
        if self.name is not None:
            self.grid.domain.field_manager.register(self)
        halo = self.grid.halo
        self.values = self.all_values[halo:-halo, halo:-halo]

    @staticmethod
    def create(grid, fill: Optional[numpy.typing.ArrayLike]=None, dtype: numpy.typing.DTypeLike=None, copy: bool=True, **kwargs) -> 'Array':
        ar = Array(**kwargs)
        if fill is not None:
            fill = numpy.asarray(fill)
        if dtype is None:
            dtype = float if fill is None else fill.dtype
        shape = (grid.ny, grid.nx)
        if copy or fill is None:
            data = numpy.empty(shape, dtype=dtype)
            if fill is not None:
                data[...] = fill
        else:
            data = numpy.broadcast_to(fill, shape)
        ar.wrap_ndarray(grid, data)
        return ar

    @property
    def ma(self) -> numpy.ma.MaskedArray:
        return numpy.ma.array(self.values, mask=self.grid.mask.values==0)

    @property
    def xarray(self) -> xarray.DataArray:
        if self._xarray is None:
            attrs = {}
            for key in ('units', 'long_name'):
                value = getattr(self, key)
                if value is not None:
                    attrs[key] = value
            dom = self.grid.domain
            coords = {}
            if self.name not in ('x' + self.grid.postfix, 'y' + self.grid.postfix, 'lon' + self.grid.postfix, 'lat' + self.grid.postfix):
                if dom.x_is_1d:
                    coords['x%s' % self.grid.postfix] = self.grid.x.xarray[0, :]
                if dom.y_is_1d:
                    coords['y%s' % self.grid.postfix] = self.grid.y.xarray[:, 0]
                coords['x%s2' % self.grid.postfix] = self.grid.x.xarray
                coords['y%s2' % self.grid.postfix] = self.grid.y.xarray
                if dom.lon is not None:
                    coords['lon%s' % self.grid.postfix] = self.grid.lon.xarray
                if dom.lat is not None:
                    coords['lat%s' % self.grid.postfix] = self.grid.lat.xarray
            self._xarray = xarray.DataArray(self.values, coords=coords, dims=('y' + self.grid.postfix, 'x' + self.grid.postfix), attrs=attrs, name=self.name)
        return self._xarray

    @property
    def plot(self):
        return self.xarray.plot

    def interp(self, target_grid, out: Optional['Array'] = None):
        if out is None:
            out = Array.create(target_grid, dtype=self.dtype)
        key = (self.grid.type, target_grid.type)
        if key == (_pygetm.UGRID, _pygetm.TGRID):
            self.grid.interp_x(self, out, offset=1)
        elif key == (_pygetm.VGRID, _pygetm.TGRID):
            self.grid.interp_y(self, out, offset=1)
        elif key in ((_pygetm.UGRID, _pygetm.UUGRID), (_pygetm.VGRID, _pygetm.UVGRID)):
            self.grid.interp_x(self, out, offset=0)
        elif key in ((_pygetm.UGRID, _pygetm.VUGRID), (_pygetm.VGRID, _pygetm.VVGRID)):
            self.grid.interp_y(self, out, offset=0)
        else:
            raise NotImplementedError('Map not implemented for grid type %i -> grid type %i' % (self.grid.type, target_grid.type))
        return out

    def __array__(self, dtype=None):
        return numpy.asarray(self.values, dtype=dtype)

    def __getitem__(self, key):
        return self.values[key]

    def __setitem__(self, key, values):
        self.values[key] = values

    @property
    def shape(self):
        return self.values.shape

    @property
    def ndim(self):
        return self.values.ndim

    @property
    def size(self):
        return self.values.size

    @property
    def dtype(self):
        return self.values.dtype

    @property
    def name(self) -> Optional[str]:
        return self._name

    @property
    def units(self) -> Optional[str]:
        return self._units

    @property
    def long_name(self) -> Optional[str]:
        return self._long_name

    # Below based on https://numpy.org/devdocs/reference/generated/numpy.lib.mixins.NDArrayOperatorsMixin.html#numpy.lib.mixins.NDArrayOperatorsMixin
    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        if method != '__call__':
            return NotImplemented

        out = kwargs.get('out', ())
        for x in inputs + out:
            # Only support operations with instances of _HANDLED_TYPES.
            # Use ArrayLike instead of type(self) for isinstance to
            # allow subclasses that don't override __array_ufunc__ to
            # handle ArrayLike objects.
            if not isinstance(x, (numpy.ndarray, numbers.Number, Array)):
                return NotImplemented
            if isinstance(x, Array) and not x.grid is self.grid:
                return NotImplemented

        # Defer to the implementation of the ufunc on unwrapped values.
        inputs = tuple(x.all_values if isinstance(x, Array) else x for x in inputs)
        if out:
            kwargs['out'] = tuple(x.all_values if isinstance(x, Array) else x for x in out)
        result = getattr(ufunc, method)(*inputs, **kwargs)

        if type(result) is tuple:
            # multiple return values
            return tuple(self.create(self.grid, fill=x) for x in result)
        elif method == 'at':
            # no return value
            return None
        else:
            # one return value
            return self.create(self.grid, fill=result)
