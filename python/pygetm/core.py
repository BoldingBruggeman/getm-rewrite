import numbers
from typing import Optional, Union, Tuple, Literal, Mapping, Any
import logging

import numpy, numpy.lib.mixins, numpy.typing
import xarray

from . import _pygetm
from . import parallel
from .constants import *

class Array(_pygetm.Array, numpy.lib.mixins.NDArrayOperatorsMixin):
    __slots__ = ('_xarray', '_scatter', '_gather', '_dist', '_name', 'attrs', '_fill_value', '_ma', 'mapped_field', 'saved', '_shape', '_ndim', '_size', '_dtype')

    def __init__(self, name: Optional[str]=None, units: Optional[str]=None, long_name: Optional[str]=None, fill_value: Optional[Union[float, int]]=None, shape: Optional[Tuple[int]]=None, dtype: Optional[numpy.typing.DTypeLike]=None, grid=None, fabm_standard_name: Optional[str]=None, constant: bool=False, attrs: Mapping[str, Any]={}):
        _pygetm.Array.__init__(self, grid)
        self._xarray: Optional[xarray.DataArray] = None
        self._scatter: Optional[parallel.Scatter] = None
        self._gather: Optional[parallel.Gather] = None
        self._dist: Optional[parallel.DistributedArray] = None
        assert fill_value is None or numpy.ndim(fill_value) == 0, 'fill_value must be a scalar value'
        self._name = name
        self.attrs = attrs.copy()
        if units:
            self.attrs['units'] = units
        if long_name:
            self.attrs['long_name'] = long_name
        self._fill_value = fill_value if fill_value is None or dtype is None else numpy.array(fill_value, dtype=dtype)
        self._ma = None
        self.mapped_field: Optional[xarray.DataArray] = None
        self.saved = False   # to be set by pygetm.output.FieldManager if this variable is requested for output
        self._shape = shape
        self._ndim = None if shape is None else len(shape)
        self._size = None if shape is None else numpy.prod(shape)
        self._dtype = dtype
        self.fabm_standard_name = fabm_standard_name
        self.constant = constant
        self.values = None

    def finish_initialization(self):
        """This is called by the underlying cython implementation after the array receives a value (self.all_values is valid)"""
        assert self.grid is not None
        self._dtype = self.all_values.dtype
        self._ndim = self.all_values.ndim
        if self._fill_value is not None:
            # Cast fill value to dtype of the array
            self._fill_value = numpy.array(self._fill_value, dtype=self._dtype)
        if self.on_boundary or self._ndim == 0:
            # boundary array or scalar
            self.values = self.all_values[...]
        else:
            self.values = self.all_values[..., self.grid.domain.haloy:-self.grid.domain.haloy, self.grid.domain.halox:-self.grid.domain.halox]
        self._shape = self.values.shape
        self._size = self.values.size

    def register(self):
        assert self.grid is not None
        if self._name is not None:
            self.grid.domain.field_manager.register(self)

    def __repr__(self) -> str:
        return super().__repr__() + self.grid.postfix

    def update_halos(self, *args, **kwargs):
        if not self.grid.domain.tiling:
            return
        if self._dist is None:
            self._dist = parallel.DistributedArray(self.grid.domain.tiling, self.all_values, self.grid.halo, overlap=self.grid.overlap)
        return self._dist.update_halos(*args, **kwargs)

    def compare_halos(self, *args, **kwargs):
        if not self.grid.domain.tiling:
            return True
        if self._dist is None:
            self._dist = parallel.DistributedArray(self.grid.domain.tiling, self.all_values, self.grid.halo, overlap=self.grid.overlap)
        return self._dist.compare_halos(*args, **kwargs)

    def scatter(self, global_data: Optional['Array']):
        if self.grid.domain.tiling.n == 1:
            self.values[...] = global_data
            return
        if self._scatter is None:
            self._scatter = parallel.Scatter(self.grid.domain.tiling, self.all_values, halo=self.grid.halo, fill_value=self._fill_value)
        self._scatter(None if global_data is None else global_data.all_values)

    def gather(self, out: Optional['Array']=None, slice_spec=()):
        if self.grid.domain.tiling.n == 1:
            if out is not None:
                out[slice_spec + (Ellipsis,)] = self.values
            return self
        if self._gather is None:
            self._gather = parallel.Gather(self.grid.domain.tiling, self.values, fill_value=self._fill_value)
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

    @staticmethod
    def create(grid, fill: Optional[numpy.typing.ArrayLike]=None, z: Literal[None, True, False, CENTERS, INTERFACES]=None, dtype: numpy.typing.DTypeLike=None, copy: bool=True, on_boundary: bool=False, **kwargs) -> 'Array':
        ar = Array(grid=grid, **kwargs)
        if fill is None and ar.fill_value is not None:
            fill = ar.fill_value
        if fill is not None:
            fill = numpy.asarray(fill)
            if z is None and not on_boundary:
                z = False if fill.ndim != 3 else (INTERFACES if fill.shape[0] == grid.nz_ + 1 else CENTERS)
        if dtype is None:
            dtype = float if fill is None else fill.dtype
        shape = [grid.nbdyp] if on_boundary else [grid.ny_, grid.nx_]
        if z:
            shape.insert(1 if on_boundary else 0, grid.nz_ + 1 if z == INTERFACES else grid.nz_,)
        if copy or fill is None:
            data = numpy.empty(shape, dtype=dtype)
            if fill is not None:
                data[...] = fill
        else:
            data = numpy.broadcast_to(fill, shape)
        ar.wrap_ndarray(data, on_boundary=on_boundary)
        ar.register()
        return ar

    def fill(self, value):
        """Set array to specified value, while respecting the mask: masked points are set to fill_value"""
        self.all_values[...] = value
        if self.fill_value is not None:
            self.all_values[..., self.grid.mask.all_values == 0] = self.fill_value

    @property
    def ma(self) -> numpy.ma.MaskedArray:
        if self._ma is None:
            mask = self.grid.mask.values == 0
            self._ma = numpy.ma.array(self.values, mask=numpy.broadcast_to(mask, self._shape))
        return self._ma

    def plot(self, **kwargs):
        if 'x' not in kwargs and 'y' not in kwargs:
            kwargs['x'] = ('lon' if self.grid.domain.spherical else 'x') + self.grid.postfix
            kwargs['y'] = ('lat' if self.grid.domain.spherical else 'y') + self.grid.postfix
        if 'shading' not in kwargs:
            kwargs['shading'] = 'auto'
        return self.xarray.plot(**kwargs)

    def interp(self, target, z: Literal[None, True, False, CENTERS, INTERFACES]=None):
        if not isinstance(target, Array):
            # Target must be a grid; we need to create the array
            target = Array.create(target, dtype=self.dtype, z=z if z is not None else self.z)
        key = (self.grid.type, target.grid.type)
        if key == (_pygetm.UGRID, _pygetm.TGRID):
            self.grid.interp_x(self, target, offset=1)
        elif key == (_pygetm.TGRID, _pygetm.UGRID):
            self.grid.interp_x(self, target, offset=0)
        elif key == (_pygetm.VGRID, _pygetm.TGRID):
            self.grid.interp_y(self, target, offset=1)
        elif key == (_pygetm.TGRID, _pygetm.VGRID):
            self.grid.interp_y(self, target, offset=0)
        elif key == (_pygetm.XGRID, _pygetm.TGRID):
            self.grid.interp_xy(self, target, ioffset=0, joffset=0)
        elif key == (_pygetm.TGRID, _pygetm.XGRID):
            self.grid.interp_xy(self, target, ioffset=1, joffset=1)
        elif key in ((_pygetm.UGRID, _pygetm.UUGRID), (_pygetm.VGRID, _pygetm.UVGRID)):
            self.grid.interp_x(self, target, offset=0)
        elif key in ((_pygetm.UGRID, _pygetm.VUGRID), (_pygetm.VGRID, _pygetm.VVGRID)):
            self.grid.interp_y(self, target, offset=0)
        elif key[0] == key[1] and self.z == INTERFACES and target.z == CENTERS:
            # vertical interpolation from layer interfaces to layer centers
            self.grid.interp_z(self, target, offset=0)
        elif key[0] == key[1] and self.z == CENTERS and target.z == INTERFACES:
            # vertical interpolation from layer centers to layer interfaces (top and bottom interfaces will be left untouched)
            self.grid.interp_z(self, target, offset=1)
        else:
            raise NotImplementedError('interp does not know how to interpolate from grid type %i to grid type %i' % (self.grid.type, target.grid.type))
        return target

    def __array__(self, dtype=None):
        return numpy.asarray(self.values, dtype=dtype)

    def isel(self, z: int, **kwargs):
        if self._ndim != 3:
            raise NotImplementedError
        if self.units is not None:
            kwargs.setdefault('units', self.units)
        if self.long_name is not None:
            kwargs.setdefault('long_name', '%s @ k=%i' % (self.long_name, z))
        ar = Array(grid=self.grid, **kwargs)
        ar.wrap_ndarray(self.all_values[z, ...])
        ar.register()
        return ar

    def __getitem__(self, key):
        return self.values[key]

    def __setitem__(self, key, values):
        self.values[key] = values

    @property
    def shape(self):
        return self._shape

    @property
    def ndim(self):
        return self._ndim

    @property
    def z(self):
        if self._ndim != (2 if self.on_boundary else 3):
            return False
        nz = self._shape[1 if self.on_boundary else 0] 
        return INTERFACES if nz == self.grid.nz_ + 1 else CENTERS

    @property
    def size(self):
        return self._size

    @property
    def dtype(self):
        return self._dtype

    @property
    def name(self) -> Optional[str]:
        return self._name

    @property
    def units(self) -> Optional[str]:
        return self.attrs.get('units')

    @property
    def long_name(self) -> Optional[str]:
        return self.attrs.get('long_name')

    @property
    def fill_value(self) -> Optional[Union[int, float]]:
        return self._fill_value

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
            return tuple(self.create(self.grid, x) for x in result)
        elif method == 'at':
            # no return value
            return None
        else:
            # one return value
            return self.create(self.grid, result)

    def set(self, value: Union[float, numpy.ndarray, xarray.DataArray], periodic_lon: bool=True, on_grid: bool=False, include_halos: Optional[bool]=None):
        """Link this array to a field or value managed by the input manager. This will perform temporal and spatial interpolation as required."""
        self.grid.domain.input_manager.add(self, value, periodic_lon=periodic_lon, on_grid=on_grid, include_halos=include_halos)

    def require_set(self, logger: Optional[logging.Logger]=None):
        """Assess whether all non-masked cells of this field have been set. If not, an error message is written to the log and False is returned."""
        valid = True
        if self._fill_value is not None:
            invalid = self.ma == self._fill_value
            if invalid.any():
                (logger or logging.getLogger()).error('%s is masked (%s) in %i active grid cells.' % (self.name, self._fill_value, invalid.sum()))
                valid = False
        return valid

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
            dims = ('y' + self.grid.postfix, 'x' + self.grid.postfix)
            if self.ndim == 3:
                dims = ('zi' if self.z == INTERFACES else 'z',) + dims
            self._xarray = xarray.DataArray(self.values, coords=coords, dims=dims, attrs=attrs, name=self.name)
        return self._xarray
