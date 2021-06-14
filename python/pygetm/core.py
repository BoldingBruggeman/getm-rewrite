import numbers
from typing import Optional

import numpy, numpy.lib.mixins
import xarray

from . import _pygetm
from . import parallel

class Array(_pygetm.Array, numpy.lib.mixins.NDArrayOperatorsMixin):
    def __init__(self):
        self._x = None

    def update_halos(self):
        pass

    def finish_initialization(self):
        assert self.grid is not None
        tiling = self.grid.domain.tiling
        if tiling is not None:
            dist = parallel.DistributedArray(tiling, self.all_values, self.grid.halo)
            self.update_halos = dist.update_halos
        halo = self.grid.halo
        self.values = self.all_values[halo:-halo, halo:-halo]

    @staticmethod
    def create(grid, fill=None, dtype=None) -> 'Array':
        if fill is not None:
            fill = numpy.asarray(fill)
        if dtype is None:
            dtype = float if fill is None else fill.dtype
        ar = Array()
        ar.empty(grid, dtype)
        if fill is not None:
            ar.all_values[...] = fill
        return ar

    @property
    def ma(self) -> numpy.ma.MaskedArray:
        return numpy.ma.array(self.values, mask=self.grid.mask.values==0)

    @property
    def x(self) -> xarray.DataArray:
        if self._x is None:
            self._x = xarray.DataArray(self.values)
        return self._x

    @property
    def plot(self):
        return self.x.plot

    def interp(self, target_grid):
        data = target_grid.map(self.all_values, self.grid)
        return Array.create(target_grid, fill=data)

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
