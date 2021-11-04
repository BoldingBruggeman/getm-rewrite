from typing import Tuple

import numpy
import numpy.typing

import netCDF4

import pygetm.core

class Base:
    def __init__(self, ndim: int, dtype: numpy.typing.DTypeLike, grid, fill_value=None, atts={}):
        self.dtype = dtype
        self.ndim = ndim
        self.grid = grid
        self.fill_value = fill_value
        self.atts = {}

    def get(self, out: numpy.typing.ArrayLike, slice_spec: Tuple[int]=(), sub: bool=False):
        raise NotImplementedError

class Field(Base):
    def __init__(self, array: pygetm.core.Array):
        atts = {}
        if array.units is not None:
            atts['units'] = array.units
        if array.long_name is not None:
            atts['long_name'] = array.long_name
        Base.__init__(self, array.ndim, array.dtype, array.grid, array.fill_value, atts)
        self.array = array

    def get(self, out: numpy.typing.ArrayLike, slice_spec: Tuple[int]=(), sub: bool=False):
        if sub:
            # Get data for subdomain only
            out[slice_spec + (Ellipsis,)] = self.array.all_values
        else:
            # Get data for global array
            self.array.gather(out, slice_spec)
        return out

class Mask(Base):
    def __init__(self, source: Base):
        Base.__init__(self, source.ndim, source.dtype, source.grid, source.fill_value, source.atts)
        self.source = source
        if self.fill_value is None:
            self.fill_value = netCDF4.default_fillvals.get(self.dtype.str[1:])
        self.mask = Field(self.source.grid.mask)

    def get(self, out: numpy.typing.ArrayLike, slice_spec: Tuple[int]=(), sub: bool=False):
        data = self.source.get(numpy.empty_like(out[slice_spec + (Ellipsis,)]))

        # Obtain mask and apply it
        mask = self.mask.get(numpy.empty(data.shape, dtype=int))
        data[mask == 0] = self.fill_value
        out[slice_spec + (Ellipsis,)] = data
        return out
