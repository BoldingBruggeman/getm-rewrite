from typing import MutableMapping, Tuple, Union, Iterable, Optional, Sequence
import collections

import numpy
import numpy.typing

import netCDF4

import pygetm.core
from pygetm.constants import INTERFACES

class Base:
    def __init__(self, ndim: int, dtype: numpy.typing.DTypeLike, grid, fill_value=None, constant: bool=False, atts={}):
        self.dtype = dtype
        self.ndim = ndim
        self.grid = grid
        self.fill_value = fill_value
        self.atts = atts
        self.constant = constant
        self.coordinates = []

    def get(self, out: numpy.typing.ArrayLike, slice_spec: Tuple[int]=(), sub: bool=False):
        raise NotImplementedError

    def get_coordinates(self) -> Sequence[str]:
        raise NotImplementedError

    def get_expression(self) -> str:
        raise NotImplementedError

class FieldCollection:
    def __init__(self, field_manager):
        self.fields: MutableMapping[str, Base] = collections.OrderedDict()
        self.expression2name = {}
        self.field_manager = field_manager

    def request(self, name: Union[str, Iterable[str]], output_name: Optional[str]=None, dtype: Optional[numpy.typing.DTypeLike]=None, mask: bool=False, generate_unique_name: bool=False):
        if not isinstance(name, str):
            # Multiple names requested; call request for each individually
            assert output_name is None, 'request is called with multiple names: %s. Therefore, output_name cannot be specified.' % (name,)
            for n in name:
                self.request(n, dtype=dtype, mask=mask)
            return

        assert name in self.field_manager.fields, 'Unknown field "%s" requested. Available: %s' % (name, ', '.join(self.field_manager.fields))

        if output_name is None:
            output_name = name
            if generate_unique_name:
                i = 0
                while output_name in self.fields:
                    output_name = '%s_%i' % (name, i)
                    i += 1
        assert output_name not in self.fields, 'A variable with name "%s" has already been added to %s.' % (output_name, self)

        array = self.field_manager.fields[name]
        array.saved = True
        dtype = dtype or array.dtype
        field = Field(array, self)
        if mask:
            field = Mask(field, self)
        self.fields[output_name] = field
        self.expression2name[field.get_expression()] = output_name
        field.coordinates = field.get_coordinates()
        return output_name

    def require(self, expression: str) -> Base:
        if expression in self.expression2name:
            return self.expression2name[expression]
        return self.request(expression, generate_unique_name=True)

class Field(Base):
    def __init__(self, array: pygetm.core.Array, collection: FieldCollection):
        atts = {}
        if array.units is not None:
            atts['units'] = array.units
        if array.long_name is not None:
            atts['long_name'] = array.long_name
        self.collection = collection
        self.array = array
        self.global_array = None
        global_domain = array.grid.domain.glob
        if global_domain and array.constant:
            self.global_array = global_domain.field_manager.fields.get(array.name)
        super().__init__(array.ndim, array.dtype, array.grid, array.fill_value, array.constant, atts)

    def get(self, out: numpy.typing.ArrayLike, slice_spec: Tuple[int]=(), sub: bool=False) -> numpy.typing.ArrayLike:
        if sub:
            # Get data for subdomain only
            out[slice_spec + (Ellipsis,)] = self.array.all_values
        else:
            # Get data for global array
            if self.global_array:
                out[...] = self.global_array.values
            self.array.gather(out, slice_spec)
        return out

    def get_coordinates(self) -> Sequence['Base']:
        coords = []
        x, y = (self.array.grid.lon, self.array.grid.lat) if self.array.grid.domain.spherical else (self.array.grid.x, self.array.grid.y)
        coords.append(self.collection.require(x.name))
        coords.append(self.collection.require(y.name))
        if self.array.z:
            z = self.array.grid.zf if self.array.z == INTERFACES else self.array.grid.zc
            coords.append(self.collection.require(z.name))
        return coords

    @property
    def z(self) -> bool:
        return self.array.z

    def get_expression(self) -> str:
        return self.array.name

class Mask(Base):
    def __init__(self, source: Base, collection: FieldCollection):
        super().__init__(source.ndim, source.dtype, source.grid, source.fill_value, source.constant, source.atts)
        self.source = source
        if self.fill_value is None:
            self.fill_value = netCDF4.default_fillvals.get(self.dtype.str[1:])
        self.mask = Field(self.source.grid.mask, collection)

    def get(self, out: numpy.typing.ArrayLike, slice_spec: Tuple[int]=(), sub: bool=False) -> numpy.typing.ArrayLike:
        data = self.source.get(None if out is None else numpy.empty_like(out[slice_spec + (Ellipsis,)]))

        # Obtain mask and apply it
        mask = self.mask.get(None if out is None else numpy.empty(data.shape, dtype=int))
        if data is not None:
            data[mask == 0] = self.fill_value
            out[slice_spec + (Ellipsis,)] = data
        return out

    def get_coordinates(self) -> Sequence['Base']:
        return self.source.get_coordinates()

    def get_expression(self) -> str:
        return 'mask(%s)' % self.source.get_expression()

    @property
    def z(self) -> bool:
        return self.source.z
