from typing import MutableMapping, Tuple, Union, Iterable, Optional, Sequence
import collections

import numpy
import numpy.typing

import netCDF4

import pygetm.core
from pygetm.constants import INTERFACES

class Base:
    __slots__ = 'dtype', 'ndim', 'grid', 'fill_value', 'atts', 'constant', 'coordinates'
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

    def is_updatable(self) -> bool:
        return False

class FieldCollection:
    def __init__(self, field_manager, default_dtype: Optional[numpy.typing.DTypeLike]=None):
        self.fields: MutableMapping[str, Base] = collections.OrderedDict()
        self.expression2name = {}
        self.field_manager = field_manager
        self.default_dtype = default_dtype
        self._updatable = []

    def request(self, field: Union[str, pygetm.core.Array, Iterable[Union[str, pygetm.core.Array]]], output_name: Optional[str]=None, dtype: Optional[numpy.typing.DTypeLike]=None, mask: Optional[bool]=None, time_average: bool=False, generate_unique_name: bool=False):
        if isinstance(field, str):
            if field not in self.field_manager.fields:
                raise Exception('Unknown field "%s" requested. Available: %s' % (field, ', '.join(self.field_manager.fields)))
            array = self.field_manager.fields[field]
        elif isinstance(field, pygetm.core.Array):
            array = field
        else:
            # Multiple names requested; call request for each individually
            if output_name is not None:
                raise Exception('Trying to add multiple fields to %s. In that case, output_name cannot be specified.' % (self,))
            for f in field:
                self.request(f, dtype=dtype, mask=mask, time_average=time_average)
            return

        if output_name is None:
            if array.name is None:
                raise Exception('Trying to add an unnamed variable to %s. In this case, output_name must be provided' % self)
            output_name = array.name
            if generate_unique_name:
                i = 0
                while output_name in self.fields:
                    output_name = '%s_%i' % (array.name, i)
                    i += 1

        if mask is None:
            mask = array.attrs.get('_mask_output', False)

        if output_name in self.fields:
            raise Exception('A variable with name "%s" has already been added to %s.' % (output_name, self))

        array.saved = True
        if dtype is None and array.dtype == float:
            dtype = self.default_dtype
        field = Field(array, self, dtype=dtype)
        if time_average:
            field = TimeAverage(field)
        if time_average or mask:
            field = Mask(field)
        if field.is_updatable():
            self._updatable.append(field)
        self.fields[output_name] = field
        self.expression2name[field.get_expression()] = output_name
        field.coordinates = field.get_coordinates()
        return output_name

    def require(self, expression: str) -> Base:
        if expression in self.expression2name:
            return self.expression2name[expression]
        return self.request(expression, generate_unique_name=True)

class Field(Base):
    __slots__ = 'collection', 'array', 'global_array', '_ncvar'
    def __init__(self, array: pygetm.core.Array, collection: FieldCollection, dtype: Optional[numpy.typing.DTypeLike]=None, atts=None):
        if atts is None:
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
        super().__init__(array.ndim, dtype or array.dtype, array.grid, array.fill_value, array.constant, atts)

    def get(self, out: numpy.typing.ArrayLike, slice_spec: Tuple[int]=(), sub: bool=False) -> numpy.typing.ArrayLike:
        if sub:
            # Get data for subdomain only
            if out is not None:
                out[slice_spec + (Ellipsis,)] = self.array.all_values
            else:
                out = self.array.all_values
        else:
            # Get data for global array
            # If we have access to the full global field (root rank only), use that, as gathering from subdomains may leave gaps.
            # Nevertheless we cannot skip the gather in that case, because only all non-root ranks will call gather anyway.
            self.array.gather(out, slice_spec)
            if self.global_array:
                out[...] = self.global_array.values
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

    @property
    def on_boundary(self) -> bool:
        return self.array.on_boundary

    def get_expression(self) -> str:
        return self.array.name

class UnivariateTransform(Field):
    __slots__ = '_source',
    def __init__(self, source: Field):
        self._source = source
        assert source.fill_value is not None, 'UnivariateTransform cannot be used on variables without fill value.'
        array = pygetm.core.Array.create(source.grid, z=source.z, dtype=source.dtype, on_boundary=source.on_boundary, fill_value=source.fill_value)
        super().__init__(array, source.collection, atts=source.atts)

    def is_updatable(self) -> bool:
        return self._source.is_updatable()

    def update(self):
        return self._source.update()

class Mask(UnivariateTransform):
    def get(self, out: numpy.typing.ArrayLike, slice_spec: Tuple[int]=(), sub: bool=False) -> numpy.typing.ArrayLike:
        self._source.get(out=self.array.all_values, sub=True)
        self.array.all_values[..., self.grid.mask.all_values == 0] = self.fill_value
        super().get(out, slice_spec, sub)

    def get_expression(self) -> str:
        return 'mask(%s)' % self._source.get_expression()

class TimeAverage(UnivariateTransform):
    __slots__ = '_n',
    def __init__(self, source: Field):
        super().__init__(source)
        self._n = 0

    def is_updatable(self) -> bool:
        return True

    def update(self):
        if self._n == 0:
            self._source.get(out=self.array.all_values, sub=True)
        else:
            self.array.all_values += self._source.get(out=None, sub=True)
        self._n += 1

    def get(self, out: numpy.typing.ArrayLike, slice_spec: Tuple[int]=(), sub: bool=False) -> numpy.typing.ArrayLike:
        if self._n > 0:
            self.array.all_values /= self._n
        super().get(out, slice_spec, sub)
        self._n = 0

    def get_expression(self) -> str:
        return 'time_average(%s)' % self._source.get_expression()
