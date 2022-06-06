from typing import MutableMapping, Tuple, Union, Optional, Sequence, Mapping
import collections

from numpy.typing import DTypeLike, ArrayLike

import pygetm.core
from pygetm.constants import INTERFACES


class Base:
    __slots__ = "dtype", "ndim", "grid", "fill_value", "atts", "constant", "coordinates"

    def __init__(
        self,
        ndim: int,
        dtype: DTypeLike,
        grid,
        fill_value=None,
        constant: bool = False,
        atts={},
    ):
        self.dtype = dtype
        self.ndim = ndim
        self.grid = grid
        self.fill_value = fill_value
        self.atts = atts
        self.constant = constant
        self.coordinates = []

    def get(
        self, out: ArrayLike, slice_spec: Tuple[int] = (), sub: bool = False,
    ):
        raise NotImplementedError

    def get_coordinates(self) -> Sequence[str]:
        raise NotImplementedError

    def get_expression(self) -> str:
        raise NotImplementedError

    def is_updatable(self) -> bool:
        return False


class FieldCollection:
    def __init__(
        self,
        available_fields: Mapping[str, pygetm.core.Array],
        default_dtype: Optional[DTypeLike] = None,
    ):
        self.fields: MutableMapping[str, Base] = collections.OrderedDict()
        self.expression2name = {}
        self.available_fields = available_fields
        self.default_dtype = default_dtype
        self._updatable = []

    def request(
        self,
        *fields: Union[str, pygetm.core.Array],
        output_name: Optional[str] = None,
        dtype: Optional[DTypeLike] = None,
        mask: Optional[bool] = None,
        time_average: bool = False,
        generate_unique_name: bool = False
    ) -> Tuple[str]:
        """Add one or more arrays to this field collection.

        Args:
            *fields: names of arrays or array objects to add. When names are provided,
                they will be looked up in the field manager.
            output_name: name to use for this field. This can only be provided if a
                single field is requested.
            dtype: data type of the field to use. Array values will be cast to this data
                type whenever the field is saved
            mask: whether to explicitly set masked values to the array's fill value when
                saving. If this argument is not provided, masking behavior is determined
                by the array's _mask_output flag.
            time_average: whether to time-average the field
            generate_unique_name: whether to generate a unique output name for requested
                fields if a field with the same name has previously been added to the
                collection. If this is not set and a field with this name was added
                previously, an exception will be raised.

            Returns:
                tuple with names of the newly added fields
        """
        if not fields:
            raise Exception(
                "One or more positional arguments must be provided"
                "to specify the field(s) requested."
            )

        # For backward compatibility: a tuple of names coud be provided
        if len(fields) == 1 and isinstance(fields[0], tuple):
            fields = fields[0]

        arrays = []
        for field in fields:
            if isinstance(field, str):
                if field not in self.available_fields:
                    raise Exception(
                        'Unknown field "%s" requested. Available: %s'
                        % (field, ", ".join(self.available_fields))
                    )
                arrays.append(self.available_fields[field])
            elif isinstance(field, pygetm.core.Array):
                arrays.append(field)
            else:
                raise Exception(
                    "Incorrect field specification %r."
                    "Expected a string or an object of type pygetm.core.Array." % field
                )

        if len(arrays) > 1 and output_name is not None:
            raise Exception(
                "Trying to add multiple fields to %r."
                "In that case, output_name cannot be specified." % (self,)
            )

        names = []
        for array in arrays:
            name = output_name
            if name is None:
                if array.name is None:
                    raise Exception(
                        "Trying to add an unnamed variable to %s."
                        "In this case, output_name must be provided" % self
                    )
                name = array.name
                if generate_unique_name:
                    i = 0
                    while name in self.fields:
                        name = "%s_%i" % (array.name, i)
                        i += 1

            if mask is None:
                mask = array.attrs.get("_mask_output", False)

            if name in self.fields:
                raise Exception(
                    'A variable with name "%s" has already been added to %s.'
                    % (name, self)
                )

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
            self.fields[name] = field
            self.expression2name[field.get_expression()] = name
            field.coordinates = field.get_coordinates()
            names.append(name)

        return tuple(names)

    def require(self, expression: str) -> str:
        """Ensure that the specified variable (or expression of variables) is included
        in the field collection. This is typically used to add coordinate variables.

        Args:
            expression: variable name or expression of variable(s)
        """
        if expression in self.expression2name:
            return self.expression2name[expression]
        (output_name,) = self.request(expression, generate_unique_name=True)
        return output_name


class Field(Base):
    __slots__ = "collection", "array", "global_array", "_ncvar"

    def __init__(
        self,
        array: pygetm.core.Array,
        collection: FieldCollection,
        dtype: Optional[DTypeLike] = None,
        atts=None,
    ):
        if atts is None:
            atts = {}
        if array.units is not None:
            atts["units"] = array.units
        if array.long_name is not None:
            atts["long_name"] = array.long_name
        self.collection = collection
        self.array = array
        self.global_array = None
        global_domain = array.grid.domain.glob
        if global_domain and array.constant:
            self.global_array = global_domain.fields.get(array.name)
        super().__init__(
            array.ndim,
            dtype or array.dtype,
            array.grid,
            array.fill_value,
            array.constant,
            atts,
        )

    def get(
        self, out: ArrayLike, slice_spec: Tuple[int] = (), sub: bool = False,
    ) -> ArrayLike:
        if sub:
            # Get data for subdomain only
            if out is not None:
                out[slice_spec + (Ellipsis,)] = self.array.all_values
            else:
                out = self.array.all_values
        else:
            # Get data for global array
            # If we have access to the full global field (root rank only), use that,
            # as gathering from subdomains may leave gaps.
            # Nevertheless we cannot skip the gather in that case,
            # because all non-root ranks will call gather anyway.
            self.array.gather(out, slice_spec)
            if self.global_array:
                out[...] = self.global_array.values
        return out

    def get_coordinates(self) -> Sequence["Base"]:
        coords = []
        x, y = (
            (self.array.grid.lon, self.array.grid.lat)
            if self.array.grid.domain.spherical
            else (self.array.grid.x, self.array.grid.y)
        )
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
    __slots__ = ("_source",)

    def __init__(self, source: Field):
        self._source = source
        assert (
            source.fill_value is not None
        ), "UnivariateTransform cannot be used on variables without fill value."
        array = pygetm.core.Array.create(
            source.grid,
            z=source.z,
            dtype=source.dtype,
            on_boundary=source.on_boundary,
            fill_value=source.fill_value,
        )
        super().__init__(array, source.collection, atts=source.atts)

    def is_updatable(self) -> bool:
        return self._source.is_updatable()

    def update(self):
        return self._source.update()


class Mask(UnivariateTransform):
    def get(
        self, out: ArrayLike, slice_spec: Tuple[int] = (), sub: bool = False,
    ) -> ArrayLike:
        self._source.get(out=self.array.all_values, sub=True)
        self.array.all_values[..., self.grid.mask.all_values == 0] = self.fill_value
        super().get(out, slice_spec, sub)

    def get_expression(self) -> str:
        return "mask(%s)" % self._source.get_expression()


class TimeAverage(UnivariateTransform):
    __slots__ = ("_n",)

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

    def get(
        self, out: ArrayLike, slice_spec: Tuple[int] = (), sub: bool = False,
    ) -> ArrayLike:
        if self._n > 0:
            self.array.all_values /= self._n
        super().get(out, slice_spec, sub)
        self._n = 0

    def get_expression(self) -> str:
        return "time_average(%s)" % self._source.get_expression()
