from typing import (
    Iterable,
    MutableMapping,
    Tuple,
    Union,
    Optional,
    Sequence,
    Mapping,
    Literal,
    Callable,
    List,
    Any,
)
import collections
import enum
import functools

import numpy as np
from numpy.typing import DTypeLike, ArrayLike

import pygetm.core
import pygetm.parallel
import pygetm.domain
import pygetm._pygetm
import pygetm.util.interpolate
from pygetm.constants import CENTERS, INTERFACES, TimeVarying


class Base:
    __slots__ = (
        "expression",
        "dtype",
        "shape",
        "ndim",
        "dims",
        "fill_value",
        "attrs",
        "time_varying",
        "coordinates",
    )

    def __init__(
        self,
        expression: str,
        shape: Tuple[int, ...],
        dims: Tuple[str, ...],
        dtype: DTypeLike,
        fill_value=None,
        time_varying: Union[TimeVarying, Literal[False]] = TimeVarying.MICRO,
        attrs: Mapping[str, Any] = {},
    ):
        self.expression = expression
        self.dtype = dtype
        self.shape = shape
        self.dims = dims
        self.ndim = len(shape)
        self.fill_value = fill_value
        self.attrs = attrs
        self.time_varying = time_varying
        self.coordinates: List[str] = []

    def get(
        self, out: Optional[ArrayLike] = None, slice_spec: Tuple[int, ...] = ()
    ) -> ArrayLike:
        raise NotImplementedError

    @property
    def coords(self) -> Sequence["Base"]:
        raise NotImplementedError

    @property
    def default_name(self) -> str:
        raise NotImplementedError

    def gather(self, tiling: pygetm.parallel.Tiling) -> "Base":
        if tiling.n == 1:
            return Slice(self)
        else:
            return Gather(self, tiling)

    @property
    def updatable(self) -> bool:
        return False

    @property
    def updater(self) -> Optional[Callable]:
        raise NotImplementedError

    @property
    def grid(self) -> Optional[pygetm.domain.Grid]:
        return None

    @property
    def mask(self) -> np.ndarray:
        if self.ndim > 2 and hasattr(self.grid, "_land3d"):
            return self.grid._land3d
        return self.grid._land


class WrappedArray(Base):
    __slots__ = ("_name", "values")

    def __init__(self, values: np.ndarray, name: str, dims: Tuple[str]):
        super().__init__(name, values.shape, dims, values.dtype, time_varying=False)
        self._name = name
        self.values = values

    def gather(self, tiling: pygetm.parallel.Tiling) -> Base:
        return self

    @property
    def default_name(self) -> str:
        return self._name

    @property
    def coords(self) -> Sequence["Base"]:
        return ()

    def get(
        self, out: Optional[ArrayLike] = None, slice_spec: Tuple[int] = ()
    ) -> ArrayLike:
        if out is None:
            return self.values
        else:
            out[slice_spec] = self.values
            return out


class Updatable(enum.Enum):
    ALWAYS = 1
    MACRO_ONLY = 2


class FieldCollection:
    def __init__(
        self,
        available_fields: Mapping[str, pygetm.core.Array],
        default_dtype: Optional[DTypeLike] = None,
        sub: bool = False,
    ):
        self.fields: MutableMapping[str, Base] = collections.OrderedDict()
        self.expression2name = {}
        self.available_fields = available_fields
        self.default_dtype = default_dtype
        self.sub = sub
        self._updaters = {}

    def request(
        self,
        *fields: Union[str, pygetm.core.Array],
        output_name: Optional[str] = None,
        dtype: Optional[DTypeLike] = None,
        mask: Optional[bool] = None,
        time_average: bool = False,
        grid: Optional[pygetm.domain.Grid] = None,
        z: Optional[Literal[None, CENTERS, INTERFACES]] = None,
        generate_unique_name: bool = False,
    ) -> Tuple[str, ...]:
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

        # For backward compatibility: a tuple of names could be provided
        if len(fields) == 1 and isinstance(fields[0], tuple):
            fields = fields[0]

        arrays = []
        for field in fields:
            if isinstance(field, str):
                if field not in self.available_fields:
                    raise Exception(
                        f"Unknown field {field!r} requested."
                        f" Available: {', '.join(self.available_fields)}"
                    )
                arrays.append(self.available_fields[field])
            elif isinstance(field, pygetm.core.Array):
                arrays.append(field)
            else:
                raise Exception(
                    f"Incorrect field specification {field!r}."
                    " Expected a string or an object of type pygetm.core.Array."
                )

        if len(arrays) > 1 and output_name is not None:
            raise Exception(
                f"Trying to add multiple fields to {self!r}."
                " In this case, output_name cannot be specified."
            )

        names = []
        for array in arrays:
            name = output_name
            if name is None:
                if array.name is None:
                    raise Exception(
                        f"Trying to add an unnamed variable to {self!r}."
                        " In this case, output_name must be provided"
                    )
                name = array.name

            mask_current = mask
            if mask_current is None:
                mask_current = array.attrs.get("_mask_output", False)

            array.saved = True
            source_grid = array.grid
            if dtype is None and array.dtype == float:
                dtype = self.default_dtype
            field = Field(array, dtype=dtype)
            if time_average and field.time_varying:
                field = TimeAverage(field)
            if grid and array.grid is not grid:
                field = Regrid(field, grid=grid)
            if z is not None and array.z:
                if isinstance(z, (Iterable, float)):
                    field = InterpZ(field, z, "z1")
                elif array.z and array.z != z:
                    field = Regrid(field, z=z)
            if time_average or mask_current or grid:
                field = Mask(field)
            if not self.sub:
                field = field.gather(source_grid.domain.tiling)
            for tf in source_grid.domain.default_output_transforms:
                field = tf(field)
            names.append(self._add_field(field, name, generate_unique_name))

        return tuple(names)

    def _add_field(self, field: Base, name: str, generate_unique_name: bool):
        final_name = name
        if generate_unique_name:
            i = 0
            while final_name in self.fields:
                final_name = f"{name}_{i}"
                i += 1
        elif final_name in self.fields:
            raise Exception(
                f"A variable with name {name!r} has already been added to {self!r}."
            )
        if field.updatable:
            triggering_macro_values = [True]
            if field.updatable == Updatable.ALWAYS:
                triggering_macro_values.append(False)
            for key in triggering_macro_values:
                self._updaters.setdefault(key, []).append(field.updater)
        self.fields[final_name] = field
        self.expression2name[field.expression] = final_name
        field.coordinates = [self.require(f) for f in field.coords]
        return final_name

    def require(self, field: Base) -> str:
        """Ensure that the specified variable (or expression of variables) is included
        in the field collection. This is typically used to add coordinate variables.

        Args:
            expression: variable name or expression of variable(s)
        """
        if field.expression in self.expression2name:
            return self.expression2name[field.expression]
        return self._add_field(field, field.default_name, generate_unique_name=True)

    def update(self, macro: bool = False):
        for updater in self._updaters.get(macro, ()):
            updater()


def grid2dims(grid: pygetm.domain.Grid, z) -> Tuple[str, ...]:
    dims = (f"y{grid.postfix}", f"x{grid.postfix}")
    if z:
        dims = ("zi" if z == INTERFACES else "z",) + dims
    return dims


class Field(Base):
    __slots__ = "collection", "array"

    def __init__(self, array: pygetm.core.Array, dtype: Optional[DTypeLike] = None):
        attrs = {}
        for key, value in array.attrs.items():
            if not key.startswith("_"):
                attrs[key] = value

        self.array = array
        default_time_varying = TimeVarying.MACRO if array.z else TimeVarying.MICRO
        time_varying = array.attrs.get("_time_varying", default_time_varying)
        shape = list(self.array.shape)
        shape[-1] += 2 * array.grid.domain.halox
        shape[-2] += 2 * array.grid.domain.haloy
        super().__init__(
            array.name,
            tuple(shape),
            grid2dims(array.grid, array.z),
            dtype or array.dtype,
            array.fill_value,
            time_varying,
            attrs,
        )

    def get(
        self, out: Optional[ArrayLike] = None, slice_spec: Tuple[int, ...] = ()
    ) -> ArrayLike:
        if out is None:
            return self.array.all_values
        else:
            out[slice_spec] = self.array.all_values
            return out

    @property
    def coords(self) -> Sequence["Base"]:
        if self.grid.domain.spherical:
            yield Field(self.grid.lon)
            yield Field(self.grid.lat)
        else:
            yield Field(self.grid.x)
            yield Field(self.grid.y)
        if self.z:
            yield Field(self.grid.zf if self.z == INTERFACES else self.grid.zc)

    @property
    def z(self) -> bool:
        return self.array.z

    @property
    def default_name(self) -> str:
        return self.array.name

    @property
    def grid(self) -> pygetm.domain.Grid:
        return self.array.grid


class UnivariateTransform(Base):
    __slots__ = "_source"

    def __init__(
        self,
        source: Base,
        shape: Optional[Tuple[int, ...]] = None,
        dims: Optional[Tuple[str, ...]] = None,
        dtype: Optional[DTypeLike] = None,
        expression: Optional[str] = None,
    ):
        if shape is None:
            shape = source.shape
        super().__init__(
            expression or f"{self.__class__.__name__}({source.expression})",
            shape,
            dims or source.dims,
            dtype or source.dtype,
            source.fill_value,
            source.time_varying,
            source.attrs,
        )
        self._source = source

    @property
    def updatable(self) -> bool:
        return self._source.updatable

    @property
    def updater(self) -> Optional[Callable]:
        return self._source.updater

    @property
    def grid(self) -> pygetm.domain.Grid:
        return self._source.grid

    @property
    def coords(self) -> Sequence[Base]:
        yield from self._source.coords

    @property
    def default_name(self) -> str:
        return self._source.default_name


class UnivariateTransformWithData(UnivariateTransform):
    __slots__ = "values"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.values = np.empty(self.shape, self.dtype)

    def get(
        self, out: Optional[ArrayLike] = None, slice_spec: Tuple[int, ...] = ()
    ) -> ArrayLike:
        if out is None:
            return self.values
        else:
            out[slice_spec] = self.values
            return out


class Gather(UnivariateTransform):
    __slots__ = "global_array", "tiling", "_slice", "_gather"

    def __init__(self, source: Base, tiling: pygetm.parallel.Tiling):
        self.tiling = tiling
        nx = tiling.nx_glob + source.shape[-1] - 4 - tiling.nx_sub
        ny = tiling.ny_glob + source.shape[-2] - 4 - tiling.ny_sub
        shape = source.shape[:-2] + (ny, nx)
        super().__init__(source, shape=shape, expression=source.expression)
        self.global_array = None
        if isinstance(source, Field):
            global_domain = source.grid.domain.glob
            if global_domain and not source.time_varying:
                self.global_array = global_domain.fields.get(source.array.name)
        local_shape = source.shape[:-2] + (source.shape[-2] - 4, source.shape[-1] - 4)
        self._slice = (Ellipsis, slice(2, -2), slice(2, -2))
        self._gather = pygetm.parallel.Gather(self.tiling, local_shape, self.dtype)

    def get(
        self, out: Optional[ArrayLike] = None, slice_spec: Tuple[int, ...] = ()
    ) -> ArrayLike:
        # Get data for global array
        # If we have access to the full global field (root rank only), use that,
        # as gathering from subdomains may leave gaps.
        # Nevertheless we cannot skip the gather in that case,
        # because all non-root ranks will call gather anyway.
        out = self._gather(self._source.get()[self._slice], out, slice_spec)
        if self.global_array:
            out[slice_spec] = self.global_array.values
        return out

    @property
    def coords(self) -> Sequence[Base]:
        for c in self._source.coords:
            yield c.gather(self.tiling)

    @property
    def grid(self) -> None:
        return None


class Slice(UnivariateTransform):
    __slots__ = "_slice"

    def __init__(self, source: Field):
        shape = list(source.shape)
        shape[-1] -= 4
        shape[-2] -= 4
        super().__init__(source, shape=tuple(shape), expression=source.expression)
        self._slice = (Ellipsis, slice(2, -2), slice(2, -2))

    def get(
        self, out: Optional[ArrayLike] = None, slice_spec: Tuple[int, ...] = ()
    ) -> ArrayLike:
        values = self._source.get()[self._slice]
        if out is None:
            return values
        else:
            out[slice_spec] = values
            return out

    @property
    def coords(self) -> Sequence[Base]:
        for c in self._source.coords:
            if c.shape[-2:] == self._source.shape[-2:]:
                c = Slice(c)
            yield c

    @property
    def grid(self) -> None:
        return None


class Mask(UnivariateTransformWithData):
    def __init__(self, source: Field):
        super().__init__(source)
        self._mask = source.mask
        assert self._mask.shape == self.shape[-self._mask.ndim :]

    def get(
        self, out: Optional[ArrayLike] = None, slice_spec: Tuple[int, ...] = ()
    ) -> ArrayLike:
        self._source.get(out=self.values)
        self.values[..., self._mask] = self.fill_value
        return super().get(out, slice_spec)


class TimeAverage(UnivariateTransformWithData):
    __slots__ = ("_n",)

    def __init__(self, source: Field):
        super().__init__(source)
        self._n = 0
        if "cell_methods" in self.attrs:
            self.attrs["cell_methods"] += " time: mean"
        else:
            self.attrs["cell_methods"] = "time: mean"
        self.values.fill(self.fill_value)

    @property
    def updatable(self) -> bool:
        if self._source.time_varying == TimeVarying.MACRO:
            return Updatable.MACRO_ONLY
        return Updatable.ALWAYS

    @property
    def updater(self) -> Optional[Callable]:
        return self.update

    def update(self):
        if self._n == 0:
            self._source.get(out=self.values)
        else:
            self.values += self._source.get()
        self._n += 1

    def get(
        self, out: Optional[ArrayLike] = None, slice_spec: Tuple[int, ...] = ()
    ) -> ArrayLike:
        if self._n > 0:
            self.values *= 1.0 / self._n
        self._n = 0
        return super().get(out, slice_spec)

    @property
    def coords(self):
        for c in super().coords:
            if c.time_varying:
                c = TimeAverage(c)
            yield c

    @property
    def default_name(self) -> str:
        return self._source.default_name + "_av"


class Regrid(UnivariateTransformWithData):
    __slots__ = ("interpolate", "_grid", "_z", "_slice")

    def __init__(
        self,
        source: Base,
        grid: Optional[pygetm.domain.Grid] = None,
        z: Optional[Literal[None, CENTERS, INTERFACES]] = None,
    ):
        assert source.grid is not None
        grid = grid or source.grid
        if grid is not source.grid:
            assert z is None
            self.interpolate = source.grid.interpolator(grid)
            shape = source.shape[:-2] + (grid.ny_, grid.nx_)
            args = f", grid={grid.postfix}"
            if source.ndim > 2:
                z = CENTERS if source.shape[0] == grid.nz else INTERFACES
        else:
            assert z is not None
            if z == CENTERS:
                assert source.shape[0] == source.grid.nz + 1
                self.interpolate = functools.partial(pygetm._pygetm.interp_z, offset=0)
                shape = (source.shape[0] - 1,) + source.shape[1:]
                args = ", z=centers"
            else:
                assert source.shape[0] == source.grid.nz
                self.interpolate = functools.partial(pygetm._pygetm.interp_z, offset=1)
                shape = (source.shape[0] + 1,) + source.shape[1:]
                args = ", z=interfaces"
        self._slice = (np.newaxis, Ellipsis) if source.ndim < 3 else (Ellipsis,)
        dims = grid2dims(grid, z)
        expression = f"{self.__class__.__name__}({source.expression}{args})"
        self._grid = grid
        self._z = z
        super().__init__(source, shape=shape, dims=dims, expression=expression)

    def get(
        self, out: Optional[ArrayLike] = None, slice_spec: Tuple[int, ...] = ()
    ) -> ArrayLike:
        self.interpolate(self._source.get()[self._slice], self.values[self._slice])
        return super().get(out, slice_spec)

    @property
    def grid(self) -> pygetm.domain.Grid:
        return self._grid

    @property
    def coords(self):
        if self._grid.domain.spherical:
            yield Field(self._grid.lon)
            yield Field(self._grid.lat)
        else:
            yield Field(self._grid.x)
            yield Field(self._grid.y)
        if self._z:
            yield Field(self._grid.zf if self._z == INTERFACES else self._grid.zc)


class InterpZ(UnivariateTransformWithData):
    __slots__ = ("z_src", "z_tgt", "z_dim")

    def __init__(self, source: Base, z: ArrayLike, dim: str):
        assert source.ndim == 3
        self.z_tgt = np.asarray(z, dtype=float)
        shape = (self.z_tgt.size,) + source.shape[-2:]
        dims = (dim,) + source.dims[1:]
        self.z_dim = dim
        expression = f"{self.__class__.__name__}({source.expression}, z={self.z_tgt})"
        at_centers = source.shape[0] == source.grid.nz
        self.z_src = (source.grid.zc if at_centers else source.grid.zf).all_values
        super().__init__(source, shape=shape, dims=dims, expression=expression)

    def get(
        self, out: Optional[ArrayLike] = None, slice_spec: Tuple[int, ...] = ()
    ) -> ArrayLike:
        ip = pygetm.util.interpolate.LinearVectorized1D(
            self.z_tgt, self.z_src, 0, self.fill_value
        )
        self.values[...] = ip(self._source.get())
        return super().get(out, slice_spec)

    @property
    def coords(self):
        for c in super().coords:
            if c.attrs.get("axis") == "Z":
                c = WrappedArray(self.z_tgt, self.z_dim, (self.z_dim,))
            yield c
