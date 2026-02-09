from typing import (
    Iterable,
    MutableMapping,
    Union,
    Optional,
    Mapping,
    Literal,
    Callable,
    Any,
    NamedTuple,
    TypeVar,
)
import collections
import enum
import functools

import numpy as np
from numpy.typing import DTypeLike, ArrayLike

import pygetm.core
import pygetm._pygetm
import pygetm.util.interpolate
import pygetm.parallel
from pygetm.constants import (
    CENTERS,
    INTERFACES,
    TimeVarying,
    CoordinateType,
    EdgeTreatment,
)


class GridInfo(NamedTuple):
    grid: pygetm.core.Grid
    z: Optional[Literal[CENTERS, INTERFACES]]
    on_boundary: bool

    def _get_gather_info(
        self, source: "Base"
    ) -> tuple[Callable, tuple, Mapping[str, Any]]:
        gatherer, local_slice, global_shape = self.grid.get_gather_info(
            shape=source.shape,
            on_boundary=self.on_boundary,
            dtype=source.dtype,
            fill_value=source.fill_value,
        )
        return gatherer, local_slice, dict(shape=global_shape)

    def get_coords(self, ndim: int) -> Iterable["Base"]:
        for array in self.grid._get_coords(ndim, self.z, self.on_boundary):
            yield Field(array)
        yield from self.grid.extra_output_coordinates


T = TypeVar("T")


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
        "_global_values",
    )

    @classmethod
    def parameterize(cls: type[T], **kwargs) -> T:
        return functools.partial(cls, **kwargs)

    def __init__(
        self,
        expression: str,
        shape: Iterable[int],
        dims: Iterable[str],
        dtype: DTypeLike,
        fill_value=None,
        time_varying: Union[TimeVarying, Literal[False]] = TimeVarying.MICRO,
        attrs: Mapping[str, Any] = {},
    ):
        self.expression = expression
        self.shape = tuple(shape)
        self.ndim = len(self.shape)
        self.dims = tuple(dims)
        assert self.ndim == len(
            self.dims
        ), f"Expected {self.ndim} dimensions but got {self.dims}"
        self.dtype = np.dtype(dtype)
        self.fill_value = fill_value
        self.attrs = attrs
        self.time_varying = time_varying
        self.coordinates: list[str] = []

    def get(
        self, out: Optional[ArrayLike] = None, slice_spec: tuple[int, ...] = ()
    ) -> ArrayLike:
        raise NotImplementedError

    @property
    def default_name(self) -> str:
        raise NotImplementedError

    def gather(self) -> "Base":
        return Gather(self)

    @property
    def updatable(self) -> bool:
        return False

    @property
    def updater(self) -> Optional[Callable]:
        raise NotImplementedError

    @property
    def grid_info(self) -> Optional[GridInfo]:
        return None

    def _get_gather_info(
        self, source: "Base"
    ) -> tuple[Callable, tuple, Mapping[str, Any]]:
        return self.grid_info._get_gather_info(source)

    @property
    def mask(self) -> np.ndarray:
        return np.broadcast_to(False, self.shape)

    @property
    def coords(self) -> Iterable["Base"]:
        yield from self.grid_info.get_coords(self.ndim)


class WrappedArray(Base):
    __slots__ = ("_name", "values", "_global_field")

    def __init__(
        self,
        values: np.ndarray,
        name: str,
        dims: tuple[str, ...],
        global_field: Optional[Base] = None,
        **kwargs,
    ):
        super().__init__(
            name, values.shape, dims, values.dtype, time_varying=False, **kwargs
        )
        self._name = name
        self.values = values
        self._global_field = global_field or self

    def gather(self) -> Base:
        return self._global_field

    @property
    def default_name(self) -> str:
        return self._name

    @property
    def coords(self) -> Iterable[Base]:
        return ()

    def get(
        self, out: Optional[ArrayLike] = None, slice_spec: tuple[int, ...] = ()
    ) -> ArrayLike:
        if out is None:
            return self.values
        else:
            out[slice_spec] = self.values
            return out


class Updatable(enum.Enum):
    ALWAYS = 1
    MACRO_ONLY = 2


class FieldCollection(Mapping[str, Base]):
    def __init__(
        self,
        available_fields: Mapping[str, pygetm.core.Array],
        default_dtype: Optional[DTypeLike] = None,
        sub: bool = False,
    ):
        self.fields: MutableMapping[str, Base] = collections.OrderedDict()
        self.expression2name: MutableMapping[str, str] = {}
        self.available_fields = available_fields
        self.default_dtype = default_dtype
        self.sub = sub
        self._updaters = {}

    def __getitem__(self, key: str) -> Base:
        return self.fields[key]

    def __iter__(self):
        return iter(self.fields)

    def __len__(self) -> int:
        return len(self.fields)

    def request(
        self,
        *fields: Union[str, pygetm.core.Array],
        output_name: Optional[str] = None,
        dtype: Optional[DTypeLike] = None,
        mask: Optional[bool] = None,
        time_average: bool = False,
        grid: Optional[pygetm.core.Grid] = None,
        z: Union[
            None, Literal[CENTERS], Literal[INTERFACES], float, Iterable[float]
        ] = None,
        generate_unique_name: bool = False,
        transforms: Iterable[type["UnivariateTransform"]] = (),
    ) -> tuple[str, ...]:
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
            grid: if provided, the field will be regridded to this grid before saving
            z: if provided, the field will be interpolated to these z levels before saving.
                This can be CENTERS or INTERFACES to indicate interpolation to the center
                or interfaces of the model grid, a 1D array of custom z levels, or a single
                float indicating the z level to interpolate to. Note that z levels are
                negative downward, with the bottom being -H.
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
            for tf in transforms:
                field = tf(field)
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
                field = field.gather()
            for tf in source_grid.default_output_transforms:
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
        return final_name

    def add_coordinates(self):
        for field in list(self.fields.values()):
            field.coordinates.extend(self.require(f) for f in field.coords)

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


class Field(Base):
    __slots__ = "collection", "array"

    def __init__(self, array: pygetm.core.Array, dtype: Optional[DTypeLike] = None):
        attrs = {}
        for key, value in array.attrs.items():
            if not key.startswith("_"):
                attrs[key] = value
        if "_global_values" in array.attrs:
            self._global_values = array.attrs["_global_values"]

        self.array = array
        default_time_varying = TimeVarying.MACRO if array.z else TimeVarying.MICRO
        time_varying = array.attrs.get("_time_varying", default_time_varying)
        shape = list(self.array.shape)
        if array.ndim >= 2 and not array.on_boundary:
            shape[-1] += 2 * array.grid.halox
            shape[-2] += 2 * array.grid.haloy
        super().__init__(
            array.name,
            shape,
            array.grid._get_dims(array.ndim, array.z, array.on_boundary),
            dtype or array.dtype,
            array.fill_value,
            time_varying,
            attrs,
        )

    def get(
        self, out: Optional[ArrayLike] = None, slice_spec: tuple[int, ...] = ()
    ) -> ArrayLike:
        if out is None:
            return self.array.all_values
        else:
            out[slice_spec] = self.array.all_values
            return out

    @property
    def default_name(self) -> str:
        return self.array.name

    @property
    def grid_info(self) -> GridInfo:
        return GridInfo(self.array.grid, self.array.z, self.array.on_boundary)

    @property
    def mask(self) -> np.ndarray:
        return self.array.all_mask


class UnivariateTransform(Base):
    __slots__ = "_source", "_grid_info_provider"

    def __init__(
        self,
        source: Base,
        shape: Optional[tuple[int, ...]] = None,
        dims: Optional[Iterable[str]] = None,
        dtype: Optional[DTypeLike] = None,
        expression: Optional[str] = None,
        inherit_grid_info: bool = True,
    ):
        if shape is None:
            shape = source.shape
        super().__init__(
            expression or f"{self.__class__.__name__}({source.expression})",
            shape,
            source.dims if dims is None else dims,
            dtype or source.dtype,
            source.fill_value,
            source.time_varying,
            source.attrs,
        )
        self._source = source
        self._grid_info_provider = source if inherit_grid_info else super()

    @property
    def default_name(self) -> str:
        return self._source.default_name

    @property
    def updatable(self) -> bool:
        return self._source.updatable

    @property
    def updater(self) -> Optional[Callable]:
        return self._source.updater

    @property
    def grid_info(self) -> Optional[GridInfo]:
        return self._grid_info_provider.grid_info

    @property
    def mask(self) -> np.ndarray:
        return self._source.mask

    @property
    def coords(self) -> Iterable[Base]:
        yield from self._grid_info_provider.coords

    def _get_gather_info(
        self, source: "Base"
    ) -> tuple[Callable, tuple, Mapping[str, Any]]:
        return self._grid_info_provider._get_gather_info(source)


class UnivariateTransformWithData(UnivariateTransform):
    __slots__ = "values"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.values = np.empty(self.shape, self.dtype)

    def get(
        self, out: Optional[ArrayLike] = None, slice_spec: tuple[int, ...] = ()
    ) -> ArrayLike:
        if out is None:
            return self.values
        else:
            out[slice_spec] = self.values
            return out


class Gather(UnivariateTransform):
    __slots__ = "root_has_global_values", "_slice", "_gather"

    def __init__(self, source: Base):
        gatherer, local_slice, kwargs = source._get_gather_info(source)
        super().__init__(source, expression=source.expression, **kwargs)
        self.root_has_global_values = hasattr(source, "_global_values")
        self._slice = local_slice
        self._gather = gatherer

    def get(
        self, out: Optional[ArrayLike] = None, slice_spec: tuple[int, ...] = ()
    ) -> ArrayLike:
        if self.root_has_global_values:
            global_values = self._source._global_values
            if global_values is not None:
                assert global_values.shape == self.shape
                if out is None:
                    return global_values
                out[slice_spec] = global_values
        else:
            local_interior = self._source.get()[self._slice]
            out = self._gather(local_interior, out, slice_spec)
        return out

    @property
    def coords(self) -> Iterable[Base]:
        for c in self._source.coords:
            yield c.gather()

    @property
    def grid_info(self) -> None:
        return None


class Mask(UnivariateTransformWithData):
    def __init__(self, source: Field):
        super().__init__(source)
        self._mask = source.mask
        assert self._mask.shape == self.shape[-self._mask.ndim :]

    def get(
        self, out: Optional[ArrayLike] = None, slice_spec: tuple[int, ...] = ()
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
        self, out: Optional[ArrayLike] = None, slice_spec: tuple[int, ...] = ()
    ) -> ArrayLike:
        if self._n > 0:
            self.values *= 1.0 / self._n
        self._n = 0
        return super().get(out, slice_spec)

    @property
    def coords(self) -> Iterable[Base]:
        for c in super().coords:
            if c.time_varying:
                c = TimeAverage(c)
            yield c

    @property
    def default_name(self) -> str:
        return self._source.default_name + "_av"


class Regrid(UnivariateTransformWithData):
    __slots__ = ("interpolate", "_grid_info", "_slice")

    def __init__(
        self,
        source: Base,
        grid: Optional[pygetm.core.Grid] = None,
        z: Optional[Literal[None, CENTERS, INTERFACES]] = None,
    ):
        source_grid_info = source.grid_info
        assert source_grid_info is not None
        assert not source_grid_info.on_boundary
        grid = grid or source_grid_info.grid
        if grid is not source_grid_info.grid:
            assert z is None
            self.interpolate = source_grid_info.grid.interpolator(grid)
            shape = source.shape[:-2] + (grid.ny_, grid.nx_)
            args = f", grid={grid.postfix}"
            if source.ndim > 2:
                z = source_grid_info.z
        else:
            assert z is not None
            if z == CENTERS:
                assert source_grid_info.z == INTERFACES
                self.interpolate = functools.partial(pygetm._pygetm.interp_z, offset=0)
                shape = (source.shape[0] - 1,) + source.shape[1:]
                args = ", z=centers"
            else:
                assert source_grid_info.z == CENTERS
                self.interpolate = functools.partial(pygetm._pygetm.interp_z, offset=1)
                shape = (source.shape[0] + 1,) + source.shape[1:]
                args = ", z=interfaces"
        super().__init__(
            source,
            shape=shape,
            dims=grid._get_dims(len(shape), z),
            expression=f"{self.__class__.__name__}({source.expression}{args})",
            inherit_grid_info=False,
        )
        self._slice = (np.newaxis, Ellipsis) if source.ndim < 3 else (Ellipsis,)
        self._grid_info = GridInfo(grid, z, False)

    def get(
        self, out: Optional[ArrayLike] = None, slice_spec: tuple[int, ...] = ()
    ) -> ArrayLike:
        self.interpolate(self._source.get()[self._slice], self.values[self._slice])
        return super().get(out, slice_spec)

    @property
    def grid_info(self) -> GridInfo:
        return self._grid_info

    @property
    def mask(self) -> np.ndarray:
        source_mask = self._source.mask.astype(float, order="C")
        mask = np.zeros(self.shape, dtype=float)
        self.interpolate(source_mask[self._slice], mask[self._slice])
        return mask > 0.99


class InterpZ(UnivariateTransformWithData):
    """Interpolate in the vertical.
    The vertical dimension must be the first dimension of the source array.
    The source array must have a vertical coordinate with ``axis`` attribute equal
    to ``Z``. This coordinate must have the same shape as the source array.
    """

    __slots__ = ("z_src", "z_tgt", "z_dim", "edges")

    def __init__(
        self,
        source: Base,
        z: ArrayLike,
        dim: str,
        edges: EdgeTreatment = EdgeTreatment.MISSING,
    ):
        self.z_tgt = np.asarray(z, dtype=float)
        shape = (self.z_tgt.size,) + source.shape[1:]
        dims = (dim,) + source.dims[1:]
        self.z_dim = dim
        self.edges = edges
        expression = f"{self.__class__.__name__}({source.expression}, z={self.z_tgt})"
        super().__init__(source, shape=shape, dims=dims, expression=expression)
        for c in source.coords:
            if c.attrs.get("axis") == "Z":
                self.z_src = c.get()
                break
        else:
            raise ValueError("Could not find source z coordinate")

    def get(
        self, out: Optional[ArrayLike] = None, slice_spec: tuple[int, ...] = ()
    ) -> ArrayLike:
        ip = pygetm.util.interpolate.LinearVectorized1D(
            self.z_tgt, self.z_src, 0, self.fill_value, edges=self.edges
        )
        self.values[...] = ip(self._source.get())
        return super().get(out, slice_spec)

    @property
    def coords(self) -> Iterable[Base]:
        for c in super().coords:
            if c.attrs.get("axis") == "Z":
                c = WrappedArray(self.z_tgt, self.z_dim, (self.z_dim,), attrs=c.attrs)
            yield c

    @property
    def mask(self) -> np.ndarray:
        return self._source.mask.all(axis=0)


class IndexXY(UnivariateTransform):
    __slots__ = "_i", "_j", "_slice", "_index_dims", "_local2global", "_index_coords"

    def __init__(
        self,
        source: Base,
        x: ArrayLike,
        y: ArrayLike,
        coordinate_type: CoordinateType = CoordinateType.IJ,
        dims: tuple[str, ...] = (),
        coords: Mapping[str, Union[Base, ArrayLike]] = {},
    ):
        index_shape = np.broadcast_shapes(np.shape(x), np.shape(y))
        assert len(dims) == len(index_shape)
        self._index_dims = dims

        source_grid_info = source.grid_info
        assert source_grid_info is not None and not source_grid_info.on_boundary
        grid = source_grid_info.grid
        self._i = np.empty(index_shape, dtype=np.intp)
        self._j = np.empty(index_shape, dtype=np.intp)
        if grid.tiling.rank == 0:
            all_mask = grid.mask.attrs["_global_values"]
            all_x = None if grid.x is None else grid.x.attrs["_global_values"]
            all_y = None if grid.y is None else grid.y.attrs["_global_values"]
            all_lon = None if grid.lon is None else grid.lon.attrs["_global_values"]
            all_lat = None if grid.lat is None else grid.lat.attrs["_global_values"]
            loc = pygetm.core.Locator(all_mask, all_x, all_y, all_lon, all_lat)
            self._i[...], self._j[...] = loc(x, y, coordinate_type=coordinate_type)
        grid.tiling.comm.Bcast(self._i)
        grid.tiling.comm.Bcast(self._j)

        if hasattr(source, "_global_values"):
            # The root (rank=0) has the global values of the source array
            if source._global_values is not None:
                # We are the root - extract the values at the desired indices
                self._global_values = source._global_values[..., self._j, self._i]
            else:
                # We are not the root - just set the attribute to flag that the root
                # does have the global values
                self._global_values = None

        # Determine which of the desired indices fall within the local subdomain
        # (excluding halos)
        i = self._i - grid.tiling.xoffset
        j = self._j - grid.tiling.yoffset
        inside = (i >= 0) & (i < grid.nx) & (j >= 0) & (j < grid.ny)

        # Determine slice needed to extract local values of interest, and mapping from
        # the resulting local values to the global array with desired indices.
        # At this point both local and global arrays are flattened (1D).
        i_in = i[inside] + grid.halox
        j_in = j[inside] + grid.haloy
        self._slice = (Ellipsis, j_in, i_in)
        self._local2global = np.flatnonzero(inside)

        flat_shape = source.shape[:-2] + self._local2global.shape
        flat_dim_name = "_".join(dims) if dims else "station"
        flat_dims = source.dims[:-2] + (flat_dim_name,)
        super().__init__(
            source, shape=flat_shape, dims=flat_dims, inherit_grid_info=False
        )

        self._index_coords = []
        for name, values in coords.items():
            values = np.broadcast_to(values, index_shape)
            global_field = WrappedArray(values, name, dims=dims)
            local_field = WrappedArray(
                values[inside], name, dims=(flat_dim_name,), global_field=global_field
            )
            self._index_coords.append(local_field)

    def get(
        self, out: Optional[ArrayLike] = None, slice_spec: tuple[int, ...] = ()
    ) -> ArrayLike:
        values = self._source.get()[self._slice]
        if out is None:
            return values
        else:
            out[slice_spec] = values
            return out

    @property
    def coords(self) -> Iterable[Base]:
        for c in self._source.coords:
            if c.dims[-2:] == self._source.dims[-2:]:
                c = IndexXY(c, self._i, self._j, dims=self._index_dims)
            yield c
        for c in self._index_coords:
            yield c

    @property
    def mask(self) -> np.ndarray:
        return self._source.mask[self._slice]

    def _get_gather_info(
        self, source: Base
    ) -> tuple[Callable, tuple, Mapping[str, Any]]:
        """Gather across all processes and reshape the flattened
        data to the original index shape"""

        class gather_serial:
            def __init__(self, final_shape):
                self.final_shape = final_shape

            def __call__(
                self, locvalues, globvalues: Optional[np.ndarray] = None, globslice=()
            ):
                locvalues = np.reshape(locvalues, self.final_shape)
                if globvalues is not None:
                    globvalues[globslice] = locvalues
                    return globvalues
                return locvalues

        tiling = self._source.grid_info.grid.tiling
        final_shape = source.shape[:-1] + self._i.shape
        final_dims = source.dims[:-1] + self._index_dims
        if tiling.n == 1:
            gatherer = gather_serial(final_shape)
        else:
            gatherer = pygetm.parallel.GatherFromIndices(
                tiling.comm,
                self._local2global,
                self._i.shape,
                source.shape[:-1],
                source.dtype,
                fill_value=source.fill_value,
            )

        return gatherer, (Ellipsis,), dict(shape=final_shape, dims=final_dims)

    @property
    def grid_info(self) -> None:
        return None
