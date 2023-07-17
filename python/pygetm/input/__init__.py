from typing import (
    Callable,
    Iterable,
    List,
    Mapping,
    Union,
    Optional,
    Sequence,
    Tuple,
    TYPE_CHECKING,
)
import glob
import numbers
import logging
import enum
import functools

import numpy as np
import numpy.typing
import numpy.lib.mixins
import xarray
import cftime

import pygetm.util.interpolate
from pygetm.constants import CENTERS, TimeVarying

if TYPE_CHECKING:
    import pygetm.core

LATITUDE_UNITS = (
    "degrees_north",
    "degree_north",
    "degree_N",
    "degrees_N",
    "degreeN",
    "degreesN",
)

LONGITUDE_UNITS = (
    "degrees_east",
    "degree_east",
    "degree_E",
    "degrees_E",
    "degreeE",
    "degreesE",
)

Z_STANDARD_NAMES = (
    "height",
    "height_above_mean_sea_level",
    "depth",
    "depth_below_geoid",
)


@xarray.register_dataarray_accessor("getm")
class GETMAccessor:
    def __init__(self, xarray_obj: xarray.DataArray):
        self._obj = xarray_obj

    @property
    def longitude(self) -> Optional[xarray.DataArray]:
        return self.coordinates.get("longitude")

    @property
    def latitude(self) -> Optional[xarray.DataArray]:
        return self.coordinates.get("latitude")

    @property
    def z(self) -> Optional[xarray.DataArray]:
        return self.coordinates.get("z")

    @property
    def time(self) -> Optional[xarray.DataArray]:
        return self.coordinates.get("time")

    @functools.cached_property
    def coordinates(self) -> Mapping[str, xarray.DataArray]:
        _coordinates = {}
        for name, coord in self._obj.coords.items():
            units = coord.attrs.get("units")
            standard_name = coord.attrs.get("standard_name")
            if standard_name in ("latitude", "longitude"):
                _coordinates[standard_name] = coord
            elif coord.attrs.get("positive") or standard_name in Z_STANDARD_NAMES:
                _coordinates["z"] = coord
            elif units in LATITUDE_UNITS:
                _coordinates["latitude"] = coord
            elif units in LONGITUDE_UNITS:
                _coordinates["longitude"] = coord
            elif coord.size > 0 and isinstance(coord.values.flat[0], cftime.datetime):
                _coordinates["time"] = coord
            elif name == "zax":
                _coordinates["z"] = coord
        return _coordinates


open_nc_files = []


def _open(path, preprocess=None, **kwargs):
    key = (path, preprocess, kwargs.copy())
    for k, ds in open_nc_files:
        if k == key:
            return ds
    ds = xarray.open_dataset(path, **kwargs)
    if preprocess:
        ds = preprocess(ds)
    open_nc_files.append((key, ds))
    return ds


def from_nc(
    paths: Union[str, Sequence[str]],
    name: str,
    preprocess: Optional[Callable[[xarray.Dataset], xarray.Dataset]] = None,
    **kwargs,
) -> xarray.DataArray:
    """Obtain a variable from one or more NetCDF files that can be used as value
    provided to :meth:`InputManager.add` and :meth:`pygetm.core.Array.set`.

    Args:
        paths: single file path, a pathname pattern containing `*` and/or `?`, or a
            sequence of file paths. If multiple paths are provided (or the pattern
            resolves to multiple valid path names), the files will be concatenated
            along their time dimension.
        preprocess: function that transforms the :class:`xarray.Dataset` opened for
            every path provided. This can be used to modify the datasets before
            concatenation in time is attempted, for instance, to cut off time indices
            that overlap between files.
        **kwargs: additional keyword arguments to be passed to
            :func:`xarray.open_dataset`
    """
    kwargs.setdefault("decode_times", True)
    kwargs["use_cftime"] = True
    kwargs["cache"] = False
    if isinstance(paths, str):
        pattern = paths
        paths = glob.glob(pattern)
        if not paths:
            raise Exception(f"No files found matching {pattern!r}")
    arrays = []
    for path in paths:
        ds = _open(path, preprocess, **kwargs)
        array = ds[name]
        # Note: we wrap the netCDF array ourselves, in order to support lazy operators
        # (e.g., add, multiply)
        lazyvar = Wrap(array.variable, name=f"from_nc({path!r}, {name!r})")
        array = xarray.DataArray(
            lazyvar,
            dims=array.dims,
            coords=array.coords,
            attrs=array.attrs,
            name=lazyvar.name,
        )
        arrays.append(array)
    if len(arrays) == 1:
        return arrays[0]
    else:
        assert all(array.getm.time is not None for array in arrays)
        return xarray.concat(
            sorted(arrays, key=lambda a: a.getm.time.values.flat[0]),
            dim=arrays[0].getm.time.dims[0],
            coords="minimal",
            combine_attrs="drop_conflicts",
        )


class LazyArray(numpy.lib.mixins.NDArrayOperatorsMixin):
    def __init__(self, shape: Iterable[int], dtype: numpy.typing.DTypeLike, name: str):
        self.shape = tuple(shape)
        self.ndim = len(self.shape)
        self.dtype = dtype
        self.name = name

    def update(self, time: cftime.datetime, numtime: np.longdouble) -> bool:
        return False

    def astype(self, dtype, **kwargs) -> np.ndarray:
        if dtype == self.dtype:
            return self
        return self.__array__(dtype)

    def __array_function__(self, func, types, args, kwargs):
        if func == np.result_type:
            args = tuple(x.dtype if isinstance(x, LazyArray) else x for x in args)
            return np.result_type(*args)
        if func == np.concatenate:
            return Concatenate(*args, **kwargs)
        args = tuple(np.asarray(x) if isinstance(x, LazyArray) else x for x in args)
        kwargs = dict(
            (k, np.asarray(v)) if isinstance(v, LazyArray) else (k, v)
            for (k, v) in kwargs.items()
        )
        return func(*args, **kwargs)

    def __array_ufunc__(self, ufunc, method: str, *inputs, **kwargs):
        if method != "__call__":
            return NotImplemented

        if "out" in kwargs:
            return NotImplemented

        COMPATIBLE_TYPES = (np.ndarray, numbers.Number, LazyArray, xarray.Variable)
        for x in inputs:
            if not isinstance(x, COMPATIBLE_TYPES):
                return NotImplemented

        return UFunc(ufunc, method, *inputs, **kwargs)

    def __array__(self, dtype=None) -> np.ndarray:
        raise NotImplementedError

    def __getitem__(self, slices) -> np.ndarray:
        return self.__array__()[slices]

    def is_time_varying(self) -> bool:
        return False

    def _finalize_slices(self, slices: Tuple):
        assert isinstance(slices, tuple)
        for i, s in enumerate(slices):
            if s is Ellipsis:
                slices = (
                    slices[:i]
                    + (slice(None),) * (self.ndim + 1 - len(slices))
                    + slices[i + 1 :]
                )
                break
        assert len(slices) == self.ndim
        return slices


class Operator(LazyArray):
    _operator_name = None

    def __init__(
        self,
        *inputs,
        passthrough=(),
        dtype=None,
        shape=None,
        name: Optional[str] = None,
        **kwargs,
    ):
        # Unpack unnamed arguments
        self.inputs = []
        self.lazy_inputs = []
        self.input_names = []
        for inp in inputs:
            assert isinstance(
                inp, (np.ndarray, numbers.Number, LazyArray, xarray.Variable)
            ), f"Input has unknown type {type(inp)}"

            # Unpack to LazyArray if possible
            if isinstance(inp, xarray.Variable) and isinstance(inp._data, LazyArray):
                inp = inp._data

            if isinstance(inp, LazyArray):
                self.input_names.append(inp.name)
            elif isinstance(inp, (np.ndarray, numbers.Number)):
                self.input_names.append(str(inp))
            else:
                # Other datatype, typically xarray.Variable.
                # Do not call object's custom str/repr, as that will cause evaluation
                # (e.g. read from file) of the entire array
                self.input_names.append(object.__repr__(inp))

            # If this is a Wrap, unwrap
            # (the wrapping was only for ufunc support)
            if isinstance(inp, Wrap):
                inp = inp._source

            if isinstance(inp, LazyArray):
                self.lazy_inputs.append(inp)
            self.inputs.append(inp)

        # Store keyword arguments as-is (no unpacking)
        self.kwargs = kwargs

        # Infer shape from inputs if not provided
        if shape is None:
            shapes = []
            for input in inputs:
                if isinstance(
                    input, (np.ndarray, LazyArray, xarray.DataArray, xarray.Variable)
                ):
                    shapes.append(input.shape)
            shape = np.broadcast_shapes(*shapes)

        for i in range(len(self.inputs)):
            if isinstance(self.inputs[i], np.ndarray):
                self.inputs[i] = np.broadcast_to(self.inputs[i], shape)

        # Process dimensions for which we can passthrough slices to inputs
        # This can be True (= all dimensions), an iterable, or a dictionary mapping
        # sliced dimensions to input dimensions (if the current operator adds or
        # removes dimensions)
        if passthrough is True:
            passthrough = range(len(shape))
        if not isinstance(passthrough, dict):
            passthrough = dict([(i, i) for i in passthrough])
        self.passthrough = passthrough
        assert all([isinstance(dim, int) for dim in self.passthrough]), (
            f"Invalid passthrough: {self.passthrough}."
            " All entries should be of type int"
        )

        # Generate a name for the variable if not provided
        if name is None:
            operator_name = self._operator_name or self.__class__.__name__
            strargs = ", ".join(self.input_names)
            strkwargs = "".join(f", {k}={v!r}" for (k, v) in self.kwargs.items())
            name = f"{operator_name}({strargs}{strkwargs})"

        super().__init__(shape, dtype or float, name)

    def update(self, *args) -> bool:
        updated = False
        for input in self.lazy_inputs:
            updated = input.update(*args) or updated
        return updated

    def is_time_varying(self) -> bool:
        return self.lazy_inputs and any(
            input.is_time_varying() for input in self.lazy_inputs
        )

    def __getitem__(self, slices) -> np.ndarray:
        preslices, postslices = [], []
        for i, slc in enumerate(self._finalize_slices(slices)):
            if i in self.passthrough:
                preslices.append(slc)
                if not isinstance(slc, (int, np.integer)):
                    postslices.append(slice(None))
            else:
                preslices.append(slice(None))
                postslices.append(slc)
        inputs = [
            inp
            if isinstance(inp, numbers.Number)
            else np.asarray(inp[tuple(preslices)])
            for inp in self.inputs
        ]
        return self.apply(*inputs)[tuple(postslices)]

    def apply(self, *inputs, dtype=None) -> np.ndarray:
        raise NotImplementedError

    def __array__(self, dtype=None) -> np.ndarray:
        return self.apply(*[np.asarray(inp) for inp in self.inputs], dtype=dtype)


class UnaryOperator(Operator):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._source = self.inputs[0]
        self._source_name = self.input_names[0]


class UFunc(Operator):
    def __init__(self, ufunc, method: str, *inputs, **kwargs):
        self._operator_name = ufunc.__name__
        super().__init__(*inputs, passthrough=True, **kwargs)
        self.ufunc = getattr(ufunc, method)

    def apply(self, *inputs, dtype=None) -> np.ndarray:
        return self.ufunc(*inputs, **self.kwargs)


class Wrap(UnaryOperator):
    def __init__(self, source: xarray.Variable, name: str, **kwargs):
        assert isinstance(source, xarray.Variable)
        super().__init__(
            source, passthrough=True, dtype=source.dtype, name=name, **kwargs
        )

    def apply(self, source, dtype=None) -> np.ndarray:
        return source

    def update(self, *args) -> bool:
        return False


class Slice(UnaryOperator):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._slices = []

    def __array__(self, dtype=None) -> np.ndarray:
        data = np.empty(self.shape, dtype or self.dtype)
        for src_slice, tgt_slice in self._slices:
            data[tgt_slice] = self._source[src_slice]
        return data

    def __getitem__(self, slices) -> np.ndarray:
        slices = self._finalize_slices(slices)
        shape = []
        for i, (l, s) in enumerate(zip(self.shape, slices)):
            if i in self.passthrough and isinstance(s, (int, np.integer)):
                # This dimension will be sliced out
                continue
            assert isinstance(s, slice), (
                f"Dimension {i} has unsupported slice type {type(s)} with value {s!r}."
                f" Passthrough: {list(self.passthrough)}"
            )
            start, stop, step = s.indices(l)
            assert i in self.passthrough or (start == 0 and stop == l and step == 1), (
                f"invalid slice for dimension {i} with length {l}:"
                f" {start}:{stop}:{step}"
            )
            shape.append(len(range(start, stop, step)))
        data = np.empty(shape, self.dtype)
        for src_slice, tgt_slice in self._slices:
            src_slice = list(src_slice)
            for iout, iin in self.passthrough.items():
                src_slice[iin] = slices[iout]
            tgt_slice = tuple(
                [
                    ori
                    for i, (cust, ori) in enumerate(zip(slices, tgt_slice))
                    if i not in self.passthrough
                    or not isinstance(cust, (int, np.integer))
                ]
            )
            data[tgt_slice] = self._source[tuple(src_slice)]
        return data


class Concatenate(UnaryOperator):
    def __init__(self, arrays, axis: int = 0, *args, **kwargs):
        shape = list(arrays[0].shape)
        for array in arrays[1:]:
            shape[axis] += array.shape[axis]
            assert all(array.shape[i] == l for i, l in enumerate(shape) if i != axis)
        self.axis = axis
        super().__init__(*arrays, shape=shape, **kwargs)

    def __array__(self, dtype=None) -> np.ndarray:
        return np.concatenate(self.inputs, axis=self.axis, dtype=dtype)

    def __getitem__(self, slices) -> np.ndarray:
        slices = list(self._finalize_slices(slices))
        if (
            isinstance(slices[self.axis], slice)
            and slices[self.axis].start is None
            and slices[self.axis].stop is None
        ):
            inputs = [input[tuple(slices)] for input in self.inputs]
            return np.concatenate(inputs, axis=self.axis)
        assert isinstance(
            slices[self.axis], (int, np.integer)
        ), f"Unsupported slice for concatenated dimension: {slices[self.axis]!r}"
        if slices[self.axis] < 0:
            slices[self.axis] += self.shape[self.axis]
        assert slices[self.axis] >= 0 and slices[self.axis] < self.shape[self.axis]
        for input in self.inputs:
            if slices[self.axis] < input.shape[self.axis]:
                return np.asarray(input[tuple(slices)])
            slices[self.axis] -= input.shape[self.axis]
        assert False, "Index out of bounds?"


def limit_region(
    source: xarray.DataArray,
    minlon: float,
    maxlon: float,
    minlat: float,
    maxlat: float,
    periodic_lon: bool = False,
    verbose: bool = False,
    require_2d: bool = True,
) -> xarray.DataArray:
    if not np.isfinite(minlon) or not np.isfinite(maxlon):
        raise Exception(f"Longitude range {minlon} - {maxlon} is not valid")
    if not np.isfinite(minlat) or not np.isfinite(maxlat):
        raise Exception(f"Latitude range {minlat} - {maxlat} is not valid")
    if minlon > maxlon:
        raise Exception(
            f"Invalid longitude range: maximum {maxlon} must be >= minimum {minlon}."
        )
    if minlat > maxlat:
        raise Exception(
            f"Invalid latitude range: maximum {maxlat} must be >= minimum {minlat}."
        )
    source_lon, source_lat = source.getm.longitude, source.getm.latitude
    if source_lon.ndim != 1:
        raise Exception(f"Source longitude must be 1D but has shape {source_lon.shape}")
    if source_lat.ndim != 1:
        raise Exception(f"Source latitude must be 1D but has shape {source_lat.shape}")
    imin = source_lon.values.searchsorted(minlon, side="right") - 1
    imax = source_lon.values.searchsorted(maxlon, side="left") + 1
    if source_lat.values[1] < source_lat.values[0]:
        jmin = (
            source_lat.size
            - source_lat.values[::-1].searchsorted(maxlat, side="left")
            - 1
        )
        jmax = (
            source_lat.size
            - source_lat.values[::-1].searchsorted(minlat, side="right")
            + 1
        )
    else:
        jmin = source_lat.values.searchsorted(minlat, side="right") - 1
        jmax = source_lat.values.searchsorted(maxlat, side="left") + 1
    if verbose:
        print(imin, imax, source_lon.values.size, jmin, jmax, source_lat.values.size)
    assert (imin >= 0 and imax <= source_lon.values.size) or periodic_lon, (
        f"Requested longitude section {minlon} - {maxlon} is not fully covered"
        f" by available range {source_lon.values[0]} - {source_lon.values[-1]}"
    )
    assert jmin >= 0 and jmax <= source_lat.values.size, (
        f"Requested latitude section {minlat} - {maxlat} is not fully covered"
        f" by available range {source_lat.values[0]} - {source_lat.values[-1]}"
    )
    if require_2d and jmax - jmin == 1:
        jmin, jmax = (jmin, jmax + 1) if jmin == 0 else (jmin - 1, jmax)
    if require_2d and imax - imin == 1:
        imin, imax = (imin, imax + 1) if imin == 0 else (imin - 1, imax)
    add_left = imin < 0
    add_right = imax >= source_lon.values.size
    imin = max(imin, 0)
    imax = min(imax, source_lon.values.size)
    ilondim = source.dims.index(source_lon.dims[0])
    ilatdim = source.dims.index(source_lat.dims[0])
    shape = list(source.shape)
    shape[ilondim] = imax - imin
    shape[ilatdim] = jmax - jmin
    center_source = tuple(
        [
            {ilondim: slice(imin, imax), ilatdim: slice(jmin, jmax)}.get(i, slice(None))
            for i in range(len(shape))
        ]
    )
    center_target = [
        {ilondim: slice(0, imax - imin), ilatdim: slice(0, jmax - jmin)}.get(
            i, slice(None)
        )
        for i in range(len(shape))
    ]
    target_lon = source_lon[center_source[ilondim]]
    target_lat = source_lat[center_source[ilatdim]]
    overlap = abs(source_lon.values[-1] - source_lon.values[0] - 360.0) < 1e-5
    if verbose:
        print(
            f"periodic longitude? {periodic_lon} Overlap?"
            f" {abs(source_lon.values[-1] - source_lon.values[0] - 360.0)} = {overlap}"
        )
    left_target = None
    right_target = None
    if add_left:
        # Periodic domain and we need to read beyond left boundary
        imin_left = source_lon.values.searchsorted(minlon + 360.0, side="right") - 1
        left_source = tuple(
            [
                {ilondim: slice(imin_left, -1 if overlap else None)}.get(i, s)
                for i, s in enumerate(center_source)
            ]
        )
        nleft = source_lon.values.size - imin_left + (-1 if overlap else 0)
        if verbose:
            print(f"adding {nleft} values on the left")
        shape[ilondim] += nleft
        left_target = tuple(
            [{ilondim: slice(0, nleft)}.get(i, s) for i, s in enumerate(center_target)]
        )
        center_target[ilondim] = slice(nleft, nleft + imax - imin)
        target_lon = xarray.concat(
            (source_lon[left_source[ilondim]] - 360.0, target_lon),
            source_lon.dims[0],
            combine_attrs="no_conflicts",
        )
    if add_right:
        # Periodic domain and we need to read beyond right boundary
        imax_right = source_lon.values.searchsorted(maxlon - 360.0, side="left") + 1
        right_source = tuple(
            [
                {ilondim: slice(1 if overlap else 0, imax_right)}.get(i, s)
                for i, s in enumerate(center_source)
            ]
        )
        nright = imax_right + (-1 if overlap else 0)
        if verbose:
            print(f"adding {nright} values on the right")
        shape[ilondim] += nright
        right_target = tuple(
            [
                {ilondim: slice(s.stop, None)}.get(i, s)
                for i, s in enumerate(center_target)
            ]
        )
        target_lon = xarray.concat(
            (target_lon, source_lon[right_source[ilondim]] + 360.0),
            source_lon.dims[0],
            combine_attrs="no_conflicts",
        )
    center_target = tuple(center_target)
    shape = tuple(shape)
    if verbose:
        print(f"final shape: {shape}")

    lazyvar = Slice(
        _as_lazyarray(source),
        shape=shape,
        passthrough=[i for i in range(len(shape)) if i not in (ilondim, ilatdim)],
    )
    lazyvar._slices.append((center_source, center_target))
    if left_target:
        lazyvar._slices.append((left_source, left_target))
    if right_target:
        lazyvar._slices.append((right_source, right_target))

    coords = dict(source.coords.items())
    coords[source_lon.name] = target_lon
    coords[source_lat.name] = target_lat
    return xarray.DataArray(
        lazyvar, dims=source.dims, coords=coords, attrs=source.attrs, name=lazyvar.name
    )


def concatenate_slices(
    source: xarray.DataArray, idim: int, slices: Tuple[slice, ...], verbose=False
) -> xarray.DataArray:
    assert idim < source.ndim
    assert all([isinstance(s, slice) for s in slices])
    shape = list(source.shape)
    shape[idim] = sum([s.stop - s.start for s in slices])
    shape = tuple(shape)
    if verbose:
        print(f"final shape: {shape}")

    istart = 0
    strslices = ""
    final_slices = []
    for s in slices:
        n = s.stop - s.start
        source_slice = [slice(None)] * source.ndim
        target_slice = [slice(None)] * source.ndim
        source_slice[idim] = s
        target_slice[idim] = slice(istart, istart + n)
        strslices += f"[{s.start}:{s.stop}],"
        final_slices.append((tuple(source_slice), tuple(target_slice)))
        istart += n
    assert istart == shape[idim]

    lazyvar = Slice(
        _as_lazyarray(source),
        shape=shape,
        passthrough=[i for i in range(len(shape)) if i != idim],
        dtype=source.dtype,
    )
    lazyvar._slices.extend(final_slices)

    coords = {}
    for name, c in source.coords.items():
        if source.dims[idim] not in c.dims:
            coords[name] = c
    return xarray.DataArray(
        lazyvar, dims=source.dims, coords=coords, attrs=source.attrs, name=lazyvar.name
    )


class Transpose(UnaryOperator):
    def __init__(self, a, axes: Iterable[int], **kwargs):
        super().__init__(a, **kwargs)
        self.axes = axes
        self.oldaxes = list(axes)
        for inew, iold in enumerate(self.axes):
            self.oldaxes[iold] = inew

    def __array__(self, dtype=None) -> np.ndarray:
        return np.asarray(self._source).transpose(self.axes)

    def __getitem__(self, slices) -> np.ndarray:
        newslices = list(self._finalize_slices(slices))
        finalslices = [slice(None)] * len(newslices)
        for inew, s in enumerate(slices):
            if isinstance(s, (int, np.integer)):
                newslices[inew] = slice(s, s + 1)
                finalslices[inew] = 0
        oldslices = tuple(newslices[inew] for inew in self.oldaxes)
        return self._source[oldslices].transpose(self.axes)[tuple(finalslices)]


def transpose(
    source: xarray.DataArray, axes: Optional[Iterable[int]] = None
) -> xarray.DataArray:
    if axes is None:
        axes = range(source.ndim)[::-1]
    dims = [source.dims[i] for i in axes]
    shape = [source.shape[i] for i in axes]
    lazyvar = Transpose(
        _as_lazyarray(source),
        axes,
        shape=shape,
        passthrough=list(range(source.ndim)),
        dtype=source.dtype,
    )
    coords = {}
    for name, c in source.coords.items():
        if c.ndim > 1:
            newcdims = []
            for d in dims:
                if d in c.dims:
                    newcdims.append(d)
            caxes = [c.dims.index(d) for d in newcdims]
            coords[name] = transpose(c, caxes)
        coords[name] = c
    return xarray.DataArray(
        lazyvar,
        dims=dims,
        coords=coords,
        attrs=source.attrs,
        name=lazyvar.name,
    )


def isel(source: xarray.DataArray, **indices) -> xarray.DataArray:
    """Index named dimensions with integers, slice objects or integer arrays"""
    advanced_indices = []
    for dim in list(indices):
        assert dim in source.dims, (
            f"indexed dimension {dim} not used by source,"
            f" which has dimensions {source.dims}"
        )
        if not isinstance(indices[dim], (int, slice)):
            advanced_indices.append(source.dims.index(dim))
            #            indices[dim] = xarray.Variable([source.dims[advanced_indices[0]]], np.asarray(indices[dim], dtype=np.intp))
            indices[dim] = xarray.Variable(
                ["__newdim"], np.asarray(indices[dim], dtype=np.intp)
            )

    # Final slices per dimension
    slices = tuple([indices.get(dim, slice(None)) for dim in source.dims])

    # Determine final shape
    shape = []
    dims = []
    passthrough = {}
    advanced_added = False
    for i, (dim, slc, l) in enumerate(zip(source.dims, slices, source.shape)):
        if i not in advanced_indices:
            # Slice is integer or slice object. if integer, it will be sliced out
            # so it does not contribute to the final shape
            if isinstance(slc, slice):
                start, stop, stride = slc.indices(l)
                dims.append(dim)
                passthrough[len(shape)] = i
                shape.append((stop - start + stride - 1) // stride)
        elif not advanced_added:
            # First advanced slice. Add the shape produced by the broadcast combination
            # of advanced indices
            assert max(advanced_indices) - min(advanced_indices) + 1 == len(
                advanced_indices
            ), "advanced indices must be side-by-side for now"
            advanced_shapes = [indices[source.dims[i]].shape for i in advanced_indices]
            for length in np.broadcast_shapes(*advanced_shapes):
                dims.append(f"dim_{len(shape)}")
                shape.append(length)
            advanced_added = True

    lazyvar = Slice(
        _as_lazyarray(source), shape=shape, passthrough=passthrough, dtype=source.dtype
    )
    lazyvar._slices.append((slices, (slice(None),) * len(shape)))

    coords = {}
    for name, c in source.coords.items():
        if all(dim not in c.dims for dim in indices):
            coords[name] = c
    return xarray.DataArray(
        lazyvar, dims=dims, coords=coords, attrs=source.attrs, name=lazyvar.name
    )


def horizontal_interpolation(
    source: xarray.DataArray,
    lon: xarray.DataArray,
    lat: xarray.DataArray,
    dtype: numpy.typing.DTypeLike = float,
    mask=None,
) -> xarray.DataArray:
    source_lon, source_lat = source.getm.longitude, source.getm.latitude
    if source_lon is None:
        raise Exception(
            f"Variable {source.name} does not have a valid longitude coordinate."
        )
    if source_lat is None:
        raise Exception(
            f"Variable {source.name} does not have a valid latitude coordinate."
        )
    assert source_lon.ndim == 1
    assert source_lat.ndim == 1
    assert np.isfinite(lon).all(), f"Some longitudes are non-finite: {lon}"
    assert np.isfinite(lat).all(), f"Some latitudes are non-finite: {lat}"
    lon, lat = np.broadcast_arrays(lon, lat)
    ilondim = source.dims.index(source_lon.dims[0])
    ilatdim = source.dims.index(source_lat.dims[0])
    assert (
        abs(ilondim - ilatdim) == 1
    ), "Longitude and latitude dimensions must be distinct and adjacent"
    dimensions = {
        0: (),
        1: (source_lon.dims[0],),
        2: (source_lat.dims[0], source_lon.dims[-1]),
    }[lon.ndim]
    shape = (
        source.shape[: min(ilondim, ilatdim)]
        + lon.shape
        + source.shape[max(ilondim, ilatdim) + 1 :]
    )
    kwargs = {"ndim_trailing": source.ndim - max(ilondim, ilatdim) - 1, "mask": mask}
    if ilondim > ilatdim:
        # Dimension order: latitude first, then longitude
        ip = pygetm.util.interpolate.Linear2DGridInterpolator(
            lat, lon, source_lat, source_lon, **kwargs
        )
    else:
        # Dimension order: longitude first, then latitude
        dimensions = dimensions[::-1]
        ip = pygetm.util.interpolate.Linear2DGridInterpolator(
            lon, lat, source_lon, source_lat, **kwargs
        )
    lon_name, lat_name = source_lon.name, source_lat.name
    if lon_name in dimensions and lon.ndim > 1:
        lon_name = lon_name + "_"
    if lat_name in dimensions and lat.ndim > 1:
        lat_name = lat_name + "_"
    lon = xarray.DataArray(lon, dims=dimensions, name=lon_name, attrs=source_lon.attrs)
    lat = xarray.DataArray(lat, dims=dimensions, name=lat_name, attrs=source_lat.attrs)
    coords = dict(
        [
            (k, v)
            for k, v in source.coords.items()
            if k not in {source_lon.name, source_lat.name}
        ]
    )
    coords[lon.name] = lon
    coords[lat.name] = lat
    dims = (
        source.dims[: min(ilondim, ilatdim)]
        + dimensions
        + source.dims[max(ilondim, ilatdim) + 1 :]
    )
    lazyvar = HorizontalInterpolation(
        ip,
        _as_lazyarray(source),
        shape,
        min(ilondim, ilatdim),
        source.ndim - max(ilondim, ilatdim) - 1,
    )
    return xarray.DataArray(
        lazyvar, dims=dims, coords=coords, attrs=source.attrs, name=lazyvar.name
    )


class HorizontalInterpolation(UnaryOperator):
    def __init__(
        self,
        ip: pygetm.util.interpolate.Linear2DGridInterpolator,
        source: LazyArray,
        shape: Iterable[int],
        npre: int,
        npost: int,
        **kwargs,
    ):
        UnaryOperator.__init__(self, source, shape=shape, **kwargs)
        self._ip = ip
        self.npre = npre
        self.npost = npost

    def __array__(self, dtype=None) -> np.ndarray:
        return self._ip(np.asarray(self._source))

    def __getitem__(self, slices) -> np.ndarray:
        src_slice, tgt_slice = [Ellipsis], [Ellipsis]
        ntrailing_dim_removed = 0
        for i, s in enumerate(self._finalize_slices(slices)):
            if i < self.npre:
                # prefixed dimension
                src_slice.insert(i, s)
            elif i >= self.ndim - self.npost:
                # trailing dimension
                if isinstance(s, (int, np.integer)):
                    ntrailing_dim_removed += 1
                    tgt_slice.append(0)
                src_slice.append(s)
            else:
                assert (
                    isinstance(s, slice)
                    and s.start is None
                    and s.stop is None
                    and s.step is None
                ), repr(s)
        source = np.asarray(self._source[tuple(src_slice)])
        source.shape = source.shape + (1,) * ntrailing_dim_removed
        result = self._ip(source)
        return result[tuple(tgt_slice)]


def vertical_interpolation(
    source: xarray.DataArray, target_z: numpy.typing.ArrayLike, itargetdim: int = 0
) -> xarray.DataArray:
    source_z = source.getm.z
    target_z = np.asarray(target_z)
    if source_z is None:
        raise Exception(
            f"Variable {source.name} does not have a valid depth coordinate."
        )
    assert source_z.ndim == 1
    izdim = source.dims.index(source_z.dims[0])
    # assert source.ndim - izdim == target_z.ndim
    # assert source.shape[izdim + 1:izdim + 3] == target_z.shape[1:], f'{source.shape[izdim + 1:izdim + 3]} vs {target_z.shape[1:]}'
    target2sourcedim = {}
    isourcedim = 0
    for i, l in enumerate(target_z.shape):
        if i == itargetdim:
            isourcedim = izdim
        else:
            while (
                isourcedim != izdim
                and isourcedim < source.ndim
                and l != source.shape[isourcedim]
            ):
                isourcedim += 1
            assert isourcedim != izdim, (
                f"Dimension with length {l} should precede depth dimension {izdim}"
                f" in {source.name}, which has shape {source.shape}"
            )
            assert isourcedim < source.ndim, (
                f"Dimension with length {l} expected after depth dimension {izdim}"
                f" in {source.name}, which has shape {source.shape}"
            )
        target2sourcedim[i] = isourcedim
        isourcedim += 1
    coords = {}
    for n, c in source.coords.items():
        if n == source.dims[izdim]:
            coords[n + "_"] = (
                [source.dims[target2sourcedim[i]] for i in range(target_z.ndim)],
                target_z,
            )
        else:
            coords[n] = c
    lazyvar = VerticalInterpolation(
        _as_lazyarray(source), target_z, izdim, source_z.values, itargetdim
    )
    return xarray.DataArray(
        lazyvar, dims=source.dims, coords=coords, attrs=source.attrs, name=lazyvar.name
    )


class VerticalInterpolation(UnaryOperator):
    def __init__(
        self,
        source: LazyArray,
        z: np.ndarray,
        izdim: int,
        source_z: np.ndarray,
        axis: int = 0,
        **kwargs,
    ):
        self.izdim = izdim
        passthrough = [idim for idim in range(source.ndim) if idim != self.izdim]
        shape = list(source.shape)
        shape[self.izdim] = z.shape[axis]
        super().__init__(source, shape=shape, passthrough=passthrough, **kwargs)
        self.z = z
        self.axis = axis
        self.source_z = source_z
        if (self.source_z >= 0.0).all():
            self.source_z = -self.source_z

    def apply(self, source, dtype=None) -> np.ndarray:
        return pygetm.util.interpolate.interp_1d(
            self.z, self.source_z, source, axis=self.axis
        )


def temporal_interpolation(
    source: xarray.DataArray, climatology: bool = False
) -> xarray.DataArray:
    time_coord = source.getm.time
    assert time_coord is not None, "No time coordinate found"
    itimedim = source.dims.index(time_coord.dims[0])
    lazyvar = TemporalInterpolation(
        _as_lazyarray(source), itimedim, time_coord.values, climatology
    )
    dims = [d for i, d in enumerate(source.dims) if i != lazyvar._itimedim]
    coords = dict(source.coords.items())
    coords[time_coord.dims[0]] = lazyvar._timecoord
    return xarray.DataArray(
        lazyvar, dims=dims, coords=coords, attrs=source.attrs, name=lazyvar.name
    )


def _as_lazyarray(array: xarray.DataArray) -> LazyArray:
    variable = array.variable
    if isinstance(variable._data, LazyArray):
        return variable._data
    else:
        name = array.name
        if name is None:
            name = object.__repr__(variable._data)
        if "source" in array.encoding:
            name += " from " + array.encoding["source"]
        return Wrap(variable, name=name)


class TemporalInterpolation(UnaryOperator):
    __slots__ = (
        "_current",
        "_itimedim",
        "_numnow",
        "_numnext",
        "_slope",
        "_inext",
        "_next",
        "_slices",
        "climatology",
        "_year",
        "_timevalues",
    )

    def __init__(
        self,
        source: LazyArray,
        itimedim: int,
        times: np.ndarray,
        climatology: bool,
        dtype=float,
        **kwargs,
    ):
        shape = list(source.shape)
        self._itimedim = itimedim
        ntime = shape[self._itimedim]
        if ntime <= 1:
            raise Exception(
                f"Cannot interpolate {source.name} in time because"
                f" its time dimension has length {ntime}."
            )
        shape.pop(self._itimedim)

        super().__init__(source, shape=shape, dtype=dtype, **kwargs)

        self._current = np.empty(shape, dtype=self.dtype)

        self.times = times
        self._timecoord = xarray.DataArray(self.times[0])
        self._timevalues = self._timecoord.values

        self._numnow = None
        self._numnext = 0.0
        self._slope = 0.0
        self._inext = -1
        self._next = 0.0
        self._slices: List[Union[int, slice]] = [slice(None)] * source.ndim

        self.climatology = climatology
        self._year = self.times[0].year
        if climatology and not all(time.year == self._year for time in self.times):
            raise Exception(
                f"{self._source_name} cannot be used as climatology because"
                " it spans more than one calendar year"
            )

    def __array__(self, dtype=None) -> np.ndarray:
        return self._current

    def __getitem__(self, slices) -> np.ndarray:
        return self._current[slices]

    def is_time_varying(self) -> bool:
        return True

    def update(
        self, time: cftime.datetime, numtime: Optional[np.longdouble] = None
    ) -> bool:
        if numtime is None:
            numtime = time.toordinal(fractional=True)

        if self._numnow is None:
            # First call to update
            self._start(time)
        elif numtime <= self._numnow:
            # Subsequent call to update, but time has not increased
            # If equal to previous time, we are done. If smaller, rewind
            if numtime == self._numnow:
                return False
            self._start(time)

        while self._numnext < numtime:
            self._move_to_next(time)

        # Do linear interpolation
        np.multiply(self._slope, numtime - self._numnext, out=self._current)
        self._current += self._next

        # Save current time
        self._numnow = numtime
        self._timevalues[...] = time
        return True

    def _start(self, time: cftime.datetime):
        if time.calendar != self.times[0].calendar:
            raise Exception(
                f"Simulation calendar {time.calendar} does not match calendar"
                f" {self.times[0].calendar} used by {self._source_name}."
            )

        # Find the time index (_inext) just before the left bound of the window we need.
        # We will subsequently call _move_to_next twice to load the actual window.
        # If a time in the series exactly matches the requested time, this will be the
        # right bound of the time window (hence side="left" below), except if that
        # is the very first time in the time series. Then it will be the left bound.
        if self.climatology:
            clim_time = time.replace(year=self._year)
            self._inext = self.times.searchsorted(clim_time, side="left") - 2
            self._year = time.year
            if self._inext < -1:
                self._inext += self.times.size
                self._year -= 1
        else:
            # Make sure the time series does not start after the requested time.
            if time < self.times[0]:
                raise Exception(
                    f"Cannot interpolate {self._source_name} to value at {time},"
                    f" because time series starts only at {self.times[0]}."
                )
            self._inext = max(-1, self.times.searchsorted(time, side="left") - 2)

        # Load left and right bound of the window encompassing the current time.
        # After that, all information for linear interpolation (_next, _slope)
        # is available.
        self._move_to_next(time)
        self._move_to_next(time)

    def _move_to_next(self, time: cftime.datetime):
        # Move to next record
        self._inext += 1
        if self._inext == self.times.size:
            if self.climatology:
                self._inext = 0
                self._year += 1
            else:
                raise Exception(
                    f"Cannot interpolate {self._source_name} to value at {time}"
                    f" because end of time series was reached ({self.times[-1]})."
                )
        old, numold = self._next, self._numnext
        self._slices[self._itimedim] = self._inext
        self._next = np.asarray(self._source[tuple(self._slices)], dtype=self.dtype)
        next_time = self.times[self._inext]
        if self.climatology:
            next_time = next_time.replace(year=self._year)
        self._numnext = next_time.toordinal(fractional=True)
        self._slope = (self._next - old) / (self._numnext - numold)


def slicespec2string(s: Union[tuple, slice, int]) -> str:
    if isinstance(s, slice):
        start = "" if s.start is None else f"{s.start}"
        stop = "" if s.stop is None else f"{s.stop}"
        if s.step in (None, 1):
            return f"{start}:{stop}"
        return f"{start}:{stop}:{s.step}"
    elif isinstance(s, tuple):
        return ",".join([slicespec2string(item) for item in s])
    return f"{s!r}"


def debug_nc_reads(logger: Optional[logging.Logger] = None):
    """Hook into :mod:`xarray` so that every read from a NetCDF file is
    written to the log.
    """
    import xarray.backends.netCDF4_

    if logger is None:
        logger = logging.getLogger("pygetm.input")
        logger.setLevel(logging.DEBUG)

    class NetCDF4ArrayWrapper2(xarray.backends.netCDF4_.NetCDF4ArrayWrapper):
        __slots__ = ()

        def _getitem(self, key):
            logger.debug(
                f"Reading {self.variable_name}[{slicespec2string(key)}]"
                f" from {self.datastore._filename}"
            )
            return super()._getitem(key)

    xarray.backends.netCDF4_.NetCDF4ArrayWrapper = NetCDF4ArrayWrapper2


class OnGrid(enum.Enum):
    #: Grids do not match. Spatially explicit data will require horizontal
    #: and - if vertically resolved - vertical interpolation.
    NONE = enum.auto()

    #: Horizontal grid matches, but vertical does not.
    #: Vertically resolved data will require vertical interpolation.
    HORIZONTAL = enum.auto()

    #: Horizontal and vertical grids match
    ALL = enum.auto()


class InputManager:
    def __init__(self):
        self._all_fields: List[Tuple[str, LazyArray, np.ndarray]] = []
        self._micro_fields: List[Tuple[str, LazyArray, np.ndarray]] = []
        self._logger = logging.getLogger()

    def debug_nc_reads(self):
        """Hook into :mod:`xarray` so that every read from a NetCDF file is written
        to the log.
        """
        _logger = self._logger.getChild("nc")
        _logger.setLevel(logging.DEBUG)
        debug_nc_reads(_logger)

    def add(
        self,
        array: pygetm.core.Array,
        value: Union[numbers.Number, np.ndarray, xarray.DataArray, LazyArray],
        periodic_lon: bool = True,
        on_grid: Union[bool, OnGrid] = False,
        include_halos: Optional[bool] = None,
        climatology: bool = False,
        mask: bool = False,
    ):
        """Link an array to the provided input. If this input is constant in time,
        the value of the array will be set immediately.

        Args:
            array: array to assign a value to
            value: input to assign. If this is time-dependent, the combination of the
                array and its linked input will be registered; the array will then be
                updated to the current time whenever :meth:`update` is called.
            periodic_lon: whether this input covers all longitudes (i.e., the entire
                globe in the horizontal) and therefore has a periodic boundary. This
                enables efficient spatial interpolation across longitude bounds of the
                input, for instance, accessing read 10 degrees West to 5 degrees East
                for an input that spans 0 to 360 degrees East.
            on_grid: whether the input is defined on the same grid (horizontal-only,
                or both horizontal and vertical) as the array that is being assigned to.
                If this is ``False``, the value will be spatially interpolated
                to the array grid. ``True`` is equivalent to :attr:`OnGrid.HORIZONTAL`.
            include_halos: whether to also update the halos of the array. If not
                provided, this default to ``True`` if the array has attributes
                ``_require_halos`` or ``_part_of_state``; otherwise it defaults to
                ``False``.
            climatology: whether the input describes a single climatological year
                (at any temporal resolution, e.g., monthly, daily) that is
                representative for any true year. This argument is relevant only if the
                provided input is time-varying. It also requires that the input does
                not span more than one year.
            mask: whether to set the array to its :attr:`pygetm.core.Array.fill_value`
                in all masked points. If not provided, only missing values in the input
                (NaNs) will be set to the fill value. This currently only has an effect
                when the input is non time-varying.
        """
        if array.all_values is None or array.all_values.size == 0:
            # The target variable does not contain data. Typically this is because
            # it specifies information on the open boundaries,
            # of which the current (sub)domain does not have any.
            self._logger.warning(
                f"Ignoring asssignment to array {array.name}"
                " because it has no associated data."
            )
            return

        if isinstance(value, (numbers.Number, np.ndarray)):
            # Constant-in-time fill value. Set it, then forget about the array
            # as it will not require further updating.
            array.fill(value)
            return

        assert isinstance(value, xarray.DataArray), (
            "If value is not numeric, it should be an xarray.DataArray,"
            f" but it is {value!r}."
        )

        if include_halos is None:
            include_halos = array.attrs.get("_require_halos", False) or array.attrs.get(
                "_part_of_state", False
            )
        if not isinstance(on_grid, OnGrid):
            on_grid = OnGrid.HORIZONTAL if on_grid else OnGrid.NONE

        grid = array.grid

        # Obtain active area of local subdomain (including halos if
        # include_halos is True) and the corresponding slice in the global domain
        # (always excluding halos)
        (
            local_slice,
            global_slice,
            local_shape,
            global_shape,
        ) = grid.domain.tiling.subdomain2slices(
            exclude_halos=not include_halos, halo_sub=2
        )

        target_slice = (Ellipsis,)
        source_lon, source_lat = value.getm.longitude, value.getm.latitude
        if array.on_boundary:
            # Open boundary information. This can either be specified for the global
            # domain (e.g., when read from NetCDF), or for only the open boundary
            # points that fall within the local subdomain. Determine which of these.
            if value.ndim >= 2 and value.shape[-2:] == global_shape:
                # on-grid data for the global domain:
                # extract data at open boundary points
                i_bnd = grid.domain.open_boundaries.i_glob
                j_bnd = grid.domain.open_boundaries.j_glob
                value = isel(value, **{value.dims[-1]: i_bnd, value.dims[-2]: j_bnd})
            elif (
                source_lon is not None
                and source_lat is not None
                and source_lon.ndim > 0
                and source_lat.ndim > 0
            ):
                # Spatially explicit input:
                # interpolate horizontally to open boundary coordinates
                if source_lon.ndim != 1:
                    raise Exception(
                        f"Unsuitable shape {source_lon.shape} of longitude coordinate"
                        f" {source_lon.name}. Off-grid boundary information can be used"
                        f" only if its longitude is 1D."
                    )
                if source_lat.ndim != 1:
                    raise Exception(
                        f"Unsuitable shape {source_lat.shape} of latitude coordinate"
                        f" {source_lat.name}. Off-grid boundary information can be used"
                        f" only if its latitude is 1D."
                    )
                ilondim = value.dims.index(source_lon.dims[0])
                ilatdim = value.dims.index(source_lat.dims[0])
                if ilondim != ilatdim:
                    lon_bnd = grid.domain.open_boundaries.lon.all_values
                    lat_bnd = grid.domain.open_boundaries.lat.all_values
                    value = limit_region(
                        value,
                        lon_bnd.min(),
                        lon_bnd.max(),
                        lat_bnd.min(),
                        lat_bnd.max(),
                        periodic_lon=periodic_lon,
                    )
                    value = pygetm.input.horizontal_interpolation(
                        value, lon_bnd, lat_bnd
                    )

            if array.z and value.getm.z is not None and value.getm.z.ndim == 1:
                # Source and target arrays are depth-explicit and the source depth is 1D
                # Ensure it is the last (fastest varying) dimension
                izdim = value.dims.index(value.getm.z.dims[0])
                if izdim != value.ndim - 1:
                    axes = list(range(value.ndim))
                    axes.append(axes.pop(izdim))
                    value = transpose(value, axes)

            idim = value.ndim - (2 if array.z else 1)
            if value.shape[idim] == grid.domain.open_boundaries.np_glob:
                # The source array covers all open boundaries (global domain).
                # If the subdomain only has a subset of those, slice out only the points
                # that fall within the current subdomain
                local_to_global = grid.domain.open_boundaries.local_to_global
                if local_to_global:
                    value = concatenate_slices(
                        value,
                        idim,
                        [slice(start, stop) for (start, stop) in local_to_global],
                    )
            elif value.shape[idim] != grid.domain.open_boundaries.np:
                raise Exception(
                    f"Extent of dimension {idim} of {value.name} is not compatible with"
                    f" open boundaries. It should have length"
                    f" {grid.domain.open_boundaries.np_glob} (number of open boundary"
                    f" points in the global domain) or {grid.domain.open_boundaries.np}"
                    f" (number of open boundary points in the current subdomain)."
                    f" Its actual extent is {value.shape[idim]}."
                )
        elif array.ndim != 0:
            # The target is a normal 2D (horizontal-only) or 3D (depth-explicit) array
            # The source data can either be on the native model grid, or at an
            # arbitrary lon, lat grid. In the latter case, we interpolate in space.
            assert array.all_values.shape[-2:] == local_shape
            target_slice = local_slice
            if source_lon is None and source_lat is None and on_grid == OnGrid.NONE:
                # time series for single location
                value = value.expand_dims(("y", "x"), (value.ndim, value.ndim + 1))
            elif on_grid == OnGrid.NONE:
                # interpolate horizontally to local array INCLUDING halos
                lon = grid.lon.all_values[target_slice]
                lat = grid.lat.all_values[target_slice]
                assert not np.isnan(lon).any()
                assert not np.isnan(lat).any()
                value = limit_region(
                    value,
                    lon.min(),
                    lon.max(),
                    lat.min(),
                    lat.max(),
                    periodic_lon=periodic_lon,
                )
                value = horizontal_interpolation(value, lon, lat)
            else:
                # the input is already on-grid, but we need to map
                # from global domain to subdomain
                assert value.shape[-2:] == global_shape, (
                    f"{array.name}: shape of values {value.shape[-2:]}"
                    f" should match that of global domain {global_shape}"
                )
                value = value[global_slice]

        if value.getm.time is not None:
            # The source data is time-dependent; during the simulation it will be
            # interpolated in time.
            if value.getm.time.size > 1:
                value = temporal_interpolation(value, climatology=climatology)
            elif value.getm.time.dims:
                time = value.getm.time.values.flat[0]
                self._logger.warning(
                    f"{array.name} is set to {value.name}, which has only one time"
                    f" point {time}. The value from this time will be used now."
                    f" {array.name} will not be further updated by the input manager"
                    " at runtime."
                )
                itimedim = value.dims.index(value.getm.time.dims[0])
                slc = [slice(None)] * value.ndim
                slc[itimedim] = 0
                value = value[tuple(slc)]

        if array.z and on_grid != OnGrid.ALL:
            # The target is a depth-explicit array.
            # The source must be defined on z coordinates
            # and interpolated to our [time-varying] depths
            coord_source = grid.domain.open_boundaries if array.on_boundary else grid
            z_coordinate = coord_source.zc if array.z == CENTERS else coord_source.zf
            z_coordinate.saved = True
            value = vertical_interpolation(
                value,
                z_coordinate.all_values[target_slice],
                itargetdim=1 if array.on_boundary else 0,
            )

        target = array.all_values[target_slice]
        try:
            np.broadcast_shapes(value.shape, target.shape)
        except ValueError:
            assert (
                False
            ), f"Source shape {value.shape} does not match target shape {target.shape}"
        data = value.variable._data
        if isinstance(data, LazyArray) and data.is_time_varying():
            time_varying = array.attrs.get("_time_varying", TimeVarying.MICRO)
            suffix = " on macrotimestep" if time_varying == TimeVarying.MACRO else ""
            self._logger.info(
                f"{array.name} will be updated dynamically from {data.name}{suffix}"
            )
            info = (array.name, data, target)
            self._all_fields.append(info)
            if time_varying == TimeVarying.MICRO:
                self._micro_fields.append(info)
        else:
            target[...] = value
            finite = np.isfinite(target)
            if array.ndim == 0 or array.on_boundary:
                unmasked = np.broadcast_to(True, target.shape)
            else:
                unmasked = np.broadcast_to(grid._water[target_slice], target.shape)
                if array.fill_value is not None:
                    keep_mask = unmasked if mask else unmasked | finite
                    target[~keep_mask] = array.fill_value
            if not finite.all(where=unmasked):
                n_unmasked = unmasked.sum()
                n_bad = n_unmasked - finite.sum(where=unmasked)
                self._logger.warning(
                    f"{array.name} is set to {value.name}, which is not finite"
                    f" (e.g., NaN) in {n_bad} of {n_unmasked} unmasked points."
                )
            minval = target.min(where=unmasked, initial=np.inf)
            maxval = target.max(where=unmasked, initial=-np.inf)
            self._logger.info(
                f"{array.name} is set to time-invariant {value.name}"
                f" (minimum: {minval}, maximum: {maxval})"
            )

    def update(self, time: cftime.datetime, macro: bool = True):
        """Update all arrays linked to time-dependent inputs to the current time.

        Args:
            time: current time
            macro: whether to also update arrays that were marked as only relevant for
                the macro (3D) time step
        """
        numtime = time.toordinal(fractional=True)
        fields = self._all_fields if macro else self._micro_fields
        for name, source, target in fields:
            self._logger.debug(f"updating {name}")
            source.update(time, numtime)
            target[...] = source

    @property
    def logger(self) -> logging.Logger:
        return self._logger

    def set_logger(self, logger: logging.Logger):
        self._logger = logger
