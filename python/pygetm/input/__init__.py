from typing import Callable, Iterable, List, Mapping, Union, Optional, Sequence, Tuple, TYPE_CHECKING
import glob
import numbers
import logging
import enum

import numpy
import numpy.typing
import numpy.lib.mixins
import xarray
import cftime

import pygetm.util.interpolate
from pygetm.constants import CENTERS

if TYPE_CHECKING:
    import pygetm.core


@xarray.register_dataarray_accessor('getm')
class GETMAccessor:
    def __init__(self, xarray_obj: xarray.DataArray):
        self._obj = xarray_obj

        self._coordinates: Optional[Mapping[str, xarray.DataArray]] = None

    @property
    def longitude(self) -> Optional[xarray.DataArray]:
        return self._interpret_coordinates().get('longitude')

    @property
    def latitude(self) -> Optional[xarray.DataArray]:
        return self._interpret_coordinates().get('latitude')

    @property
    def z(self) -> Optional[xarray.DataArray]:
        return self._interpret_coordinates().get('z')

    @property
    def time(self) -> Optional[xarray.DataArray]:
        return self._interpret_coordinates().get('time')

    def _interpret_coordinates(self) -> Mapping[str, xarray.DataArray]:
        if self._coordinates is None:
            self._coordinates = {}
            for name, coord in self._obj.coords.items():
                units = coord.attrs.get('units')
                standard_name = coord.attrs.get('standard_name')
                if units in ('degrees_north', 'degree_north', 'degree_N', 'degrees_N', 'degreeN', 'degreesN') or standard_name == 'latitude':
                    self._coordinates['latitude'] = coord
                elif units in ('degrees_east', 'degree_east', 'degree_E', 'degrees_E', 'degreeE', 'degreesE') or standard_name == 'longitude':
                    self._coordinates['longitude'] = coord
                elif coord.size > 0 and isinstance(coord.values.flat[0], cftime.datetime):
                    self._coordinates['time'] = coord
                elif name == 'zax':
                    self._coordinates['z'] = coord
        return self._coordinates


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


def from_nc(paths: Union[str, Sequence[str]], name: str, preprocess: Optional[Callable[[xarray.Dataset], xarray.Dataset]]=None, **kwargs) -> xarray.DataArray:
    """Obtain a variable from one or more NetCDF files that can be used as value provided to
    :meth:`InputManager.add` and :meth:`pygetm.core.Array.set`.

    Args:
        paths: single file path, a pathname pattern containing `*` and/or `?`, or a sequence of file paths.
            If multiple paths are provided (or the pattern resolves to multiple valid path names), the files
            will be concatenated along their time dimension.
        preprocess: function that transforms the :class:`xarray.Dataset` opened for every path provided.
            This can be used to modify the datasets before concatenation in time is attempted,
            for instance, to cut off time indices that overlap between files.
        **kwargs: additional keyword arguments to be passed to :func:`xarray.open_dataset`
    """
    kwargs.setdefault('decode_times', True)
    kwargs['use_cftime'] = True
    kwargs['cache'] = False
    if isinstance(paths, str):
        pattern = paths
        paths = glob.glob(pattern)
        if not paths:
            raise Exception('No files found matching %s' % pattern)
    arrays = []
    for path in paths:
        ds = _open(path, preprocess, **kwargs)
        array = ds[name]
        # Note: we wrap the netCDF array ourselves, in order to support lazy operators (e.g., add, multiply)
        lazyvar = WrappedArray(array.variable, name='from_nc(%s, %s)' % (path, name))
        array = xarray.DataArray(lazyvar, dims=array.dims, coords=array.coords, attrs=array.attrs, name=lazyvar.name)
        arrays.append(array)
    if len(arrays) == 1:
        return arrays[0]
    else:
        assert all(array.getm.time is not None for array in arrays)
        return xarray.concat(sorted(arrays, key=lambda a: a.getm.time.values.flat[0]), dim=arrays[0].getm.time.dims[0])


class LazyArray(numpy.lib.mixins.NDArrayOperatorsMixin):
    def __init__(self, shape: Iterable[int], dtype: numpy.typing.DTypeLike, name: str):
        self.shape = tuple(shape)
        self.ndim = len(self.shape)
        self.dtype = dtype
        self.name = name

    def update(self, time: cftime.datetime, numtime: numpy.longdouble) -> bool:
        return False

    def astype(self, dtype, **kwargs) -> numpy.ndarray:
        if dtype == self.dtype:
            return self
        return self.__array__(dtype)

    def __array_function__(self, func, types, args, kwargs):
        if func == numpy.result_type:
            args = tuple(x.dtype if isinstance(x, LazyArray) else x for x in args)
            return numpy.result_type(*args)
        if func == numpy.concatenate:
            return ConcatenatedArray(*args, **kwargs)
        args = tuple(numpy.asarray(x) if isinstance(x, LazyArray) else x for x in args)
        kwargs = dict((k, numpy.asarray(v)) if isinstance(v, LazyArray) else (k, v) for (k, v) in kwargs.items())
        return func(*args, **kwargs)

    def __array_ufunc__(self, ufunc, method: str, *inputs, **kwargs):
        if method != '__call__':
            return NotImplemented

        if 'out' in kwargs:
            return NotImplemented

        for x in inputs:
            if not isinstance(x, (numpy.ndarray, numbers.Number, LazyArray, xarray.Variable)):
                return NotImplemented

        return UFuncResult(getattr(ufunc, method), *inputs, **kwargs)

    def __array__(self, dtype=None) -> numpy.ndarray:
        raise NotImplementedError

    def __getitem__(self, slices) -> numpy.ndarray:
        return self.__array__()[slices]

    def is_time_varying(self) -> bool:
        return False

    def _finalize_slices(self, slices: Tuple):
        assert isinstance(slices, tuple)
        for i, s in enumerate(slices):
            if s is Ellipsis:
                slices = slices[:i] + (slice(None),) * (self.ndim + 1 - len(slices)) + slices[i + 1:]
                break
        assert len(slices) == self.ndim
        return slices


class OperatorResult(LazyArray):
    def __init__(self, *inputs, passthrough=(), dtype=None, shape=None, name: Optional[str]=None, **kwargs):
        # Unpack unnamed arguments
        self.inputs = []
        self.lazy_inputs = []
        self.input_names = []
        for inp in inputs:
            assert isinstance(inp, (numpy.ndarray, numbers.Number, LazyArray, xarray.Variable)), 'Input has unknown type %s' % type(inp)

            # Unpack to LazyArray if possible
            if isinstance(inp, xarray.Variable) and isinstance(inp._data, LazyArray):
                inp = inp._data

            if isinstance(inp, xarray.Variable) and isinstance(inp._data, LazyArray):
                self.input_names.append(inp._data.name)
            elif isinstance(inp, (numpy.ndarray, numbers.Number)):
                self.input_names.append(str(inp))
            else:
                # Other datatype, typically xarray.Variable.
                # Do not call str/repr, as that will cause evaluation (e..g read from file) of the entire array
                self.input_names.append(str(type(inp)))

            # If this is a WrappedArray, unwrap (the wrapping was only for ufunc support)
            if isinstance(inp, WrappedArray):
                inp = inp.inputs[0]

            if isinstance(inp, LazyArray):
                self.lazy_inputs.append(inp)
            self.inputs.append(inp)

        # Store keyword arguments as-is (no unpacking)
        self.kwargs = kwargs

        # Infer shape from inputs if not provided
        if shape is None:
            shapes = []
            for input in inputs:
                if isinstance(input, (numpy.ndarray, LazyArray, xarray.DataArray, xarray.Variable)):
                    shapes.append(input.shape)
            shape = numpy.broadcast_shapes(*shapes)

        # Process dimensions for which we can passthrough slices to inputs
        # This can be True (= all dimensions), an iterable, or a dictionary mapping sliced dimensions
        # to input dimensions (if the current operator adds or removes dimensions)
        if passthrough is True:
            passthrough = range(len(shape))
        if not isinstance(passthrough, dict):
            passthrough = dict([(i, i) for i in passthrough])
        self.passthrough = passthrough
        assert all([isinstance(dim, int) for dim in self.passthrough]), 'Invalid passthrough: %s. All entries should be of type int' % (self.passthrough,)

        # Generate a name for the variable if not provided
        if name is None:
            name = '%s(%s%s)' % (self.__class__.__name__, ', '.join(self.input_names), ''.join(', %s=%s' % item for item in self.kwargs.items()))

        super().__init__(shape, dtype or float, name)

    def update(self, *args) -> bool:
        updated = False
        for input in self.lazy_inputs:
            updated = input.update(*args) or updated
        return updated

    def is_time_varying(self) -> bool:
        return self.lazy_inputs and any(input.is_time_varying() for input in self.lazy_inputs)

    def __getitem__(self, slices) -> numpy.ndarray:
        preslices, postslices = [], []
        for i, slc in enumerate(self._finalize_slices(slices)):
            if i in self.passthrough:
                preslices.append(slc)
                if not isinstance(slc, (int, numpy.integer)): postslices.append(slice(None))
            else:
                preslices.append(slice(None))
                postslices.append(slc)
        inputs = [inp if isinstance(inp, numbers.Number) else numpy.asarray(inp[tuple(preslices)]) for inp in self.inputs]
        return self.apply(*inputs)[tuple(postslices)]

    def apply(self, *inputs, dtype=None) -> numpy.ndarray:
        raise NotImplementedError

    def __array__(self, dtype=None) -> numpy.ndarray:
        return self.apply(*[numpy.asarray(inp) for inp in self.inputs], dtype=dtype)


class UnaryOperatorResult(OperatorResult):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._source = self.inputs[0]
        self._source_name = self.input_names[0]


class UFuncResult(OperatorResult):
    def __init__(self, ufunc, *inputs, **kwargs):
        super().__init__(*inputs, passthrough=True, **kwargs)
        self.ufunc = ufunc

    def apply(self, *inputs, dtype=None) -> numpy.ndarray:
        return self.ufunc(*inputs, **self.kwargs)


class WrappedArray(UnaryOperatorResult):
    def __init__(self, source: xarray.Variable, **kwargs):
        assert isinstance(source, xarray.Variable)
        super().__init__(source, passthrough=True, dtype=source.dtype, **kwargs)

    def apply(self, source, dtype=None) -> numpy.ndarray:
        return source

    def update(self, *args) -> bool:
        return False


class SliceArray(UnaryOperatorResult):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._slices = []

    def __array__(self, dtype=None) -> numpy.ndarray:
        data = numpy.empty(self.shape, dtype or self.dtype)
        for src_slice, tgt_slice in self._slices:
            data[tgt_slice] = self._source[src_slice]
        return data

    def __getitem__(self, slices) -> numpy.ndarray:
        slices = self._finalize_slices(slices)
        shape = []
        for i, (l, s) in enumerate(zip(self.shape, slices)):
            if i in self.passthrough and isinstance(s, (int, numpy.integer)):
                # This dimension will be sliced out
                continue
            assert isinstance(s, slice), 'Dimension %i has unsupported slice type %s with value %s. Passthrough: %s' % (i, type(s), s, list(self.passthrough))
            start, stop, step = s.indices(l)
            assert i in self.passthrough or (start == 0 and stop == l and step == 1), 'invalid slice for dimension %i with length %i: %i:%i:%i' % (i, l, start, stop, step)
            shape.append(len(range(start, stop, step)))
        data = numpy.empty(shape, self.dtype)
        for src_slice, tgt_slice in self._slices:
            src_slice = list(src_slice)
            for iout, iin in self.passthrough.items():
                src_slice[iin] = slices[iout]
            tgt_slice = tuple([ori for i, (cust, ori) in enumerate(zip(slices, tgt_slice)) if i not in self.passthrough or not isinstance(cust, (int, numpy.integer))])
            data[tgt_slice] = self._source[tuple(src_slice)]
        return data


class ConcatenatedArray(UnaryOperatorResult):
    def __init__(self, arrays, axis: int=0, *args, **kwargs):
        shape = list(arrays[0].shape)
        for array in arrays[1:]:
            shape[axis] += array.shape[axis]
            assert all(array.shape[i] == l for i, l in enumerate(shape) if i != axis)
        self.axis = axis
        super().__init__(*arrays, shape=shape, **kwargs)

    def __array__(self, dtype=None) -> numpy.ndarray:
        return numpy.concatenate(self.inputs, axis=self.axis, dtype=dtype)

    def __getitem__(self, slices) -> numpy.ndarray:
        slices = list(self._finalize_slices(slices))
        assert isinstance(slices[self.axis], (int, numpy.integer)), 'Unsupported slice for concatenated dimension: %s' % (slices[self.axis],)
        if slices[self.axis] < 0:
            slices[self.axis] += self.shape[self.axis]
        assert slices[self.axis] >= 0 and slices[self.axis] < self.shape[self.axis]
        for input in self.inputs:
            if slices[self.axis] < input.shape[self.axis]:
                return numpy.asarray(input[tuple(slices)])
            slices[self.axis] -= input.shape[self.axis]
        assert False, 'Index out of bounds?'


def limit_region(source: xarray.DataArray, minlon: float, maxlon: float, minlat: float, maxlat: float, periodic_lon: bool=False, verbose: bool=False, require_2d: bool=True) -> xarray.DataArray:
    assert numpy.isfinite(minlon) and numpy.isfinite(maxlon), 'Longitude range %s - %s is not valid' % (minlon, maxlon)
    assert numpy.isfinite(minlat) and numpy.isfinite(maxlat), 'Latitude range %s - %s is not valid' % (minlat, maxlat)
    assert minlon <= maxlon, 'Minimum longitude %s must be smaller than, or equal to, maximum longitude %s.' % (minlon, maxlon)
    assert minlat <= maxlat, 'Minimum latitude %s must be smaller than, or equal to, maximum latitude %s.' % (minlat, maxlat)
    source_lon, source_lat = source.getm.longitude, source.getm.latitude
    assert source_lon.ndim == 1, 'Source longitude must be 1D but has shape %s' % (source_lon.shape,)
    assert source_lat.ndim == 1, 'Source latitude must be 1D but has shape %s' % (source_lat.shape,)
    imin = source_lon.values.searchsorted(minlon, side='right') - 1
    imax = source_lon.values.searchsorted(maxlon, side='left') + 1
    if source_lat.values[1] < source_lat.values[0]:
        jmin = source_lat.size - source_lat.values[::-1].searchsorted(maxlat, side='left') - 1
        jmax = source_lat.size - source_lat.values[::-1].searchsorted(minlat, side='right') + 1
    else:
        jmin = source_lat.values.searchsorted(minlat, side='right') - 1
        jmax = source_lat.values.searchsorted(maxlat, side='left') + 1
    if verbose:
        print(imin, imax, source_lon.values.size, jmin, jmax, source_lat.values.size)
    assert (imin >= 0 and imax <= source_lon.values.size) or periodic_lon, 'Requested longitude section %s - %s is not fully covered by available range %s - %s' % (minlon, maxlon, source_lon.values[0], source_lon.values[-1])
    assert jmin >= 0 and jmax <= source_lat.values.size, 'Requested latitude section %s - %s is not fully covered by available range %s - %s' % (minlat, maxlat, source_lat.values[0], source_lat.values[-1])
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
    center_source = tuple([{ilondim: slice(imin, imax), ilatdim: slice(jmin, jmax)}.get(i, slice(None)) for i in range(len(shape))])
    center_target = [{ilondim: slice(0, imax - imin), ilatdim: slice(0, jmax - jmin)}.get(i, slice(None)) for i in range(len(shape))]
    target_lon = source_lon[center_source[ilondim]]
    target_lat = source_lat[center_source[ilatdim]]
    overlap = abs(source_lon.values[-1] - source_lon.values[0] - 360.) < 1e-5
    if verbose:
        print('periodic longitude? %s Overlap? %s = %s' % (periodic_lon, abs(source_lon.values[-1] - source_lon.values[0] - 360.), overlap))
    left_target = None
    right_target = None
    if add_left:
        # Periodic domain and we need to read beyond left boundary
        imin_left = source_lon.values.searchsorted(minlon + 360., side='right') - 1
        left_source = tuple([{ilondim: slice(imin_left, -1 if overlap else None)}.get(i, s) for i, s in enumerate(center_source)])
        nleft = source_lon.values.size - imin_left + (-1 if overlap else 0)
        if verbose:
            print('adding %i values on the left' % (nleft,))
        shape[ilondim] += nleft
        left_target = tuple([{ilondim: slice(0, nleft)}.get(i, s) for i, s in enumerate(center_target)])
        center_target[ilondim] = slice(nleft, nleft + imax - imin)
        target_lon = xarray.concat((source_lon[left_source[ilondim]] - 360., target_lon), source_lon.dims[0], combine_attrs='no_conflicts')
    if add_right:
        # Periodic domain and we need to read beyond right boundary
        imax_right = source_lon.values.searchsorted(maxlon - 360., side='left') + 1
        right_source = tuple([{ilondim: slice(1 if overlap else 0, imax_right)}.get(i, s) for i, s in enumerate(center_source)])
        if verbose:
            print('adding %i values on the right' % (imax_right + (-1 if overlap else 0),))
        shape[ilondim] += imax_right + (-1 if overlap else 0)
        right_target = tuple([{ilondim: slice(s.stop, None)}.get(i, s) for i, s in enumerate(center_target)])
        target_lon = xarray.concat((target_lon, source_lon[right_source[ilondim]] + 360.), source_lon.dims[0], combine_attrs='no_conflicts')
    center_target = tuple(center_target)
    shape = tuple(shape)
    if verbose:
        print('final shape: %s' % (shape,))

    lazyvar = SliceArray(source.variable, shape=shape, passthrough=[i for i in range(len(shape)) if i not in (ilondim, ilatdim)], name='limit_region(%s, minlon=%s, maxlon=%s, minlat=%s, maxlat=%s)' % (source.name, minlon, maxlon, minlat, maxlat))
    lazyvar._slices.append((center_source, center_target))
    if left_target:
        lazyvar._slices.append((left_source, left_target))
    if right_target:
        lazyvar._slices.append((right_source, right_target))

    coords = dict(source.coords.items())
    coords[source_lon.name] = target_lon
    coords[source_lat.name] = target_lat
    return xarray.DataArray(lazyvar, dims=source.dims, coords=coords, attrs=source.attrs, name=lazyvar.name)


def concatenate_slices(source: xarray.DataArray, idim: int, slices: Tuple[slice], verbose=False) -> xarray.DataArray:
    assert idim < source.ndim
    assert all([isinstance(s, slice) for s in slices])
    shape = list(source.shape)
    shape[idim] = sum([s.stop - s.start for s in slices])
    shape = tuple(shape)
    if verbose:
        print('final shape: %s' % (shape,))

    istart = 0
    strslices = ''
    final_slices = []
    for s in slices:
        n = s.stop - s.start
        source_slice = [slice(None)] * source.ndim
        target_slice = [slice(None)] * source.ndim
        source_slice[idim] = s
        target_slice[idim] = slice(istart, istart + n)
        strslices += '[%i:%i],' % (s.start, s.stop)
        final_slices.append((tuple(source_slice), tuple(target_slice)))
        istart += n
    assert istart == shape[idim]

    lazyvar = SliceArray(source.variable, shape=shape, passthrough=[i for i in range(len(shape)) if i != idim], dtype=source.dtype, name='concatenate_slices(%s, slices=(%s))' % (source.name, strslices))
    lazyvar._slices.extend(final_slices)

    coords = {}
    for name, c in source.coords.items():
        if source.dims[idim] not in c.dims:
            coords[name] = c
    return xarray.DataArray(lazyvar, dims=source.dims, coords=coords, attrs=source.attrs, name=lazyvar.name)


class Transpose(UnaryOperatorResult):
    def __array__(self, dtype=None) -> numpy.ndarray:
        return numpy.asarray(self._source).transpose()

    def __getitem__(self, slices) -> numpy.ndarray:
        return self._source[slices[::-1]].transpose()


def transpose(source: xarray.DataArray) -> xarray.DataArray:
    lazyvar = Transpose(source.variable, shape=source.shape[::-1], passthrough=list(range(source.ndim)), dtype=source.dtype, name='transpose(%s)' % (source.name,))
    coords = {}
    for name, c in source.coords.items():
        coords[name] = c.transpose()
    return xarray.DataArray(lazyvar, dims=source.dims[::-1], coords=coords, attrs=source.attrs, name=lazyvar.name)


def isel(source: xarray.DataArray, **indices) -> xarray.DataArray:
    """Index named dimensions with integers, slice objects or integer arrays"""
    advanced_indices = []
    for dim in list(indices):
        assert dim in source.dims, 'indexed dimension %s not used by source, which has dimensions %s' % (dim, source.dims)
        if not isinstance(indices[dim], (int, slice)):
            advanced_indices.append(source.dims.index(dim))
#            indices[dim] = xarray.Variable([source.dims[advanced_indices[0]]], numpy.asarray(indices[dim], dtype=numpy.intp))
            indices[dim] = xarray.Variable(['__newdim'], numpy.asarray(indices[dim], dtype=numpy.intp))

    # Final slices per dimension
    slices = tuple([indices.get(dim, slice(None)) for dim in source.dims])

    # Determine final shape
    shape = []
    dims = []
    passthrough = {}
    advanced_added = False
    for i, (dim, slc, l) in enumerate(zip(source.dims, slices, source.shape)):
        if i not in advanced_indices:
            # Slice is integer or slice object. if integer, it will be sliced out so it does not contribute to the final shape
            if isinstance(slc, slice):
                start, stop, stride = slc.indices(l)
                dims.append(dim)
                passthrough[len(shape)] = i
                shape.append((stop - start + stride - 1) // stride)
        elif not advanced_added:
            # First advanced slice. Add the shape produced by the broadcast combination of advanced indices
            assert max(advanced_indices) - min(advanced_indices) + 1 == len(advanced_indices), 'advanced indices must be side-by-side for now'
            for l in numpy.broadcast_shapes(*[indices[source.dims[i]].shape for i in advanced_indices]):
                dims.append('dim_%i' % len(shape))
                shape.append(l)
            advanced_added = True

    lazyvar = SliceArray(source.variable, shape=shape, passthrough=passthrough, dtype=source.dtype, name='isel(%s, %s)' % (source.name, ''.join([', %s=%s' % (name, value) for name, value in indices.items()])))
    lazyvar._slices.append((slices, (slice(None),) * len(shape)))

    coords = {}
    for name, c in source.coords.items():
        if all(dim not in c.dims for dim in indices):
            coords[name] = c
    return xarray.DataArray(lazyvar, dims=dims, coords=coords, attrs=source.attrs, name=lazyvar.name)


def horizontal_interpolation(source: xarray.DataArray, lon: xarray.DataArray, lat: xarray.DataArray, dtype: numpy.typing.DTypeLike=float, mask=None) -> xarray.DataArray:
    assert source.getm.longitude is not None, 'Variable %s does not have a valid longitude coordinate.' % source.name
    assert source.getm.latitude is not None, 'Variable %s does not have a valid latitude coordinate.' % source.name
    source_lon, source_lat = source.getm.longitude, source.getm.latitude
    assert source_lon.ndim == 1
    assert source_lat.ndim == 1
    assert numpy.isfinite(lon).all(), 'Some longitudes are non-finite: %s' % (lon,)
    assert numpy.isfinite(lat).all(), 'Some latitudes are non-finite: %s' % (lat,)
    lon, lat = numpy.broadcast_arrays(lon, lat)
    ilondim = source.dims.index(source_lon.dims[0])
    ilatdim = source.dims.index(source_lat.dims[0])
    assert abs(ilondim - ilatdim) == 1, 'Longitude and latitude dimensions must be distinct and adjacent'
    dimensions = {0: (), 1: (source_lon.dims[0],), 2: (source_lat.dims[0],  source_lon.dims[-1])}[lon.ndim]
    shape = source.shape[:min(ilondim, ilatdim)] + lon.shape + source.shape[max(ilondim, ilatdim) + 1:]
    kwargs = {'ndim_trailing': source.ndim - max(ilondim, ilatdim) - 1, 'mask': mask}
    if ilondim > ilatdim:
        # Dimension order: latitude first, then longitude
        ip = pygetm.util.interpolate.Linear2DGridInterpolator(lat, lon, source_lat, source_lon, **kwargs)
    else:
        # Dimension order: longitude first, then latitude
        dimensions = dimensions[::-1]
        ip = pygetm.util.interpolate.Linear2DGridInterpolator(lon, lat, source_lon, source_lat, **kwargs)
    lon_name, lat_name = source_lon.name, source_lat.name
    if lon_name in dimensions and lon.ndim > 1:
        lon_name = lon_name + '_'
    if lat_name in dimensions and lat.ndim > 1:
        lat_name = lat_name + '_'
    lon = xarray.DataArray(lon, dims=dimensions, name=lon_name, attrs=source_lon.attrs)
    lat = xarray.DataArray(lat, dims=dimensions, name=lat_name, attrs=source_lat.attrs)
    coords = dict([(k, v) for k, v in source.coords.items() if k not in {source_lon.name, source_lat.name}])
    coords[lon.name] = lon
    coords[lat.name] = lat
    dims = source.dims[:min(ilondim, ilatdim)] + dimensions + source.dims[max(ilondim, ilatdim) + 1:]
    lazyvar = SpatialInterpolation(ip, source.variable, shape, min(ilondim, ilatdim), source.ndim - max(ilondim, ilatdim) - 1, name='horizontal_interpolation(%s)' % source.name)
    return xarray.DataArray(lazyvar, dims=dims, coords=coords, attrs=source.attrs, name=lazyvar.name)


class SpatialInterpolation(UnaryOperatorResult):
    def __init__(self, ip: pygetm.util.interpolate.Linear2DGridInterpolator, source: xarray.Variable, shape: Iterable[int], npre: int, npost: int, **kwargs):
        UnaryOperatorResult.__init__(self, source, shape=shape, **kwargs)
        self._ip = ip
        self.npre = npre
        self.npost = npost

    def __array__(self, dtype=None) -> numpy.ndarray:
        return self._ip(numpy.asarray(self._source))

    def __getitem__(self, slices) -> numpy.ndarray:
        src_slice, tgt_slice = [Ellipsis], [Ellipsis]
        ntrailing_dim_removed = 0
        for i, s in enumerate(self._finalize_slices(slices)):
            if i < self.npre:
                # prefixed dimension
                src_slice.insert(i, s)
            elif i >= self.ndim - self.npost:
                # trailing dimension
                if isinstance(s, (int, numpy.integer)):
                    ntrailing_dim_removed += 1
                    tgt_slice.append(0)
                src_slice.append(s)
            else:
                assert isinstance(s, slice) and s.start is None and s.stop is None and s.step is None, '%s' % s
        source = numpy.asarray(self._source[tuple(src_slice)])
        source.shape = source.shape + (1,) * ntrailing_dim_removed
        result = self._ip(source)
        return result[tuple(tgt_slice)]


def vertical_interpolation(source: xarray.DataArray, target_z: numpy.ndarray, itargetdim: int=0) -> xarray.DataArray:
    source_z = source.getm.z
    assert source_z is not None, 'Variable %s does not have a valid depth coordinate.' % source.name
    assert source_z.ndim == 1
    izdim = source.dims.index(source_z.dims[0])
    #assert source.ndim - izdim == target_z.ndim
    #assert source.shape[izdim + 1:izdim + 3] == target_z.shape[1:], '%s vs %s' % (source.shape[izdim + 1:izdim + 3], target_z.shape[1:])
    target2sourcedim = {}
    isourcedim = 0
    for i, l in enumerate(target_z.shape):
        if i == itargetdim:
            isourcedim = izdim
        else:
            while isourcedim != izdim and isourcedim < source.ndim and l != source.shape[isourcedim]:
                isourcedim += 1
            assert isourcedim != izdim, 'Dimension with length %i should precede depth dimension %i in %s, which has shape %s' % (l, izdim, source.name, source.shape)
            assert isourcedim < source.ndim, 'Dimension with length %i expected after depth dimension %i in %s, which has shape %s' % (l, izdim, source.name, source.shape)
        target2sourcedim[i] = isourcedim
        isourcedim += 1
    coords = {}
    for n, c in source.coords.items():
        if n == source.dims[izdim]:
            coords[n + '_'] = ([source.dims[target2sourcedim[i]] for i in range(target_z.ndim)], target_z)
        else:
            coords[n] = c
    lazyvar = VerticalInterpolation(source.variable, target_z, izdim, source_z.values, itargetdim, name='horizontal_interpolation(%s)' % source.name)
    return xarray.DataArray(lazyvar, dims=source.dims, coords=coords, attrs=source.attrs, name=lazyvar.name)


class VerticalInterpolation(UnaryOperatorResult):
    def __init__(self, source: xarray.Variable, z: numpy.ndarray, izdim: int, source_z: numpy.ndarray, axis: int=0, **kwargs):
        self.izdim = izdim
        passthrough = [idim for idim in range(source.ndim) if idim != self.izdim]
        shape = list(source.shape)
        shape[self.izdim] = z.shape[axis]
        super().__init__(source, shape=shape, passthrough=passthrough, **kwargs)
        self.z = z
        self.axis = axis
        self.source_z = source_z
        if (self.source_z >= 0.).all():
            self.source_z = -self.source_z

    def apply(self, source, dtype=None) -> numpy.ndarray:
        return pygetm.util.interpolate.interp_1d(self.z, self.source_z, source, axis=self.axis)


def temporal_interpolation(source: xarray.DataArray, climatology: bool=False) -> xarray.DataArray:
    time_coord = source.getm.time
    assert time_coord is not None, 'No time coordinate found'
    itimedim = source.dims.index(time_coord.dims[0])
    lazyvar = TemporalInterpolationResult(source.variable, itimedim, time_coord.values, climatology, name='temporal_interpolation(%s)' % source.name)
    dims = [d for i, d in enumerate(source.dims) if i != lazyvar._itimedim]
    coords = dict(source.coords.items())
    coords[time_coord.dims[0]] = lazyvar._timecoord
    return xarray.DataArray(lazyvar, dims=dims, coords=coords, attrs=source.attrs, name=lazyvar.name)


class TemporalInterpolationResult(UnaryOperatorResult):
    __slots__ = '_current', '_itimedim', '_numnow', '_numnext', '_slope', '_inext', '_next', '_slices', 'climatology', '_year', '_timevalues'

    def __init__(self, source: xarray.Variable, itimedim: int, times: numpy.ndarray, climatology: bool, dtype=float, **kwargs):
        shape = list(source.shape)
        self._itimedim = itimedim
        assert shape[self._itimedim] > 1, 'Cannot interpolate %s in time because its time dimension has length %i.' % (source.name, shape[self._itimedim])
        shape.pop(self._itimedim)

        super().__init__(source, shape=shape, dtype=dtype, **kwargs)

        self._current = numpy.empty(shape, dtype=self.dtype)

        self.times = times
        self._timecoord = xarray.DataArray(self.times[0])
        self._timevalues = self._timecoord.values

        self._numnow = None
        self._numnext = 0.
        self._slope = 0.
        self._inext = -1
        self._next = 0.
        self._slices: List[Union[int, slice]] = [slice(None)] * source.ndim

        self.climatology = climatology
        assert not climatology or all(time.year == self.times[0].year for time in self.times)
        self._year = self.times[0].year

    def __array__(self, dtype=None) -> numpy.ndarray:
        return self._current

    def __getitem__(self, slices) -> numpy.ndarray:
        return self._current[slices]

    def is_time_varying(self) -> bool:
        return True

    def update(self, time: cftime.datetime, numtime: Optional[numpy.longdouble]=None) -> bool:
        if numtime is None:
            numtime = time.toordinal(fractional=True)

        if self._numnow is None:
            # First call to update - make sure the time series does not start after the requested time.
            if time.calendar != self.times[0].calendar:
                raise Exception('Simulation calendar %s does not match calendar %s used by %s.' % (time.calendar, self.times[0].calendar, self._source_name))
            if self.climatology:
                self._inext = self.times.searchsorted(time.replace(year=self._year), side='right') - 2
                self._year = time.year
                if self._inext < -1:
                    self._inext += self.times.size
                    self._year -= 1
            else:
                if time < self.times[0]:
                    raise Exception('Cannot interpolate %s to value at %s, because time series starts only at %s.' % (self._source_name, time.strftime(), self.times[0].strftime()))
                self._inext = self.times.searchsorted(time, side='right') - 2
            while self._inext < 1:
                self._move_to_next()
        elif numtime <= self._numnow:
            # Subsequent call to update, but time has not increased
            # If equal to previous time, we are done. If smaller, raise an exception
            if numtime == self._numnow:
                return False
            raise Exception('Time can only increase, but previous time was %s, new time %s' % (self._timevalues.flat[0].strftime(), time.strftime()))

        while self._numnext < numtime:
            self._move_to_next()

        # Do linear interpolation
        numpy.multiply(self._slope, numtime - self._numnext, out=self._current)
        self._current += self._next

        # Save current time
        self._numnow = numtime
        self._timevalues[...] = time
        return True

    def _move_to_next(self):
        # Move to next record
        self._inext += 1
        if self._inext == self.times.size:
            if self.climatology:
                self._inext = 0
                self._year += 1
            else:
                raise Exception('Cannot interpolate %s to value at %s because end of time series was reached (%s).' % (self._source_name, time.strftime(), self.times[-1].strftime()))
        old, numold = self._next, self._numnext
        self._slices[self._itimedim] = self._inext
        self._next = numpy.asarray(self._source[tuple(self._slices)], dtype=self.dtype)
        next_time = self.times[self._inext]
        if self.climatology:
            next_time = next_time.replace(year=self._year)
        self._numnext = next_time.toordinal(fractional=True)
        self._slope = (self._next - old) / (self._numnext - numold)

def debug_nc_reads(logger: Optional[logging.Logger]=None):
    """Hook into :mod:`xarray` so that every read from a NetCDF file is written to the log."""
    import xarray.backends.netCDF4_
    if logger is None:
        logger = logging.getLogger('pygetm.input')
        logger.setLevel(logging.DEBUG)
    class NetCDF4ArrayWrapper2(xarray.backends.netCDF4_.NetCDF4ArrayWrapper):
        __slots__ = ()
        def _getitem(self, key):
            logger.debug('Reading %s[%s] from %s' % (self.variable_name, key, self.datastore._filename))
            return super()._getitem(key)
    xarray.backends.netCDF4_.NetCDF4ArrayWrapper = NetCDF4ArrayWrapper2


class OnGrid(enum.Enum):
    NONE = enum.auto()        #: grids do not match. Spatially explicit data will require horizontal and - if vertically resolved - vertical interpolation.
    HORIZONTAL = enum.auto()  #: horizontal grid matches, but vertical does not. Vertically resolved data will require vertical interpolation.
    ALL = enum.auto()         #: horizontal and vertical grids match


class InputManager:
    def __init__(self):
        self.fields = []
        self._logger = logging.getLogger()

    def debug_nc_reads(self):
        """Hook into :mod:`xarray` so that every read from a NetCDF file is written to the log."""
        _logger = self._logger.getChild('nc')
        _logger.setLevel(logging.DEBUG)
        debug_nc_reads(_logger)

    def add(self, array: pygetm.core.Array, value: Union[numbers.Number, numpy.ndarray, xarray.DataArray, LazyArray], periodic_lon: bool=True, on_grid: Union[bool, OnGrid]=False, include_halos: Optional[bool]=None, climatology: bool=False, mask: bool=False):
        """Link an array to the provided input. If this input is constant in time, the value of the array will be set immediately.
        
        Args:
            array: array to assign a value to
            value: input to assign. If this is time-dependent, the combination of the array and its linked input will be
                registered; the array will then be updated to the current time whenever :meth:`update` is called.
            periodic_lon: whether this input covers all longitudes (i.e., the entire globe in the horizontal) and therefore
                has a periodic boundary. This enables efficient spatial interpolation across longitude bounds of the input,
                for instance, accessing read 10 degrees West to 5 degrees East for an input that spans 0 to 360 degrees East.
            on_grid: whether the input is defined on the same grid (horizontal-only, or both horizontal and vertical)
                as the array that is being assigned to. if this is ``False``, the value will be spatially interpolated
                to the array grid. ``True`` is equivalent to :attr:`OnGrid.HORIZONTAL`.
            include_halos: whether to also update the halos of the array. If not provided, this default to ``True`` if
                the array has attributes ``_require_halos`` or ``_part_of_state``; otherwise it defaults to ``False``.
            climatology: whether the input describes a single climatological year (at any temporal resolution, e.g.,
                monthly, daily) that is representative for any true year. This argument is relevant only if the provided
                input is time-varying. It also requires that the input does not span more than one year.
            mask: whether to set the array to its :attr:`pygetm.core.Array.fill_value` in all masked points.
                If not provided, only missing values in the input (NaNs) will be set to the fill value.
                This currently only has an effect when the input is non time-varying.
        """
        if array.all_values is None or array.all_values.size == 0:
            # The target variable does not contain data. Typically this is because it specifies information on the open boundaries,
            # of which the current (sub)domain does not have any.
            self._logger.warning('Ignoring asssignment to array %s because it has no associated data.' % array.name)
            return

        if isinstance(value, (numbers.Number, numpy.ndarray)):
            # Constant-in-time fill value. Set it, then forget about the array as it will not require further updating.
            array.fill(value)
            return

        assert isinstance(value, xarray.DataArray), 'If value is not numeric, it should be an xarray.DataArray, but it is %s (type %s).' % (value, type(value))

        if include_halos is None:
            include_halos = array.attrs.get('_require_halos', False) or array.attrs.get('_part_of_state', False)
        if not isinstance(on_grid, OnGrid):
            on_grid = OnGrid.HORIZONTAL if on_grid else OnGrid.NONE

        grid = array.grid

        # Obtain active area of local subdomain (including halos if include_halos is True)
        # and the corresponding slice in the global domain (always excluding halos)
        local_slice, global_slice, local_shape, global_shape = grid.domain.tiling.subdomain2slices(exclude_halos=not include_halos, halo_sub=2)

        target_slice = (Ellipsis,)
        if array.on_boundary:
            # Open boundary information. This can either be specified for the global domain (e.g., when read from netCDF),
            # or for only the open boundary points that fall within the local subdomain. Determine which of these.
            source_lon, source_lat = value.getm.longitude, value.getm.latitude
            if value.ndim >= 2 and value.shape[-2:] == global_shape:
                # on-grid data for the global domain - extract data at open boundary points
                value = isel(value, **{value.dims[-1]: grid.domain.open_boundaries.i_glob, value.dims[-2]: grid.domain.open_boundaries.j_glob})
                if array.z:
                    # open boundary arrays have z dimension last (fastest varying), but gridded 3D data have z first
                    value = transpose(value)
            elif source_lon is not None and source_lat is not None and source_lon.ndim > 0 and source_lat.ndim > 0:
                # Spatially explicit input: interpolate horizontally to open boundary coordinates
                if source_lon.ndim != 1:
                    raise Exception('Unsuitable shape %s of longitude coordinate %s. Off-grid boundary information can be used only if its longitude is 1D.' % (source_lon.shape, source_lon.name))
                if source_lat.ndim != 1:
                    raise Exception('Unsuitable shape %s of latitude coordinate %s. Off-grid boundary information can be used only if its latitude is 1D.' % (source_lat.shape, source_lat.name))
                ilondim = value.dims.index(source_lon.dims[0])
                ilatdim = value.dims.index(source_lat.dims[0])
                if ilondim != ilatdim:
                    lon, lat = grid.domain.open_boundaries.lon.all_values, grid.domain.open_boundaries.lat.all_values
                    value = limit_region(value, lon.min(), lon.max(), lat.min(), lat.max(), periodic_lon=periodic_lon)
                    value = pygetm.input.horizontal_interpolation(value, lon, lat)
            idim = value.ndim - (2 if array.z else 1)
            if value.shape[idim] == grid.domain.open_boundaries.np_glob:
                # The source array covers all open boundaries (global domain).
                # If the subdomain only has a subset of those, slice out only the points
                # that fall within the current subdomain
                if grid.domain.open_boundaries.local_to_global:
                    value = concatenate_slices(value, idim, [slice(start, stop) for (start, stop) in grid.domain.open_boundaries.local_to_global])
            elif value.shape[idim] != grid.domain.open_boundaries.np:
                raise Exception('Extent of dimension %i of %s is not compatible with open boundaries. It should have length %i (number of open boundary points in the global domain) or %i (number of open boundary points in the subdomain). Its actual extent is %i.' % (idim, value.name, grid.domain.open_boundaries.np_glob, grid.domain.open_boundaries.np, value.shape[idim]))
        elif array.ndim != 0:
            # The target is a normal 2D (horizontal-only) or 3D (depth-explicit) array
            # The source data can either be on the native model grid, or at an arbitrary lon, lat grid.
            # In the latter case, we interpolate in space.
            assert array.all_values.shape[-2:] == local_shape
            target_slice = local_slice
            if on_grid == OnGrid.NONE:
                # interpolate horizontally to local array INCLUDING halos
                lon, lat = grid.lon.all_values[target_slice], grid.lat.all_values[target_slice]
                assert not numpy.isnan(lon).any()
                assert not numpy.isnan(lat).any()
                value = limit_region(value, lon.min(), lon.max(), lat.min(), lat.max(), periodic_lon=periodic_lon)
                value = horizontal_interpolation(value, lon, lat)
            else:
                # the input is already on-grid, but we need to map from global domain to subdomain
                assert value.shape[-2:] == global_shape, "%s: shape of values %s should match that of global domain %s" % (array.name, value.shape[-2:], global_shape)
                value = value[global_slice]

        if value.getm.time is not None:
            # The source data is time-dependent; during the simulation it will be interpolated in time.
            if value.getm.time.size > 1:
                value = temporal_interpolation(value, climatology=climatology)
            elif value.getm.time.dims:
                self._logger.warning('%s is set to %s, which has only one time point %s. The value from this time will be used now. %s will not be further updated by the input manager at runtime.' % (array.name, value.name, value.getm.time.values.flat[0].strftime(), array.name))
                itimedim = value.dims.index(value.getm.time.dims[0])
                value = value[tuple([0 if idim == itimedim else slice(None) for idim in range(value.ndim)])]

        if array.z and on_grid != OnGrid.ALL:
            # The target is a depth-explicit array.
            # The source must be defined on z coordinates and interpolated to our [time-varying] depths
            if array.on_boundary:
                z_coordinate = grid.domain.open_boundaries.zc if array.z == CENTERS else grid.domain.open_boundaries.zf
            else:
                z_coordinate = grid.zc if array.z == CENTERS else grid.zf
            z_coordinate.saved = True
            value = vertical_interpolation(value, z_coordinate.all_values[target_slice], itargetdim=1 if array.on_boundary else 0)

        target = array.all_values[target_slice]
        assert value.shape == target.shape, 'Source shape %s does not match target shape %s' % (value.shape, target.shape)
        if isinstance(value.variable._data, LazyArray) and value.variable._data.is_time_varying():
            _3d_only = array.attrs.get('_3d_only', False)
            self._logger.info('%s will be updated dynamically from %s%s' % (array.name, value.name, ' on macrotimestep' if _3d_only else ''))
            self.fields.append((array.name, value.variable._data, target, not _3d_only))
        else:
            target[...] = value
            finite = numpy.isfinite(target)
            if array.ndim == 0 or array.on_boundary:
                unmasked = numpy.broadcast_to(True, target.shape)
            else:
                unmasked = numpy.broadcast_to(grid.mask.all_values[target_slice] != 0, target.shape)
                if array.fill_value is not None:
                    keep_mask = unmasked if mask else numpy.logical_or(unmasked, finite)
                    target[~keep_mask] = array.fill_value
            if not finite.all(where=unmasked):
                n_unmasked = unmasked.sum()
                self._logger.warning('%s is set to %s, which is not finite (e.g., NaN) in %i of %i unmasked points.' % (array.name, value.name, n_unmasked - finite.sum(where=unmasked), n_unmasked))
            self._logger.info('%s is set to time-invariant %s (minimum: %s, maximum: %s)' % (array.name, value.name, target.min(where=unmasked, initial=numpy.inf), target.max(where=unmasked, initial=-numpy.inf)))

    def update(self, time: cftime.datetime, include_3d: bool=True):
        """Update all arrays linked to time-dependent inputs to the current time.
        
        Args:
            time: current time
            include_3d: whether to also update arrays that were marked as only relevant for the macro (3D) time step
        """
        numtime = time.toordinal(fractional=True)
        for name, source, target, update_always in self.fields:
            if include_3d or update_always:
                self._logger.debug('updating %s' % name)
                source.update(time, numtime)
                target[...] = source

    @property
    def logger(self) -> logging.Logger:
        return self._logger

    def set_logger(self, logger: logging.Logger):
        self._logger = logger
