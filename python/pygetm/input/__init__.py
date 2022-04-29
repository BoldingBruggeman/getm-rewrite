from typing import Iterable, List, Mapping, Union, Optional, Mapping, Sequence, Tuple
import glob
import numbers
import logging

import numpy
import numpy.typing
import numpy.lib.mixins
import xarray
import cftime

import pygetm.util.interpolate
from pygetm.constants import CENTERS

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
def from_nc(paths: Union[str, Sequence[str]], name: str, preprocess=None, cache=False, **kwargs) -> xarray.DataArray:
    key = (paths, preprocess, cache, kwargs.copy())
    for k, ds in open_nc_files:
        if k == key:
            break
    else:
        kwargs['decode_times'] = True
        kwargs['use_cftime'] = True
        kwargs['cache'] = cache
        if isinstance(paths, str):
            paths = glob.glob(paths)
        if len(paths) == 1:
            ds = xarray.open_dataset(paths[0], **kwargs)
            if preprocess:
                ds = preprocess(ds)
        else:
            ds = xarray.open_mfdataset(paths, preprocess=preprocess, **kwargs)
        open_nc_files.append((key, ds))
    array = ds[name]
    # Note: we wrap the netCDF array ourselves, in order to support lazy operators (e.g., add, multiply)
    return xarray.DataArray(WrappedArray(array), dims=array.dims, coords=array.coords, attrs=array.attrs, name='from_nc(%s, %s)' % (paths, name))

class LazyArray(numpy.lib.mixins.NDArrayOperatorsMixin):
    def __init__(self, shape: Iterable[int], dtype):
        self.shape = tuple(shape)
        self.ndim = len(self.shape)
        self.dtype = dtype

    def update(self, time: cftime.datetime, numtime: numpy.longdouble) -> bool:
        return False

    def astype(self, dtype, **kwargs) -> numpy.ndarray:
        return self.__array__(dtype)

    def __array_function__(self, func, types, args, kwargs):
        args = tuple(numpy.asarray(x) if isinstance(x, LazyArray) else x for x in args)
        kwargs = dict((k, numpy.asarray(v)) if isinstance(v, LazyArray) else (k, v) for (k, v) in kwargs.items())
        return func(*args, **kwargs)

    def __array_ufunc__(self, ufunc, method: str, *inputs, **kwargs):
        if method != '__call__':
            return NotImplemented

        if 'out' in kwargs:
            return NotImplemented

        for x in inputs:
            if not isinstance(x, (numpy.ndarray, numbers.Number, LazyArray, xarray.DataArray)):
                return NotImplemented

        return UFuncResult(getattr(ufunc, method), *inputs, **kwargs)

    def __array__(self, dtype=None) -> numpy.ndarray:
        raise NotImplementedError

    def __getitem__(self, slices) -> numpy.ndarray:
        return self.__array__()[slices]

    def is_time_varying(self) -> bool:
        return False

class OperatorResult(LazyArray):
    def __init__(self, *inputs, passthrough=(), dtype=None, shape=None, **kwargs):
        self.inputs = []
        for inp in inputs:
            if isinstance(inp, xarray.DataArray):
                inp = inp.variable
            if isinstance(inp, xarray.Variable) and isinstance(inp._data, WrappedArray):
                inp = inp._data.inputs[0]
            self.inputs.append(inp)
        self.input_names = [getattr(inp, 'name', None) for inp in inputs]
        self.lazy_inputs = []
        for input in inputs:
            if isinstance(input, xarray.DataArray) and isinstance(input.variable._data, LazyArray):
                self.lazy_inputs.append(input.variable._data)
            elif isinstance(input, LazyArray):
                self.lazy_inputs.append(input)
        self.kwargs = kwargs
        if shape is None:
            # Infer shape from inputs
            shapes = []
            for input in inputs:
                if isinstance(input, (numpy.ndarray, LazyArray, xarray.DataArray)):
                    shapes.append(input.shape)
            shape = numpy.broadcast_shapes(*shapes)
        if passthrough is True:
            passthrough = range(len(shape))
        if not isinstance(passthrough, dict):
            passthrough = dict([(i, i) for i in passthrough])
        self.passthrough = passthrough
        assert all([isinstance(dim, int) for dim in self.passthrough]), 'Invalid passthrough: %s. All entries should be of type int' % (self.passthrough,)
        super().__init__(shape, dtype or float)

    def update(self, *args) -> bool:
        updated = False
        for input in self.lazy_inputs:
            updated = input.update(*args) or updated
        return updated

    def is_time_varying(self) -> bool:
        return self.lazy_inputs and any(input.is_time_varying() for input in self.lazy_inputs)

    def __getitem__(self, slices) -> numpy.ndarray:
        assert isinstance(slices, tuple)
        for i, s in enumerate(slices):
            if s is Ellipsis:
                slices = slices[:i] + (slice(None),) * (self.ndim + 1 - len(slices)) + slices[i + 1:]
                break
        assert len(slices) == self.ndim
        preslices, postslices = [], []
        for i, slc in enumerate(slices):
            if i in self.passthrough:
                preslices.append(slc)
                if not isinstance(slc, int): postslices.append(slice(None))
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

class WrappedArray(OperatorResult):
    def __init__(self, source: xarray.DataArray):
        assert not isinstance(source.variable._data, LazyArray)
        super().__init__(source, passthrough=True, dtype=source.dtype)

    def apply(self, source, dtype=None) -> numpy.ndarray:
        return source

    def update(self, *args) -> bool:
        return False

class ConcatenatedSliceArray(UnaryOperatorResult):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._slices = []

    def __array__(self, dtype=None) -> numpy.ndarray:
        data = numpy.empty(self.shape, dtype or self.dtype)
        for src_slice, tgt_slice in self._slices:
            data[tgt_slice] = self._source[src_slice]
        return data

    def __getitem__(self, slices) -> numpy.ndarray:
        assert isinstance(slices, tuple)
        if Ellipsis in slices:
            i = slices.index(Ellipsis)
            slices = slices[:i] + (slice(None),) * (self.ndim + 1 - len(slices)) + slices[i + 1:]
        assert len(slices) == self.ndim
        shape = []
        for i, (l, s) in enumerate(zip(self.shape, slices)):
            if i in self.passthrough and isinstance(s, int):
                # This dimension will be sliced out
                continue
            start, stop, stride = s.indices(l)
            assert i in self.passthrough or (start == 0 and stop == l and stride == 1), 'invalid slice for dimension %i with length %i: %i:%i:%i' % (i, l, start, stop, stride)
            shape.append((stop - start + stride - 1) // stride)
        data = numpy.empty(shape, self.dtype)
        for src_slice, tgt_slice in self._slices:
            src_slice = list(src_slice)
            for iout, iin in self.passthrough.items():
                src_slice[iin] = slices[iout]
            tgt_slice = tuple([ori for i, (cust, ori) in enumerate(zip(slices, tgt_slice)) if i not in self.passthrough or not isinstance(cust, int)])
            data[tgt_slice] = self._source[tuple(src_slice)]
        return data

def limit_region(source: xarray.DataArray, minlon: float, maxlon: float, minlat: float, maxlat: float, periodic_lon=False, verbose=False) -> xarray.DataArray:
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

    data = ConcatenatedSliceArray(source, shape=shape, passthrough=[i for i in range(len(shape)) if i not in (ilondim, ilatdim)])
    data._slices.append((center_source, center_target))
    if left_target:
        data._slices.append((left_source, left_target))
    if right_target:
        data._slices.append((right_source, right_target))

    coords = dict(source.coords.items())
    coords[source_lon.name] = target_lon
    coords[source_lat.name] = target_lat
    return xarray.DataArray(data, dims=source.dims, coords=coords, attrs=source.attrs, name='limit_region(%s, minlon=%s, maxlon=%s, minlat=%s, maxlat=%s)' % (source.name, minlon, maxlon, minlat, maxlat))


def concatenate_slices(source: xarray.DataArray, idim: int, slices: Tuple[slice], verbose=False) -> xarray.DataArray:
    assert idim < source.ndim
    assert all([isinstance(s, slice) for s in slices])
    shape = list(source.shape)
    shape[idim] = sum([s.stop - s.start for s in slices])
    shape = tuple(shape)
    if verbose:
        print('final shape: %s' % (shape,))

    data = ConcatenatedSliceArray(source, shape=shape, passthrough=[i for i in range(len(shape)) if i != idim], dtype=source.dtype)

    istart = 0
    strslices = ''
    for s in slices:
        n = s.stop - s.start
        source_slice = [slice(None)] * source.ndim
        target_slice = [slice(None)] * source.ndim
        source_slice[idim] = s
        target_slice[idim] = slice(istart, istart + n)
        strslices += '[%i:%i],' % (s.start, s.stop)
        data._slices.append((tuple(source_slice), tuple(target_slice)))
        istart += n
    assert istart == shape[idim]

    coords = {}
    for name, c in source.coords.items():
        if source.dims[idim] not in c.dims:
            coords[name] = c
    return xarray.DataArray(data, dims=source.dims, coords=coords, attrs=source.attrs, name='concatenate_slices(%s, slices=(%s))' % (source.name, strslices))


class Transpose(UnaryOperatorResult):
    def __array__(self, dtype=None) -> numpy.ndarray:
        return numpy.asarray(self._source).transpose()

    def __getitem__(self, slices) -> numpy.ndarray:
        return self._source[slices[::-1]].transpose()

def transpose(source: xarray.DataArray) -> xarray.DataArray:
    data = Transpose(source, shape=source.shape[::-1], passthrough=list(range(source.ndim)), dtype=source.dtype)
    coords = {}
    for name, c in source.coords.items():
        coords[name] = c.transpose()
    return xarray.DataArray(data, dims=source.dims[::-1], coords=coords, attrs=source.attrs, name='transpose(%s)' % (source.name,))

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

    data = ConcatenatedSliceArray(source, shape=shape, passthrough=passthrough, dtype=source.dtype)
    data._slices.append((slices, (slice(None),) * len(shape)))

    coords = {}
    for name, c in source.coords.items():
        if all(dim not in c.dims for dim in indices):
            coords[name] = c
    return xarray.DataArray(data, dims=dims, coords=coords, attrs=source.attrs, name='isel(%s, %s)' % (source.name, ''.join([', %s=%s' % (name, value) for name, value in indices.items()])))

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
    return xarray.DataArray(SpatialInterpolation(ip, source, shape, min(ilondim, ilatdim), source.ndim - max(ilondim, ilatdim) - 1), dims=dims, coords=coords, attrs=source.attrs, name='horizontal_interpolation(%s)' % source.name)

class SpatialInterpolation(UnaryOperatorResult):
    def __init__(self, ip: pygetm.util.interpolate.Linear2DGridInterpolator, source: xarray.DataArray, shape: Iterable[int], npre: int, npost: int):
        UnaryOperatorResult.__init__(self, source, shape=shape)
        self._ip = ip
        self.npre = npre
        self.npost = npost

    def __array__(self, dtype=None) -> numpy.ndarray:
        return self._ip(self._source.values)

    def __getitem__(self, slices) -> numpy.ndarray:
        assert isinstance(slices, tuple)
        if Ellipsis in slices:
            i = slices.index(Ellipsis)
            slices = slices[:i] + (slice(None),) * (self.ndim + 1 - len(slices)) + slices[i + 1:]
        assert len(slices) == self.ndim
        src_slice, tgt_slice = [Ellipsis], [Ellipsis]
        for i, s in enumerate(slices):
            if i < self.npre:
                # prefixed dimension
                src_slice.insert(i, s)
            elif i >= self.ndim - self.npost:
                # trailing dimension
                if isinstance(s, int):
                    s = slice(s, s + 1)
                    tgt_slice.append(0)
                src_slice.append(s)
            else:
                assert isinstance(s, slice) and s.start is None and s.stop is None and s.step is None, '%s' % s
        source = self._source[tuple(src_slice)]
        result = self._ip(source.values)
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
    return xarray.DataArray(VerticalInterpolation(source, target_z, itargetdim), dims=source.dims, coords=coords, attrs=source.attrs, name='vertical_interpolation(%s)' % source.name)

class VerticalInterpolation(UnaryOperatorResult):
    def __init__(self, source: xarray.DataArray, z: numpy.ndarray, axis: int=0):
        source_z = source.getm.z
        self.izdim = source.dims.index(source_z.dims[0])
        passthrough = [idim for idim in range(source.ndim) if idim != self.izdim]
        shape = list(source.shape)
        shape[self.izdim] = z.shape[axis]
        UnaryOperatorResult.__init__(self, source, shape=shape, passthrough=passthrough)
        self.z = z
        self.axis = axis
        self.source_z = source_z.values
        if (self.source_z >= 0.).all():
            self.source_z = -self.source_z

    def apply(self, source, dtype=None) -> numpy.ndarray:
        return pygetm.util.interpolate.interp_1d(self.z, self.source_z, source, axis=self.axis)

def temporal_interpolation(source: xarray.DataArray, climatology: bool=False) -> xarray.DataArray:
    time_coord = source.getm.time
    assert time_coord is not None, 'No time coordinate found'
    result = TemporalInterpolationResult(source, climatology)
    dims = [d for i, d in enumerate(source.dims) if i != result._itimedim]
    coords = dict(source.coords.items())
    coords[time_coord.dims[0]] = result._timecoord
    return xarray.DataArray(result, dims=dims, coords=coords, attrs=source.attrs, name='temporal_interpolation(%s)' % source.name)

class TemporalInterpolationResult(UnaryOperatorResult):
    __slots__ = '_current', '_itimedim', '_numnow', '_numnext', '_slope', '_inext', '_next', '_slices', 'climatology', '_year'

    def __init__(self, source: xarray.DataArray, climatology: bool, dtype=float):
        time_coord = source.getm.time
        shape = list(source.shape)
        self._itimedim = source.dims.index(time_coord.dims[0])
        assert shape[self._itimedim] > 1, 'Cannot interpolate %s in time because its time dimension has length %i.' % (source.name, shape[self._itimedim])
        shape.pop(self._itimedim)

        super().__init__(source, shape=shape, dtype=dtype)

        self._current = numpy.empty(shape, dtype=self.dtype)

        self.times = time_coord.values
        self._timecoord = xarray.DataArray(self.times[0])

        self._numnow = None
        self._numnext = 0.
        self._slope = 0.
        self._inext = -1
        self._next = 0.
        self._slices: List[Union[int, slice]] = [slice(None)] * source.ndim
        self.name = source.name

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

        if numtime == self._numnow:
            return False

        if self._numnow is None:
            # First call to update - make sure the time series does not start after the requested time.
            if time.calendar != self.times[0].calendar:
                raise Exception('Simulation calendar %s does not match calendar %s used by %s.' % (time.calendar, self.times[0].calendar, self.name))
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
        elif numtime < self._numnow:
            # Subsequent call to update - make sure the requested time equals or exceeds the previously requested value
            raise Exception('Time can only increase, but previous time was %s, new time %s' % (self._timecoord.values.flat[0].strftime(), time.strftime()))

        while self._inext < 1 or self._numnext < numtime:
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
            self._next = numpy.asarray(self._source[tuple(self._slices)].values, dtype=self.dtype)
            next_time = self.times[self._inext]
            if self.climatology:
                next_time = next_time.replace(year=self._year)
            self._numnext = next_time.toordinal(fractional=True)
            self._slope = (self._next - old) / (self._numnext - numold)

        # Do linear interpolation
        self._current[...] = self._next + (numtime - self._numnext) * self._slope

        # Save current time
        self._numnow = numtime
        self._timecoord.values[...] = time
        return True

def debug_nc_reads(logger: Optional[logging.Logger]=None):
    """Hook into xarray so that every read from a NetCDF file is written to the log."""
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

class InputManager:
    def __init__(self):
        self.fields = []
        self._logger = logging.getLogger()

    def debug_nc_reads(self):
        """Hook into xarray so that every read from a NetCDF file is written to the log."""
        _logger = self._logger.getChild('nc')
        _logger.setLevel(logging.DEBUG)
        debug_nc_reads(_logger)

    def add(self, array, value: Union[numbers.Number, numpy.ndarray, xarray.DataArray, LazyArray], periodic_lon: bool=True, on_grid: bool=False, include_halos: Optional[bool]=None, climatology: bool=False):
        """Link an array to the provided input. If this input is constant in time, the value of the array will be set immediately.
        If the input is time-dependent, the array and its linked input will be registered with the input manager; the array
        will then be updated to the current time when InputManager.update is called."""
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

        grid = array.grid

        target_slice = (Ellipsis,)
        if array.on_boundary:
            # Open boundary information. This can either be specified for the global domain (e.g., when read from netCDF),
            # or for only the open boundary points that fall within the local subdomain. Determine which of these.
            if value.ndim >= 2 and value.shape[-1] == grid.domain.nx and value.shape[-2] == grid.domain.ny:
                # on-grid data for the global domain - extract data at open boundary points
                value = isel(value, **{value.dims[-1]: grid.domain.open_boundaries.i_glob, value.dims[-2]: grid.domain.open_boundaries.j_glob})
                if array.z:
                    # open boundary arrays have z dimension last (fastest varying), but gridded 3D data have z first
                    value = transpose(value)
            idim = value.ndim - (2 if array.z else 1)
            if value.shape[idim] != grid.domain.open_boundaries.np:
                # The source array covers all open boundaries (global domain).
                # Slice out only the points that fall within the current subdomain
                if value.shape[idim] != grid.domain.open_boundaries.np_glob:
                    self._logger.error('Dimension %i of %s does not have expected extent %i (number of open boundary points in the global domain). Its actual extent is %i' % (idim, value.name, grid.domain.open_boundaries.np_glob, value.shape[idim]))
                if grid.domain.open_boundaries.local_to_global:
                    value = concatenate_slices(value, idim, [slice(start, stop) for (start, stop) in grid.domain.open_boundaries.local_to_global])
        elif array.ndim != 0:
            # The target is a normal 2D (horizontal-only) or 3D (depth-explicit) array
            # The source data can either be on the native model grid, or at an arbitrary lon, lat grid.
            # In the latter case, we interpolate in space.
            if not on_grid:
                # interpolate horizontally to local array INCLUDING halos
                target_slice, _, _, _ = grid.domain.tiling.subdomain2slices(exclude_halos=not include_halos, halo_sub=2, halo_glob=2, exclude_global_halos=True)
                lon, lat = grid.lon.all_values[target_slice], grid.lat.all_values[target_slice]
                assert not numpy.isnan(lon).any()
                assert not numpy.isnan(lat).any()
                value = limit_region(value, lon.min(), lon.max(), lat.min(), lat.max(), periodic_lon=periodic_lon)
                value = horizontal_interpolation(value, lon, lat)
            else:
                # the input is already on-grid, but we may need to map from global domain to subdomain
                # source is global array EXCLUDING halos, but subdomain that we assign to does include halos (i.e., no halo exchange needed after assignment)
                target_slice, global_slice, _, _ = grid.domain.tiling.subdomain2slices(exclude_halos=not include_halos, halo_sub=2)
                value = value[global_slice]

        if value.getm.time is not None:
            # The source data is time-dependent; during the simulation it will be interpolated in time.
            if value.getm.time.size > 1:
                value = temporal_interpolation(value, climatology=climatology)
            elif value.getm.time.dims:
                self._logger.warning('%s is set to %s, which has only one time point %s. The value from this time will be used now. %s will not be further updated by the input manager at runtime.' % (array.name, value.name, value.getm.time.values.flat[0].strftime(), array.name))
                itimedim = value.dims.index(value.getm.time.dims[0])
                value = value[tuple([0 if idim == itimedim else slice(None) for idim in range(value.ndim)])]

        if array.z:
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
            unmasked = True if (array.ndim == 0 or array.on_boundary) else numpy.broadcast_to(grid.mask.all_values[target_slice] != 0, target.shape)
            finite = numpy.isfinite(target)
            if not finite.all(where=unmasked):
                n_unmasked = unmasked.sum()
                self._logger.warning('%s is set to %s, which is not finite (e.g., NaN) in %i of %i unmasked points.' % (array.name, value.name, n_unmasked - finite.sum(where=unmasked), n_unmasked))
            self._logger.info('%s is set to time-invariant %s (minimum: %s, maximum: %s)' % (array.name, value.name, target.min(where=unmasked, initial=numpy.inf), target.max(where=unmasked, initial=-numpy.inf)))

    def update(self, time: cftime.datetime, include_3d: bool=True):
        """Update all arrays linked to time-dependent inputs to the current time."""
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
