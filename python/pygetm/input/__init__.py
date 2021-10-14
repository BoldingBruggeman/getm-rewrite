from typing import Iterable, List, Mapping, Union, Optional, Mapping, Sequence, Tuple
import glob
import numbers

import numpy
import numpy.typing
import numpy.lib.mixins
import xarray
import cftime

import pygetm.util.interpolate

register = []

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
    def time(self) -> Optional[xarray.DataArray]:
        return self._interpret_coordinates().get('time')

    def _interpret_coordinates(self) -> Mapping[str, xarray.DataArray]:
        if self._coordinates is None:
            self._coordinates = {}
            #print(self._obj)
            for name, coord in self._obj.coords.items():
                units = coord.attrs.get('units')
                standard_name = coord.attrs.get('standard_name')
                #print(name, units, standard_name)
                if units in ('degrees_north', 'degree_north', 'degree_N', 'degrees_N', 'degreeN', 'degreesN') or standard_name == 'latitude':
                    self._coordinates['latitude'] = coord
                elif units in ('degrees_east', 'degree_east', 'degree_E', 'degrees_E', 'degreeE', 'degreesE') or standard_name == 'longitude':
                    self._coordinates['longitude'] = coord
                elif units is not None and ' since ' in units:
                    self._coordinates['time'] = coord
        return self._coordinates

    def update(self, time: cftime.datetime) -> bool:
        if isinstance(self._obj.data, LazyArray):
            return self._obj.data.update(time)
        return False

def get_from_nc(paths: Union[str, Sequence[str]], name: str, preprocess=None, cache=False, **kwargs) -> xarray.DataArray:
    kwargs['decode_times'] = False
    kwargs['cache'] = cache
    if isinstance(paths, str):
        paths = glob.glob(paths)
    if len(paths) == 1:
        ds = xarray.open_dataset(paths[0], **kwargs)
        if preprocess:
            ds = preprocess(ds)
    else:
        ds = xarray.open_mfdataset(paths, preprocess=preprocess, **kwargs)
    return ds[name]

class LazyArray(numpy.lib.mixins.NDArrayOperatorsMixin):
    def __init__(self, shape: Iterable[int], dtype, passthrough=()):
        self.shape = tuple(shape)
        self.ndim = len(self.shape)
        self.dtype = dtype
        self._slices = []
        self.passthrough = frozenset(passthrough)

    def update(self, time: cftime.datetime) -> bool:
        return False

    def astype(self, dtype, **kwargs):
        return self.__array__(dtype)

    def __array_function__(self, func, types, args, kwargs):
        args = tuple(numpy.asarray(x) if isinstance(x, LazyArray) else x for x in args)
        kwargs = dict((k, numpy.asarray(v)) if isinstance(v, LazyArray) else (k, v) for (k, v) in kwargs.items())
        return func(*args, **kwargs)

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        if method != '__call__':
            return NotImplemented

        out = kwargs.get('out', ())
        for x in inputs + out:
            if not isinstance(x, (numpy.ndarray, numbers.Number, LazyArray)):
                return NotImplemented

        inputs = tuple(numpy.asarray(x) if isinstance(x, LazyArray) else x for x in inputs)
        return getattr(ufunc, method)(*inputs, **kwargs)

    def __array__(self, dtype=None):
        raise NotImplementedError

    def __getitem__(self, slices):
        return numpy.asarray(self)[slices]

class UnaryOperatorResult(LazyArray):
    def __init__(self, source: xarray.DataArray, shape: Iterable[int]=(), dtype=None, passthrough=()):
        LazyArray.__init__(self, shape, dtype or source.dtype, passthrough)
        self._source = source

    def __array__(self, dtype=None):
        raise NotImplementedError

    def update(self, time: cftime.datetime) -> bool:
        return self._source.getm.update(time)

class LimitedRegionArray(UnaryOperatorResult):
    def __array__(self, dtype=None):
        data = numpy.empty(self.shape, dtype or self.dtype)
        for src_slice, tgt_slice in self._slices:
            data[tgt_slice] = self._source[src_slice]
        return data

    def __getitem__(self, slices):
        if Ellipsis in slices:
            i = slices.index(Ellipsis)
            slices = slices[:i] + (slice(None),) * (self.ndim + 1 - len(slices)) + slices[i + 1:]
        assert isinstance(slices, tuple)
        assert len(slices) == self.ndim
        shape = [l for i, (l, s) in enumerate(zip(self.shape, slices)) if i not in self.passthrough or not isinstance(s, int)]
        data = numpy.empty(shape, self.dtype)
        for src_slice, tgt_slice in self._slices:
            src_slice = tuple([(cust if i in self.passthrough else ori) for i, (cust, ori) in enumerate(zip(slices, src_slice))])
            tgt_slice = tuple([ori for i, (cust, ori) in enumerate(zip(slices, tgt_slice)) if i not in self.passthrough or not isinstance(cust, int)])
            data[tgt_slice] = self._source.variable[src_slice]
        return data

def limit_region(source: xarray.DataArray, minlon: float, maxlon: float, minlat: float, maxlat: float, periodic_lon=False, verbose=False):
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

    data = LimitedRegionArray(source, shape, passthrough=[i for i in range(len(shape)) if i not in (ilondim, ilatdim)])
    data._slices.append((center_source, center_target))
    if left_target:
        data._slices.append((left_source, left_target))
    if right_target:
        data._slices.append((right_source, right_target))

    coords = dict(source.coords.items())
    coords[source_lon.name] = target_lon
    coords[source_lat.name] = target_lat
    return xarray.DataArray(data, dims=source.dims, coords=coords, attrs=source.attrs)

def spatial_interpolation(source: xarray.DataArray, lon: xarray.DataArray, lat: xarray.DataArray, dtype: numpy.typing.DTypeLike=float, mask=None):
    assert source.getm.longitude is not None, 'Variable %s does not have a valid longitude coordinate.' % source.name
    assert source.getm.latitude is not None, 'Variable %s does not have a valid latitude coordinate.' % source.name
    source_lon, source_lat = source.getm.longitude, source.getm.latitude
    assert source_lon.ndim == 1
    assert source_lat.ndim == 1
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
    return xarray.DataArray(SpatialInterpolation(ip, source, shape), dims=dims, coords=coords, attrs=source.attrs)

class SpatialInterpolation(UnaryOperatorResult):
    def __init__(self, ip: pygetm.util.interpolate.Linear2DGridInterpolator, source: xarray.DataArray, shape: Iterable[int]):
        UnaryOperatorResult.__init__(self, source, shape)
        self._ip = ip

    def __array__(self, dtype=None):
        return self._ip(self._source.values)

def temporal_interpolation(source: xarray.DataArray):
    time_coord = source.getm.time
    if time_coord is None:
        return source
    result = TemporalInterpolationResult(source)
    dims = [d for i, d in enumerate(source.dims) if i != result._itimedim]
    coords = dict(source.coords.items())
    coords[time_coord.dims[0]] = result._timecoord
    return xarray.DataArray(result, dims=dims, coords=coords, attrs=source.attrs)

class TemporalInterpolationResult(UnaryOperatorResult):
    def __init__(self, source: xarray.DataArray):
        time_coord = source.getm.time
        shape = list(source.shape)
        self._itimedim = source.dims.index(time_coord.dims[0])
        assert shape[self._itimedim] > 1, 'Cannot interpolate %s in time because its time dimension has length %i.' % (source.name, shape[self._itimedim])
        shape.pop(self._itimedim)

        UnaryOperatorResult.__init__(self, source, shape)

        self._current = numpy.empty(shape, dtype=source.dtype)

        self._numtimes = time_coord.values
        self._time_units = time_coord.attrs['units']
        self._time_calendar = time_coord.attrs.get('calendar')
        self._timecoord = xarray.DataArray(cftime.num2date(self._numtimes[0], self._time_units, self._time_calendar))

        self._numnow = None
        self._numnext = 0.
        self._slope = 0.
        self._inext = -1
        self._next = 0.
        self.slices = [slice(None)] * source.ndim

    def __array__(self, dtype=None):
        return self._current

    def update(self, time: cftime.datetime) -> bool:
        # Convert the time into a numerical value matching units and calendar of the variable's time coordinate
        numtime = cftime.date2num(time, self._time_units, self._time_calendar)

        if numtime == self._numnow:
            return False

        if self._numnow is None:
            # First call to update - make sure the time series does not start after the requested time.
            if self._numtimes[0] > numtime:
                raise Exception('Cannot interpolate %s to value at %s, because time series starts after (%s).' % (self._source.name, numtime, self._numtimes[0]))
        elif numtime < self._numnow:
            # Subsequent call to update - make sure the requested time equals or exceeds the previously requested value
            raise Exception('Time can only increase, but previous time was %s, new time %s' % (self._numnow, numtime))

        while self._inext < 1 or self._numnext < numtime:
            # Move to next record
            self._inext += 1
            if self._inext == self._numtimes.size:
                raise Exception('Cannot interpolate %s to value at %s because end of time series was reached (%s).' % (self._source.name, numtime, self._numtimes[self._inext - 1]))
            old, numold = self._next, self._numnext
            self.slices[self._itimedim] = self._inext
            self._next = self._source.data[tuple(self.slices)]
            self._numnext = self._numtimes[self._inext]
            self._slope = (self._next - old) / (self._numnext - numold)

        # Do linear interpolation
        self._current[...] = self._next + (numtime - self._numnext) * self._slope

        # Save current time
        self._numnow = numtime
        self._timecoord.values[...] = time
        return True
