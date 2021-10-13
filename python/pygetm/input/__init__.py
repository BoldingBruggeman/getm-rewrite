from typing import List, Mapping, Union, Optional, Mapping, Sequence
import glob

import cftime
import xarray
import numpy
import numpy.typing
import scipy.interpolate
from xarray.core.dataarray import DataArray

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

    def interp(self, lon: xarray.DataArray, lat: xarray.DataArray, transpose: Optional[bool]=None) -> xarray.DataArray:
        assert self.longitude is not None, 'Variable %s does not have a valid longitude coordinate.' % self._obj.name
        assert self.latitude is not None, 'Variable %s does not have a valid latitude coordinate.' % self._obj.name
        lon, lat = numpy.broadcast_arrays(lon, lat)
        dimensions = {0: (), 1: (self.longitude.dims[0],), 2: (self.latitude.dims[0],  self.longitude.dims[-1])}[lon.ndim]
        lon_name, lat_name = self.longitude.name, self.latitude.name
        if lon_name in dimensions:
            lon_name = lon_name + '_'
        if lat_name in dimensions:
            lat_name = lat_name + '_'
        lon = xarray.DataArray(lon, dims=dimensions, name=lon_name)
        lat = xarray.DataArray(lat, dims=dimensions, name=lat_name)
        data = InterpolatedData(self._obj, self.longitude.values, self.latitude.values, lon.values, lat.values, transpose=transpose)
        return xarray.DataArray(data, dims=dimensions, coords={lon.name: lon, lat.name: lat})

class Variable:
    @classmethod
    def get(cls, *args, **kwargs) -> 'Variable':
        return cls(*args, **kwargs)

    def __init__(self, x: xarray.DataArray):
        self._x = x

    def update(self, time: cftime.datetime) -> bool:
        return False

    @property
    def x(self) -> xarray.DataArray:
        return self._x

class NetCDFVariable(Variable):
    @classmethod
    def get(cls, paths: Union[str, Sequence[str]], name: str, preprocess=None, **kwargs) -> 'NetCDFVariable':
        kwargs['decode_times'] = False
        if isinstance(paths, str):
            paths = glob.glob(paths)
        if len(paths) == 1:
            ds = xarray.open_dataset(paths[0], **kwargs)
            if preprocess:
                ds = preprocess(ds)
        else:
            ds = xarray.open_mfdataset(paths, preprocess=preprocess, **kwargs)
        return NetCDFVariable(ds[name])

class UnaryOperator(Variable):
    def __init__(self, source: Variable, result: xarray.DataArray):
        self.source = source
        Variable.__init__(self, result)
        self.apply()

    def update(self, time: cftime.datetime) -> bool:
        if not self.source.update(time):
            return False
        self.apply()
        return True

    def apply(self):
        pass

class LimitRegion(UnaryOperator):
    def __init__(self, source: Variable, minlon: float, maxlon: float, minlat: float, maxlat: float, periodic_lon=False, verbose=False):
        assert minlon <= maxlon, 'Minimum longitude %s must be smaller than, or equal to, maximum longitude %s.' % (minlon, maxlon)
        assert minlat <= maxlat, 'Minimum latitude %s must be smaller than, or equal to, maximum latitude %s.' % (minlat, maxlat)
        source_lon, source_lat = source.x.getm.longitude, source.x.getm.latitude
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
        ilondim = source.x.dims.index(source_lon.dims[0])
        ilatdim = source.x.dims.index(source_lat.dims[0])
        shape = list(source.x.shape)
        shape[ilondim] = imax - imin
        shape[ilatdim] = jmax - jmin
        self._center_source = tuple([{ilondim: slice(imin, imax), ilatdim: slice(jmin, jmax)}.get(i, slice(None)) for i in range(len(shape))])
        self._center_target = [{ilondim: slice(0, imax - imin), ilatdim: slice(0, jmax - jmin)}.get(i, slice(None)) for i in range(len(shape))]
        target_lon = source_lon[self._center_source[ilondim]]
        target_lat = source_lat[self._center_source[ilatdim]]
        overlap = abs(source_lon.values[-1] - source_lon.values[0] - 360.) < 1e-5
        if verbose:
            print('periodic longitude? %s Overlap? %s = %s' % (periodic_lon, abs(source_lon.values[-1] - source_lon.values[0] - 360.), overlap))
        self._left_target = None
        self._right_target = None
        if add_left:
            # Periodic domain and we need to read beyond left boundary
            imin_left = source_lon.values.searchsorted(minlon + 360., side='right') - 1
            self._left_source = tuple([{ilondim: slice(imin_left, -1 if overlap else None)}.get(i, s) for i, s in enumerate(self._center_source)])
            nleft = source_lon.values.size - imin_left + (-1 if overlap else 0)
            if verbose:
                print('adding %i values on the left' % (nleft,))
            shape[ilondim] += nleft
            self._left_target = tuple([{ilondim: slice(0, nleft)}.get(i, s) for i, s in enumerate(self._center_target)])
            self._center_target[ilondim] = slice(nleft, nleft + imax - imin)
            target_lon = xarray.concat((source_lon[self._left_source[ilondim]] - 360., target_lon), source_lon.dims[0], combine_attrs='no_conflicts')
        if add_right:
            # Periodic domain and we need to read beyond right boundary
            imax_right = source_lon.values.searchsorted(maxlon - 360., side='left') + 1
            self._right_source = tuple([{ilondim: slice(1 if overlap else 0, imax_right)}.get(i, s) for i, s in enumerate(self._center_source)])
            if verbose:
                print('adding %i values on the right' % (imax_right + (-1 if overlap else 0),))
            shape[ilondim] += imax_right + (-1 if overlap else 0)
            self._right_target = tuple([{ilondim: slice(s.stop, None)}.get(i, s) for i, s in enumerate(self._center_target)])
            target_lon = xarray.concat((target_lon, source_lon[self._right_source[ilondim]] + 360.), source_lon.dims[0], combine_attrs='no_conflicts')
        self._center_target = tuple(self._center_target)
        shape = tuple(shape)
        if verbose:
            print('final shape: %s' % (shape,))
        self._data = numpy.empty(shape, dtype=source.x.dtype)
        coords = dict(source.x.coords.items())
        coords[source_lon.name] = target_lon
        coords[source_lat.name] = target_lat
        result = xarray.DataArray(self._data, dims=source.x.dims, coords=coords, attrs=source.x.attrs)
        UnaryOperator.__init__(self, source, result)

    def apply(self):
        if self._left_target:
            self._data[self._left_target] = self.source.x.data[self._left_source]
        self._data[self._center_target] = self.source.x.data[self._center_source]
        if self._right_target:
            self._data[self._right_target] = self.source.x.data[self._right_source]

class SpatialInterpolation(UnaryOperator):
    def __init__(self, source: Variable, lon: xarray.DataArray, lat: xarray.DataArray, dtype: numpy.typing.DTypeLike=float):
        assert source.x.getm.longitude is not None, 'Variable %s does not have a valid longitude coordinate.' % source.x.name
        assert source.x.getm.latitude is not None, 'Variable %s does not have a valid latitude coordinate.' % source.x.name
        source_lon, source_lat = source.x.getm.longitude, source.x.getm.latitude
        assert source_lon.ndim == 1
        assert source_lat.ndim == 1
        lon, lat = numpy.broadcast_arrays(lon, lat)
        dimensions = {0: (), 1: (source_lon.dims[0],), 2: (source_lat.dims[0],  source_lon.dims[-1])}[lon.ndim]
        if source.x.dims.index(source_lon.dims[0]) > source.x.dims.index(source_lat.dims[0]):
            # Dimension order: latitude first, then longitude
            self._ip = pygetm.util.interpolate.Linear2DGridInterpolator(lat, lon, source_lat, source_lon)
        else:
            # Dimension order: longitude first, then latitude
            dimensions = dimensions[::-1]
            self._ip = pygetm.util.interpolate.Linear2DGridInterpolator(lon, lat, source_lon, source_lat)
        lon_name, lat_name = source_lon.name, source_lat.name
        if lon_name in dimensions and lon.ndim > 1:
            lon_name = lon_name + '_'
        if lat_name in dimensions and lat.ndim > 1:
            lat_name = lat_name + '_'
        lon = xarray.DataArray(lon, dims=dimensions, name=lon_name, attrs=source_lon.attrs)
        lat = xarray.DataArray(lat, dims=dimensions, name=lat_name, attrs=source_lat.attrs)
        self._data = numpy.empty(lon.shape, dtype=dtype)
        coords = dict([(k, v) for k, v in source.x.coords.items() if k not in {source_lon.name, source_lat.name}])
        coords[lon.name] = lon
        coords[lat.name] = lat
        result = xarray.DataArray(self._data, dims=dimensions, coords=coords, attrs=source.x.attrs)
        UnaryOperator.__init__(self, source, result)

    def apply(self):
        self._data[...] = self._ip(self.source.x.values)

class TemporalInterpolation(UnaryOperator):
    @classmethod
    def get(cls, source: Variable) -> Variable:
        time_coord = source.x.getm.time
        if time_coord is None:
            return source
        return TemporalInterpolation(source)

    def __init__(self, source: Variable):
        time_coord = source.x.getm.time
        assert time_coord is not None
        shape = list(source.x.shape)
        self._itimedim = source.x.dims.index(time_coord.dims[0])
        assert shape[self._itimedim] > 1, 'Cannot interpolate %s in time because its time dimension has length %i.' % (self.source.x.name, shape[self._itimedim])
        shape.pop(self._itimedim)
        self._current = numpy.empty(shape, dtype=source.x.dtype)

        self._numtimes = time_coord.values
        self._time_units = time_coord.attrs['units']
        self._time_calendar = time_coord.attrs.get('calendar')
        self._timecoord = xarray.DataArray(cftime.num2date(self._numtimes[0], self._time_units, self._time_calendar))

        dims = [d for i, d in enumerate(source.x.dims) if i != self._itimedim]
        coords = dict(source.x.coords.items())
        coords[time_coord.dims[0]] = self._timecoord

        self._numnow = None
        self._numnext = 0.
        self._slope = 0.
        self._inext = -1
        self._next = 0.

        result = xarray.DataArray(self._current, dims=dims, coords=coords, attrs=source.x.attrs)
        UnaryOperator.__init__(self, source, result)

    def update(self, time: cftime.datetime) -> bool:
        # Convert the time into a numerical value matching units and calendar of the variable's time coordinate
        numtime = cftime.date2num(time, self._time_units, self._time_calendar)

        if numtime == self._numnow:
            return False

        if self._numnow is None:
            # First call to update - make sure the time series does not start after the requested time.
            if self._numtimes[0] > numtime:
                raise Exception('Cannot interpolate %s to value at %s, because time series starts after (%s).' % (self.source.x.name, numtime, self._numtimes[0]))
        elif numtime < self._numnow:
            # Subsequent call to update - make sure the requested time equals or exceeds the previously requested value
            raise Exception('Time can only increase, but previous time was %s, new time %s' % (self._numnow, numtime))

        while self._inext < 1 or self._numnext < numtime:
            # Move to next record
            self._inext += 1
            if self._inext == self._numtimes.size:
                raise Exception('Cannot interpolate %s to value at %s because end of time series was reached (%s).' % (self.source.x.name, numtime, self._numtimes[self._inext - 1]))
            old, numold = self._next, self._numnext
            slices: List[Union[int, slice]] = [slice(None) for _ in self._current.shape]
            slices.insert(self._itimedim, self._inext)
            self._next = self.source.x.data[tuple(slices)]
            self._numnext = self._numtimes[self._inext]
            self._slope = (self._next - old) / (self._numnext - numold)

        # Do linear interpolation
        self._current[...] = self._next + (numtime - self._numnext) * self._slope

        # Save current time
        self._numnow = numtime
        self._timecoord.values[...] = time
        return True

class InterpolatedData:
    def __init__(self, xarray_obj: xarray.DataArray, source_lon, source_lat, target_lon, target_lat, transpose: Optional[bool]=None):
        self._obj = xarray_obj
        assert target_lon.shape == target_lat.shape
        assert source_lon.ndim == 1, 'Longitude of source grid must be one-dimensional, but has shape %s' % (source_lon.shape,)
        assert source_lat.ndim == 1, 'Latitude of source grid must be one-dimensional, but has shape %s' % (source_lat.shape,)
        if transpose is None:
            transpose = xarray_obj.shape == (source_lon.size, source_lat.size)
        minlon, maxlon = target_lon.min(), target_lon.max()
        minlat, maxlat = target_lat.min(), target_lat.max()
        assert source_lon[0] <= minlon and source_lon[-1] >= maxlon, 'Requested longitude range (%s - %s) does not fall completely within the source grid (%s - %s).' % (minlon, maxlon, source_lon[0], source_lon[-1])
        assert source_lat[0] <= minlat and source_lat[-1] >= maxlat, 'Requested latitude range (%s - %s) does not fall completely within the source grid (%s - %s).' % (minlat, maxlat, source_lat[0], source_lat[-1])
        self.imin = source_lon.searchsorted(minlon, side='right') - 1
        self.imax = source_lon.searchsorted(maxlon, side='left') + 1
        self.jmin = source_lat.searchsorted(minlat, side='right') - 1
        self.jmax = source_lat.searchsorted(maxlat, side='left') + 1
        self.source_lon = source_lon[self.imin:self.imax]
        self.source_lat = source_lat[self.jmin:self.jmax]
        self.target_lon = target_lon
        self.target_lat = target_lat
        self.shape = target_lon.shape
        self.dtype = self._obj.dtype
        self.transpose = transpose

    def __array_function__(self, func, types, args, kwargs):
        #print('__array_function__', func)
        return func(self.__array__(), *args, **kwargs)

    def __array__(self):
        #print('__array__')
        if self.transpose:
            source_data = self._obj[self.imin:self.imax, self.jmin:self.jmax].values.T
        else:
            source_data = self._obj[self.jmin:self.jmax, self.imin:self.imax].values
        #print(source_data)
        ip = scipy.interpolate.RectBivariateSpline(self.source_lat, self.source_lon, source_data, kx=1, ky=1)
        data = ip(self.target_lat, self.target_lon, grid=False)
        #print(data)
        data.shape = self.target_lat.shape
        return data

def request_from_netcdf(path: str, name: str) -> xarray.DataArray:
    ds = xarray.open_dataset(path)
    return ds[name]