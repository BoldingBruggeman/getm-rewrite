from typing import List, Mapping, Union, Optional, Mapping

import cftime
import xarray
import numpy
import scipy.interpolate
from xarray.core.dataarray import DataArray

register = []

@xarray.register_dataarray_accessor('getm')
class GETMAccessor:
    def __init__(self, xarray_obj: xarray.DataArray):
        self._obj = xarray_obj

        self._coordinates: Optional[Mapping[str, xarray.DataArray]] = None

        # Static time descriptors used by time interpolation
        self._itimedim = -1
        self._time_units = None
        self._time_calendar = None
        self._numtimes = None

        # Mutable helper variable for time interpolation:
        self._numnow = None
        self._current = None
        self._slope = None
        self._inext = -1
        self._inumnext = 0.
        self._next = 0.

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

    def now(self, time: cftime.datetime) -> numpy.ndarray:
        if self._current is None:
            time_coord = self.time
            if time_coord is not None:
                shape = list(self._obj.shape)
                self._itimedim = self._obj.dims.index(time_coord.dims[0])
                assert shape[self._itimedim] > 1, 'Cannot interpolate %s in time because its time dimension has length %i.' % (self._obj.name, shape[self._itimedim])
                self._numtimes = time_coord.values
                self._time_units = time_coord.attrs['units']
                self._time_calendar = time_coord.attrs.get('calendar')
                self._current = numpy.empty(shape, dtype=self._obj.dtype)
            else:
                self._current = self._obj.values

        # If this variable is time-independent, or already at the requested time, return directly
        if self._itimedim == -1:
            return self._current

        # Convert the time into a numerical value mathcing units and calendar of the variable's time coordinate
        numtime = cftime.date2num(time, self._time_units, self._time_calendar)

        if numtime == self._numtime:
            return self._current

        if self.now is None:
            # First call to update - make sure the time series does not start after the requested time.
            if self._numtimes[0] > numtime:
                raise Exception('Cannot interpolate %s to value at %s, because time series starts after (%s).' % (self._obj.name, numtime, self._numtimes[0]))
        elif numtime < self.now:
            # Subsequent call to update - make sure the requested time equals or exceeds the previously requested value
            raise Exception('Time can only increase, but previous time was %s, new time %s' % (self.now, numtime))

        while self._inext < 1 or self._numnext < numtime:
            # Move to next record
            self._inext += 1
            if self._inext == self._numtimes.size:
                raise Exception('Cannot interpolate %s to value at %s because end of time series was reached (%s).' % (self._obj.name, numtime, self._numtimes[self._inext - 1]))
            old, numold = self._next, self._numnext
            slices: List[Union[int, slice]] = [slice(None) for _ in self._current.shape]
            slices.insert(self._itimedim, self._inext)
            self._next = self._obj[tuple(slices)]
            self._slope = (self._next - old) / (self._next - numold)
            self._numnext = self._numtimes[self._inext]

        # Do linear interpolation
        self._current[...] = self._next + (numtime - self._numnext) * self._slope

        # Save current time
        self._numtime = numtime
        return self._current

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