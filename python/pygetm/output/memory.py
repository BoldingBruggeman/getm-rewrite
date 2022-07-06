from typing import Optional

import numpy as np
import cftime
import xarray

from . import File


class MemoryFile(File):
    def start_now(
        self, seconds_passed: float, time: Optional[cftime.datetime], *args, **kwargs
    ):
        self._times = None if time is None else []
        self._seconds = []
        self._recorded_fields = {}
        for name, field in self.fields.items():
            if field.time_varying:
                self._recorded_fields[name] = (field, [])
        self._x = None

    def save_now(self, seconds_passed: float, time: Optional[cftime.datetime]):
        self._seconds.append(seconds_passed)
        if self._times is not None:
            self._times.append(time)
        for field, values in self._recorded_fields.values():
            values.append(np.array(field.get(), copy=True))

    @property
    def x(self) -> xarray.Dataset:
        # If the dataset is already cached, return that
        if self._x is not None:
            return self._x

        # Time coordinates
        sec = xarray.Variable(("time",), self._seconds, {"units": "s"})
        timecoords = {"seconds": sec}
        if self._times is not None:
            timecoords["time"] = xarray.Variable(("time",), np.array(self._times))

        # Create variables
        name2var = timecoords.copy()
        renames = {}
        for name, field in self.fields.items():
            if field.time_varying:
                dims = ("time",) + field.dims
                values = self._recorded_fields[name][1]
            else:
                dims = field.dims
                values = field.get()
            values = np.ma.masked_equal(values, field.fill_value)
            name2var[name] = xarray.Variable(dims, values, field.atts)
            if name in dims and len(dims) > 1:
                renames[name] = name + "_"

        # Create dataarrays (variables + coordinates)
        arrays = {}
        for name, var in name2var.items():
            coords = {}
            if "time" in var.dims:
                coords.update(timecoords)
            if name in self.fields:
                for c in self.fields[name].coordinates:
                    coords[renames.get(c, c)] = name2var[c]
            name = renames.get(name, name)
            arrays[name] = xarray.DataArray(var, coords=coords, name=name)

        # Return dataset with all arrays
        return xarray.Dataset(arrays)

    def close_now(self, seconds_passed: float, time: Optional[cftime.datetime]):
        self._x = self.x
        self._times = self._seconds = self._recorded_fields = None
