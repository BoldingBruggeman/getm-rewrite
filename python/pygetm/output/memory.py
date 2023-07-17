from typing import Optional

import numpy as np
import cftime
import xarray as xr

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
    def x(self) -> xr.Dataset:
        # If the dataset is already cached, return that
        if self._x is not None:
            return self._x

        # Time coordinates
        sec = xr.Variable(("time",), self._seconds, {"units": "s"})
        timecoords = {"seconds": sec}
        if self._times is not None:
            timecoords["time"] = xr.Variable(("time",), np.array(self._times))

        renames = {}
        for name, field in self.fields.items():
            if name in field.dims and field.ndim > 1:
                renames[name] = name + "_"

        # Create variables
        name2var = timecoords.copy()
        for name, field in self.fields.items():
            coord_names = [renames.get(c, c) for c in field.coordinates]
            if field.time_varying:
                dims = ("time",) + field.dims
                values = self._recorded_fields[name][1]
                coord_names.append("time" if "time" in timecoords else "seconds")
            else:
                dims = field.dims
                values = field.get()
            attrs = field.attrs.copy()
            if coord_names:
                attrs["coordinates"] = " ".join(coord_names)
            name = renames.get(name, name)
            values = np.ma.masked_equal(values, field.fill_value)
            name2var[name] = xr.Variable(dims, values, attrs)

        # Create dataarrays (variables + coordinates)
        arrays = {}
        for name, var in name2var.items():
            coords = {}
            for c in var.attrs.get("coordinates", "").split():
                coords[c] = name2var[c]
            arrays[name] = xr.DataArray(var, coords=coords, name=name)

        # Return dataset with all arrays
        return xr.Dataset(arrays)

    def close_now(self, *args, **kwargs):
        self._x = self.x
        self._times = self._seconds = self._recorded_fields = None

    def __getattr__(self, name: str):
        return getattr(self.x, name)
