from typing import Optional

import numpy as np
import cftime
import xarray

from . import File


class MemoryFile(File):
    def start_now(self, *args, **kwargs):
        self._times = []
        self._seconds = []
        self._recorded_fields = {}
        for name, field in self.fields.items():
            if field.time_varying:
                self._recorded_fields[name] = (field, [])

    def save_now(self, seconds_passed: float, time: Optional[cftime.datetime]):
        self._times.append(time)
        self._seconds.append(seconds_passed)
        for field, values in self._recorded_fields.values():
            values.append(np.array(field.get(), copy=True))

    @property
    def x(self) -> xarray.Dataset:
        sec = xarray.Variable(("time",), self._seconds, {"units": "s"})
        timecoords = {"seconds": sec}
        if None not in self._times:
            timecoords["time"] = xarray.Variable(("time",), np.array(self._times))
        name2var = timecoords.copy()
        for name, field in self.fields.items():
            if not field.time_varying:
                values = np.ma.masked_equal(field.get(), field.fill_value)
                name2var[name] = xarray.Variable(field.dims, values, field.atts)
        for name, (field, values) in self._recorded_fields.items():
            dims = ("time",) + field.dims
            values = np.ma.masked_equal(values, field.fill_value)
            name2var[name] = xarray.Variable(dims, values, field.atts)
        renames = {}
        for name, var in name2var.items():
            if name in var.dims and var.ndim > 1:
                renames[name] = name + "_"
        arrays = {}
        for name, var in name2var.items():
            coords = {}
            if "time" in var.dims:
                coords.update(timecoords)
            field = self.fields.get(name)
            if field:
                for c in field.coordinates:
                    coords[renames.get(c, c)] = name2var[c]
            name = renames.get(name, name)
            arrays[name] = xarray.DataArray(var, coords=coords, name=name)
        return xarray.Dataset(arrays)
