from typing import Optional, Mapping
import os.path
import logging
import datetime
import sys

import cftime
import netCDF4

from . import File
import pygetm.core
import pygetm._pygetm


class NetCDFFile(File):
    def __init__(
        self,
        available_fields: Mapping[str, pygetm.core.Array],
        logger: logging.Logger,
        path: str,
        rank: int,
        sync_interval: Optional[int] = 1,
        time_reference: Optional[cftime.datetime] = None,
        format: str = "NETCDF4",
        **kwargs,
    ):
        """Create a NetCDF file for output

        Args:
            available_fields: collection of model fields that may be added
            logger: target for log messages
            path: file to create. If it exists it will be clobbered.
            rank: rank of this subdomain. This will be used to determine whether
                we are the root (rank 0) all output is gathered and written to a single
                file. Otherwise the rank will be used as suffix for the
                subdomain-specific files.
            sync_interval: frequency to call NetCDF sync, which forces all output to
                be written to disk. If set to None, syncronization will happen only
                when the file is closed as the end of a simulation.
            time_reference: time reference (epoch) to use as offset for the time
                coordinate.
            format: underlying file format (see :class:`netCDF4.Dataset` documentation)
            **kwargs: additional keyword arguments passed to :class:`pygetm.output.File`
        """
        super().__init__(available_fields, logger, path=path, **kwargs)
        name, ext = os.path.splitext(path)
        if self.sub:
            name += f"_{rank:05}"
        self.path = name + ext
        self.nc = None
        self.itime = 0
        self.is_root = rank == 0
        self.time_offset = 0.0
        self.time_reference = time_reference
        self.sync_interval = sync_interval
        self.format = format
        self._field2nc = {}
        self._varying_fields = []
        self.nctime_bnds = None

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.path!r})"

    def start_now(
        self,
        seconds_passed: float,
        time: Optional[cftime.datetime],
        default_time_reference: Optional[cftime.datetime],
    ):
        if self.is_root or self.sub:
            # Create the NetCDF file
            self.nc = netCDF4.Dataset(self.path, "w", format=self.format)
            now = datetime.datetime.now()
            cmdline = " ".join(sys.argv)
            self.nc.history = f"{now:%Y-%m-%d %H:%M:%S} {cmdline}"
            self.nc.source = f"pygetm {pygetm._pygetm.get_version()}"
            for field in self.fields.values():
                for dim, length in zip(field.dims, field.shape):
                    if dim not in self.nc.dimensions:
                        self.nc.createDimension(dim, length)
                    elif length != self.nc.dimensions[dim].size:
                        raise Exception(
                            f"Existing dimension {dim} has incompatible length"
                            f" {self.nc.dimensions[dim].size} (need {length})"
                        )
            self.nc.createDimension("time")
            self.nctime = self.nc.createVariable("time", float, ("time",))
            self.nctime.axis = "T"
            if time is not None:
                time_reference = self.time_reference or default_time_reference
                self.nctime.units = f"seconds since {time_reference:%Y-%m-%d %H:%M:%S}"
                self.nctime.calendar = time.calendar
                self.time_offset = (
                    time - time_reference
                ).total_seconds() - seconds_passed
            else:
                self.nctime.units = "s"
                self.nctime.standard_name = "time"
            self.ncvars = []
            needs_time_bounds = False
            for output_name, field in self.fields.items():
                dims = field.dims
                if field.time_varying:
                    dims = ("time",) + dims
                ncvar = self.nc.createVariable(
                    output_name, field.dtype, dims, fill_value=field.fill_value
                )
                ncvar.set_auto_maskandscale(False)
                ncvar.expression = field.expression
                for att, value in field.attrs.items():
                    setattr(ncvar, att, value)
                if field.coordinates:
                    ncvar.coordinates = " ".join(field.coordinates)
                needs_time_bounds |= "time: mean" in field.attrs.get("cell_methods", "")
                self._field2nc[field] = ncvar
            if needs_time_bounds:
                self.nc.createDimension("nv", 2)
                self.nctime_bnds = self.nc.createVariable(
                    "time_bnds", float, ("time", "nv")
                )
                self.nctime.bounds = "time_bnds"
                self.previous_time_coord = self.time_offset + seconds_passed

        for field in self.fields.values():
            if field.time_varying:
                self._varying_fields.append(field)
            else:
                field.get(self._field2nc.get(field))

    def save_now(self, seconds_passed: float, time: Optional[cftime.datetime]):
        if self.nc is not None:
            time_coord = self.time_offset + seconds_passed
            self.nctime[self.itime] = time_coord
            if self.nctime_bnds is not None:
                self.nctime_bnds[self.itime, :] = [self.previous_time_coord, time_coord]
                self.previous_time_coord = time_coord
        for field in self._varying_fields:
            field.get(self._field2nc.get(field), slice_spec=(self.itime,))
        self.itime += 1
        if (
            self.nc is not None
            and self.sync_interval is not None
            and self.itime % self.sync_interval == 0
        ):
            self.nc.sync()

    def close_now(self, seconds_passed: float, time: Optional[cftime.datetime]):
        if self.nc is not None:
            self.nc.close()
