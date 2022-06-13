from typing import Optional, Mapping
import os.path
import logging

import cftime
import netCDF4

from . import File
from ..constants import INTERFACES
from .operators import TimeVarying
import pygetm.core


class NetCDFFile(File):
    def __init__(
        self,
        available_fields: Mapping[str, pygetm.core.Array],
        logger: logging.Logger,
        path: str,
        rank: int,
        sub: bool = False,
        sync_interval: Optional[int] = 1,
        time_reference: Optional[cftime.datetime] = None,
        **kwargs
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
            sub: whether to write to separate files per subdomain
            sync_interval: frequency to call NetCDF sync, which forces all output to
                be written to disk. If set to None, syncronization will happen only
                when the file is closed as the end of a simulation.
            time_reference: time reference (epoch) to use as offset for the time
                coordinate.
            **kwargs: additional keyword arguments passed to :class:`pygetm.output.File`
        """
        super().__init__(available_fields, logger, path=path, **kwargs)
        name, ext = os.path.splitext(path)
        if sub:
            name += "_%05i" % rank
        self.path = name + ext
        self.nc = None
        self.itime = 0
        self.is_root = rank == 0
        self.sub = sub
        self.created = False
        self.time_offset = 0.0
        self.time_reference = time_reference
        self.sync_interval = sync_interval
        self._field2nc = {}
        self._varying_fields = []

    def __repr__(self) -> str:
        return "%s(%r)" % (self.__class__.__name__, self.path)

    def _create(self, seconds_passed: float, time: Optional[cftime.datetime]):
        self.created = True

        if self.is_root or self.sub:
            # Create the NetCDF file
            self.nc = netCDF4.Dataset(self.path, "w")
            for grid in frozenset(field.grid for field in self.fields.values()):
                if self.sub:
                    # Output data (including halos) for the current subdomain
                    nx, ny, nz = grid.nx_, grid.ny_, grid.nz_
                else:
                    # Output data for entire (global) domain
                    nx = grid.domain.tiling.nx_glob + grid.nx - grid.domain.T.nx
                    ny = grid.domain.tiling.ny_glob + grid.ny - grid.domain.T.ny
                    nz = grid.nz
                xname, yname, zname, ziname = (
                    "x%s" % grid.postfix,
                    "y%s" % grid.postfix,
                    "z",
                    "zi",
                )
                if xname not in self.nc.dimensions:
                    self.nc.createDimension(xname, nx)
                if yname not in self.nc.dimensions:
                    self.nc.createDimension(yname, ny)
                if zname not in self.nc.dimensions:
                    self.nc.createDimension(zname, nz)
                if ziname not in self.nc.dimensions:
                    self.nc.createDimension(ziname, nz + 1)
            self.nc.createDimension("time",)
            self.nctime = self.nc.createVariable("time", float, ("time",))
            if time is not None:
                self.nctime.units = "seconds since %s" % self.time_reference.strftime(
                    "%Y-%m-%d %H:%M:%S"
                )
                self.nctime.calendar = time.calendar
                self.time_offset = (
                    time - self.time_reference
                ).total_seconds() - seconds_passed
            else:
                self.nctime.units = "s"
            self.ncvars = []
            for output_name, field in self.fields.items():
                dims = ("y%s" % field.grid.postfix, "x%s" % field.grid.postfix)
                if field.z:
                    dims = ("zi" if field.z == INTERFACES else "z",) + dims
                if field.time_varying != TimeVarying.NO:
                    dims = ("time",) + dims
                ncvar = self.nc.createVariable(
                    output_name, field.dtype, dims, fill_value=field.fill_value
                )
                ncvar.set_auto_maskandscale(False)
                for att, value in field.atts.items():
                    setattr(ncvar, att, value)
                if field.coordinates:
                    setattr(ncvar, "coordinates", " ".join(field.coordinates))
                self._field2nc[field] = ncvar

        for field in self.fields.values():
            if field.time_varying == TimeVarying.NO:
                field.get(self._field2nc.get(field), sub=self.sub)
            else:
                self._varying_fields.append(field)

    def start_now(
        self,
        itimestep: int,
        time: Optional[cftime.datetime],
        default_time_reference: Optional[cftime.datetime],
    ):
        if self.time_reference is None:
            self.time_reference = default_time_reference

    def save_now(self, seconds_passed: float, time: Optional[cftime.datetime]):
        if not self.created:
            self._create(seconds_passed, time)
        if self.nc is not None:
            self.nctime[self.itime] = self.time_offset + seconds_passed
        for field in self._varying_fields:
            field.get(self._field2nc.get(field), slice_spec=(self.itime,), sub=self.sub)
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
