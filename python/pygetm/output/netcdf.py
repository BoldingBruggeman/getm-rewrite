from typing import Optional, Mapping, Dict, List, Tuple, Union
import os
import logging
import datetime
import sys
from pathlib import Path

import cftime
import netCDF4

from . import File
from . import operators
import pygetm.core
import pygetm._pygetm


def _create_file(path: Union[os.PathLike, str], **kwargs) -> netCDF4.Dataset:
    nc = netCDF4.Dataset(path, "w", **kwargs)
    now = datetime.datetime.now()
    cmdline = " ".join(sys.argv)
    nc.history = f"{now:%Y-%m-%d %H:%M:%S} {cmdline}"
    nc.source = f"pygetm {pygetm._pygetm.get_version()}"
    return nc


def _create_dimensions(
    nc: netCDF4.Dataset, fields: Mapping[str, operators.Base]
) -> Tuple[bool, bool]:
    needs_time = False
    needs_time_bounds = False
    for output_name, field in fields.items():
        for dim, length in zip(field.dims, field.shape):
            if dim not in nc.dimensions:
                assert length > 0
                nc.createDimension(dim, length)
            elif length != nc.dimensions[dim].size:
                raise Exception(
                    f"Error adding {output_name} with shape {field.shape}:"
                    f" existing dimension {dim} has incompatible length"
                    f" {nc.dimensions[dim].size} (needs {length})"
                )
        needs_time |= bool(field.time_varying)
        needs_time_bounds |= "time: mean" in field.attrs.get("cell_methods", "")

    if needs_time_bounds:
        nc.createDimension("nv", 2)
    if needs_time:
        nc.createDimension("time", None)
    return needs_time, needs_time_bounds


def _add_time_coordinate(
    nc: netCDF4.Dataset,
    time: Optional[cftime.datetime],
    seconds_passed: float,
    time_reference: Optional[cftime.datetime],
) -> Tuple[netCDF4.Variable, float]:
    nctime = nc.createVariable("time", float, ("time",))
    nctime.axis = "T"
    if time is not None:
        nctime.units = f"seconds since {time_reference:%Y-%m-%d %H:%M:%S}"
        nctime.calendar = time.calendar
        time_offset = (time - time_reference).total_seconds() - seconds_passed
    else:
        nctime.units = "s"
        nctime.standard_name = "time"
        time_offset = 0.0
    return nctime, time_offset


def _add_time_bounds(nc: netCDF4.Dataset, nctime: netCDF4.Variable) -> netCDF4.Variable:
    nctime_bnds = nc.createVariable("time_bnds", float, ("time", "nv"))
    nctime.bounds = "time_bnds"
    for att in ("units", "calendar"):
        if hasattr(nctime, att):
            setattr(nctime_bnds, att, getattr(nctime, att))
    return nctime_bnds


def _add_variable(
    nc: netCDF4.Dataset, name: str, field: operators.Base, **kwargs
) -> netCDF4.Variable:
    dims = field.dims
    if field.time_varying:
        dims = ("time",) + dims
    ncvar = nc.createVariable(
        name, field.dtype, dims, fill_value=field.fill_value, **kwargs
    )

    # Use try-except for set_auto_maskandscale because only some NetCDF
    # engines support it (netCDF4 does, h5netcdf.legacyapi does not)
    try:
        ncvar.set_auto_maskandscale(False)
    except AttributeError:
        pass

    # Variable attributes
    ncvar.expression = field.expression
    for att, value in field.attrs.items():
        setattr(ncvar, att, value)
    if field.coordinates:
        ncvar.coordinates = " ".join(field.coordinates)

    return ncvar


class NetCDFFile(File):
    def __init__(
        self,
        available_fields: Mapping[str, pygetm.core.Array],
        logger: logging.Logger,
        path: Union[os.PathLike[str], str],
        rank: int,
        sync_interval: Optional[int] = 1,
        time_reference: Optional[cftime.datetime] = None,
        format: str = "NETCDF4",
        compression: Optional[str] = None,
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
                be written to disk. If set to None, synchronization will happen only
                when the file is closed as the end of a simulation.
            time_reference: time reference (epoch) to use as offset for the time
                coordinate.
            format: underlying file format (see :class:`netCDF4.Dataset` documentation)
            compression: compression algorithm to apply to all variables
                (see :meth:`netCDF4.Dataset.createVariable` documentation)
            **kwargs: additional keyword arguments passed to :class:`pygetm.output.File`
        """
        super().__init__(available_fields, logger, **kwargs)
        path = Path(path)
        if self.sub:
            path = path.with_stem(f"{path.stem}_{rank:05}")
        self.path = path
        self.nc: Optional[netCDF4.Dataset] = None
        self.itime = 0
        self.is_root = rank == 0
        self.time_offset = 0.0
        self.time_reference = time_reference
        self.sync_interval = sync_interval
        self.format = format
        self.compression = compression
        self._field2nc: Dict[operators.Base, netCDF4.Variable] = {}
        self._varying_fields: List[operators.Base] = []
        self.nctime: Optional[netCDF4.Variable] = None
        self.nctime_bnds: Optional[netCDF4.Variable] = None

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}('{self.path}')"

    def start_now(
        self,
        seconds_passed: float,
        time: Optional[cftime.datetime],
        default_time_reference: Optional[cftime.datetime],
    ) -> bool:
        if self.is_root or self.sub:
            included_fields = {}
            for output_name, field in self.fields.items():
                if 0 in field.shape:
                    self._logger.warning(
                        f"Skipping {output_name} because it contains no data"
                        f" (shape={field.shape})"
                    )
                else:
                    included_fields[output_name] = field

            if included_fields:
                self.nc = _create_file(self.path, format=self.format)
                has_time, has_time_bounds = _create_dimensions(self.nc, included_fields)
                if has_time:
                    time_reference = self.time_reference or default_time_reference
                    self.nctime, self.time_offset = _add_time_coordinate(
                        self.nc, time, seconds_passed, time_reference
                    )
                    if has_time_bounds:
                        self.nctime_bnds = _add_time_bounds(self.nc, self.nctime)
                        self.previous_time_coord = self.time_offset + seconds_passed
                for output_name, field in included_fields.items():
                    self._field2nc[field] = _add_variable(
                        self.nc, output_name, field, compression=self.compression
                    )

        for field in self.fields.values():
            if field.time_varying:
                # Store field for update at each time step
                self._varying_fields.append(field)
            else:
                # Write static field now
                field.get(self._field2nc.get(field))

        return len(self._varying_fields) > 0

    def save_now(self, seconds_passed: float, time: Optional[cftime.datetime]):
        # Update time coordinate(s), if used
        if self.nctime is not None:
            time_coord = self.time_offset + seconds_passed
            self.nctime[self.itime] = time_coord
            if self.nctime_bnds is not None:
                self.nctime_bnds[self.itime, :] = [self.previous_time_coord, time_coord]
                self.previous_time_coord = time_coord

        # Update time-varying fields
        for field in self._varying_fields:
            field.get(self._field2nc.get(field), slice_spec=(self.itime,))

        # Increment time index and sync to disk if needed
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
