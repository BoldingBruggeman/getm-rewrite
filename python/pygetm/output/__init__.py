import logging
from typing import Mapping, Optional, List, Union
import datetime
import enum

from numpy.typing import DTypeLike
import cftime

from .. import core
from . import operators


class TimeUnit(enum.Enum):
    TIMESTEPS = enum.auto()
    SECONDS = enum.auto()
    MONTHS = enum.auto()
    YEARS = enum.auto()


class File(operators.FieldCollection):
    def __init__(
        self,
        available_fields: Mapping[str, core.Array],
        logger: logging.Logger,
        interval: Union[datetime.timedelta, int] = 1,
        interval_units: TimeUnit = TimeUnit.TIMESTEPS,
        path: Optional[str] = None,
        default_dtype: Optional[DTypeLike] = None,
    ):
        super().__init__(available_fields, default_dtype=default_dtype)
        self._logger = logger
        self.next = None
        if isinstance(interval, datetime.timedelta):
            interval = interval.total_seconds()
            interval_units = TimeUnit.SECONDS
        self.interval = interval
        self.interval_units = interval_units
        self.save_on_close_only = self.interval_units == TimeUnit.TIMESTEPS and self.interval == -1
        self.path = path

    def _next_time(
        self, seconds_passed: float, itimestep: int, time: Optional[cftime.datetime],
    ) -> Union[float, int]:
        if self.interval_units == TimeUnit.SECONDS:
            return seconds_passed + self.interval
        elif self.interval_units == TimeUnit.TIMESTEPS:
            return itimestep + self.interval
        elif self.interval_units == TimeUnit.MONTHS:
            next_month = time.month + self.interval
            next_year = time.year + (next_month - 1) // 12
            next_month = 1 + (next_month - 1) % 12
            next_time = time.replace(month=next_month, year=next_year)
            return seconds_passed + (next_time - time).total_seconds()
        elif self.interval_units == TimeUnit.YEARS:
            next_time = time.replace(year=time.year + self.interval)
            return seconds_passed + (next_time - time).total_seconds()

    def start(
        self,
        itimestep: int,
        time: Optional[cftime.datetime],
        save: bool,
        default_time_reference: Optional[cftime.datetime],
    ):
        self.start_now(itimestep, time, default_time_reference)
        if save and not self.save_on_close_only:
            self._logger.debug("Saving initial state")
            self.save_now(0.0, time)
        self.next = self._next_time(0.0, itimestep, time)

    def save(
        self,
        seconds_passed: float,
        itimestep: int,
        time: Optional[cftime.datetime],
        macro: bool,
    ):
        self.update(macro=macro)
        if self.save_on_close_only:
            return
        now = itimestep if self.interval == TimeUnit.TIMESTEPS else seconds_passed
        if now >= self.next:
            self._logger.debug("Saving")
            self.save_now(seconds_passed, time)
            self.next = self._next_time(seconds_passed, itimestep, time)

    def close(self, seconds_passed: float, time: Optional[cftime.datetime]):
        if self.save_on_close_only:
            self.save_now(seconds_passed, time)
        self._logger.debug("Closing")
        self.close_now(seconds_passed, time)

    def start_now(
        self,
        itimestep: int,
        time: Optional[cftime.datetime],
        default_time_reference: Optional[cftime.datetime],
    ):
        pass

    def save_now(self, seconds_passed: float, time: Optional[cftime.datetime]):
        raise NotImplementedError

    def close_now(self, seconds_passed: float, time: Optional[cftime.datetime]):
        pass


# Note: netcdf cannot be imported before File is defined,
# because that is imported by netcdf itself.
# If this import statement came earlier, it would cause a circular dependency error.
from . import netcdf


class OutputManager:
    def __init__(
        self,
        fields: Mapping[str, core.Array],
        rank: int,
        logger: Optional[logging.Logger] = None,
    ):
        self.fields = fields
        self.rank = rank
        self.files: List[File] = []
        self._logger = logger or logging.getLogger()

    def add_netcdf_file(self, path: str, **kwargs) -> netcdf.NetCDFFile:
        """Add a NetCDF file for output.

        Args:
            path: NetCDF file to write to.
            **kwargs: additional keyword arguments passed to
                :class:`pygetm.output.netcdf.NetCDFFile`
        """
        self._logger.debug("Adding NetCDF file %s" % path)
        file = netcdf.NetCDFFile(
            self.fields, self._logger.getChild(path), path, rank=self.rank, **kwargs
        )
        self.files.append(file)
        return file

    def add_restart(self, path: str, **kwargs) -> netcdf.NetCDFFile:
        """Add a restart file to write to.

        Args:
            path: NetCDF file to write to.
            **kwargs: additional keyword arguments passed to
                :class:`pygetm.output.netcdf.NetCDFFile`

        This is a wrapper around :meth:`add_netcdf_file` that automatically adds all
        arrays with the ``_part_of_state`` flag. By default, the restart file will be
        configured to be written at the end of the simulation, but this can be
        customized by providing argument ``interval``.
        """
        kwargs.setdefault("interval", -1)
        file = self.add_netcdf_file(path, **kwargs)
        for field in self.fields.values():
            if field.attrs.get("_part_of_state"):
                file.request(field)
        return file

    def start(
        self,
        itimestep: int = 0,
        time: Optional[cftime.datetime] = None,
        save: bool = True,
        default_time_reference: Optional[cftime.datetime] = None,
    ):
        for file in self.files:
            file.start(itimestep, time, save, default_time_reference or time)

    def save(
        self,
        seconds_passed: float,
        itimestep: int,
        time: Optional[cftime.datetime] = None,
        macro: bool = True,
    ):
        for file in self.files:
            file.save(seconds_passed, itimestep, time, macro)

    def close(self, seconds_passed: float, time: Optional[cftime.datetime] = None):
        for file in self.files:
            self._logger.debug("Closing %s" % file.path)
            file.close(seconds_passed, time)
