from typing import Optional
import os.path
import logging

import numpy.typing
import cftime
import netCDF4

from .. import core
from . import FieldManager, File
from .. import _pygetm
from ..constants import INTERFACES

class NetCDFFile(File):
    def __init__(self, field_manager: FieldManager, logger: logging.Logger, path: str, rank: int, sub: bool=False, sync_interval: Optional[int]=1, time_reference: Optional[cftime.datetime]=None, **kwargs):
        super().__init__(field_manager, logger, path=path, **kwargs)
        name, ext = os.path.splitext(path)
        if sub:
            name += '_%05i' % rank
        self.path = name + ext
        self.nc = None
        self.itime = 0
        self.is_root = rank == 0
        self.sub = sub
        self.created = False
        self.time_offset = 0.
        self.time_reference = time_reference
        self.sync_interval = sync_interval

    def __repr__(self) -> str:
        return 'NetCDFFile(\'%s\')' % self.path

    def _create(self, seconds_passed: float, time: Optional[cftime.datetime]):
        self.created = True

        if self.is_root or self.sub:
            # Create the NetCDF file
            self.nc = netCDF4.Dataset(self.path, 'w')
            for grid in frozenset(field.grid for field in self.fields.values()):
                if self.sub:
                    # Output data (including halos) for the current subdomain
                    nx, ny, nz = grid.nx_, grid.ny_, grid.nz_
                else:
                    # Output data for entire (global) domain
                    nx = grid.domain.tiling.nx_glob + grid.nx - grid.domain.T.nx
                    ny = grid.domain.tiling.ny_glob + grid.ny - grid.domain.T.ny
                    nz = grid.nz
                xname, yname, zname, ziname = 'x%s' % grid.postfix, 'y%s' % grid.postfix, 'z', 'zi'
                if (xname not in self.nc.dimensions): self.nc.createDimension(xname, nx)
                if (yname not in self.nc.dimensions): self.nc.createDimension(yname, ny)
                if (zname not in self.nc.dimensions): self.nc.createDimension(zname, nz)
                if (ziname not in self.nc.dimensions): self.nc.createDimension(ziname, nz + 1)
            self.nc.createDimension('time',)
            self.nctime = self.nc.createVariable('time', float, ('time',))
            if time is not None:
                self.nctime.units = 'seconds since %s' % self.time_reference.strftime('%Y-%m-%d %H:%M:%S')
                self.nctime.calendar = time.calendar
                self.time_offset = (time - self.time_reference).total_seconds() - seconds_passed
            else:
                self.nctime.units = 's'
            self.ncvars = []
            for output_name, field in self.fields.items():
                dims = ('y%s' % field.grid.postfix, 'x%s' % field.grid.postfix)
                if field.z:
                    dims = ('zi' if field.z == INTERFACES else 'z',) + dims
                if not field.constant:
                    dims = ('time',) + dims
                ncvar = self.nc.createVariable(output_name, field.dtype, dims, fill_value=field.fill_value)
                for att, value in field.atts.items():
                    setattr(ncvar, att, value)
                if field.coordinates:
                    setattr(ncvar, 'coordinates', ' '.join(field.coordinates))
                field._ncvar = ncvar

        # Save time-invariant fields
        for field in self.fields.values():
            if field.constant:
                field.get(getattr(field, '_ncvar', None), sub=self.sub)

    def start(self, itimestep: int, time: Optional[cftime.datetime], save: bool):
        if self.time_reference is None:
            self.time_reference = time
        super().start(itimestep, time, save)

    def save_now(self, seconds_passed: float, time: Optional[cftime.datetime]):
        if not self.created:
            self._create(seconds_passed, time)
        if self.nc is not None:
            self.nctime[self.itime] = self.time_offset + seconds_passed
        for field in self.fields.values():
            if not field.constant:
                field.get(getattr(field, '_ncvar', None), slice_spec=(self.itime,), sub=self.sub)
        self.itime += 1
        if self.nc is not None and self.sync_interval is not None and self.itime % self.sync_interval == 0:
            self.nc.sync()

    def close_now(self, seconds_passed: float, time: Optional[cftime.datetime]):
        if self.nc is not None:
            self.nc.close()