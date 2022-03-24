from typing import  MutableSequence, Optional
import os.path

import numpy.typing
import netCDF4

from .. import core
from . import FieldManager, File
from .. import _pygetm
from ..constants import INTERFACES

class NetCDFFile(File):
    def __init__(self, field_manager: FieldManager, path: str, rank: int, sub: bool=False, sync_interval: int=1, **kwargs):
        super().__init__(field_manager, path=path, **kwargs)
        name, ext = os.path.splitext(path)
        if sub:
            name += '_%05i' % rank
        self.path = name + ext
        self.nc = None
        self.itime = 0
        self.is_root = rank == 0
        self.sub = sub
        self.created = False
        self.sync_interval = sync_interval

    def _create(self):
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
                field.ncvar = ncvar

        # Save time-invariant fields
        for field in self.fields.values():
            if field.constant:
                field.get(getattr(field, 'ncvar', None), sub=self.sub)

    def save_now(self):
        if not self.created:
            self._create()
        if self.nc is not None:
            self.nctime[self.itime] = self.itime
        for field in self.fields.values():
            if not field.constant:
                field.get(getattr(field, 'ncvar', None), slice_spec=(self.itime,), sub=self.sub)
        self.itime += 1
        if self.nc is not None and self.itime % self.sync_interval == 0:
            self.nc.sync()

    def close(self):
        if self.nc is not None:
            self.nc.close()