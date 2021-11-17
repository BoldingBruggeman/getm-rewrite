from typing import  MutableSequence, Optional
import os.path

import numpy.typing
import netCDF4

from .. import core
from . import FieldManager, File
from .. import _pygetm

class NetCDFFile(File):
    def __init__(self, field_manager: FieldManager, path: str, rank: int, sub: bool=False, **kwargs):
        super().__init__(field_manager, **kwargs)
        name, ext = os.path.splitext(path)
        if sub:
            name += '_%05i' % rank
        self.path = name + ext
        self.nc = None
        self.itime = 0
        self.is_root = rank == 0
        self.sub = sub
        self.created = False

    def _create(self):
        self.created = True
        if not (self.is_root or self.sub):
            return
        self.nc = netCDF4.Dataset(self.path, 'w')
        for grid in frozenset(field.grid for field in self.fields.values()):
            if self.sub:
                # Output data (including halos) for the current subdomain
                nx, ny, nz = grid.nx_, grid.ny_, grid.nz_
            else:
                # Output data for entire (global) domain
                nx, ny, nz = grid.nx * grid.domain.tiling.ncol, grid.ny * grid.domain.tiling.nrow, grid.nz
            xname, yname, zname = 'x%s' % grid.xypostfix, 'y%s' % grid.xypostfix, 'z%s' % grid.zpostfix
            if (xname not in self.nc.dimensions): self.nc.createDimension(xname, nx)
            if (yname not in self.nc.dimensions): self.nc.createDimension(yname, ny)
            if (zname not in self.nc.dimensions): self.nc.createDimension(zname, nz)
        self.nc.createDimension('time',)
        self.ncvars = []
        for output_name, field in self.fields.items():
            dims = ('y%s' % field.grid.xypostfix, 'x%s' % field.grid.xypostfix)
            if field.ndim == 3:
                dims = ('z%s' % field.grid.zpostfix,) + dims
            dims = ('time',) + dims
            ncvar = self.nc.createVariable(output_name, field.dtype, dims, fill_value=field.fill_value)
            for att, value in field.atts.items():
                setattr(ncvar, att, value)
            field.ncvar = ncvar

    def save_now(self):
        if not self.created:
            self._create()
        for field in self.fields.values():
            field.get(getattr(field, 'ncvar', None), slice_spec=(self.itime,), sub=self.sub)
        if self.nc is not None:
            self.nc.sync()
        self.itime += 1

    def close(self):
        if self.nc is not None:
            self.nc.close()