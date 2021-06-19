from typing import  MutableSequence, Optional
import os.path

import numpy.typing
import netCDF4

from .. import core
from . import FieldManager, File
from .. import _pygetm

class NetCDFFile(File):
    def __init__(self, field_manager: FieldManager, path: str, rank: int, interval: int=1, sub: bool=False):
        File.__init__(self, field_manager)
        name, ext = os.path.splitext(path)
        if sub:
            name += '_%05i' % rank
        self.path = name + ext
        self.nc = None
        self.ncvars: MutableSequence[Optional[netCDF4.Variable]] = []
        self.itime = 0
        self.is_root = rank == 0
        self.interval = interval
        self.wait = 0
        self.sub = sub
        self.created = False

    def request(self, name: str, output_name: Optional[str]=None, dtype: Optional[numpy.typing.DTypeLike]=None):
        File.request(self, name, output_name, dtype)
        self.ncvars.append(None)

    def _create(self):
        self.created = True
        if not (self.is_root or self.sub):
            return
        self.nc = netCDF4.Dataset(self.path, 'w')
        for grid in frozenset(array.grid for array in self.fields):
            nx, ny = grid.nx, grid.ny
            if not self.sub:
                nx, ny = (nx - 2 * grid.halo) * grid.domain.tiling.ncol, (ny - 2 * grid.halo) * grid.domain.tiling.nrow
            self.nc.createDimension('x%s' % grid.postfix, nx)
            self.nc.createDimension('y%s' % grid.postfix, ny)
        self.nc.createDimension('time',)
        self.ncvars = []
        for output_name, array in zip(self.order, self.fields):
            dims = ('y%s' % array.grid.postfix, 'x%s' % array.grid.postfix)
            if array.ndim == 3: dims = ('z',) + dims
            dims = ('time',) + dims
            ncvar = self.nc.createVariable(output_name, array.dtype, dims)
            atts = {}
            if array.units is not None:
                atts['units'] = array.units
            if array.long_name is not None:
                atts['long_name'] = array.long_name
            for att, value in atts.items():
                setattr(ncvar, att, value)
            self.ncvars.append(ncvar)

    def save(self):
        if self.wait == 0:
            if not self.created:
                self._create()
            for array, ncvar in zip(self.fields, self.ncvars):
                if self.sub:
                    ncvar[self.itime, ...] = array.all_values
                else:
                    array.gather(ncvar, slice_spec=(self.itime,))
            if self.nc is not None:
                self.nc.sync()
            self.itime += 1
            self.wait = self.interval
        self.wait -= 1

    def close(self):
        if self.nc is not None:
            self.nc.close()