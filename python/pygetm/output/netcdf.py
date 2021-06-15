from typing import  MutableSequence, Optional

import numpy.typing
import netCDF4

from .. import core
from . import FieldManager, File
from .. import _pygetm

class NetCDFFile(File):
    def __init__(self, field_manager: FieldManager, path: str, is_root: bool, interval: int=1):
        File.__init__(self, field_manager)
        self.path = path
        self.nc = None
        self.ncvars: MutableSequence[Optional[netCDF4.Variable]] = []
        self.itime = 0
        self.is_root = is_root
        self.interval = interval
        self.wait = 0

    def request(self, name: str, output_name: Optional[str]=None, dtype: Optional[numpy.typing.DTypeLike]=None):
        File.request(self, name, output_name, dtype)
        self.ncvars.append(None)

    def _create(self):
        self.nc = netCDF4.Dataset(self.path, 'w')
        postfixes = {_pygetm.TGRID: 't', _pygetm.UGRID: 'u', _pygetm.VGRID: 'v', _pygetm.XGRID: 'x'}
        for grid in frozenset(array.grid for array in self.fields):
            nx, ny = (grid.nx - 2 * grid.halo) * grid.domain.tiling.ncol, (grid.ny - 2 * grid.halo) * grid.domain.tiling.nrow
            postfix = postfixes[grid.type]
            self.nc.createDimension('x%s' % postfix, nx)
            self.nc.createDimension('y%s' % postfix, ny)
            self.nc.createDimension('time',)
        self.ncvars = []
        for output_name, array in zip(self.order, self.fields):
            postfix = postfixes[array.grid.type]
            dims = ('y%s' % postfix, 'x%s' % postfix)
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
            if self.is_root and self.nc is None:
                self._create()
            for array, ncvar in zip(self.fields, self.ncvars):
                array.gather(ncvar, slice_spec=(self.itime,))
            if self.nc is not None:
                self.nc.sync()
            self.itime += 1
            self.wait = self.interval
        self.wait -= 1
