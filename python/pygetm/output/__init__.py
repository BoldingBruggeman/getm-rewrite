import logging
from typing import MutableMapping, Optional, Union, Iterable, List
import collections

import numpy.typing

from .. import core
from . import operators

class FieldManager:
    def __init__(self):
        self.fields: MutableMapping[str, core.Array] = {}

    def register(self, array: core.Array):
        assert array.name is not None, 'Cannot register field without name.'
        assert array.name not in self.fields, 'A field with name "%s" has already been registered.' % array.name
        self.fields[array.name] = array

class File:
    def __init__(self, field_manager: FieldManager, interval: int=1):
        self.field_manager = field_manager
        self.fields: MutableMapping[str, operators.Base] = collections.OrderedDict()
        self.wait = interval
        self.interval = interval

    def close(self):
        pass

    def save(self):
        self.wait -= 1
        if self.wait == 0:
            self.save_now()
            self.wait = self.interval

    def save_now(self):
        raise NotImplementedError

    def request(self, name: Union[str, Iterable[str]], output_name: Optional[str]=None, dtype: Optional[numpy.typing.DTypeLike]=None, mask: bool=False):
        if not isinstance(name, str):
            assert output_name is None
            for n in name:
                self.request(n, dtype=dtype, mask=mask)
            return
        assert name in self.field_manager.fields, 'Unknown field "%s" requested. Available: %s' % (name, ', '.join(self.field_manager.fields))
        if output_name is None:
            output_name = name
        assert output_name not in self.fields, 'A variable with name "%s" has already been added to %s.' % (output_name, self)
        array = self.field_manager.fields[name]
        array.saved = True
        dtype = dtype or array.dtype
        field = operators.Field(array)
        if mask:
            field = operators.Mask(field)
        self.fields[output_name] = field

class OutputManager(FieldManager):
    def __init__(self, rank: int, logger: Optional[logging.Logger]=None):
        FieldManager.__init__(self)
        self.rank = rank
        self.files: List[File] = []
        self._logger = logger or logging.getLogger()

    def add_netcdf_file(self, path: str, **kwargs):
        self._logger.debug('Adding NetCDF file %s' % path)
        from . import netcdf
        file = netcdf.NetCDFFile(self, path, rank=self.rank, **kwargs)
        self.files.append(file)
        return file

    def start(self, save: bool=True):
        if save:
            self._logger.info('Saving initial state')
            for file in self.files:
                file.save_now()

    def save(self):
        for file in self.files:
            self._logger.debug('Saving values to %s' % file)
            file.save()

    def close(self):
        for file in self.files:
            self._logger.debug('Closing %s' % file)
            file.close()
