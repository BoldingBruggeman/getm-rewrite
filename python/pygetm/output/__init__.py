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

class File(operators.FieldCollection):
    def __init__(self, field_manager: FieldManager, interval: int=1, path: Optional[str]=None):
        super().__init__(field_manager)
        self.wait = interval
        self.interval = interval
        self.path = path

    def close(self):
        pass

    def save(self):
        self.wait -= 1
        if self.wait == 0:
            self.save_now()
            self.wait = self.interval

    def save_now(self):
        raise NotImplementedError

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
            self._logger.debug('Saving values to %s' % file.path)
            file.save()

    def close(self):
        for file in self.files:
            self._logger.debug('Closing %s' % file.path)
            file.close()
