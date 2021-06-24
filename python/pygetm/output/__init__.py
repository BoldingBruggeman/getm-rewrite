from typing import MutableMapping, MutableSequence, Optional

import numpy.typing

from .. import core

class FieldManager:
    def __init__(self):
        self.fields: MutableMapping[str, core.Array] = {}

    def register(self, array: core.Array):
        assert array.name is not None, 'Cannot register field without name.'
        assert array.name not in self.fields, 'A field with name "%s" has already been registered.' % array.name
        self.fields[array.name] = array

class File:
    def __init__(self, field_manager: FieldManager):
        self.field_manager = field_manager
        self.order: MutableSequence[str] = []
        self.fields: MutableSequence[core.Array] = []

    def close(self):
        pass

    def request(self, name: str, output_name: Optional[str]=None, dtype: Optional[numpy.typing.DTypeLike]=None):
        assert name in self.field_manager.fields, 'Unknown field "%s" requested. Available: %s' % (name, ', '.join(self.field_manager.fields))
        if output_name is None:
            output_name = name
        assert output_name not in self.order, 'A variable with name "%s" has already been added to %s.' % (output_name, self)
        array = self.field_manager.fields[name]
        dtype = dtype or array.dtype
        self.order.append(output_name)
        self.fields.append(array)

class OutputManager(FieldManager):
    def __init__(self, rank: int):
        FieldManager.__init__(self)
        self.rank = rank
        self.files = []

    def add_netcdf_file(self, path: str, **kwargs):
        from . import netcdf
        file = netcdf.NetCDFFile(self, path, rank=self.rank, **kwargs)
        self.files.append(file)
        return file

    def save(self):
        for file in self.files:
            file.save()

    def close(self):
        for file in self.files:
            file.close()
