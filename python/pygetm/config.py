from typing import Any, Mapping, Iterator, Optional
import collections.abc
import logging

import yaml

from . import input
from . import domain

class Node(collections.abc.Mapping):
    def __init__(self, dictionary: Mapping[str, Any], prefix: str=''):
        self.prefix = prefix
        assert isinstance(dictionary, Mapping)
        self.dictionary = {}
        for name, value in dictionary.items():
            if isinstance(value, Mapping):
                value = Node(value, prefix = '%s%s/' % (self.prefix, name))
            self.dictionary[name] = value
        self.retrieved = set()
    
    def __getitem__(self, path: str):
        components = path.split('/', 1)
        value = self.dictionary[components[0]]
        self.retrieved.add(components[0])
        if len(components) > 1:
            assert isinstance(value, Node)
            value = value[components[1]]
        return value

    def __iter__(self) -> Iterator:
        return self.dictionary.__iter__()

    def __len__(self) -> int:
        return self.dictionary.__len__()

    def check(self):
        unused = []
        for name, value in self.dictionary.items():
            if name not in self.retrieved:
                unused.append('%s%s' % (self.prefix, name))
            if isinstance(value, Node):
                unused += value.check()
        return unused

    def get_input(self, name: str, domain: domain.Domain, is_3d: bool=False, default: Optional[float]=None, logger: Optional[logging.Logger]=None):
        info = self.get(name)
        if info is None:
            if default is None:
                raise Exception('%s%s is required' % (self.prefix, name))
            info = default
        if isinstance(info, Node):
            preprocess = None
            time_slice = info.get('time_slice')
            if time_slice is not None:
                start, stop = time_slice.split(':')
                preprocess = lambda ds: ds.isel(time=slice(int(start), int(stop)))
            path, variable = info['path'], info['variable']
            if logger is not None:
                logger.info('%s%s = variable %s from %s' % (self.prefix, name, variable, path))
            return domain.T.map(input.get_from_nc(path, variable, preprocess=preprocess), name=name)
        else:
            if logger is not None:
                logger.info('%s%s = %s' % (self.prefix, name, info))
            return domain.T.array(fill=info, is_3d=is_3d, name=name)        

def configure(path: str):
    with open(path) as f:
        settings = yaml.safe_load(f)
    if not isinstance(settings, Mapping):
        raise Exception('%s should contain a mapping with configuration information, but instead contains %s' % (path, settings))
    return Node(settings)
