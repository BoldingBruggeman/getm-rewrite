import os.path
from typing import Tuple, Mapping

import numpy
import cftime
import xarray
import otps2

import pygetm.input

ROOT = '../../../igotm/data/TPXO9'

def get(lon, lat, variable='h', verbose: bool=False, root=ROOT, scale_factor: float=1.) -> xarray.DataArray:
    assert variable in ('h', 'u', 'v', 'hz', 'hu', 'hv')

    lon = numpy.asarray(lon)
    lat = numpy.asarray(lat)
    def select(ncvar) -> xarray.DataArray:
        out = pygetm.input.limit_region(ncvar, lon.min(), lon.max(), lat.min(), lat.max(), periodic_lon=True)
        out = pygetm.input.horizontal_interpolation(out, lon, lat)
        return out

    if variable in ('hz', 'hu', 'hv'):
        axis = variable[1]
        with xarray.open_dataset(os.path.join(root, 'grid_tpxo9_atlas.nc')) as ds:
            ds = ds.set_coords(('lat_%s' % axis, 'lon_%s' % axis))
            return select(ds[variable])

    scale_factor *= {'h': 1e-3, 'u': 1e-4, 'v': 1e-4}.get(variable, 1.)
    axis = {'h': 'z'}.get(variable, variable)
    file_prefix = {'v': 'u'}.get(variable, variable)
    components: Mapping[str, Tuple[numpy.ndarray, numpy.ndarray]] = {}
    for component in ('m2', 's2', 'n2', 'k2', 'k1', 'o1', 'p1', 'q1', 'm4', 'ms4', 'mn4', '2n2'):
        if verbose:
            print('TPXO: reading %s constituent of %s...' % (component, variable))
        name = '%s_%s_tpxo9_atlas_30.nc' % (file_prefix, component)
        path = os.path.join(root, name)
        with xarray.open_dataset(path) as ds:
            ds = ds.set_coords(('lat_%s' % axis, 'lon_%s' % axis))
            x = select(ds['%sIm' % variable])
            components[component] = scale_factor * select(ds['%sRe' % variable]).values, scale_factor * x.values
    return xarray.DataArray(Data(components, lat), dims=x.dims, coords=x.coords, name='tpxo(%s, %s)' % (root, variable))

class Data(pygetm.input.LazyArray):
    def __init__(self, components: Mapping[str, Tuple[numpy.ndarray, numpy.ndarray]], lat: numpy.ndarray):
        pygetm.input.LazyArray.__init__(self, lat.shape, numpy.float64)
        self.components = components
        self.lat = lat
        self.time = None

    def update(self, time: cftime.datetime) -> bool:
        self.time = time
        return True

    def __array__(self, dtype=None) -> numpy.ndarray:
        assert self.time is not None, 'update has not yet been called'
        return otps2.predict_tide_2d(self.components, self.lat, self.time, ntime=1, delta_time=0)[0, ...]

    def is_time_varying(self) -> bool:
        return True