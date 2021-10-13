import os.path
import xarray
import otps2

import pygetm.input

ROOT = '../../../igotm/data/TPXO9'

class Dataset:
    def __init__(self, lon, lat, variable='h', verbose: bool=False, root=ROOT):
        assert variable in ('h', 'u', 'u')
        self.lat = lat
        self.components = {}
        for component in ('m2', 's2', 'n2', 'k2', 'k1', 'o1', 'p1', 'q1', 'm4', 'ms4', 'mn4', '2n2'):
            if verbose:
                print('TPXO: reading %s constituent of %s...' % (component, variable))
            name = '%s_%s_tpxo9_atlas_30.nc' % ('h' if variable == 'h' else 'u', component)
            path = os.path.join(root, name)
            with xarray.open_dataset(path) as ds:
                ds = ds.set_coords(('lat_z', 'lon_z'))
                re = pygetm.input.Variable(ds['%sRe' % variable])
                im = pygetm.input.Variable(ds['%sIm' % variable])
                re = pygetm.input.LimitRegion(re, lon.min(), lon.max(), lat.min(), lat.max(), periodic_lon=True)
                im = pygetm.input.LimitRegion(im, lon.min(), lon.max(), lat.min(), lat.max(), periodic_lon=True)
                re = pygetm.input.SpatialInterpolation(re, lon, lat)
                im = pygetm.input.SpatialInterpolation(im, lon, lat)
                self.components[component] = re.x.values, im.x.values

    def update(self, date):
        return otps2.predict_tide_2d(self.components, self.lat, date, ntime=1, delta_time=0)[0, ...]
