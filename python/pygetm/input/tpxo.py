import os.path
import xarray
import otps2

ROOT = '../../../igotm/data/TPXO9'

class Dataset:
    def __init__(self, lon, lat, variable='h', verbose: bool=False, root=ROOT):
        assert variable in ('h', 'u', 'u')
        self.lon = lon % 360.
        self.lat = lat
        self.components = {}
        for component in ('m2', 's2', 'n2', 'k2', 'k1', 'o1', 'p1', 'q1', 'm4', 'ms4', 'mn4', '2n2'):
            if verbose:
                print('TPXO: reading %s constituent of %s...' % (component, variable))
            name = '%s_%s_tpxo9_atlas_30.nc' % ('h' if variable == 'h' else 'u', component)
            path = os.path.join(root, name)
            with xarray.open_dataset(path) as ds:
                ds = ds.set_coords(('lat_z', 'lon_z'))
                re = ds['%sRe' % variable]
                im = ds['%sIm' % variable]
                self.components[component] = re.getm.interp(self.lon, self.lat).values, im.getm.interp(self.lon, self.lat).values

    def update(self, date):
        return otps2.predict_tide_2d(self.components, self.lat, date, ntime=1, delta_time=0)[0, ...]
