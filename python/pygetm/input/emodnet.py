import urllib.request
import rioxarray

BASE_URL = 'https://ows.emodnet-bathymetry.eu/wcs?service=wcs&version=1.0.0&'

def get(minlon, maxlon, minlat, maxlat, resolution=1. / 16 / 60):
    assert maxlon > minlon, 'Maximum longitude %s must exceed minimum longitude %s' % (maxlon, minlon)
    assert maxlat > minlat, 'Maximum longitude %s must exceed minimum longitude %s' % (maxlat, minlat)
    url = BASE_URL + 'request=GetCoverage&coverage=emodnet:mean&crs=EPSG:4326&BBOX=%s,%s,%s,%s&format=GeoTIFF&interpolation=nearest&resx=%s&resy=%s' % (minlon, minlat, maxlon, maxlat, resolution, resolution)
    with urllib.request.urlopen(url) as f, open('map.tif', 'wb') as fout:
        fout.write(f.read())
    da = rioxarray.open_rasterio('map.tif')
    da['x'].attrs['units'] = 'degrees_east'
    da['y'].attrs['units'] = 'degrees_north'
    return da

if __name__ == '__main__':
    url = BASE_URL + '&request=GetCapabilities'
    url = BASE_URL + '&request=DescribeCoverage&coverage=emodnet:mean'
    with urllib.request.urlopen(url) as f, open('dump.xml', 'wb') as fout:
        fout.write(f.read())
    bathy = get(-5, -3, 49.5, 51.5)
    from matplotlib import pyplot
    bathy.plot()
    pyplot.show()