import urllib.request
import tempfile
import os
import shutil
import xarray
import rioxarray

BASE_URL = 'https://ows.emodnet-bathymetry.eu/wcs?service=wcs&version=1.0.0&'

def get(minlon: float, maxlon: float, minlat: float, maxlat: float, reslon: float=1. / 60 / 16, reslat: float=1. / 60 / 16) -> xarray.DataArray:
    assert maxlon > minlon, 'Maximum longitude %s must exceed minimum longitude %s' % (maxlon, minlon)
    assert maxlat > minlat, 'Maximum longitude %s must exceed minimum longitude %s' % (maxlat, minlat)
    url = BASE_URL + 'request=GetCoverage&coverage=emodnet:mean&crs=EPSG:4326&BBOX=%s,%s,%s,%s&format=GeoTIFF&interpolation=nearest&resx=%s&resy=%s' % (minlon, minlat, maxlon, maxlat, reslon, reslat)
#KB - to allow view the tiff file    with urllib.request.urlopen(url) as f, tempfile.NamedTemporaryFile(delete=False) as fout:
    with urllib.request.urlopen(url) as f, open('map.tif', 'wb') as fout:
        fout.write(f.read())
    da = rioxarray.open_rasterio(fout.name)
    #os.remove(fout.name)
    da['x'].attrs['units'] = 'degrees_east'
    da['y'].attrs['units'] = 'degrees_north'
    return da

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--area', nargs=4, type=float, help="bounding box: minlon, maxlon, minlat, maxlat", default=[-5.0, -3.0, 49.5, 51.5])
    parser.add_argument('-r', '--res', nargs=2, type=float, help="resolution: reslon, reslat", default=[1./60/16, 1./60/16])
    parser.add_argument('--tiff', action='store_true', help="display the raw TIFF image")
    parser.add_argument('--plot', action='store_true', help="plot the bathymetry")
    parser.add_argument('--ncfile', nargs="?", type=argparse.FileType('wb'), help="save NetCDF file")
    parser.add_argument('--verbose', action='store_true', help="be verbose about input")

    args=parser.parse_args()

    if args.verbose:
       print("area:")
       print(args.area)
       print("resolution:")
       print(args.res)

    url = BASE_URL + '&request=GetCapabilities'
    url = BASE_URL + '&request=DescribeCoverage&coverage=emodnet:mean'
    with urllib.request.urlopen(url) as f, open('dump.xml', 'wb') as fout:
        shutil.copyfileobj(f, fout)
    bathy = get(args.area[0], args.area[1], args.area[2], args.area[3], reslon=args.res[0], reslat=args.res[1])
    from matplotlib import pyplot
    if args.tiff:
        img = pyplot.imread('map.tif')
        pyplot.imshow(img[:, :, 0], cmap=pyplot.cm.coolwarm)
        pyplot.show()
    if args.plot:
        bathy.where(bathy <= 0).plot()
        pyplot.show()
    if args.ncfile:
        bathy.to_netcdf(args.ncfile)
