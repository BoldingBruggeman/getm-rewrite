import urllib.request
import tempfile
from typing import Optional
import argparse

import xarray as xr

try:
    import rioxarray
except ImportError:
    raise Exception("You need rioxarray. See https://corteva.github.io/rioxarray")

BASE_URL = "https://ows.emodnet-bathymetry.eu/wcs?service=wcs&version=1.0.0&"
DEFAULT_RESOLUTION = 1.0 / 60 / 16


def get(
    minlon: float,
    maxlon: float,
    minlat: float,
    maxlat: float,
    reslon: float = DEFAULT_RESOLUTION,
    reslat: float = DEFAULT_RESOLUTION,
    tiff_path: Optional[str] = None,
    verbose: bool = False,
) -> xr.DataArray:
    if maxlon < minlon:
        raise Exception(
            f"Maximum longitude {maxlon} must exceed minimum longitude {minlon}"
        )
    if maxlat <= minlat:
        raise Exception(
            f"Maximum latitude {maxlat} must exceed minimum latitude {minlat}"
        )
    url = BASE_URL + (
        f"request=GetCoverage&coverage=emodnet:mean&crs=EPSG:4326"
        f"&BBOX={minlon-reslon},{minlat-reslat},{maxlon+reslon},{maxlat+reslat}"
        f"&format=GeoTIFF&interpolation=nearest"
        f"&resx={reslon}&resy={reslat}"
    )
    if tiff_path is None:
        tiff = tempfile.NamedTemporaryFile(delete=False)
    else:
        tiff = open(tiff_path, "wb")
    if verbose:
        print("Downloading TIFF...")
    with urllib.request.urlopen(url) as f, tiff as fout:
        fout.write(f.read())
    if verbose:
        print("Opening TIFF as xarray...")
    da = rioxarray.open_rasterio(fout.name, default_name="bath")
    # os.remove(fout.name)
    da["x"].attrs["units"] = "degrees_east"
    da["y"].attrs["units"] = "degrees_north"
    if da.ndim == 3 and da.shape[0] == 1:
        # Squeeze out band dimension
        da = da[0, :, :]
    return da


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("minlon", help="minimum longitude (degrees East)", type=float)
    parser.add_argument("maxlon", help="maximum longitude (degrees East)", type=float)
    parser.add_argument("minlat", help="minimum latitude (degrees North)", type=float)
    parser.add_argument("maxlat", help="maximum latitude (degrees North)", type=float)
    parser.add_argument("ncfile", nargs="?", help="NetCDF file to save bathymetry to")
    parser.add_argument(
        "--resolution",
        "-r",
        type=float,
        help="resolution (degrees)",
        default=DEFAULT_RESOLUTION,
    )
    parser.add_argument("--tiff", help="TIFF file to save raw bathymetry to")
    parser.add_argument("--plot", action="store_true", help="plot the bathymetry")
    parser.add_argument("--verbose", "-v", action="store_true", help="debug output")
    args = parser.parse_args()

    if args.verbose:
        print(f"Longitude: {args.minlon} - {args.maxlon} degrees East")
        print(f"Latitude: {args.minlat} - {args.maxlat} degrees North")
        print(f"Resolution: {args.resolution} degrees")

    # url = BASE_URL + "&request=GetCapabilities"
    # url = BASE_URL + "&request=DescribeCoverage&coverage=emodnet:mean"
    # with urllib.request.urlopen(url) as f, open("dump.xml", "wb") as fout:
    #     shutil.copyfileobj(f, fout)
    bathy = get(
        args.minlon,
        args.maxlon,
        args.minlat,
        args.maxlat,
        reslon=args.resolution,
        reslat=args.resolution,
        tiff_path=args.tiff,
        verbose=args.verbose,
    )
    from matplotlib import pyplot

    # if args.tiff:
    #     img = pyplot.imread("map.tif")
    #     pyplot.imshow(img[:, :, 0], cmap=pyplot.cm.coolwarm)
    #     pyplot.show()
    if args.plot:
        bathy.where(bathy <= 0).plot()
        pyplot.show()
    if args.ncfile is not None:
        bathy.to_netcdf(args.ncfile)
