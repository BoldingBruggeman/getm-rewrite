import urllib.request
import tempfile
from typing import Optional
import shutil
import argparse

import xarray
import rioxarray

BASE_URL = "https://ows.emodnet-bathymetry.eu/wcs?service=wcs&version=1.0.0&"


def get(
    minlon: float,
    maxlon: float,
    minlat: float,
    maxlat: float,
    reslon: float = 1.0 / 60 / 16,
    reslat: float = 1.0 / 60 / 16,
    tiff_path: Optional[str] = None,
    verbose: bool = False,
) -> xarray.DataArray:
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
        f"&BBOX={minlon},{minlat},{maxlon},{maxlat}"
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
    da = rioxarray.open_rasterio(fout.name)
    # os.remove(fout.name)
    da["x"].attrs["units"] = "degrees_east"
    da["y"].attrs["units"] = "degrees_north"
    return da


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-a",
        "--area",
        nargs=4,
        type=float,
        help="bounding box: minlon, maxlon, minlat, maxlat",
        default=[-5.0, -3.0, 49.5, 51.5],
    )
    parser.add_argument(
        "-r",
        "--res",
        nargs=2,
        type=float,
        help="resolution in degrees: reslon, reslat",
        default=[1.0 / 60 / 16, 1.0 / 60 / 16],
    )
    parser.add_argument(
        "--tiff", type=str, help="TIFF file to save raw bathymetry to",
    )
    parser.add_argument("--plot", action="store_true", help="plot the bathymetry")
    parser.add_argument(
        "ncfile", nargs="?", type=str, help="NetCDF file to save bathymetry to",
    )
    parser.add_argument("--verbose", "-v", action="store_true", help="debug output")

    args = parser.parse_args()

    if args.verbose:
        print(f"area: {args.area}")
        print(f"resolution: {args.res}")

    # url = BASE_URL + "&request=GetCapabilities"
    # url = BASE_URL + "&request=DescribeCoverage&coverage=emodnet:mean"
    # with urllib.request.urlopen(url) as f, open("dump.xml", "wb") as fout:
    #     shutil.copyfileobj(f, fout)
    bathy = get(
        args.area[0],
        args.area[1],
        args.area[2],
        args.area[3],
        reslon=args.res[0],
        reslat=args.res[1],
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
