from typing import Iterable, Optional
import multiprocessing
import os
import argparse

import yaml

try:
    import cdsapi
except ImportError:
    raise Exception("You need cdsapi. See https://cds.climate.copernicus.eu/api-how-to")

VARIABLES = {
    "u10": "10m_u_component_of_wind",
    "v10": "10m_v_component_of_wind",
    "t2m": "2m_temperature",
    "d2m": "2m_dewpoint_temperature",
    "sp": "surface_pressure",
    "tcc": "total_cloud_cover",
    "tp": "total_precipitation",
    "ssr": "surface_net_solar_radiation",
}
DEFAULT_VARIABLES = ("u10", "v10", "t2m", "d2m", "sp", "tcc", "tp")


def _download_year(
    year: int, area: list[float], variable: str, path: str, **cds_settings
):
    c = cdsapi.Client(verify=1, progress=False, **cds_settings)
    request = {
        "variable": [variable],
        "product_type": "reanalysis",
        "format": "netcdf",
        "year": "%04i" % year,
        "month": ["%02i" % m for m in range(1, 13)],
        "day": ["%02i" % d for d in range(1, 32)],
        "time": ["%02i:00" % h for h in range(0, 24)],
        "grid": ["0.25/0.25"],
        "area": area,
    }
    r = c.retrieve("reanalysis-era5-single-levels", request)
    r.download(path)
    return path


def get(
    minlon: float,
    maxlon: float,
    minlat: float,
    maxlat: float,
    start_year: int,
    stop_year: int,
    selected_variables: Iterable[str] = DEFAULT_VARIABLES,
    target_dir: str = ".",
    cdsapirc: Optional[str] = None,
):
    area = [maxlat, minlon, minlat, maxlon]
    cds_settings = {}
    if cdsapirc:
        with open(cdsapirc, "r") as f:
            cds_settings.update(yaml.safe_load(f))

    pool = multiprocessing.Pool(
        processes=(stop_year - start_year + 1) * len(selected_variables)
    )
    results = []
    for key in selected_variables:
        for year in range(start_year, stop_year + 1):
            path = os.path.join(target_dir, f"era5_{key}_{year}.nc")
            print("Queuing download of %s..." % path)
            results.append(
                pool.apply_async(
                    _download_year,
                    args=(year, area, VARIABLES[key], path),
                    kwds=cds_settings,
                )
            )
    for res in results:
        path = res.get()
        print("%s done." % path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("minlon", help="minimum longitude (degrees East)", type=float)
    parser.add_argument("maxlon", help="maximum longitude (degrees East)", type=float)
    parser.add_argument("minlat", help="minimum latitude (degrees North)", type=float)
    parser.add_argument("maxlat", help="maximum latitude (degrees North)", type=float)
    parser.add_argument("start_year", help="start year", type=int)
    parser.add_argument("stop_year", help="stop year", type=int)
    parser.add_argument(
        "--cdsapirc",
        help="path to CDS configuration file (see https://cds.climate.copernicus.eu/api-how-to)",
    )
    args = parser.parse_args()
    get(
        args.minlon,
        args.maxlon,
        args.minlat,
        args.maxlat,
        args.start_year,
        args.stop_year,
        cdsapirc=args.cdsapirc,
    )
