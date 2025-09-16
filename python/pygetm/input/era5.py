from typing import Iterable, Optional, List, Mapping
import multiprocessing
import os
import argparse
import logging

import yaml

try:
    import cdsapi
except ImportError:
    raise Exception("You need cdsapi. See https://cds.climate.copernicus.eu/how-to-api")

VARIABLES = {
    "u10": "10m_u_component_of_wind",
    "v10": "10m_v_component_of_wind",
    "t2m": "2m_temperature",
    "d2m": "2m_dewpoint_temperature",
    "sp": "surface_pressure",
    "tcc": "total_cloud_cover",
    "tp": "total_precipitation",
    "ssr": "surface_net_solar_radiation",
    "ssrd": "surface_solar_radiation_downwards",
    "str": "surface_net_thermal_radiation",
    "strd": "surface_thermal_radiation_downwards",
    "tco3": "total_column_ozone",
    "tcwv": "total_column_water_vapour",
    "tclw": "total_column_cloud_liquid_water",
    "siconc": "sea_ice_cover",
}
DEFAULT_VARIABLES = ("u10", "v10", "t2m", "d2m", "sp", "tcc", "tp")


def _download_year(
    year: int,
    area: List[float],
    variables: List[str],
    fmt: str,
    path: str,
    **cds_settings,
):
    c = cdsapi.Client(verify=1, progress=False, **cds_settings)
    request = {
        "product_type": ["reanalysis"],
        "variable": variables,
        "year": [f"{year:04}"],
        "month": [f"{m:02}" for m in range(1, 13)],
        "day": [f"{d:02}" for d in range(1, 32)],
        "time": [f"{h:02}:00" for h in range(0, 24)],
        "grid": ["0.25/0.25"],
        "data_format": fmt,
        "download_format": "unarchived",
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
    stop_year: Optional[int] = None,
    variables: Iterable[str] = DEFAULT_VARIABLES,
    fmt: str = "netcdf",
    target_dir: str = ".",
    cdsapirc: Optional[str] = None,
    logger: Optional[logging.Logger] = None,
) -> Mapping[int, str]:
    logging.basicConfig(level=logging.INFO)
    logger = logger or logging.getLogger()

    assert minlon >= -360.0 and maxlon <= 360.0, (
        "Longitude must be between -360 and 360"
    )
    assert minlat >= -90.0 and maxlat <= 90.0, "Latitude must be between -360 and 360"

    minlon -= minlon % 0.25
    maxlon += -maxlon % 0.25
    minlat -= minlat % 0.25
    maxlat += -maxlat % 0.25
    minlon = max(-360.0, minlon)
    maxlon = min(360.0, maxlon)
    logger.info(
        f"Final area: longitude = {minlon} - {maxlon}, latitude = {minlat} - {maxlat}"
    )

    if stop_year is None:
        stop_year = start_year
    logger.info(f"Downloading {', '.join(variables)} for {start_year} - {stop_year}")

    area = [maxlat, minlon, minlat, maxlon]
    cds_settings = {}
    if cdsapirc:
        with open(cdsapirc, "r") as f:
            cds_settings.update(yaml.safe_load(f))

    results = []
    years = list(range(start_year, stop_year + 1))
    pool = multiprocessing.Pool(processes=len(years) * len(variables))
    os.makedirs(target_dir, exist_ok=True)
    ext = "nc" if fmt == "netcdf" else "grib"
    keys = []
    for year in years:
        logger.info(f"  {year}:")
        for variable in variables:
            path = os.path.join(target_dir, f"era5_{variable}_{year}.{ext}")
            full_name = VARIABLES[variable]
            logger.info(f"    {full_name}: {path}")
            results.append(
                pool.apply_async(
                    _download_year,
                    args=(year, area, (full_name,), fmt, path),
                    kwds=cds_settings,
                )
            )
            keys.append((year, variable))
    key2path = {}
    for key, res in zip(keys, results):
        key2path[key] = res.get()
    return key2path


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("minlon", help="minimum longitude (degrees East)", type=float)
    parser.add_argument("maxlon", help="maximum longitude (degrees East)", type=float)
    parser.add_argument("minlat", help="minimum latitude (degrees North)", type=float)
    parser.add_argument("maxlat", help="maximum latitude (degrees North)", type=float)
    parser.add_argument("start_year", help="start year", type=int)
    parser.add_argument("stop_year", help="stop year", type=int)
    parser.add_argument(
        "-v",
        help=f"extra variable to download ({', '.join(VARIABLES)})",
        action="append",
        dest="variables",
        metavar="VARIABLE",
        default=[],
    )
    parser.add_argument(
        "--no_default_variables",
        action="store_false",
        dest="default_variables",
        help=(
            f"do not include default variables ({', '.join(DEFAULT_VARIABLES)})"
            " unless explicitly specified"
        ),
    )
    parser.add_argument(
        "--cdsapirc",
        help=(
            "path to CDS configuration file"
            " (see https://cds.climate.copernicus.eu/how-to-api)"
        ),
    )
    parser.add_argument(
        "--grib", action="store_const", const="grib", dest="fmt", default="netcdf"
    )
    args = parser.parse_args()
    vars = set(args.variables)
    if args.default_variables:
        vars.update(DEFAULT_VARIABLES)

    get(
        args.minlon,
        args.maxlon,
        args.minlat,
        args.maxlat,
        args.start_year,
        args.stop_year,
        variables=vars,
        cdsapirc=args.cdsapirc,
        fmt=args.fmt,
    )
