from typing import Iterable, Optional, Mapping, Union
import multiprocessing
import os
import argparse
import logging
from pathlib import Path
import contextlib


import yaml
import xarray as xr

try:
    import cdsapi
except ImportError:
    raise Exception("You need cdsapi. See https://cds.climate.copernicus.eu/how-to-api")

from . import util
from . import limit_region

VARIABLES = {
    "u10": "10m_u_component_of_wind",
    "v10": "10m_v_component_of_wind",
    "t2m": "2m_temperature",
    "d2m": "2m_dewpoint_temperature",
    "sp": "surface_pressure",
    "mslp": "mean_sea_level_pressure",
    "tcc": "total_cloud_cover",
    "tp": "total_precipitation",
    "e": "evaporation",
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
    area: list[float],
    variables: list[str],
    fmt: str,
    path: Path,
    complevel: int = 0,
    **cds_settings,
) -> Path:
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
    if fmt == "netcdf":
        tmpname = path.with_suffix(".tmp")
        with xr.open_dataset(path) as ds:
            util.configure_chunking_and_compression(ds, complevel=complevel)
            ds.to_netcdf(tmpname)
        tmpname.replace(path)
    return path


def get(
    minlon: float,
    maxlon: float,
    minlat: float,
    maxlat: float,
    start_year: int,
    stop_year: Optional[int] = None,
    variables: Iterable[str] = DEFAULT_VARIABLES,
    target_dir: Union[os.PathLike[str], str] = ".",
    logger: Optional[logging.Logger] = None,
    source: str = "cds",
    **kwargs,
) -> Mapping[tuple[int, str], Path]:
    logging.basicConfig(level=logging.INFO)
    logger = logger or logging.getLogger()

    if minlon >= maxlon:
        raise ValueError(
            f"Maximum longitude {maxlon} must exceed minimum longitude {minlon}"
        )
    if minlat >= maxlat:
        raise ValueError(
            f"Maximum latitude {maxlat} must exceed minimum latitude {minlat}"
        )
    if minlon < -360.0 or maxlon > 360.0:
        raise ValueError(
            f"Longitude range {minlon} - {maxlon} must fall within -360 - 360°"
        )
    if minlat < -90.0 or maxlat > 90.0:
        raise ValueError(
            f"Latitude range {minlat} - {maxlat} must fall within -90 - 90°"
        )

    minlon -= minlon % 0.25
    maxlon += -maxlon % 0.25
    minlat -= minlat % 0.25
    maxlat += -maxlat % 0.25
    minlon = max(-360.0, minlon)
    maxlon = min(360.0, maxlon)
    logger.info(f"Final area:")
    logger.info(f"  longitude: {minlon} - {maxlon} °East")
    logger.info(f"  latitude: {minlat} - {maxlat} °North")

    if stop_year is None:
        stop_year = start_year
    logger.info(f"Period: {start_year} - {stop_year}")
    logger.info(f"Selected variables: {', '.join(variables)}")

    area = [maxlat, minlon, minlat, maxlon]

    target_dir = Path(target_dir)
    target_dir.mkdir(parents=True, exist_ok=True)

    years = range(start_year, stop_year + 1)

    getter = {"cds": _get_cds, "arco": _get_arco}[source.lower()]
    return getter(
        variables=variables,
        years=years,
        area=area,
        target_dir=target_dir,
        logger=logger,
        **kwargs,
    )


def _get_arco(
    variables: Iterable[str],
    years: Iterable[str],
    area: list[float] = [90.0, -180.0, -90.0, 180.0],
    target_dir: Path = Path("."),
    logger: Optional[logging.Logger] = None,
    url: str = "gs://gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3",
    **kwargs,
):
    logger = logger or logging.getLogger()

    import zarr.storage
    import zarr.experimental.cache_store
    from dask.diagnostics import ProgressBar

    remote_store = zarr.storage.FsspecStore.from_url(
        url, read_only=True, storage_options=dict(token="anon")
    )
    cache_store = zarr.storage.MemoryStore()
    cached_store = zarr.experimental.cache_store.CacheStore(
        store=remote_store, cache_store=cache_store, max_size=1024**3
    )

    maxlat, minlon, minlat, maxlon = area
    results = {}
    with xr.open_zarr(cached_store, chunks=None) as ds:
        for year in years:
            ds_current = ds.sel(time=slice(f"{year}-01-01", f"{year}-12-31T23:30:00"))
            logger.info(f"  {year}:")
            for variable in variables:
                path = target_dir / f"era5_{variable}_{year}.nc"
                full_name = VARIABLES[variable]
                da = ds_current[full_name]
                da = limit_region(da, minlon, maxlon, minlat, maxlat, periodic_lon=True)
                shape = tuple(map(int, da.shape))
                logger.info(f"    {full_name}: {path} ({da.dtype}, shape {shape})")
                with ProgressBar():
                    da.chunk(time=24).rename(variable).to_netcdf(path)
                results[(year, variable)] = path

    # To avoid hang on Windows
    with contextlib.suppress(Exception):
        zarr.core.sync.cleanup_resources()

    return results


def _get_cds(
    variables: Iterable[str],
    years: Iterable[str],
    area: list[float] = [90.0, -180.0, -90.0, 180.0],
    target_dir: Path = Path("."),
    logger: Optional[logging.Logger] = None,
    fmt: str = "netcdf",
    cdsapirc: Union[os.PathLike, str, bytes, None] = None,
    **kwargs,
):
    logger = logger or logging.getLogger()
    if cdsapirc:
        with open(cdsapirc, "r") as f:
            kwargs.update(yaml.safe_load(f))

    pool = multiprocessing.Pool(processes=len(years) * len(variables))
    ext = "nc" if fmt == "netcdf" else "grib"
    tasks: list[tuple[int, str, multiprocessing.pool.AsyncResult]] = []
    for year in years:
        logger.info(f"  {year}:")
        for variable in variables:
            path = target_dir / f"era5_{variable}_{year}.{ext}"
            full_name = VARIABLES[variable]
            logger.info(f"    {full_name}: {path}")
            result = pool.apply_async(
                _download_year,
                args=(year, area, (full_name,), fmt, path),
                kwds=kwargs,
            )
            tasks.append((year, variable, result))
    return {(year, variable): res.get() for year, variable, res in tasks}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("minlon", help="minimum longitude (°East)", type=float)
    parser.add_argument("maxlon", help="maximum longitude (°East)", type=float)
    parser.add_argument("minlat", help="minimum latitude (°North)", type=float)
    parser.add_argument("maxlat", help="maximum latitude (°North)", type=float)
    parser.add_argument("start_year", help="first year to download", type=int)
    parser.add_argument(
        "stop_year", help="last year to download", type=int, nargs="?", default=None
    )
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
    parser.add_argument("--source", choices=["cds", "arco"], default="cds")
    parser.add_argument(
        "--complevel",
        help="NetCDF compression level (0-9). Default = 0 = no compression",
        type=int,
        default=0,
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
        complevel=args.complevel,
        source=args.source,
    )
