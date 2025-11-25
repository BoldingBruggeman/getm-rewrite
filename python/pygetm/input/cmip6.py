import argparse
import logging
import os
from typing import Optional, Union
import tempfile
import shutil
from pathlib import Path
from collections.abc import Iterable, Iterator
import contextlib

import xarray as xr

# import numpy as np
# import pandas as pd
# import gcsfs
# import zarr  # indirectly needed for xarray to open zarr stores

try:
    import intake_esgf
except ImportError:
    raise Exception("You need intake_esgf. See https://intake-esgf.readthedocs.io/")


import pygetm.input
import pygetm.input.util

# url for the CSV file that contains the data catalog
# URL_CATALOG = "https://storage.googleapis.com/cmip6/cmip6-zarr-consolidated-stores.csv"


METEO_VARS_INSTANT = (
    "tas",  # Near-Surface Air Temperature
    "uas",  # Eastward Near-Surface Wind
    "vas",  # Northward Near-Surface Wind
    "huss",  # Near-Surface Specific Humidity
    "ps",  # Surface Air Pressure
)

METEO_VARS_ACCUM = (
    "rlds",  # Surface Downwelling Longwave Radiation
    "rsds",  # Surface Downwelling Shortwave Radiation
    "rlus",  # Surface Upwelling Longwave Radiation
    "rsus",  # Surface Upwelling Shortwave Radiation
    "pr",  # Precipitation
)

METEO_VARS = METEO_VARS_INSTANT + METEO_VARS_ACCUM


def list_available_models():
    intake_esgf.conf.set(all_indices=True)
    cat = intake_esgf.ESGFCatalog()
    cat = cat.search(
        # experiment_id="ssp245",
        variable_id=METEO_VARS[0],
        frequency="3hr",
    )
    print(f"Available experiment_ids: {', '.join(cat.df.experiment_id.unique())}")
    print(f"Available source_ids: {', '.join(cat.df.source_id.unique())}")


def _create_logger() -> logging.Logger:
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    return logger


@contextlib.contextmanager
def get_global_meteo(
    source_id: str = "GFDL-ESM4",
    experiment_id: str = "ssp245",
    logger: Optional[logging.Logger] = None,
    cache_dir: Union[os.PathLike[str], str, None] = None,
    variables: Iterable[str] = METEO_VARS,
    prefer_streaming: bool = True,
) -> Iterator[xr.Dataset]:
    logger = logger or _create_logger()

    # # open the data catalog with pandas, and take a peek at how it's formatted
    # logger.info(f"Reading catalog from {URL_CATALOG}")
    # df_catalog = pd.read_csv(URL_CATALOG)

    # # prepare the search criteria as a string
    # search_string = f"table_id == '3hr' & source_id == '{source_id}' & experiment_id == '{experiment_id}' & variable_id=='{varname}'"
    # df_search = df_catalog.query(search_string)

    # # authenticate access to Google Cloud
    # fs = gcsfs.GCSFileSystem(token="anon", access="read_only")

    # create a MutableMapping from a store URL
    # mapper = fs.get_mapper(df_search.zstore.values[-1])

    # open using xarray
    # ds = xr.open_zarr(mapper, consolidated=True)

    purge_cache_dir = cache_dir is None
    if cache_dir is None:
        cache_dir = tempfile.mkdtemp(prefix="esgf_cache_")
    cache_dir = Path(cache_dir)
    logger.info(f"Using temporary directory {cache_dir} for ESGF cache")
    intake_esgf.conf.set(all_indices=True, local_cache=cache_dir)

    def on_rm_error(function, path, excinfo):
        logger.warning(f"Failed to remove {path}: {excinfo}")

    try:
        cat = intake_esgf.ESGFCatalog()
        result = []
        for varname in variables:
            logger.info(f"Processing variable {varname}")

            cat = cat.search(
                experiment_id=experiment_id,
                source_id=source_id,
                variable_id=varname,
                frequency="3hr",
            )
            cat.remove_ensembles()
            paths = cat.to_path_dict(prefer_streaming=prefer_streaming)[varname]

            logger.info(f"  opening {len(paths)} files")
            da = pygetm.input.from_nc(paths, varname)

            result.append(da.rename(varname))

        yield xr.merge(result, compat="identical")
    finally:
        for _, ds in pygetm.input.open_nc_files:
            ds.close()

        if purge_cache_dir and cache_dir.exists():
            logger.info(f"Clearing temporary directory {cache_dir}")
            shutil.rmtree(cache_dir, onexc=on_rm_error)


def get_meteo(
    minlon: float,
    maxlon: float,
    minlat: float,
    maxlat: float,
    source_id: str = "GFDL-ESM4",
    experiment_id: str = "ssp245",
    target_dir: Union[os.PathLike[str], str] = ".",
    logger: Optional[logging.Logger] = None,
    complevel: int = 0,
    **kwargs,
):
    logger = logger or _create_logger()
    target_dir = Path(target_dir)
    for name in METEO_VARS:
        with get_global_meteo(
            source_id, experiment_id, logger, variables=[name], **kwargs
        ) as ds:
            logger.info(f"  subsetting {name}")
            da_subset = pygetm.input.limit_region(
                ds[name], minlon, maxlon, minlat, maxlat, periodic_lon=True
            )
            logger.info(f"  saving")
            ds = da_subset.to_dataset(name=name)
            pygetm.input.util.configure_chunking_and_compression(
                ds, complevel=complevel
            )
            ds.chunk(time=1000).to_netcdf(target_dir / f"{name}.nc")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("minlon", help="minimum longitude (°East)", type=float)
    parser.add_argument("maxlon", help="maximum longitude (°East)", type=float)
    parser.add_argument("minlat", help="minimum latitude (°North)", type=float)
    parser.add_argument("maxlat", help="maximum latitude (°North)", type=float)
    parser.add_argument("--source_id", help="source_id", default="GFDL-ESM4")
    parser.add_argument("--experiment_id", help="experiment_id", default="ssp245")
    parser.add_argument(
        "--no-streaming",
        help="Disable streaming",
        action="store_false",
        dest="prefer_streaming",
    )
    parser.add_argument(
        "--complevel",
        help="NetCDF compression level (0-9). Default = 0 = no compression",
        type=int,
        default=0,
    )
    parser.add_argument(
        "--list",
        help="List available experiment_ids and source_ids",
        action="store_true",
    )
    parser.add_argument(
        "--cache_dir", help="directory for caching CMIP6 downloads", default=None
    )
    args = parser.parse_args()
    if args.list:
        list_available_models()
    else:
        get_meteo(
            args.minlon,
            args.maxlon,
            args.minlat,
            args.maxlat,
            source_id=args.source_id,
            experiment_id=args.experiment_id,
            cache_dir=args.cache_dir,
            complevel=args.complevel,
            prefer_streaming=args.prefer_streaming,
        )
