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
import cftime

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

# CMIP6 atmospheric pCO2:
# * historical: https://doi.org/10.5194/gmd-10-2057-2017
#   https://metagrid.esgf-west.org/search/input4MIPs/?project=input4MIPs&activeFacets=%7B%22mip_era%22%3A%22CMIP6%22%2C%22institution_id%22%3A%22UoM%22%2C%22grid_label%22%3A%22gr-0p5x360deg%22%2C%22data_node%22%3A%22esgf-node.ornl.gov%22%2C%22target_mip_list%22%3A%22CMIP%22%2C%22variable_id%22%3A%22mole-fraction-of-carbon-dioxide-in-air%22%7D
# * scenarios: https://doi.org/10.5194/gmd-9-3461-2016
#   https://metagrid.esgf-west.org/search/input4MIPs/?project=input4MIPs&activeFacets=%7B%22mip_era%22%3A%22CMIP6%22%2C%22institution_id%22%3A%22UoM%22%2C%22grid_label%22%3A%22gr-0p5x360deg%22%2C%22target_mip_list%22%3A%22ScenarioMIP%22%2C%22variable_id%22%3A%22mole_fraction_of_carbon_dioxide_in_air%22%2C%22data_node%22%3A%22esgf-node.ornl.gov%22%2C%22source_version%22%3A%221.2.1%22%7D
EXPERIMENT2CO2ATM = dict(
    historical="mole-fraction-of-carbon-dioxide-in-air_input4MIPs_GHGConcentrations_CMIP_UoM-CMIP-1-2-0_gr-0p5x360deg_000001-201412.nc",
    ssp119="mole-fraction-of-carbon-dioxide-in-air_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-IMAGE-ssp119-1-2-1_gr-0p5x360deg_201501-250012.nc",
    ssp126="mole-fraction-of-carbon-dioxide-in-air_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-IMAGE-ssp126-1-2-1_gr-0p5x360deg_201501-250012.nc",
    ssp245="mole-fraction-of-carbon-dioxide-in-air_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-MESSAGE-GLOBIOM-ssp245-1-2-1_gr-0p5x360deg_201501-250012.nc",
    ssp370="mole-fraction-of-carbon-dioxide-in-air_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-AIM-ssp370-1-2-1_gr-0p5x360deg_201501-250012.nc",
    ssp585="mole-fraction-of-carbon-dioxide-in-air_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-REMIND-MAGPIE-ssp585-1-2-1_gr-0p5x360deg_201501-250012.nc",
)

# Extra usful variables:
# co2 Amon
# siconc SIday

METEO_VARS = METEO_VARS_INSTANT + METEO_VARS_ACCUM


def list_available_models():
    intake_esgf.conf.set()  # all_indices=True)
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
def get_global_cmip6(
    source_id: str = "GFDL-ESM4",
    experiment_id: str = "ssp245",
    logger: Optional[logging.Logger] = None,
    cache_dir: Union[os.PathLike[str], str, None] = None,
    variables: Iterable[str] = METEO_VARS,
    prefer_streaming: bool = True,
    frequency: str = "3hr",
    **kwargs,
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

    with _get_catalog(logger, cache_dir=cache_dir) as cat:
        result = []
        for varname in variables:
            logger.info(f"Processing variable {varname}")

            cat = cat.search(
                experiment_id=experiment_id,
                source_id=source_id,
                variable_id=varname,
                frequency=frequency,
                **kwargs,
            )
            cat.remove_ensembles()
            paths = cat.to_path_dict(prefer_streaming=prefer_streaming)[varname]

            logger.info(f"  opening {len(paths)} files")
            da = pygetm.input.from_nc(paths, varname)

            result.append(da.rename(varname))

        yield xr.merge(result, compat="identical")


@contextlib.contextmanager
def _get_catalog(
    logger: logging.Logger, cache_dir: Union[os.PathLike[str], str, None] = None
) -> Iterator[xr.Dataset]:
    purge_cache_dir = cache_dir is None
    if cache_dir is None:
        cache_dir = tempfile.mkdtemp(prefix="esgf_cache_")
    cache_dir = Path(cache_dir)
    logger.info(f"Using temporary directory {cache_dir} for ESGF cache")
    intake_esgf.conf.set(
        local_cache=cache_dir,
        all_indices=True,
        # indices={"esgf.ceda.ac.uk": True},
        # no_indices=True,
    )

    def on_rm_error(function, path, excinfo):
        logger.warning(f"Failed to remove {path}: {excinfo}")

    try:
        yield intake_esgf.ESGFCatalog()
    finally:
        for _, ds in pygetm.input.open_nc_files:
            ds.close()

        if purge_cache_dir and cache_dir.exists():
            logger.info(f"Clearing temporary directory {cache_dir}")
            shutil.rmtree(cache_dir, onexc=on_rm_error)


@contextlib.contextmanager
def get_global_pco2(
    experiment_id: str = "ssp245",
    logger: Optional[logging.Logger] = None,
    cache_dir: Union[os.PathLike[str], str, None] = None,
    prefer_streaming: bool = True,
    add_longitude: bool = True,
    **kwargs,
) -> Iterator[xr.Dataset]:
    """Get global time- and latitude-dependent atmospheric pCO2.
    This was used to force CMIP6 models, and accordingly it comes from the input4MIPs project.
    It is available for both historical and future scenarios.
    The data is on a 0.5° latitude grid at monthly resolution."""
    logger = logger or _create_logger()

    def _handle_year_zero(ds: xr.Dataset) -> xr.Dataset:
        # Historical forcing uses standard/gregorian calendar and an offset of year 0,
        # which the CF convention does not support:
        # https://cfconventions.org/Data/cf-conventions/cf-conventions-1.13/cf-conventions.html#calendar
        # Therefore, we override the units string to use year -1 as the offset
        # (year -1 is followed by year 1: https://unidata.github.io/cftime/api.html#cftime.date2num).
        def _convert(name: str) -> xr.DataArray:
            values = cftime.num2date(
                ds[name].values,
                units=ds["time"].attrs["units"].replace("0-1-1", "-1-1-1"),
                calendar=ds["time"].attrs["calendar"],
            )
            return (ds[name].dims, values)

        ds.update({"time": _convert("time"), "time_bnds": _convert("time_bnds")})
        return ds

    fn = EXPERIMENT2CO2ATM[experiment_id]
    kwargs["source_id"] = (
        fn.split("_GHGConcentrations_", 1)[1].split("_", 1)[1].split("_gr-", 1)[0]
    )

    with _get_catalog(logger, cache_dir=cache_dir) as cat:
        cat = cat.search(
            project="input4MIPs",
            grid_label="gr-0p5x360deg",
            variable="mole_fraction_of_carbon_dioxide_in_air",
            **kwargs,
        )
        paths = cat.to_path_dict(prefer_streaming=prefer_streaming)
        assert len(paths) == 1, f"Expected exactly one file, but found {len(paths)}"
        path = list(paths.values())[0][0]
        if experiment_id == "historical":
            ds = _handle_year_zero(xr.open_dataset(path, decode_times=False))
        else:
            ds = xr.open_dataset(path)
        pco2 = ds["mole_fraction_of_carbon_dioxide_in_air"]
        if add_longitude:
            pco2 = (
                pco2.expand_dims("lon", axis=-1)
                .pad(lon=(0, 1), mode="edge")
                .assign_coords(lon=[0.0, 360.0])
            )
            pco2.coords["lon"].attrs["units"] = "degrees_east"
        yield pygetm.input.wrap(pco2)
        ds.close()


def get(
    minlon: float,
    maxlon: float,
    minlat: float,
    maxlat: float,
    source_id: str = "GFDL-ESM4",
    target_dir: Union[os.PathLike[str], str] = ".",
    logger: Optional[logging.Logger] = None,
    complevel: int = 0,
    target: str = "meteo",
    **kwargs,
):
    logger = logger or _create_logger()
    target_dir = Path(target_dir)

    def _save(da, name):
        logger.info(f"  subsetting {name}")
        da_subset = pygetm.input.limit_region(
            da, minlon, maxlon, minlat, maxlat, periodic_lon=True
        )
        path = target_dir / f"{name}.nc"
        logger.info(f"  saving to {path}")
        ds = da_subset.to_dataset(name=name)
        pygetm.input.util.configure_chunking_and_compression(ds, complevel=complevel)
        ds.chunk(time=1000).to_netcdf(path)

    if target == "pco2":
        with get_global_pco2(logger=logger, **kwargs) as da:
            _save(da, "pco2")
    if target == "ts":
        for name in ("thetao", "so"):
            with get_global_cmip6(
                variables=[name],
                source_id=source_id,
                frequency="mon",
                grid_label="gr",
                logger=logger,
                **kwargs,
            ) as ds:
                _save(ds[name], name)
    else:
        for name in METEO_VARS:
            with get_global_cmip6(
                source_id, logger=logger, variables=[name], **kwargs
            ) as ds:
                _save(ds[name], name)


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
    parser.add_argument(
        "-t",
        "--target",
        choices=["meteo", "pco2", "ts"],
        default="meteo",
        help="which data to download",
    )
    args = parser.parse_args()
    if args.list:
        list_available_models()
    else:
        get(
            args.minlon,
            args.maxlon,
            args.minlat,
            args.maxlat,
            source_id=args.source_id,
            experiment_id=args.experiment_id,
            cache_dir=args.cache_dir,
            complevel=args.complevel,
            prefer_streaming=args.prefer_streaming,
            target=args.target,
        )
