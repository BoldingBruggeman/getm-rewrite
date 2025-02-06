import tarfile
import shutil
import glob
import tempfile
import os.path
import netCDF4
import argparse
import logging
import urllib.request
from typing import Optional

import numpy as np

from pygetm import pygsw
from pygetm.util.nctools import copy_variable
from pygetm.util.fill import Filler
from pygetm.input import from_nc, horizontal_interpolation, limit_region

URL = "https://www.nodc.noaa.gov/archive/arc0107/0162565/1.1/data/0-data/mapped/GLODAPv2.2016b_MappedClimatologies.tar.gz"


def download(outfile: str, logger: logging.Logger, source: str = URL):
    filename = os.path.basename(source)
    with tempfile.TemporaryDirectory() as root:
        target = os.path.join(root, filename)

        if source.startswith("http"):
            logger.info(f"Downloading {source}...")
            with urllib.request.urlopen(source) as response, open(target, "wb") as fout:
                shutil.copyfileobj(response, fout)
            source = target

        logger.info(f"Extracting {filename}...")
        with tarfile.open(source) as tar:
            tar.extractall(root)

        logger.info("Collecting variables...")
        with netCDF4.Dataset(outfile, "w", format="NETCDF4") as ncout:
            ncout.set_fill_off()
            for path in glob.glob(os.path.join(root, "**/*.nc"), recursive=True):
                with netCDF4.Dataset(path) as nc:
                    nc.set_auto_maskandscale(False)
                    varname = os.path.basename(path).rsplit(".", 2)[-2]
                    ncvar = nc[varname]
                    logger.info(f"  - {varname}: {ncvar.long_name} ({ncvar.units})...")
                    if not ncout.variables:
                        copy_variable(nc["Depth"], ncout)
                        copy_variable(nc["lon"], ncout).units = "degrees_east"
                        copy_variable(nc["lat"], ncout).units = "degrees_north"
                    ncvarout = copy_variable(ncvar, ncout)
                    ncvarout.coordinates = "lon lat Depth"


def fill(path: str, logger: logging.Logger):
    with netCDF4.Dataset(path, "r+") as nc:
        nc.set_auto_maskandscale(False)
        for name in nc.variables:
            if name not in ("lon", "lat", "Depth"):
                ncvar = nc[name]
                data = ncvar[:, :, :]
                mask = data == ncvar._FillValue
                logger.info(f"Filling {mask.sum()} masked points for {name}...")
                filler = Filler(mask, dim_weights=(0.9, 1.0, 1.0))
                filler(data)
                ncvar[:, :, :] = data


def add_density(path: str, logger: logging.Logger):
    with netCDF4.Dataset(path, "r+") as nc:
        lon, lat, z = nc["lon"][:] % 360, nc["lat"][:], nc["Depth"][:]
        lon, lat, z = np.broadcast_arrays(
            lon[np.newaxis, np.newaxis, :],
            lat[np.newaxis, :, np.newaxis],
            z[:, np.newaxis, np.newaxis],
        )
        salt = nc["salinity"][...]
        temp = nc["temperature"][...]
        unmasked = ~(np.ma.getmaskarray(salt) | np.ma.getmaskarray(temp))

        z = z[unmasked]

        logger.info("Calculating absolute salinity from practical salinity...")
        sa = pygsw.sa_from_sp(lon[unmasked], lat[unmasked], z, salt[unmasked])

        logger.info("Calculating conservative temperature...")
        pt = pygsw.pt0_from_t(sa, temp[unmasked], z)
        ct = pygsw.ct_from_pt(sa, pt)

        logger.info("Calculating density...")
        rho = pygsw.rho(sa, ct, z)

        res = np.full_like(nc["salinity"], nc["salinity"]._FillValue)

        res[unmasked] = sa
        ncvar = copy_variable(nc["salinity"], nc, name="sa", copy_data=False)
        ncvar.long_name = "absolute salinity"
        ncvar.units = "g kg-1"
        ncvar[:, :, :] = res

        res[unmasked] = ct
        ncvar = copy_variable(nc["temperature"], nc, name="ct", copy_data=False)
        ncvar.long_name = "conservative temperature"
        ncvar.units = "degrees_Celsius"
        ncvar[:, :, :] = res

        res[unmasked] = rho
        ncvar = copy_variable(nc["salinity"], nc, name="density", copy_data=False)
        ncvar.long_name = "in situ density"
        ncvar.units = "kg m-3"
        ncvar[:, :, :] = res


def regrid(
    infile: str,
    lon,
    lat,
    outfile: str,
    *,
    clamp: bool = False,
    logger: Optional[logging.Logger] = None,
):
    logging.basicConfig(level=logging.INFO)
    logger = logger or logging.getLogger()
    lat_ip = lat
    with netCDF4.Dataset(infile) as nc:
        target_variables = [n for n, v in nc.variables.items() if v.ndim == 3]
        if clamp:
            lat_glodap = nc.variables["lat"][:]
            lat_ip = np.clip(lat, lat_glodap.min(), lat_glodap.max())
            if (lat != lat_ip).any():
                logger.warning(
                    f"Clipping requested latitude ({lat.min():.5f} - {lat.max():.5f})"
                    f" to available GLODAP range ({lat_glodap.min()} - {lat_glodap.max()})"
                )
    logger.info(
        f"Interpolating {len(target_variables)} variables"
        f" and writing them to {outfile}..."
    )
    with netCDF4.Dataset(outfile, "w") as ncbath:
        for name in target_variables:
            values = from_nc(infile, name)
            values = limit_region(values, -180.0, 180.0, -80.0, 89.5, periodic_lon=True)
            if not ncbath.variables:
                ncbath.createDimension("x", lon.shape[1])
                ncbath.createDimension("y", lon.shape[0])
                ncbath.createVariable("lon", float, ("y", "x"))[:, :] = lon
                ncbath.createVariable("lat", float, ("y", "x"))[:, :] = lat
                glodap_depth = values.coords["Depth"]
                ncbath.createDimension("depth", len(glodap_depth))
                ncdepth = ncbath.createVariable("depth", float, ("depth",))
                ncdepth[:] = glodap_depth
                ncdepth.positive = glodap_depth.attrs["positive"]
            values_ip = horizontal_interpolation(values, lon, lat_ip)
            values_ip = np.asarray(values_ip)
            logger.info(f"  {name}: {values.long_name} ({values.units})...")
            ncvar = ncbath.createVariable(
                name, float, ("depth", "y", "x"), fill_value=-1
            )
            ncvar.long_name = values.long_name
            ncvar.units = values.units
            ncvar[:, :, :] = values_ip


def download_and_process(
    outfile: str, source: str = URL, logger: Optional[logging.Logger] = None
):
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger()
    download(outfile, logger=logger, source=source)
    add_density(outfile, logger=logger)
    fill(outfile, logger=logger)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "outfile", nargs="?", default="glodap.nc", help="NetCDF file to write output to"
    )
    parser.add_argument(
        "--source", help="Source of the original GLODAP file (tar.gz).", default=URL
    )
    args = parser.parse_args()

    download_and_process(args.outfile, args.source)
