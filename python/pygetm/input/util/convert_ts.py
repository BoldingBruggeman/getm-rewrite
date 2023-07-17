from typing import Mapping, Optional
import logging
import argparse
import datetime
import sys

import numpy as np
import netCDF4

import pygsw
from ..woa import copyNcVariable


def convert_ts(
    sfile: str,
    tfile: str,
    outfile: str,
    sname: Optional[str] = None,
    tname: Optional[str] = None,
    saname: str = "sa",
    ctname: str = "ct",
    logger: Optional[logging.Logger] = None,
):
    def _replace_nc_attrs(ncvar: netCDF4.Variable, mapping: Mapping[str, str]):
        for k in ncvar.ncattrs():
            v = getattr(ncvar, k)
            if isinstance(v, str):
                for m1, m2 in mapping.items():
                    v = v.replace(m1, m2)
                setattr(ncvar, k, v)

    def _find_standard_variable(ncfile: netCDF4.Dataset, standard_name: str) -> str:
        for n, v in ncfile.variables.items():
            if getattr(v, "standard_name", None) == standard_name:
                logger.info(f"Found {standard_name} = {n} in {ncfile.filepath()}")
                return n
        else:
            raise Exception(
                f"No variable with standard_name {standard_name} found in"
                f" {ncfile.filepath()}. Please specify the variable name as argument."
            )

    logger = logger or logging.getLogger()
    nc_s_file = netCDF4.Dataset(sfile, "r")
    nc_t_file = netCDF4.Dataset(tfile, "r")
    ncout = netCDF4.Dataset(outfile, "w")
    logger.info(f"Writing {saname} and {ctname} to {ncout.filepath()}")
    with nc_t_file, nc_s_file, ncout:
        if sname is None:
            sname = _find_standard_variable(nc_s_file, "sea_water_salinity")
        if tname is None:
            tname = _find_standard_variable(nc_t_file, "sea_water_temperature")

        cmdline = " ".join(sys.argv)
        ncout.history = f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} {cmdline}"

        nc_t = nc_t_file[tname]
        nc_s = nc_s_file[sname]
        assert nc_t.ndim == 4
        assert nc_t.shape == nc_s.shape

        nc_ct = copyNcVariable(nc_t, ncout, name=ctname, copy_data=False)
        nc_sa = copyNcVariable(nc_s, ncout, name=saname, copy_data=False)
        copyNcVariable(nc_t_file["depth"], ncout)
        copyNcVariable(nc_t_file["lon"], ncout)
        copyNcVariable(nc_t_file["lat"], ncout)
        copyNcVariable(nc_t_file["time"], ncout)
        _replace_nc_attrs(nc_sa, {"sea_water_salinity": "sea_water_absolute_salinity"})
        nc_sa.units = "g kg-1"
        _replace_nc_attrs(
            nc_ct, {"sea_water_temperature": "sea_water_conservative_temperature"}
        )

        lon = nc_t_file["lon"][:].astype(float)[np.newaxis, np.newaxis, :]
        lat = nc_t_file["lat"][:].astype(float)[np.newaxis, :, np.newaxis]
        p = nc_t_file["depth"][:].astype(float)[:, np.newaxis, np.newaxis]
        lon, lat, p = np.broadcast_arrays(lon, lat, p)
        ct_slice = np.full_like(nc_ct[0, ...], nc_ct._FillValue)
        sa_slice = np.full_like(nc_sa[0, ...], nc_sa._FillValue)
        valid = ~np.ma.getmaskarray(nc_t[0, ...])
        logger.info(f"{valid.sum()} of {valid.size} cells are water")
        p = p[valid]
        lon = lon[valid]
        lat = lat[valid]
        for itime in range(nc_t.shape[0]):
            logger.info(f"  - writing time {itime}")
            s = nc_s[itime, ...][valid].astype(float)
            t = nc_t[itime, ...][valid].astype(float)
            sa = pygsw.sa_from_sp(lon, lat, p, s)
            pt0 = pygsw.pt0_from_t(sa, t, p)
            ct_slice[valid] = pygsw.ct_from_pt(sa, pt0)
            sa_slice[valid] = sa
            nc_ct[itime, ...] = ct_slice
            nc_sa[itime, ...] = sa_slice


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger()

    parser = argparse.ArgumentParser()
    parser.add_argument("sfile", help="file with practical salinity")
    parser.add_argument("tfile", help="file with in-situ temperature")
    parser.add_argument(
        "outfile",
        help="output file to write absolute salinity and conservative temperature to",
    )
    parser.add_argument("--sname", help="name of practical salinity variable (input)")
    parser.add_argument("--tname", help="name of in-situ temperature variable (input)")
    parser.add_argument(
        "--saname", help="name of absolute salinity variable (output)", default="sa"
    )
    parser.add_argument(
        "--ctname",
        help="name of conservative temperature variable (output)",
        default="ct",
    )
    args = parser.parse_args()
    convert_ts(
        args.sfile,
        args.tfile,
        args.outfile,
        sname=args.sname,
        tname=args.tname,
        saname=args.saname,
        ctname=args.ctname,
        logger=logger,
    )
