import argparse
from typing import Optional, Literal
import datetime
import logging

import netCDF4
import scipy.spatial
import numpy as np


def copyNcVariable(
    ncvar,
    nctarget,
    dimensions=None,
    copy_data=True,
    chunksizes=None,
    name=None,
    zlib=False,
):
    if name is None:
        name = ncvar.name
    if dimensions is None:
        dimensions = ncvar.dimensions
    for dim in dimensions:
        if dim not in nctarget.dimensions:
            length = ncvar.shape[ncvar.dimensions.index(dim)]
            nctarget.createDimension(dim, length)
    fill_value = None if not hasattr(ncvar, "_FillValue") else ncvar._FillValue
    ncvarnew = nctarget.createVariable(
        name,
        ncvar.dtype,
        dimensions,
        fill_value=fill_value,
        chunksizes=chunksizes,
        zlib=zlib,
    )
    ncvarnew.setncatts(
        {att: getattr(ncvar, att) for att in ncvar.ncattrs() if att != "_FillValue"}
    )
    ncvarnew.set_auto_maskandscale(False)
    if copy_data:
        ncvarnew[...] = ncvar[...]
    return ncvarnew


def convert(
    name: Literal["t", "s", "n", "p", "i", "o"],
    outfile: str,
    decade: Optional[
        Literal[
            "decav", "5564", "6574", "7584", "8594", "95A4", "A5B7", "decav81B0", "all"
        ]
    ] = None,
    grid: Optional[Literal["01", "04", "5d"]] = None,
    template: Optional[str] = None,
    logger: Optional[logging.Logger] = None,
) -> str:
    logger = logger or logging.getLogger()
    if grid is None:
        grid = "04" if name in "ts" else "01"

    long_grid = {"04": "0.25", "01": "1.00", "5d": "5deg"}[grid]
    long_name = {
        "t": "temperature",
        "s": "salinity",
        "n": "nitrate",
        "p": "phosphate",
        "i": "silicate",
        "o": "oxygen",
    }[name]
    if decade is None:
        decade = "decav" if name in "ts" else "all"
    if template is None:
        template = f"https://www.ncei.noaa.gov/thredds-ocean/dodsC/ncei/woa/{long_name}/{decade}/{long_grid}/woa18_{decade}_{name}%02i_{grid}.nc"
    varname = name + "_an"
    logger.info(f"Combining annual mean and monthly values for {name}...")
    with netCDF4.Dataset(outfile, "w", format="NETCDF3_64BIT_OFFSET") as ncout:
        ncout.set_fill_off()
        ncout.createDimension("time", 12)
        logger.info("  - mean climatology")
        with netCDF4.Dataset(template % 0) as ncin:
            ncin.set_auto_maskandscale(False)
            copyNcVariable(ncin["depth"], ncout)
            copyNcVariable(ncin["lon"], ncout)
            copyNcVariable(ncin["lat"], ncout)
            ncvar = ncin[varname]
            andata = ncvar[...]
            dim = ncvar.dimensions[1]
        nctime = ncout.createVariable("time", float, ("time",))
        time_ref = datetime.datetime(2000, 1, 1)
        nctime.units = f"days since {time_ref}"
        nctime.calendar = "standard"
        logger.info("  - monthly data:")
        for imonth in range(1, 13):
            logger.info("    - %i" % imonth)
            dt = datetime.datetime(2000, imonth, 15, 12)
            with netCDF4.Dataset(template % imonth) as ncin:
                nctime[imonth - 1] = (dt - time_ref).days
                ncin.set_auto_maskandscale(False)
                ncvar = ncin[varname]
                if imonth == 1:
                    ncvar_mo = copyNcVariable(
                        ncvar,
                        ncout,
                        dimensions=("time", dim, "lat", "lon"),
                        copy_data=False,
                        name=varname,
                    )
                    ncvar_mo[:, ncvar.shape[1] :, :, :] = andata[
                        :, ncvar.shape[1] :, :, :
                    ]
                ncvar_mo[imonth - 1, : ncvar.shape[1], :, :] = ncvar[0, ...]
    return varname


def fill(path: str, varname: str, logger: Optional[logging.Logger] = None):
    logger = logger or logging.getLogger()
    logger.info("Filling missing values...")

    with netCDF4.Dataset(path, "r+") as nc:
        ncvar = nc[varname]
        coords = np.moveaxis(np.indices(ncvar.shape[1:], dtype=float), 0, -1)
        coords[:, 0] *= 0.9  # favour nearby depths over nearby lon/lat
        masked = np.ma.getmaskarray(ncvar[0, :, :, :])
        unmasked = ~masked
        logger.info("  - building KDTree")
        tree = scipy.spatial.cKDTree(coords[unmasked, :])
        logger.info("  - finding nearest neighbors")
        dist, inearest = tree.query(coords[masked, :], workers=-1)
        for itime in range(ncvar.shape[0]):
            logger.info(f"  - writing time {itime}")
            data = ncvar[itime, ...]
            data[masked] = data[unmasked][inearest]
            ncvar[itime, :, :, :] = data


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger()

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "name",
        choices="tsnpio",
        help="character identifying the variable to operate on (t=temperature, s=salinity, n=nitrate, p=phosphate, i=silicate, o=oxygen)",
    )
    parser.add_argument(
        "outfile", nargs="?", help="output file, defaults to <variable character>.nc"
    )
    parser.add_argument(
        "--template",
        default=None,
        help="python format string (%% based) specifying the NetCDF files to operate on. It must contain %%02i as placeholder for the month",
    )
    parser.add_argument(
        "--grid",
        choices=("01", "04", "5d"),
        help="Resolution to use (only used if --template is not provided), defaults to highest available",
    )
    parser.add_argument(
        "--decade",
        default=None,
        help="Decade to use (only used if --template is not provided), defaults to average over entire time period",
    )
    parser.add_argument(
        "--fill", action="store_true", help="whether to fill missing values"
    )
    args = parser.parse_args()
    if args.outfile is None:
        args.outfile = f"{args.name}.nc"
    varname = convert(
        args.name,
        args.outfile,
        grid=args.grid,
        decade=args.decade,
        template=args.template,
        logger=logger,
    )
    if args.fill:
        fill(args.outfile, varname, logger=logger)
    logger.info(f"Results saved to {args.outfile}")
