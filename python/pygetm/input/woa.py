import argparse
from typing import Optional
import datetime

import netCDF4
import xarray


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


def convert(woapaths: str, outfile: str, name: Optional[str]) -> str:
    print('Combining annual mean and monthly values...')
    with netCDF4.Dataset(outfile, "w", format="NETCDF3_64BIT_OFFSET") as ncout:
        ncout.set_fill_off()
        ncout.createDimension("time", 12)
        print("  - mean climatology")
        with netCDF4.Dataset(woapaths % 0) as ncin:
            if name is None:
                for name in ncin.variables:
                    if name.endswith("_an"):
                        break
            ncin.set_auto_maskandscale(False)
            copyNcVariable(ncin["depth"], ncout)
            copyNcVariable(ncin["lon"], ncout)
            copyNcVariable(ncin["lat"], ncout)
            ncvar = ncin[name]
            andata = ncvar[...]
            dim = ncvar.dimensions[1]
        nctime = ncout.createVariable("time", float, ("time",))
        time_ref = datetime.datetime(2000, 1, 1)
        nctime.units = f"days since {time_ref}"
        print("  - monthly data:")
        for imonth in range(1, 13):
            print("    - %i" % imonth)
            dt = datetime.datetime(2000, imonth, 15, 12)
            with netCDF4.Dataset(woapaths % imonth) as ncin:
                nctime[imonth - 1] = (dt - time_ref).days
                ncin.set_auto_maskandscale(False)
                ncvar = ncin[name]
                if imonth == 1:
                    ncvar_mo = copyNcVariable(
                        ncvar,
                        ncout,
                        dimensions=("time", dim, "lat", "lon"),
                        copy_data=False,
                        name=name,
                    )
                    ncvar_mo[:, ncvar.shape[1] :, :, :] = andata[
                        :, ncvar.shape[1] :, :, :
                    ]
                ncvar_mo[imonth - 1, : ncvar.shape[1], :, :] = ncvar[0, ...]
    return name

def fill(path: str, name: str):
    print('Filling missing values...')
    ds = xarray.open_dataset(path)
    da = ds[name].interpolate_na(dim='lat', method='nearest')
    da.to_netcdf(path)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("source")
    parser.add_argument("outfile")
    parser.add_argument("--name", default=None)
    parser.add_argument("--fill", action='store_true')
    args = parser.parse_args()
    name = convert(args.source, args.outfile, name=args.name)
    if args.fill:
        fill(args.outfile, name)
