from typing import Optional
import xarray as xr
import cftime


def replace_calendar(da: xr.DataArray, calendar: str) -> xr.DataArray:
    tmcoord = da.getm.time
    new_time = [
        cftime.datetime(
            dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second, calendar=calendar
        )
        for dt in tmcoord.values
    ]
    return da.assign_coords({tmcoord.name: new_time})


def configure_chunking_and_compression(
    ds: xr.Dataset, complevel: Optional[int] = None
) -> None:
    """Configure chunking and compression for all DataArrays in a Dataset.
    This is done by manipulating the `encoding` attribute of each DataArray in-place.
    This attribute controls how the data will be written to a NetCDF file.
    Data variables are compressed if `complevel` is non-zero, otherwise decompressed.
    Compressed variables have chunk size 1 along the time dimension,
    and chunk sizes equal to the full dimension size along other dimensions.
    Coordinate variables are never compressed.

    Args:
        ds: Dataset to configure
        complevel: Compression level (0-9) for data variables, or None to use
            the original compression level (if any). If 0, no compression
            is applied.
    """

    def _configure(da: xr.DataArray, complevel: Optional[int]):
        if complevel is not None:
            da.encoding["complevel"] = complevel
        if da.encoding.get("complevel", 0) == 0:
            # Decompress (no chunking at all)
            da.encoding.pop("chunksizes", None)
        else:
            # Set chunk size for time dimension to 1 to allow
            # efficient access to values @ single time points
            da.encoding["chunksizes"] = list(da.shape)
            da.encoding["chunksizes"][0] = 1
            da.encoding["shuffle"] = True
            da.encoding["zlib"] = True
        da.encoding.pop("contiguous", None)

    for da in ds.coords.values():
        _configure(da, 0)
    for da in ds.data_vars.values():
        _configure(da, complevel)
