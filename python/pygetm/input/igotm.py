import http.client
import json
import logging
from typing import Optional
import argparse

import urllib.parse
import numpy as np
import numpy.typing as npt
import cftime
import xarray as xr


ERA5_LONG_NAMES = dict(
    u10="eastward wind speed",
    v10="northward wind speed",
    sp="surface pressure",
    t2m="air temperature @ 2 m",
    d2m="dew point temperature @ 2 m",
    tcc="total cloud cover",
    ssr="surface net solar radiation",
    tp="total precipitation",
)
ERA5_UNITS = dict(
    u10="m s-1",
    v10="m s-1",
    sp="Pa",
    t2m="degrees_Celsius",
    d2m="degrees_Celsius",
    tcc="1",
    ssr="W m-2",
    tp="m s-1",
)


def download_era5(
    lng: npt.ArrayLike,
    lat: npt.ArrayLike,
    start_year: int,
    stop_year: Optional[int] = None,
    dims: Optional[tuple[str]] = None,
    logger: Optional[logging.Logger] = None,
) -> xr.Dataset:
    lng = np.asarray(lng)
    lat = np.asarray(lat)
    assert lng.shape == lat.shape
    if dims is None:
        dims = tuple(f"dim_{i}" for i in range(lng.ndim))
    assert lng.ndim == len(dims), (
        f"Number of coordinate dimensions {lng.ndim}"
        f" does not match number of dimension names in {dims}"
    )
    if stop_year is None:
        stop_year = start_year
    logger = logger or logging.getLogger()

    data_vars = dict(
        lng=xr.DataArray(
            lng, dims=dims, attrs=dict(long_name="longitude", units="degrees_east")
        ),
        lat=xr.DataArray(
            lat, dims=dims, attrs=dict(long_name="latitude", units="degrees_north")
        ),
    )
    it = np.nditer([lng, lat], flags=["multi_index"])
    conn = http.client.HTTPSConnection("igotm.bolding-bruggeman.com")
    query = {"target": "era5", "minyear": start_year, "maxyear": stop_year}
    for lo, la in it:
        logger.info(f"Retrieving ERA5 for {lo:.6f} °East, {la:.6f} °North")
        url = urllib.parse.urlencode(query | {"lng": lo, "lat": la})
        conn.request("GET", "/get?" + url)
        response = conn.getresponse()
        assert response.status == 200, f"Server returned error code {response.status}"
        data = json.load(response)
        for k, v in data.items():
            if k not in ("start", "step"):
                if "time" not in data_vars:
                    numtime = np.arange(len(v)) * data["step"]
                    time_units = f"seconds since {data['start'].replace('T', ' ')}"
                    data_vars["time"] = xr.DataArray(
                        cftime.num2date(numtime, time_units), dims=("time",)
                    )
                if k not in data_vars:
                    values = np.empty((len(v),) + lng.shape, dtype=np.float32)
                    coords = {k: data_vars[k] for k in ("time", "lat", "lng")}
                    attrs = dict(long_name=ERA5_LONG_NAMES[k], units=ERA5_UNITS[k])
                    data_vars[k] = xr.DataArray(
                        values, coords=coords, dims=("time",) + dims, attrs=attrs
                    )
                data_vars[k].values[(slice(None),) + it.multi_index] = v
    conn.close()
    return xr.Dataset(data_vars)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("lng", type=float, help="longitude")
    parser.add_argument("lat", type=float, help="latitude")
    parser.add_argument("year", type=int, help="year")
    parser.add_argument("outfile", help="NetCDF file to write forcing to")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    ds = download_era5(args.lng, args.lat, args.year)
    ds.to_netcdf(args.outfile)
