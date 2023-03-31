import argparse
import datetime
from typing import Iterable, Mapping

import cftime
import numpy as np
import xarray


METEO_COLUMNS = ("u10", "v10", "sp", "t2m", "hum", "tcc")
METEO_LONG_NAMES = dict(
    u10="Eastward wind speed",
    v10="Northward wind speed",
    sp="surface pressure",
    t2m="air temperature @ 2 m",
    hum="humidity @ 2 m",
    tcc="total cloud cover",
)
METEO_UNITS = dict(u10="m s-1", v10="m s-1", tcc="1")


def get_timeseries(
    path: str,
    columns: Iterable[str],
    units: Mapping[str, str] = {},
    long_names: Mapping[str, str] = {},
) -> xarray.Dataset:
    dts = []
    all_values = []
    with open(path) as f:
        last_dt = None
        for line in f:
            if line.startswith("#"):
                continue
            dt = datetime.datetime.strptime(line[:19], "%Y-%m-%d %H:%M:%S")
            assert last_dt is None or last_dt <= dt
            items = line[20:].rstrip("\n").split()
            assert len(items) == len(columns)
            dts.append(
                cftime.datetime(
                    dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second
                )
            )
            all_values.append(list(map(float, items)))

    timevar = xarray.DataArray(np.array(dts), dims=("time",))
    name2var = dict(time=timevar)
    for name, values in zip(columns, np.transpose(all_values)):
        attrs = {}
        if name in units:
            attrs["units"] = units[name]
        if name in long_names:
            attrs["long_name"] = long_names[name]
        name2var[name] = xarray.DataArray(
            values, coords={"time": timevar}, dims=("time",), name=name, attrs=attrs
        )
    return xarray.Dataset(name2var)


def get_meteo(path: str) -> xarray.Dataset:
    return get_timeseries(path, METEO_COLUMNS, METEO_UNITS, METEO_LONG_NAMES)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("file")
    args = parser.parse_args()
    ds = get_meteo(args.file)
    print(ds)
