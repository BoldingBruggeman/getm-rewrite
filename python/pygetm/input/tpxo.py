import os.path
from typing import Tuple, Mapping

import numpy
import cftime
import xarray
import otps2
from numpy.typing import ArrayLike

import pygetm.input

ROOT = "../../../igotm/data/TPXO9"

COMPONENTS = ("m2", "s2", "n2", "k2", "k1", "o1", "p1", "q1", "m4", "ms4", "mn4", "2n2")


def get(
    lon: ArrayLike,
    lat: ArrayLike,
    variable: str = "h",
    verbose: bool = False,
    root: str = ROOT,
    scale_factor: float = 1.0,
) -> xarray.DataArray:
    assert variable in ("h", "u", "v", "hz", "hu", "hv")

    lon = numpy.asarray(lon)
    lat = numpy.asarray(lat)
    if lon.size == 0:
        return xarray.DataArray(numpy.empty_like(lon))

    def select(ncvar) -> xarray.DataArray:
        out = pygetm.input.limit_region(
            ncvar, lon.min(), lon.max(), lat.min(), lat.max(), periodic_lon=True
        )
        out = pygetm.input.horizontal_interpolation(out, lon, lat)
        return out

    postfix = ""
    if os.path.isfile(os.path.join(root, "grid_tpxo9_atlas_30_v5.nc")):
        postfix = "_v5"

    if variable in ("hz", "hu", "hv"):
        # Water depth at z, u, or v points.
        # These are static variables definied in the grid file
        axis = variable[1]
        with xarray.open_dataset(
            os.path.join(root, f"grid_tpxo9_atlas{postfix}.nc")
        ) as ds:
            ds = ds.set_coords((f"lat_{axis}", f"lon_{axis}"))
            return select(ds[variable])

    scale_factor *= {"h": 1e-3, "u": 1e-4, "v": 1e-4}.get(variable, 1.0)
    axis = {"h": "z"}.get(variable, variable)
    file_prefix = {"v": "u"}.get(variable, variable)
    components: Mapping[str, Tuple[numpy.ndarray, numpy.ndarray]] = {}
    for component in COMPONENTS:
        if verbose:
            print(f"TPXO: reading {component} constituent of {variable}...")
        name = f"{file_prefix}_{component}_tpxo9_atlas_30{postfix}.nc"
        path = os.path.join(root, name)
        with xarray.open_dataset(path) as ds:
            ds = ds.set_coords((f"lat_{axis}", f"lon_{axis}"))
            x = select(ds[f"{variable}Im"])
            components[component] = (
                scale_factor * select(ds[f"{variable}Re"]).values,
                scale_factor * x.values,
            )
    lazyvar = Data(components, lat, name=f"tpxo({root!r}, {variable!r})")
    return xarray.DataArray(lazyvar, dims=x.dims, coords=x.coords, name=lazyvar.name)


class Data(pygetm.input.LazyArray):
    def __init__(
        self,
        components: Mapping[str, Tuple[numpy.ndarray, numpy.ndarray]],
        lat: numpy.ndarray,
        **kwargs,
    ):
        super().__init__(lat.shape, numpy.float64, **kwargs)
        self.components = components
        self.lat = lat
        self.time = None

    def update(self, time: cftime.datetime, numtime: numpy.longdouble) -> bool:
        self.time = time
        return True

    def __array__(self, dtype=None) -> numpy.ndarray:
        assert self.time is not None, "update has not yet been called"
        return otps2.predict_tide_2d(
            self.components, self.lat, self.time, ntime=1, delta_time=0
        )[0, ...]

    def is_time_varying(self) -> bool:
        return True
