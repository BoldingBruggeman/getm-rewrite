import os.path
from typing import Mapping

import numpy as np
import numpy.typing as npt
import cftime
import xarray as xr
from pygetm import otps2

import pygetm.input

ROOT = "../../../igotm/data/TPXO9"

COMPONENTS = ("m2", "s2", "n2", "k2", "k1", "o1", "p1", "q1", "m4", "ms4", "mn4", "2n2")


def get(
    lon: npt.ArrayLike,
    lat: npt.ArrayLike,
    variable: str = "h",
    verbose: bool = False,
    root: str = ROOT,
    scale_factor: float = 1.0,
) -> xr.DataArray:
    assert variable in ("h", "u", "v", "hz", "hu", "hv")

    lon = np.asarray(lon)
    lat = np.asarray(lat)
    if lon.size == 0:
        return xr.DataArray(np.empty_like(lon))

    def select(ncvar) -> xr.DataArray:
        out = pygetm.input.limit_region(
            ncvar, lon.min(), lon.max(), lat.min(), lat.max(), periodic_lon=True
        )
        out = pygetm.input.horizontal_interpolation(out, lon, lat)
        return out

    # Detect version
    if os.path.isfile(os.path.join(root, "grid_tpxo10atlas_v2.nc")):
        version, postfix = "10", "_v2"
    elif os.path.isfile(os.path.join(root, "grid_tpxo9_atlas_30_v5.nc")):
        version, postfix = "9", "_v5"
    else:
        version, postfix = "9", ""

    if variable in ("hz", "hu", "hv"):
        grid_file = os.path.join(root, f"grid_tpxo{version}_atlas{postfix}.nc")
        if not os.path.isfile(grid_file):
            # The template for the grid file is unfortunately different in TPXO10
            grid_file = os.path.join(root, f"grid_tpxo{version}atlas{postfix}.nc")

        # Water depth at z, u, or v points.
        # These are static variables defined in the grid file
        axis = variable[1]
        with xr.open_dataset(grid_file) as ds:
            ds = ds.set_coords((f"lat_{axis}", f"lon_{axis}"))
            return select(ds[variable])

    scale_factor *= {"h": 1e-3, "u": 1e-4, "v": 1e-4}.get(variable, 1.0)
    axis = {"h": "z"}.get(variable, variable)
    file_prefix = {"v": "u"}.get(variable, variable)
    components: Mapping[str, tuple[np.ndarray, np.ndarray]] = {}
    for component in COMPONENTS:
        if verbose:
            print(f"TPXO: reading {component} constituent of {variable}...")
        name = f"{file_prefix}_{component}_tpxo{version}_atlas_30{postfix}.nc"
        path = os.path.join(root, name)
        with xr.open_dataset(path) as ds:
            ds = ds.set_coords((f"lat_{axis}", f"lon_{axis}"))
            re = select(ds[f"{variable}Re"])
            im = select(ds[f"{variable}Im"])
            components[component] = (scale_factor * re.values, scale_factor * im.values)
    lazyvar = Data(components, lat, name=f"tpxo({root!r}, {variable!r})")
    return xr.DataArray(lazyvar, dims=re.dims, coords=re.coords, name=lazyvar.name)


class Data(pygetm.input.LazyArray):
    def __init__(
        self,
        components: Mapping[str, tuple[np.ndarray, np.ndarray]],
        lat: np.ndarray,
        **kwargs,
    ):
        super().__init__(lat.shape, np.float64, **kwargs)
        self.components = components
        self.lat = lat
        self.time = None

    def update(self, time: cftime.datetime, numtime: np.longdouble) -> bool:
        self.time = time
        return True

    def __array__(self, dtype=None, copy=None) -> np.ndarray:
        assert self.time is not None, "update has not yet been called"
        return otps2.predict_tide_2d(
            self.components, self.lat, self.time, ntime=1, delta_time=0
        )[0, ...]

    def is_time_varying(self) -> bool:
        return True
