from typing import Optional, Union
import os

import numpy as np
import netCDF4

import pygetm.domain
import pygetm.open_boundaries


def domain_from_topo(
    path: Union[str, os.PathLike[str]], **kwargs
) -> pygetm.domain.Domain:
    """Create a domain object from a topo.nc file used by legacy GETM.

    Args:
        path: NetCDF file with legacy domain topography
        **kwargs: keyword arguments that are ultimately passed to
            :class:`pygetm.domain.Domain`

    Bottom roughness can be prescribed with keyword argument ``z0``,
    which will be passed to :class:`pygetm.domain.Domain`. This argument `must`
    be provided if the topo file does not contain bottom roughness ``z0``.
    If ``z0`` is present in the argument list as well as the topo file, the
    argument takes priority.
    """

    def _get_metrics(nc: netCDF4.Dataset):
        H = nc.variables["bathymetry"]
        yname, xname = H.dimensions
        if hasattr(H, "missing_value"):
            H = np.ma.masked_equal(H, H.missing_value)
        mask = np.where(np.ma.getmaskarray(H), 0, 1)

        # Follow legacy GETM in inferring regularly spaced grids from
        # first and last value. This could improve accuracy if values in the topo
        # file are not stored in double precision.
        ncx = nc.variables[xname]
        ncy = nc.variables[yname]
        x = np.linspace(ncx[0], ncx[-1], ncx.size, dtype=float)
        y = np.linspace(ncy[0], ncy[-1], ncy.size, dtype=float)

        return x, y, H, mask

    with netCDF4.Dataset(path) as nc:
        grid_type = int(np.reshape(nc["grid_type"], ()))
        if grid_type == 1:
            # Cartesian
            x, y, H, mask = _get_metrics(nc)
            kwargs.setdefault("lon", nc.variables.get("lonc", None))
            kwargs.setdefault("lat", nc.variables.get("latc", None))
            domain = pygetm.domain.create_cartesian(x, y, H=H, mask=mask, **kwargs)
        elif grid_type == 2:
            # spherical
            lon, lat, H, mask = _get_metrics(nc)
            domain = pygetm.domain.create_spherical(lon, lat, H=H, mask=mask, **kwargs)
        elif grid_type == 3:
            # planar curvilinear
            raise NotImplementedError(
                "No support yet for planar curvilinear coordinates"
            )
        elif grid_type == 4:
            # spherical curvilinear
            raise NotImplementedError(
                "No support yet for spherical curvilinear coordinates"
            )
        else:
            raise NotImplementedError(f"Unknown grid_type {grid_type} found")

        if "z0" not in kwargs:
            if "z0" not in nc.variables:
                raise Exception(
                    f"Bottom roughness z0 is not present in {path}; you need to"
                    " provide it as keyword argument to domain_from_topo instead."
                )
            domain.z0 = nc["z0"]

    return domain


class DatFile:
    """Support for reading GETM dat files with comments indicated by ! or #.
    Whitespace-only lines are skipped."""

    def __init__(self, path: Union[str, os.PathLike[str]]):
        self.path = path
        self.f = open(path)

    def get_line(self) -> str:
        """Return next non-empty line"""
        line = None
        while not line:
            line = self.f.readline()
            assert (
                line != ""
            ), f"End-of-file reached in {self.path} while trying to read next line."
            line = line.split("#", 1)[0].split("!", 1)[0].strip()
        return line

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.f.close()


def load_bdyinfo(
    domain: pygetm.domain.Domain,
    path: Union[str, os.PathLike[str]],
    type_2d: Optional[int] = None,
    type_3d: Optional[int] = None,
):
    """Add open boundaries from bdyinfo.dat to domain.

    Args:
        domain: domain to add open boundaries to
        path: data file with open boundary information
        type_2d: type of 2D open boundary condition to use.
            If provided, this overrides the type configured in the file.
        type_3d: type of 3D open boundary condition to use.
            If provided, this overrides the type configured in the file.
    """
    with DatFile(path) as f:
        for side in (
            pygetm.open_boundaries.Side.WEST,
            pygetm.open_boundaries.Side.NORTH,
            pygetm.open_boundaries.Side.EAST,
            pygetm.open_boundaries.Side.SOUTH,
        ):
            n = int(f.get_line())
            for _ in range(n):
                # Note: for Western and Eastern boundaries, l and m are indices in x
                # and y dimensions, respectively, but that is the other way around
                # (y and x, respectively) for Northern and Southern boundaries.
                # Note that indices are 1-based as in Fortran. We convert to the Python
                # convention: 0-based indices, with the upper bound being the first
                # index that is EXcluded.
                l, mstart, mstop, type_2d_, type_3d_ = map(int, f.get_line().split())
                domain.open_boundaries.add_by_index(
                    side,
                    l - 1,
                    mstart - 1,
                    mstop,
                    type_2d_ if type_2d is None else type_2d,
                    type_3d_ if type_3d is None else type_3d,
                )


def load_riverinfo(domain: pygetm.domain.Domain, path: Union[str, os.PathLike[str]]):
    """Add rivers from riverinfo.dat to domain

    Args:
        domain: domain to add rivers to
        path: data file with river information
    """
    # First count how many times each river appears.
    # Rivers that appear multiple times are split over multiple cells
    # and will have an index appended to their name in the final list of rivers.
    name2split = {}
    with DatFile(path) as f:
        n = int(f.get_line())
        for _ in range(n):
            items = f.get_line().split()
            assert len(items) in (3, 5)
            name = items[2]
            name2split[name] = name2split.get(name, 0) + 1

    name2count = {}
    with DatFile(path) as f:
        n = int(f.get_line())
        for _ in range(n):
            items = f.get_line().split()
            i, j = int(items[0]), int(items[1])  # Indices of the river cell (1-based!)
            mouth_name = name = items[2]

            # Depth extent: zl and zu are depths i.e. measured from the surface.
            # Negative values have the following meaning - if zl < 0 use bottom
            # and if zu < 0 use surface. zl and zu are optional but either none
            # or both must be specified
            zl, zu = np.inf, 0.0
            if len(items) == 5:
                zl, zu = float(items[3]), float(items[4])
                zu = max(0.0, zu)
                if zl < 0:
                    zl = np.inf

            if name2split[name] > 1:
                # This river is split over multiple cells; append an index to its name
                imouth = name2count.get(name, 0)
                mouth_name = f"{name}[{imouth}]"
                name2count[name] = imouth + 1

            # Note: we convert from 1-based indices to 0-based indices!
            river = domain.rivers.add_by_index(mouth_name, i - 1, j - 1, zl=zl, zu=zu)

            river.split = name2split[name]  # number of cells this river is split over
            river.original_name = name  # the original river name (without index)
