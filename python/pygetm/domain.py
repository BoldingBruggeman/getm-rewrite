from typing import Mapping, Optional, Tuple, Union, TYPE_CHECKING
import enum
import functools
import logging

import numpy as np
import numpy.typing as npt

from . import core
from . import parallel
from . import rivers
from . import open_boundaries
from .constants import CoordinateType, GRAVITY

if TYPE_CHECKING:
    import matplotlib.figure
    import matplotlib.colors


class EdgeTreatment(enum.Enum):
    MISSING = enum.auto()
    CLAMP = enum.auto()
    PERIODIC = enum.auto()
    EXTRAPOLATE = enum.auto()
    EXTRAPOLATE_PERIODIC = enum.auto()


def _get_rectangle_overlap(
    istart1: int,
    istop1: int,
    jstart1: int,
    jstop1: int,
    istart2: int,
    istop2: int,
    jstart2: int,
    jstop2: int,
):
    istart = max(istart1, istart2)
    ni = min(istop1, istop2) - istart
    jstart = max(jstart1, jstart2)
    nj = min(jstop1, jstop2) - jstart
    iskip = istart - istart1
    jskip = jstart - jstart1
    slice1 = (slice(jskip, jskip + nj), slice(iskip, iskip + ni))
    iskip = istart - istart2
    jskip = jstart - jstart2
    slice2 = (slice(jskip, jskip + nj), slice(iskip, iskip + ni))
    return slice1, slice2


DEG2RAD = np.pi / 180  # degree to radian conversion
RAD2DEG = 180 / np.pi  # radian to degree conversion
R_EARTH = 6378815.0  # radius of the earth (m)
OMEGA = (
    2.0 * np.pi / 86164.0
)  # rotation rate of the earth (rad/s), 86164 is number of seconds in a sidereal day


def coriolis(lat: npt.ArrayLike) -> np.ndarray:
    """Calculate Coriolis parameter f for the given latitude.

    Args:
        lat: latitude (°North)

    Returns:
        Coriolis parameter f
    """
    return (2.0 * OMEGA) * np.sin(DEG2RAD * np.asarray(lat, dtype=float))


def centers_to_supergrid_1d(
    source: npt.ArrayLike,
    *,
    dtype: Optional[npt.DTypeLike] = None,
    edges: EdgeTreatment = EdgeTreatment.MISSING,
    missing_value=np.nan,
) -> np.ndarray:
    source = np.asarray(source)
    assert source.ndim == 1, "source must be one-dimensional"
    assert source.size > 1, "source must have at least 2 elements"
    out_shape = (source.size * 2 + 1,)
    out = np.empty_like(source, shape=out_shape, dtype=dtype)
    out[1::2] = source
    out[2:-2:2] = 0.5 * (source[1:] + source[:-1])
    edges = EdgeTreatment(edges)
    if edges == EdgeTreatment.PERIODIC:
        # Reconstruct first and last interface by averaging very first and last values
        out[0] = out[-1] = 0.5 * (out[1] + out[-2])
    elif edges == EdgeTreatment.EXTRAPOLATE:
        # Reconstruct first and last interface through linear extrapolation,
        # using nearest interior difference
        out[0] = 2 * out[1] - out[2]
        out[-1] = 2 * out[-2] - out[-3]
    elif edges == EdgeTreatment.EXTRAPOLATE_PERIODIC:
        # Reconstruct first and last interface through linear extrapolation,
        # using oppposite interior difference
        out[0] = out[1] + (out[-3] - out[-2])
        out[-1] = out[-2] + (out[2] - out[1])
    elif edges == EdgeTreatment.CLAMP:
        out[0] = out[1]
        out[-1] = out[-2]
    else:
        out[0] = out[-1] = missing_value
    return out


def expand_2d(
    source: npt.ArrayLike,
    *,
    dtype: Optional[npt.DTypeLike] = None,
    edges_x: EdgeTreatment = EdgeTreatment.MISSING,
    edges_y: EdgeTreatment = EdgeTreatment.MISSING,
    missing_value=np.nan,
) -> np.ndarray:
    # Create an array to hold data at centers (T points),
    # with strips of size 1 on all sides to support interpolation to interfaces
    source_shape = np.shape(source)
    out_shape = source_shape[:-2] + (source_shape[-2] + 2, source_shape[-1] + 2)
    out = np.empty_like(source, shape=out_shape, dtype=dtype)
    out[..., 1:-1, 1:-1] = source

    edges_x = EdgeTreatment(edges_x)
    edges_y = EdgeTreatment(edges_y)

    if edges_x == EdgeTreatment.EXTRAPOLATE:
        out[..., 0] = 2 * out[..., 1] - out[..., 2]
        out[..., -1] = 2 * out[..., -2] - out[..., -3]
    elif edges_x == EdgeTreatment.EXTRAPOLATE_PERIODIC:
        out[..., 0] = out[..., 1] + (out[..., -3] - out[..., -2])
        out[..., -1] = out[..., -2] + (out[..., 2] - out[..., 1])
    elif edges_x == EdgeTreatment.CLAMP:
        out[..., 0] = out[..., 1]
        out[..., -1] = out[..., -2]
    elif edges_x == EdgeTreatment.PERIODIC:
        out[..., 0] = out[..., -2]
        out[..., -1] = out[..., 1]
    elif edges_x == EdgeTreatment.MISSING:
        out[..., 0] = missing_value
        out[..., -1] = missing_value

    if edges_y == EdgeTreatment.EXTRAPOLATE:
        out[..., 0, :] = 2 * out[..., 1, :] - out[..., 2, :]
        out[..., -1, :] = 2 * out[..., -2, :] - out[..., -3, :]
    elif edges_y == EdgeTreatment.EXTRAPOLATE_PERIODIC:
        out[..., 0, :] = out[..., 1, :] + (out[..., -3, :] - out[..., -2, :])
        out[..., -1, :] = out[..., -2, :] + (out[..., 2, :] - out[..., 1, :])
    elif edges_y == EdgeTreatment.CLAMP:
        out[..., 0, :] = out[..., 1, :]
        out[..., -1, :] = out[..., -2, :]
    elif edges_y == EdgeTreatment.PERIODIC:
        out[..., 0, :] = out[..., -2, :]
        out[..., -1, :] = out[..., 1, :]
    elif edges_y == EdgeTreatment.MISSING:
        out[..., 0, :] = missing_value
        out[..., -1, :] = missing_value

    return out


def centers_to_supergrid_2d(
    source: npt.ArrayLike,
    *,
    dtype: Optional[npt.DTypeLike] = None,
    edges_x: EdgeTreatment = EdgeTreatment.MISSING,
    edges_y: EdgeTreatment = EdgeTreatment.MISSING,
    missing_value=np.nan,
) -> np.ndarray:
    source = np.asanyarray(source)
    ny, nx = source.shape

    source_ex = expand_2d(
        source,
        dtype=dtype,
        edges_x=edges_x,
        edges_y=edges_y,
        missing_value=missing_value,
    )

    out_shape = source_ex.shape[:-2] + (ny * 2 + 1, nx * 2 + 1)
    out = np.empty_like(source_ex, shape=out_shape)

    if_ip_shape = (4,) + source_ex.shape[:-2] + (ny + 1, nx + 1)
    data_if_ip = np.empty_like(source_ex, shape=if_ip_shape)
    data_if_ip[0, ...] = source_ex[..., :-1, :-1]
    data_if_ip[1, ...] = source_ex[..., 1:, :-1]
    data_if_ip[2, ...] = source_ex[..., :-1, 1:]
    data_if_ip[3, ...] = source_ex[..., 1:, 1:]
    out[..., 1::2, 1::2] = source_ex[1:-1, 1:-1]  # T points
    out[..., ::2, ::2] = data_if_ip.mean(axis=0)  # X points

    data_if_ip[0, ..., :-1] = source_ex[..., :-1, 1:-1]
    data_if_ip[1, ..., :-1] = source_ex[..., 1:, 1:-1]
    out[..., ::2, 1::2] = data_if_ip[:2, ..., :-1].mean(axis=0)  # V points

    data_if_ip[0, ..., :-1, :] = source_ex[..., 1:-1, :-1]
    data_if_ip[1, ..., :-1, :] = source_ex[..., 1:-1, 1:]
    out[..., 1::2, ::2] = data_if_ip[:2, ..., :-1, :].mean(axis=0)  # U points

    return np.ma.filled(out, missing_value)


def interfaces_to_supergrid_1d(
    source: npt.ArrayLike,
    *,
    dtype: Optional[npt.DTypeLike] = None,
    out: Optional[np.ndarray] = None,
) -> np.ndarray:
    source = np.asarray(source)
    assert source.ndim == 1, "data must be one-dimensional"
    assert source.size > 1, "data must have at least 2 elements"
    if out is None:
        out = np.empty_like(source, shape=(source.size * 2 - 1,), dtype=dtype)
    out[0::2] = source
    out[1::2] = 0.5 * (source[1:] + source[:-1])
    return out


def interfaces_to_supergrid_2d(
    source: npt.ArrayLike,
    *,
    dtype: Optional[npt.DTypeLike] = None,
    out: Optional[np.ndarray] = None,
) -> np.ndarray:
    source = np.asarray(source)
    assert source.ndim == 2, "data must be two-dimensional"
    assert (
        source.shape[0] > 1 and source.shape[1] > 1
    ), "dimensions must have length >= 2"
    if out is None:
        out_shape = (source.shape[0] * 2 - 1, source.shape[1] * 2 - 1)
        out = np.empty_like(source, shape=out_shape, dtype=dtype)
    out[0::2, 0::2] = source
    out[1::2, 0::2] = 0.5 * (source[:-1, :] + source[1:, :])
    out[0::2, 1::2] = 0.5 * (source[:, :-1] + source[:, 1:])
    out[1::2, 1::2] = 0.25 * (
        source[:-1, :-1] + source[:-1, 1:] + source[1:, :-1] + source[1:, 1:]
    )
    return out


def create_cartesian(
    x: npt.ArrayLike,
    y: npt.ArrayLike,
    *,
    interfaces=False,
    central_lon: Optional[float] = None,
    central_lat: Optional[float] = None,
    **kwargs,
) -> "Domain":
    """Create Cartesian domain from x and y coordinates.

    Args:
        x: array with x coordinates (1d or 2d)
            (at cell interfaces if `interfaces=True`, else at cell centers)
        y: array with y coordinates (1d or 2d)
            (at cell interfaces if `interfaces=True`, else at cell centers)
        interfaces: coordinates are given at cell interfaces, rather than cell centers.
        **kwargs: additional arguments passed to :class:`Domain`
    """
    x = np.asarray(x)
    y = np.asarray(y)
    if y.ndim == 1:
        y = y[:, np.newaxis]
    if interfaces:
        nx, ny = x.shape[-1] - 1, y.shape[0] - 1
    else:
        nx, ny = x.shape[-1], y.shape[0]
    if central_lon is not None and central_lat is not None:
        relx = x - 0.5 * (x.min() + x.max())
        rely = y - 0.5 * (y.min() + y.max())
        m_per_degree = DEG2RAD * R_EARTH
        lats = rely / m_per_degree + central_lat
        lons = relx / m_per_degree / np.cos(DEG2RAD * lats) + central_lon
        lons, lats = np.broadcast_arrays(lons, lats)
        kwargs.setdefault("lon", lons)
        kwargs.setdefault("lat", lats)
    return Domain(nx, ny, x=x, y=y, coordinate_type=CoordinateType.XY, **kwargs)


def create_spherical(
    lon: npt.ArrayLike, lat: npt.ArrayLike, *, interfaces=False, **kwargs
) -> "Domain":
    """Create spherical domain from longitudes and latitudes.

    Args:
        lon: array with longitude coordinates (1d or 2d)
            (at cell interfaces if `interfaces=True`, else at cell centers)
        lat: array with latitude coordinates (1d or 2d)
            (at cell interfaces if `interfaces=True`, else at cell centers)
        interfaces: coordinates are given at cell interfaces, rather than cell centers.
        **kwargs: additional arguments passed to :class:`Domain`
    """
    lon = np.asarray(lon)
    lat = np.asarray(lat)
    if lat.ndim == 1:
        lat = lat[:, np.newaxis]

    if interfaces:
        nx, ny = lon.shape[-1] - 1, lat.shape[0] - 1
    else:
        nx, ny = lon.shape[-1], lat.shape[0]
    return Domain(
        nx, ny, lon=lon, lat=lat, coordinate_type=CoordinateType.LONLAT, **kwargs
    )


def create_spherical_at_resolution(
    minlon: float,
    maxlon: float,
    minlat: float,
    maxlat: float,
    resolution: float,
    **kwargs,
) -> "Domain":
    """Create spherical domain encompassing the specified longitude range and latitude
    range and desired resolution in m.

    Args:
        minlon: minimum longitude (°East)
        maxlon: maximum longitude (°East)
        minlat: minimum latitude (°North)
        maxlat: maximum latitude (°North)
        resolution: maximum grid cell length and width (m)
        **kwargs: additional arguments passed to :class:`Domain`
    """
    if maxlon <= minlon:
        raise Exception(
            f"Maximum longitude {maxlon} must exceed minimum longitude {minlon}"
        )
    if maxlat <= minlat:
        raise Exception(
            f"Maximum latitude {maxlat} must exceed minimum latitude {minlat}"
        )
    if resolution <= 0.0:
        raise Exception(f"Desired resolution must exceed 0, but is {resolution} m")
    dlat = resolution / (DEG2RAD * R_EARTH)
    minabslat = min(abs(minlat), abs(maxlat))
    dlon = resolution / (DEG2RAD * R_EARTH) / np.cos(DEG2RAD * minabslat)
    nx = int(np.ceil((maxlon - minlon) / dlon)) + 1
    ny = int(np.ceil((maxlat - minlat) / dlat)) + 1
    return create_spherical(
        np.linspace(minlon, maxlon, nx),
        np.linspace(minlat, maxlat, ny),
        interfaces=True,
        **kwargs,
    )


class DomainArray:
    def __init__(self, writeable: bool = True, **kwargs):
        self.writable = writeable
        self.kwargs = kwargs

    def __set_name__(self, owner, name):
        self.private_name = f"_{name}"

    def __get__(self, domain: "Domain", objtype=None) -> Optional[np.ndarray]:
        values = getattr(domain, self.private_name, None)
        if values is not None:
            if values.flags.writeable and not self.writable:
                values = values.view()
                values.flags.writeable = False
            elif not values.flags.writeable and self.writable:
                values = values.copy()
                setattr(domain, self.private_name, values)
        return values

    def __set__(self, domain: "Domain", values: Optional[npt.ArrayLike]):
        assert self.writable
        if domain.comm.rank == 0:
            values = domain._map_array(values, **self.kwargs)
            setattr(domain, self.private_name, values)


def calculate_and_bcast(method):
    @functools.wraps(method)
    def wrapper(self: "Domain", *args, **kwargs):
        if self.comm.rank == 0:
            result = method(self, *args, **kwargs)
        else:
            result = None
        return self.comm.bcast(result)

    return wrapper


def _rotation(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    # For each point, draw lines to the nearest neighbor (1/2 a grid cell) on
    # the left, right, top and bottom.
    rot_left = np.arctan2(y[:, 1:-1] - y[:, :-2], x[:, 1:-1] - x[:, :-2])
    rot_right = np.arctan2(y[:, 2:] - y[:, 1:-1], x[:, 2:] - x[:, 1:-1])
    rot_bot = np.arctan2(y[1:-1, :] - y[:-2, :], x[1:-1, :] - x[:-2, :]) - 0.5 * np.pi
    rot_top = np.arctan2(y[2:, :] - y[1:-1, :], x[2:, :] - x[1:-1, :]) - 0.5 * np.pi
    x_dum = (
        np.cos(rot_left[1:-1, :])
        + np.cos(rot_right[1:-1, :])
        + np.cos(rot_bot[:, 1:-1])
        + np.cos(rot_top[:, 1:-1])
    )
    y_dum = (
        np.sin(rot_left[1:-1, :])
        + np.sin(rot_right[1:-1, :])
        + np.sin(rot_bot[:, 1:-1])
        + np.sin(rot_top[:, 1:-1])
    )
    return np.arctan2(y_dum, x_dum)


class Domain:
    # Grid metrics frozen at initialization
    x = DomainArray(writeable=False)
    y = DomainArray(writeable=False)
    lon = DomainArray(writeable=False)
    lat = DomainArray(writeable=False)
    f = DomainArray(writeable=False)
    dx = DomainArray(writeable=False)
    dy = DomainArray(writeable=False)
    rotation = DomainArray(writeable=False)
    area = DomainArray(writeable=False)

    # Grid metrics that can be manipulated after the domain is created
    mask = DomainArray(missing_value=0, dtype=int)
    H = DomainArray()
    z0 = DomainArray()

    def __init__(
        self,
        nx: int,
        ny: int,
        *,
        lon: Optional[npt.ArrayLike] = None,
        lat: Optional[npt.ArrayLike] = None,
        x: Optional[npt.ArrayLike] = None,
        y: Optional[npt.ArrayLike] = None,
        coordinate_type: Optional[CoordinateType] = None,
        mask: Optional[npt.ArrayLike] = 1,
        H: Optional[npt.ArrayLike] = None,
        z0: Optional[npt.ArrayLike] = 0.0,
        f: Optional[npt.ArrayLike] = None,
        periodic_x: bool = False,
        periodic_y: bool = False,
        comm: Optional[parallel.MPI.Comm] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Create new domain.

        Args:
            nx: number of tracer points in x-direction
            ny: number of tracer points in y-direction
            lon: longitude (°East)
            lat: latitude (°North)
            x: x coordinate (m)
            y: y coordinate (m)
            coordinate_type: preferred coordinate type for plots and output
            mask: initial mask (0: land, 1: water)
            H: initial bathymetric depth. This is the distance between the bottom and
                some arbitrary depth reference (m, positive if bottom lies below the
                depth reference). Typically the depth reference is mean sea level.
            z0: minimum hydrodynamic bottom roughness (m)
            f: Coriolis parameter. If not provided, it will be calculated from latitude.
                In that case argument `lat` must be provided.
            periodic_x: use periodic boundary in x-direction (left == right)
            periodic_y: use periodic boundary in y-direction (top == bottom)
            comm: MPI communicator that comprises all processes that should get access
                to the domain
            logger: logger for diagnostic messages
        """
        if nx <= 0:
            raise Exception(f"Number of x points is {nx} but must be > 0")
        if ny <= 0:
            raise Exception(f"Number of y points is {ny} but must be > 0")

        self.comm = comm or parallel.mpi4py_autofree(parallel.MPI.COMM_WORLD.Dup())
        has_xy = self.comm.bcast(x is not None and y is not None)
        has_lonlat = self.comm.bcast(lon is not None and lat is not None)
        if not (has_xy or has_lonlat):
            raise Exception("Either x and y, or lon and lat, must be provided")
        has_f = self.comm.bcast(f is not None or lat is not None)
        if not has_f:
            raise Exception(
                "Either lat of f must be provided to determine the Coriolis parameter."
            )
        if coordinate_type is None:
            coordinate_type = CoordinateType.XY if has_xy else CoordinateType.LONLAT
        assert (coordinate_type == CoordinateType.XY and has_xy) or (
            coordinate_type == CoordinateType.LONLAT
            and has_lonlat
            or (coordinate_type == CoordinateType.IJ and (has_xy or has_lonlat))
        )

        self.nx = nx
        self.ny = ny
        self.periodic_x = periodic_x
        self.periodic_y = periodic_y
        self.coordinate_type = coordinate_type
        self.root_logger = logger or parallel.get_logger()
        self.logger = self.root_logger.getChild("domain")

        self.open_boundaries = open_boundaries.OpenBoundaries(
            nx, ny, self.logger.getChild("open_boundaries")
        )
        self.rivers = rivers.Rivers(
            nx, ny, coordinate_type, self.logger.getChild("rivers")
        )
        self.default_output_transforms = []
        self.extra_output_coordinates = []
        self.input_grid_mappers = []

        if self.comm.rank != 0:
            return

        self._x = self._map_array(x, edges=EdgeTreatment.EXTRAPOLATE)
        self._y = self._map_array(y, edges=EdgeTreatment.EXTRAPOLATE)
        if lon is not None and np.shape(lon) != (1 + ny * 2, 1 + nx * 2):
            # Interpolate longitude in cos-sin space to handle periodic boundary
            # condition, but skip this if longitude is already provided on supergrid
            # (no interpolation needed) to improve accuracy of rotation tests.
            lon = np.asarray(lon).astype(float, copy=False)
            lon_rad = DEG2RAD * lon
            coslon = np.cos(lon_rad)
            sinlon = np.sin(lon_rad)
            coslon = self._map_array(coslon, edges=EdgeTreatment.EXTRAPOLATE)
            sinlon = self._map_array(sinlon, edges=EdgeTreatment.EXTRAPOLATE)
            self._lon = np.arctan2(sinlon, coslon) * RAD2DEG
            central_lon = 0.5 * (lon.min() + lon.max())
            wrap_lon = central_lon - 180.0
            self._lon = (self._lon - wrap_lon) % 360.0 + wrap_lon
        else:
            self._lon = self._map_array(lon)
        self._lat = self._map_array(lat, edges=EdgeTreatment.EXTRAPOLATE)
        self._f = self._map_array(f)
        self._mask = self._map_array(mask, missing_value=0, dtype=int)
        self._H = self._map_array(H, edges=EdgeTreatment.CLAMP)
        self._z0 = self._map_array(z0)

        kwargs_expand = {}
        if self.periodic_x:
            kwargs_expand["edges_x"] = EdgeTreatment.EXTRAPOLATE_PERIODIC
        if self.periodic_y:
            kwargs_expand["edges_y"] = EdgeTreatment.EXTRAPOLATE_PERIODIC
        if has_xy:
            # Expand x, y by 1 in each direction to calculate dx, dy, rotation
            x_ex = expand_2d(self._x, **kwargs_expand)
            y_ex = expand_2d(self._y, **kwargs_expand)
        if has_lonlat:
            # Expand lon, lat by 1 in each direction to calculate dx, dy, rotation
            lon_rad_ex = DEG2RAD * expand_2d(self._lon, **kwargs_expand)
            lat_rad_ex = DEG2RAD * expand_2d(self._lat, **kwargs_expand)

        if has_xy:
            dx_x = x_ex[1:-1, 2:] - x_ex[1:-1, :-2]
            dy_x = y_ex[1:-1, 2:] - y_ex[1:-1, :-2]
            dx_y = x_ex[2:, 1:-1] - x_ex[:-2, 1:-1]
            dy_y = y_ex[2:, 1:-1] - y_ex[:-2, 1:-1]
            scale = 1.0
        else:
            dlon_rad_x = lon_rad_ex[1:-1, 2:] - lon_rad_ex[1:-1, :-2]
            dlat_rad_x = lat_rad_ex[1:-1, 2:] - lat_rad_ex[1:-1, :-2]
            coslat_x = np.cos(0.5 * (lat_rad_ex[1:-1, 2:] + lat_rad_ex[1:-1, :-2]))
            dx_x = coslat_x * np.sin(0.5 * dlon_rad_x)
            dy_x = np.sin(0.5 * dlat_rad_x) * np.cos(0.5 * dlon_rad_x)
            dlon_rad_y = lon_rad_ex[2:, 1:-1] - lon_rad_ex[:-2, 1:-1]
            dlat_rad_y = lat_rad_ex[2:, 1:-1] - lat_rad_ex[:-2, 1:-1]
            coslat_y = np.cos(0.5 * (lat_rad_ex[2:, 1:-1] + lat_rad_ex[:-2, 1:-1]))
            dx_y = coslat_y * np.sin(0.5 * dlon_rad_y)
            dy_y = np.sin(0.5 * dlat_rad_y) * np.cos(0.5 * dlon_rad_y)
            scale = R_EARTH * 2.0
        self._dx = scale * np.hypot(dx_x, dy_x)
        self._dy = scale * np.hypot(dx_y, dy_y)

        if has_lonlat and not (
            (self._lat == self._lat[0, 0]).any()
            and (self._lon == self._lon[0, 0]).any()
            and has_xy
        ):
            # Proper rotation with respect to true North
            rotation = _rotation(lon_rad_ex, lat_rad_ex)
        else:
            # Rotation with respect to y-axis - assumes y-axis always points to
            # true North (can be valid only for infinitesimally small domain)
            rotation = _rotation(x_ex, y_ex)
        self._rotation = rotation if rotation[1:-1, 1:-1].any() else None

        self._area = self._dx * self._dy

    def _map_array(
        self,
        values: Optional[npt.ArrayLike],
        *,
        dtype: npt.DTypeLike = float,
        edges: EdgeTreatment = EdgeTreatment.MISSING,
        missing_value=np.nan,
    ) -> Optional[np.ndarray]:
        if self.comm.rank != 0 or values is None:
            return None

        source_shape = np.shape(values)
        source_shape = (1,) * (2 - len(source_shape)) + source_shape  # broadcast
        target_shape = (self.ny * 2 + 1, self.nx * 2 + 1)

        def can_cast(*target_shape: int) -> bool:
            assert len(target_shape) == len(source_shape)
            return all((l == 1 or l == lr) for l, lr in zip(source_shape, target_shape))

        if source_shape[0] == 1 and source_shape[1] == 1:
            # scalar value
            mapped_values = np.array(values, dtype=dtype)
        elif can_cast(self.ny, self.nx):
            # values provided at cell centers
            edges_x = edges_y = edges
            if self.periodic_x and edges_x != EdgeTreatment.EXTRAPOLATE:
                edges_x = EdgeTreatment.PERIODIC
            if self.periodic_y and edges_y != EdgeTreatment.EXTRAPOLATE:
                edges_y = EdgeTreatment.PERIODIC
            if source_shape[0] == 1:
                values_sup = centers_to_supergrid_1d(
                    np.ravel(values),
                    edges=edges_x,
                    missing_value=missing_value,
                    dtype=dtype,
                )
                mapped_values = values_sup[np.newaxis, :]
            elif source_shape[1] == 1:
                values_sup = centers_to_supergrid_1d(
                    np.ravel(values),
                    edges=edges_y,
                    missing_value=missing_value,
                    dtype=dtype,
                )
                mapped_values = values_sup[:, np.newaxis]
            else:
                mapped_values = centers_to_supergrid_2d(
                    values,
                    edges_x=edges_x,
                    edges_y=edges_y,
                    missing_value=missing_value,
                    dtype=dtype,
                )
        elif can_cast(self.ny + 1, self.nx + 1):
            # values provided at cell corners
            if source_shape[0] == 1:
                values_sup = interfaces_to_supergrid_1d(np.ravel(values), dtype=dtype)
                mapped_values = values_sup[np.newaxis, :]
            elif source_shape[1] == 1:
                values_sup = interfaces_to_supergrid_1d(np.ravel(values), dtype=dtype)
                mapped_values = values_sup[:, np.newaxis]
            else:
                mapped_values = interfaces_to_supergrid_2d(values)
        else:
            # values provided on supergrid
            assert can_cast(
                *target_shape
            ), f"Cannot map array with shape {values.shape} to supergrid with shape {target_shape}"
            mapped_values = np.array(values, dtype=dtype)
        if mapped_values.shape != target_shape:
            mapped_values = np.broadcast_to(mapped_values, target_shape)
        return mapped_values

    def create_tiling(self) -> parallel.Tiling:
        mask = None if self.comm.rank != 0 else self._mask[1::2, 1::2]
        mask = self.comm.bcast(mask)
        return parallel.Tiling.autodetect(
            mask,
            periodic_x=self.periodic_x,
            periodic_y=self.periodic_y,
            comm=self.comm,
            logger=self.logger.getChild("subdomain_decomposition"),
        )

    def create_grids(
        self,
        nz: int,
        halox: int,
        haloy: int,
        fields: Optional[Mapping[str, core.Array]] = None,
        tiling: Optional[parallel.Tiling] = None,
        input_manager=None,
        velocity_grids: int = 0,
        t_postfix: str = "",
    ) -> core.Grid:
        if self.comm.rank == 0:
            if self._H is None:
                raise Exception("Water depth at rest (H) has not been provided")
            # We are the root node - update the global mask
            self.infer_UVX_masks2()
            self.open_boundaries.adjust_mask(self.mask)

            # Map river coordinates to global grid indices
            self._map_rivers()

        if tiling is None:
            tiling = self.create_tiling()
        elif tiling.nx_glob is None:
            tiling.set_extent(self.nx, self.ny)

        # NB the fields argument cannot have a default of {}, as that causes
        # that global dictionary to be shared among all calls to create_grids.
        if fields is None:
            fields = {}

        def create_grid(
            postfix: str,
            ioffset: int,
            joffset: int,
            *,
            overlap: int = 0,
            **kwargs,
        ) -> core.Grid:
            grid = core.Grid(
                tiling.nx_sub + overlap,
                tiling.ny_sub + overlap,
                nz,
                halox=halox,
                haloy=haloy,
                postfix=postfix,
                fields=fields,
                tiling=tiling,
                ioffset=ioffset,
                joffset=joffset,
                overlap=overlap,
                **kwargs,
            )
            self._populate_grid(grid)
            grid.input_manager = input_manager
            grid.default_output_transforms = self.default_output_transforms
            grid.extra_output_coordinates = self.extra_output_coordinates
            grid.input_grid_mappers = [m(grid=grid) for m in self.input_grid_mappers]
            return grid

        U = V = X = UU = UV = VU = VV = None

        if velocity_grids > 1:
            UU = create_grid("_uu_adv", 3, 1)
            UV = create_grid("_uv_adv", 2, 2)
            VU = create_grid("_vu_adv", 2, 2)
            VV = create_grid("_vv_adv", 1, 3)

        if velocity_grids > 0:
            U = create_grid("u", 2, 1, ugrid=UU, vgrid=UV)
            V = create_grid("v", 1, 2, ugrid=VU, vgrid=VV)
            X = create_grid("x", 0, 0, overlap=1, istart=0, jstart=0)

        T = create_grid(t_postfix, 1, 1, ugrid=U, vgrid=V, xgrid=X)
        T.rivers = self.rivers

        if velocity_grids > 0:
            T.infer_water_contact()
            T.close_flux_interfaces()

        if velocity_grids > 1:
            U.close_flux_interfaces()
            V.close_flux_interfaces()

            # No transport between velocity points along an open boundary (just outside)
            # This is done to state that no valid values (in e.g. h and D) are required
            # in these points.
            UV.mask.all_values[:-1, :][
                (U.mask.all_values[:-1, :] == 4) & (U.mask.all_values[1:, :] == 4)
            ] = 0
            VU.mask.all_values[:, :-1][
                (V.mask.all_values[:, :-1] == 4) & (V.mask.all_values[:, 1:] == 4)
            ] = 0

        T.freeze()

        self.open_boundaries.initialize(T)
        self.rivers.initialize(T)

        return T

    def _populate_grid(self, grid: core.Grid):
        NAMEMAP = dict(z0="z0b_min", f="cor")

        edges_x = EdgeTreatment.PERIODIC if self.periodic_x else EdgeTreatment.MISSING
        edges_y = EdgeTreatment.PERIODIC if self.periodic_y else EdgeTreatment.MISSING

        ioffset = grid.ioffset % 2
        joffset = grid.joffset % 2
        istart_glob = (ioffset - grid.ioffset) // 2
        jstart_glob = (joffset - grid.joffset) // 2
        assert istart_glob <= 0, "Supergrid starts after first x"
        assert jstart_glob <= 0, "Supergrid starts after first y"
        nx_glob = self.nx + grid.overlap
        ny_glob = self.ny + grid.overlap

        def _transfer_variable(source: np.ndarray, target: core.Array):
            sendbuf = None
            if grid.tiling.rank == 0:
                # We are on the root node that has the global values
                if grid.joffset > 1 or grid.ioffset > 1:
                    # UU or VV grid that needs one more strip of 1 cell
                    # beyond the end of the supergrid. By default, that is
                    # left at missing_value, but that is inappropriate for
                    # periodic boundaries that need mirroring. Therefore we
                    # expand the grid properly, any mirroring included.
                    source = expand_2d(
                        source,
                        edges_x=edges_x,
                        edges_y=edges_y,
                        missing_value=target.fill_value,
                    )[1:, 1:]
                global_values = source[joffset::2, ioffset::2]
                istop_glob = istart_glob + global_values.shape[-1]
                jstop_glob = jstart_glob + global_values.shape[-2]
                assert istop_glob >= nx_glob, "Supergrid stops before last x"
                assert jstop_glob >= ny_glob, "Supergrid stops before last y"

                sendbuf = grid.tiling._get_work_array(
                    (grid.tiling.comm.size,) + target.all_values.shape,
                    target.dtype,
                    target.fill_value,
                )
                sendbuf.fill(target.fill_value)
                for irow, icol, rank in parallel._iterate_rankmap(grid.tiling.map):
                    if rank < 0:
                        continue

                    jslice, islice = grid.tiling.subdomain2rawslices(
                        irow,
                        icol,
                        halox_sub=grid.halox,
                        haloy_sub=grid.haloy,
                        share=grid.overlap,
                    )
                    slice_glob, slice_loc = _get_rectangle_overlap(
                        istart_glob,
                        istop_glob,
                        jstart_glob,
                        jstop_glob,
                        islice.start,
                        islice.stop,
                        jslice.start,
                        jslice.stop,
                    )
                    sendbuf[(rank,) + slice_loc] = global_values[slice_glob]

                # Keep a pointer to the full global field
                # This will be used preferentially for full-domain output
                target.attrs["_global_values"] = global_values[
                    -jstart_glob : ny_glob - jstart_glob,
                    -istart_glob : nx_glob - istart_glob,
                ]
                assert target.attrs["_global_values"].shape == (ny_glob, nx_glob)
            else:
                target.attrs["_global_values"] = None

            grid.tiling.comm.Scatter(sendbuf, target.all_values)
            if self.periodic_x or self.periodic_y:
                target.update_halos()

        retrieved_from_domain = set()
        for name in (
            "x",
            "y",
            "lon",
            "lat",
            "f",
            "dx",
            "dy",
            "rotation",
            "area",
            "mask",
            "H",
            "z0",
        ):
            source = getattr(self, f"_{name}", None)
            available = grid.tiling.comm.bcast(source is not None)
            if available:
                target = grid.create_array(NAMEMAP.get(name, name))
                if target is not None:
                    _transfer_variable(source, target)
                    retrieved_from_domain.add(name)

        # If the Coriolis parameter was not set explicitly at domain level,
        # calculate it from latitude
        if "f" not in retrieved_from_domain:
            grid.create_array("cor").all_values[...] = coriolis(grid._lat.all_values)

        # Set default horizontal coordinates (e.g., for output and online plotting)
        # based on the coordinate type set at domain level.
        if self.coordinate_type == CoordinateType.XY:
            grid.horizontal_coordinates += [grid.x, grid.y]
        elif self.coordinate_type == CoordinateType.LONLAT:
            grid.horizontal_coordinates += [grid.lon, grid.lat]

    @calculate_and_bcast
    def cfl_check(
        self, z: float = 0.0, return_location: bool = False
    ) -> Union[float, Tuple[float, int, int, float]]:
        """Determine maximum time step (s) for depth-integrated equations

        Args:
            z: maximum surface elevation (m)
            return_location: whether to also return the location
                and depth that determined the maximum step

        Note: this returns global indices for the T grid, not the supergrid
        """
        dx = self._dx[1::2, 1::2]
        dy = self._dy[1::2, 1::2]
        D = self._H[1::2, 1::2] + z
        mask = (self._mask[1::2, 1::2] > 0) & (D > 0.0)
        denom2 = (2.0 * GRAVITY) * D * (dx**2 + dy**2)
        maxdts = dx * dy / np.sqrt(denom2, where=mask, out=np.ones_like(D))
        maxdts[~mask] = np.inf
        maxdt = maxdts.min()
        if return_location:
            j, i = np.unravel_index(np.argmin(maxdts), maxdts.shape)
            return (maxdt, i, j, self._H[1 + 2 * j, 1 + 2 * i])
        return maxdt

    @property
    def maxdt(self) -> float:
        """Maximum time step (s) for depth-integrated equations"""
        return self.cfl_check()

    # @calculate_and_bcast
    def infer_UVX_masks2(self):
        tmask = self.mask[1::2, 1::2]
        umask = self.mask[1::2, ::2]
        vmask = self.mask[::2, 1::2]
        xmask = self.mask[::2, ::2]

        edges_x = EdgeTreatment.PERIODIC if self.periodic_x else EdgeTreatment.MISSING
        edges_y = EdgeTreatment.PERIODIC if self.periodic_y else EdgeTreatment.MISSING
        tmask_ex = expand_2d(tmask, edges_x=edges_x, edges_y=edges_y, missing_value=0)

        # Now mask U,V,X points unless all their T neighbors are valid - this mask will
        # be sent to Fortran and determine which points are computed
        bad = (tmask_ex[1:-1, 1:] == 0) | (tmask_ex[1:-1, :-1] == 0)
        umask[bad & (umask == 1)] = 0
        bad = (tmask_ex[1:, 1:-1] == 0) | (tmask_ex[:-1, 1:-1] == 0)
        vmask[bad & (vmask == 1)] = 0
        bad = (
            (tmask_ex[1:, 1:] == 0)
            | (tmask_ex[:-1, 1:] == 0)
            | (tmask_ex[1:, :-1] == 0)
            | (tmask_ex[:-1, :-1] == 0)
        )
        xmask[bad] = 0

    @calculate_and_bcast
    def mask_shallow(self, minimum_depth: float):
        """Mask all points shallower less the specified value.

        Args:
            minimum_depth: minimum bathmetric depth :attr:`H`; points that are
                shallower will be masked
        """
        self.mask[self._H < minimum_depth] = 0

    @calculate_and_bcast
    def mask_subbasins(self, nkeep: int = 1):
        """Identify all separate basins (each a collection of unmasked cell
        centers connected via top/bottom/left/right interfaces), and mask all
        but the largest one(s).

        Args:
            nkeep: number of basins to keep
        """
        import skimage.segmentation

        tmask = np.minimum(self.mask[1::2, 1::2], 1)
        next_id = -1
        basin2size = {}
        while True:
            indices = (tmask == 1).nonzero()
            if indices[0].size == 0:
                break
            seed_point = (indices[0][0], indices[1][0])
            skimage.segmentation.flood_fill(
                tmask, seed_point, next_id, connectivity=1, in_place=True
            )
            basin2size[next_id] = (tmask == next_id).sum()
            next_id -= 1

        ordered = sorted(basin2size.keys(), key=lambda x: basin2size[x], reverse=True)
        for v in ordered[nkeep:]:
            tmask[tmask == v] = 0
        self.mask = np.where(tmask == 0, 0, self.mask[1::2, 1::2])

    @calculate_and_bcast
    def limit_velocity_depth(self, critical_depth: float = np.inf):
        """Decrease bathymetric depth of velocity (U, V) points to the minimum of the
        bathymetric depth of both neighboring T points, wherever one of these two
        points is shallower than the specified critical depth.

        Args:
            critical_depth: neighbor depth at which the limiting starts. If either
                neighbor (T grid) is shallower than this value, the depth of velocity
                point (U or V grid) is restricted.
        """
        # NB this is guaranteed to return a writeable H,
        # even if self.H was read-only before
        H = self.H

        tdepth = H[1::2, 1::2]
        Vsel = (tdepth[1:, :] <= critical_depth) | (tdepth[:-1, :] <= critical_depth)
        np.putmask(H[2:-2:2, 1::2], Vsel, np.minimum(tdepth[1:, :], tdepth[:-1, :]))
        Usel = (tdepth[:, 1:] <= critical_depth) | (tdepth[:, :-1] <= critical_depth)
        np.putmask(H[1::2, 2:-2:2], Usel, np.minimum(tdepth[:, 1:], tdepth[:, :-1]))
        self.logger.info(
            f"limit_velocity_depth has decreased depth in {Usel.sum()} U points"
            f" ({Usel.sum(where=self._mask[1::2, 2:-2:2] > 0)} currently unmasked),"
            f" {Vsel.sum()} V points ({Vsel.sum(where=self._mask[2:-2:2, 1::2] > 0)}"
            " currently unmasked)."
        )

    @calculate_and_bcast
    def mask_rectangle(
        self,
        xmin: Optional[float] = None,
        xmax: Optional[float] = None,
        ymin: Optional[float] = None,
        ymax: Optional[float] = None,
        mask_value: int = 0,
        coordinate_type: Optional[CoordinateType] = None,
    ):
        """Mask all points that fall within the specified rectangle.

        Args:
            xmin: lower native x coordinate of the rectangle to mask
                (default: left boundary of the domain)
            xmax: upper native x coordinate of the rectangle to mask
                (default: right boundary of the domain)
            ymin: lower native y coordinate of the rectangle to mask
                (default: bottom boundary of the domain)
            ymax: upper native y coordinate of the rectangle to mask
                (default: top boundary of the domain)

        Coordinates will be interpreted as longitude, latitude if the domain is
        configured as spherical; otherwise they will be interpreted as Cartesian
        x and y (m).
        """
        selected = np.ones(self.mask.shape, dtype=bool)
        coordinate_type = coordinate_type or self.coordinate_type
        if coordinate_type == CoordinateType.LONLAT:
            x, y = (self._lon, self._lat)
        elif coordinate_type == CoordinateType.XY:
            x, y = (self._x, self._y)
        else:
            x = np.linspace(-0.5, self.nx + 0.5, 1 + 2 * self.nx)
            y = np.linspace(-0.5, self.ny + 0.5, 1 + 2 * self.ny)
            x, y = np.broadcast_arrays(x[np.newaxis, :], y[:, np.newaxis])
        if xmin is not None:
            selected &= x >= xmin
        if xmax is not None:
            selected &= x <= xmax
        if ymin is not None:
            selected &= y >= ymin
        if ymax is not None:
            selected &= y <= ymax
        self.mask[selected] = mask_value

    @calculate_and_bcast
    def mask_indices(
        self, istart: int, istop: int, jstart: int, jstop: int, mask_value: int = 0
    ):
        """Mask all points that fall within the specified rectangle.
        Indices must be provided for the T grid, and thus range between 0 and `nx`
        for i, and between 0 and `ny` for j.

        Args:
            istart: lower x index (first that is included)
            istop: upper x index (first that is EXcluded)
            jstart: lower y index (first that is included)
            jstop: upper y index (first that is EXcluded)
        """
        istart_T = 1 + 2 * istart
        istop_T = 1 + 2 * istop
        jstart_T = 1 + 2 * jstart
        jstop_T = 1 + 2 * jstop
        self.mask[jstart_T:jstop_T, istart_T:istop_T] = mask_value

    def rotate(self) -> "Domain":
        """Return a copy of the domain rotated 90° clockwise"""

        def tp(array):
            return None if array is None else np.transpose(array)[::-1, :]

        kwargs = dict(
            coordinate_type=self.coordinate_type,
            comm=self.comm,
            logger=self.root_logger,
        )
        if self.comm.rank == 0:
            kwargs.update(
                lon=tp(self._lon),
                lat=tp(self._lat),
                x=tp(self._x),
                y=tp(self._y),
                mask=np.minimum(tp(self._mask), 1),
                H=tp(self._H),
                z0=tp(self._z0),
                f=tp(self._f),
            )
        rotated_domain = Domain(self.ny, self.nx, **kwargs)
        MAP = {
            open_boundaries.Side.WEST: open_boundaries.Side.NORTH,
            open_boundaries.Side.NORTH: open_boundaries.Side.EAST,
            open_boundaries.Side.EAST: open_boundaries.Side.SOUTH,
            open_boundaries.Side.SOUTH: open_boundaries.Side.WEST,
        }
        for b in self.open_boundaries:
            mstart, mstop, l = b.mstart_glob, b.mstop_glob, b.l_glob
            if b.side in (open_boundaries.Side.NORTH, open_boundaries.Side.SOUTH):
                mstart, mstop = (self.nx - 1 - mstart, self.nx - 1 - mstop)
            else:
                l = self.nx - 1 - l
            rotated_domain.open_boundaries.add_by_index(
                MAP[b.side], l, mstart, mstop, type_2d=b.type_2d, type_3d=b.type_3d
            )
        for r in self.rivers.values():
            if r.i_glob is not None and r.j_glob is not None:
                rot_r = rotated_domain.rivers.add_by_index(
                    r.name, r.j_glob, self.nx - 1 - r.i_glob
                )
            for att in ("original_name", "split", "zl", "zu"):
                if hasattr(r, att):
                    setattr(rot_r, att, getattr(r, att))
        rotated_domain.open_boundaries.sponge.tmrlx = self.open_boundaries.sponge.tmrlx
        return rotated_domain

    def _map_rivers(self):
        assert self.comm.rank == 0
        x_ = None if self._x is None else self._x[1::2, 1::2]
        y_ = None if self._y is None else self._y[1::2, 1::2]
        lon_ = None if self._lon is None else self._lon[1::2, 1::2]
        lat_ = None if self._lat is None else self._lat[1::2, 1::2]
        self.rivers.map_to_grid(self._mask[1::2, 1::2], x_, y_, lon_, lat_)

    def plot(
        self,
        field: Optional[np.ndarray] = None,
        *,
        fig: Optional["matplotlib.figure.Figure"] = None,
        show_bathymetry: bool = True,
        show_mask: bool = False,
        show_mesh: bool = False,
        show_rivers: bool = True,
        show_open_boundaries: bool = True,
        show_subdomains: bool = False,
        editable: bool = False,
        coordinate_type: Optional[CoordinateType] = None,
        tiling: Optional[parallel.Tiling] = None,
        label: Optional[str] = None,
        cmap: Union[None, "matplotlib.colors.Colormap", str] = None,
    ) -> Optional["matplotlib.figure.Figure"]:
        """Plot the domain, optionally including bathymetric depth, mesh and
        river positions.

        Args:
            field: 2D field to plot. it must be defined on the supergrid.
            fig: :class:`matplotlib.figure.Figure` instance to plot to.
                If not provided, a new figure is created.
            show_bathymetry: show bathymetry as color map
            show_mask: show mask as color map (this disables ``show_bathymetry``)
            show_mesh: show model grid
            show_rivers: show rivers with position and name
            show_subdomains: show subdomain decompositon
            editable: allow interactive selection of rectangular regions in the domain
                plot that are subsequently masked out
            coordinate_type: coordinates to use for x and y axes (x/y, lon/lat, or i/j)
            tiling: subdomain decomposition. This must be provided if `show_subdomains`
                is set
            label: label for the colorbar of the plotted field
            cmap: colormap to use

        Returns:
            :class:`matplotlib.figure.Figure` instance for processes with rank 0,
            otherwise ``None``
        """
        if self.comm.rank != 0:
            return

        import matplotlib.pyplot
        import matplotlib.collections
        import matplotlib.widgets

        if fig is None:
            fig, ax = matplotlib.pyplot.subplots(
                figsize=(0.15 * self.nx, 0.15 * self.ny)
            )
        else:
            ax = fig.gca()

        if coordinate_type is None:
            coordinate_type = self.coordinate_type
        if coordinate_type == CoordinateType.LONLAT:
            x, y = (self._lon, self._lat)
            xlabel, ylabel = "longitude (°East)", "latitude (°North)"
        elif coordinate_type == CoordinateType.XY:
            x, y = (self._x, self._y)
            xlabel, ylabel = "x (m)", "y (m)"
        else:
            x = -0.5 + 0.5 * np.arange(1 + self.nx * 2)
            y = -0.5 + 0.5 * np.arange(1 + self.ny * 2)
            x, y = np.broadcast_arrays(x, y[:, np.newaxis])
            xlabel, ylabel = "cell index", "cell index"

        if field is None:
            if show_mask:
                field = self._mask
                label = "mask value"
            elif show_bathymetry:
                import cmocean

                cmap = cmocean.cm.deep
                cmap.set_bad("gray")
                field = np.ma.array(self._H, mask=self._mask == 0)
                label = "undisturbed water depth (m)"

        c = ax.pcolormesh(
            x,
            y,
            field,
            alpha=0.5 if show_mesh else 1,
            shading="auto",
            cmap=cmap,
        )
        # c = ax.contourf(x, y, np.ma.array(self.H, mask=self.mask==0), 20, alpha=0.5 if show_mesh else 1)
        cb = fig.colorbar(c)
        if label is not None:
            cb.set_label(label)

        if show_rivers and self.rivers:
            self._map_rivers()
            for river in self.rivers.global_rivers:
                i_sup, j_sup = 1 + river.i_glob * 2, 1 + river.j_glob * 2
                river_x, river_y = x[j_sup, i_sup], y[j_sup, i_sup]
                ax.plot([river_x], [river_y], ".r")
                ax.text(river_x, river_y, river.name, color="r")

        if show_open_boundaries:
            x_X, y_X = x[::2, ::2], y[::2, ::2]
            for b in self.open_boundaries:
                if b.side in (open_boundaries.Side.WEST, open_boundaries.Side.EAST):
                    imin, imax = b.l_glob, b.l_glob + 1
                    jmin = min(b.mstart_glob, b.mstop_glob - b.mstep)
                    jmax = max(b.mstart_glob, b.mstop_glob - b.mstep) + 1
                else:
                    jmin, jmax = b.l_glob, b.l_glob + 1
                    imin = min(b.mstart_glob, b.mstop_glob - b.mstep)
                    imax = max(b.mstart_glob, b.mstop_glob - b.mstep) + 1
                x_b = x_X[jmin : jmax + 1, imin : imax + 1]
                y_b = y_X[jmin : jmax + 1, imin : imax + 1]
                if x_b.shape[1] == 2:
                    x_b, y_b = x_b.T, y_b.T
                x_b = np.hstack((x_b[0, :], x_b[1, ::-1]))
                y_b = np.hstack((y_b[0, :], y_b[1, ::-1]))
                ax.fill(x_b, y_b, alpha=0.3, fc="r", ec="None")
                ax.fill(x_b, y_b, fc="None", ec="r")

        def plot_mesh(ax, x, y, **kwargs):
            segs1 = np.stack((x, y), axis=2)
            segs2 = segs1.transpose(1, 0, 2)
            ax.add_collection(matplotlib.collections.LineCollection(segs1, **kwargs))
            ax.add_collection(matplotlib.collections.LineCollection(segs2, **kwargs))

        if show_mesh:
            plot_mesh(
                ax, x[::2, ::2], y[::2, ::2], colors="k", linestyle="-", linewidth=0.3
            )
            # ax.pcolor(x[1::2, 1::2], y[1::2, 1::2], np.ma.array(x[1::2, 1::2], mask=True), edgecolors='k', linestyles='--', linewidth=.2)
            # pc = ax.pcolormesh(x[1::2, 1::2], y[1::2, 1::2],  np.ma.array(x[1::2, 1::2], mask=True), edgecolor='gray', linestyles='--', linewidth=.2)
            ax.plot(x[::2, ::2], y[::2, ::2], "xk", markersize=3.0)
            ax.plot(x[1::2, 1::2], y[1::2, 1::2], ".k", markersize=2.5)

        if show_subdomains:
            assert tiling is not None
            for icol in range(tiling.ncol + 1):
                i = 2 * (tiling.xoffset_global + icol * tiling.nx_sub)
                if i >= 0 and i < x.shape[-1]:
                    ax.plot(x[:, i], y[:, i], "-k", linewidth=2.0)
                    ax.plot(x[:, i], y[:, i], "--w", linewidth=2.0)
            for irow in range(tiling.nrow + 1):
                j = 2 * (tiling.yoffset_global + irow * tiling.ny_sub)
                if j >= 0 and j < x.shape[-2]:
                    ax.plot(x[j, :], y[j, :], "-k", linewidth=2.0)
                    ax.plot(x[j, :], y[j, :], "--w", linewidth=2.0)

            # Rank of each subdomain (cross for subdomains that are not used, e.g. land)
            for irow in range(tiling.nrow):
                for icol in range(tiling.ncol):
                    imin = 2 * (tiling.xoffset_global + icol * tiling.nx_sub)
                    imax = imin + 2 * tiling.nx_sub
                    jmin = 2 * (tiling.yoffset_global + irow * tiling.ny_sub)
                    jmax = jmin + 2 * tiling.ny_sub
                    imin, imax = max(0, imin), min(x.shape[-1] - 1, imax)
                    jmin, jmax = max(0, jmin), min(x.shape[-2] - 1, jmax)
                    imid = (imin + imax) // 2
                    jmid = (jmin + jmax) // 2
                    rank = tiling.map[irow, icol]
                    if rank >= 0:
                        ax.text(
                            x[jmid, imid],
                            y[jmid, imid],
                            str(rank),
                            horizontalalignment="center",
                            verticalalignment="center",
                        )
                    else:
                        ax.plot(
                            [x[jmin, imin], x[jmax, imax]],
                            [y[jmin, imin], y[jmax, imax]],
                            "-k",
                        )
                        ax.plot(
                            [x[jmin, imax], x[jmax, imin]],
                            [y[jmin, imax], y[jmax, imin]],
                            "-k",
                        )

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if coordinate_type != CoordinateType.LONLAT:
            ax.axis("equal")
        xmin, xmax = np.nanmin(x), np.nanmax(x)
        ymin, ymax = np.nanmin(y), np.nanmax(y)
        xmargin = 0.05 * (xmax - xmin)
        ymargin = 0.05 * (ymax - ymin)
        ax.set_xlim(xmin - xmargin, xmax + xmargin)
        ax.set_ylim(ymin - ymargin, ymax + ymargin)

        def on_select(eclick, erelease):
            xmin, xmax = (
                min(eclick.xdata, erelease.xdata),
                max(eclick.xdata, erelease.xdata),
            )
            ymin, ymax = (
                min(eclick.ydata, erelease.ydata),
                max(eclick.ydata, erelease.ydata),
            )
            self.mask_rectangle(xmin, xmax, ymin, ymax)
            c.set_array(np.ma.array(self._H, mask=self._mask == 0).ravel())
            fig.canvas.draw()
            # self.sel.set_active(False)
            # self.sel = None
            # ax.draw()
            # fig.clf()
            # self.plot(fig=fig, show_mesh=show_mesh)

        if editable:
            self.sel = matplotlib.widgets.RectangleSelector(
                ax, on_select, useblit=True, button=[1], interactive=False
            )
        return fig


# def load(path: str, nz: int, **kwargs) -> "Domain":
#     """Load domain from file. Typically this is a file created by :meth:`Domain.save`.

#     Args:
#         path: NetCDF file to load from
#         nz: number of vertical layers
#         **kwargs: additional keyword arguments to pass to :class:`Domain`
#     """
#     name2kwarg = dict(z0b_min="z0", cor="f")
#     with netCDF4.Dataset(path) as nc:
#         for name in ("lon", "lat", "x", "y", "H", "mask", "z0b_min", "cor", "z0", "f"):
#             if name in nc.variables and name not in kwargs:
#                 kwargs[name2kwarg.get(name, name)] = nc.variables[name][...]
#     spherical = kwargs.get("spherical", "lon" in kwargs)
#     ny, nx = kwargs["lon" if spherical else "x"].shape
#     nx = (nx - 1) // 2
#     ny = (ny - 1) // 2
#     return create(nx, ny, nz, spherical=spherical, **kwargs)

#     def save(self, path: str, full: bool = False, sub: bool = False):
#         """Save grid to a NetCDF file that can be interpreted by :func:`load`.

#         Args:
#             path: NetCDF file to save to
#             full: also save every field including its halos separately
#                 This can be useful for debugging.
#             sub: save the local subdomain, not the global domain
#         """
#         if self.glob is not self and not sub:
#             # We need to save the global domain; not the current subdomain.
#             # If we are the root, divert the plot command to the global domain.
#             # Otherwise just ignore this and return.
#             if self.glob:
#                 self.glob.save(path, full)
#             return

#         with netCDF4.Dataset(path, "w") as nc:

#             def create(
#                 name, units, long_name, values, coordinates: str, dimensions=("y", "x")
#             ):
#                 fill_value = None
#                 if np.ma.getmask(values) is not np.ma.nomask:
#                     fill_value = values.fill_value
#                 ncvar = nc.createVariable(
#                     name, values.dtype, dimensions, fill_value=fill_value
#                 )
#                 ncvar.units = units
#                 ncvar.long_name = long_name
#                 # ncvar.coordinates = coordinates
#                 ncvar[...] = values

#             def create_var(name, units, long_name, values, values_):
#                 if values is None:
#                     return
#                 create(
#                     name,
#                     units,
#                     long_name,
#                     values,
#                     coordinates="lon lat" if self.spherical else "x y",
#                 )
#                 if full:
#                     create(
#                         name + "_",
#                         units,
#                         long_name,
#                         values_,
#                         dimensions=("y_", "x_"),
#                         coordinates="lon_ lat_" if self.spherical else "x_ y_",
#                     )

#             nc.createDimension("x", self.H.shape[1])
#             nc.createDimension("y", self.H.shape[0])
#             if full:
#                 nc.createDimension("x_", self.H_.shape[1])
#                 nc.createDimension("y_", self.H_.shape[0])
#                 create_var("dx", "m", "dx", self.dx, self.dx_)
#                 create_var("dy", "m", "dy", self.dy, self.dy_)
#                 create_var("area", "m2", "area", self.dx * self.dy, self.dx_ * self.dy_)
#             create_var("lat", "degrees_north", "latitude", self.lat, self.lat_)
#             create_var("lon", "degrees_east", "longitude", self.lon, self.lon_)
#             create_var("x", "m", "x", self.x, self.x_)
#             create_var("y", "m", "y", self.y, self.y_)
#             create_var("H", "m", "undisturbed water depth", self.H, self.H_)
#             create_var("mask", "", "mask", self.mask, self.mask_)
#             create_var("z0", "m", "bottom roughness", self.z0b_min, self.z0b_min_)
#             create_var("f", "", "Coriolis parameter", self.cor, self.cor_)

#     def contains(
#         self,
#         x: float,
#         y: float,
#         include_halos: bool = False,
#         spherical: Optional[bool] = None,
#     ) -> bool:
#         """Determine whether the domain contains the specified point.

#         Args:
#             x: native x coordinate (longitude for spherical grids, Cartesian coordinate
#                 in m otherwise)
#             y: native y coordinate (latitude for spherical grids, Cartesian coordinate
#                 in m otherwise)
#             include_halos: whether to also search the halos

#         Returns:
#             True if the point falls within the domain, False otherwise
#         """
#         if spherical is None:
#             spherical = self.spherical
#         local_slice, _, _, _ = self.tiling.subdomain2slices(
#             halox_sub=2 * self.halox,
#             haloy_sub=2 * self.haloy,
#             halox_glob=2 * self.halox,
#             haloy_glob=2 * self.haloy,
#             scale=2,
#             share=1,
#             exclude_halos=not include_halos,
#             exclude_global_halos=True,
#         )
#         allx, ally = (self.lon_, self.lat_) if spherical else (self.x_, self.y_)
#         allx, ally = allx[local_slice], ally[local_slice]
#         ny, nx = allx.shape

#         # Determine whether point falls within current subdomain
#         # based on https://wrf.ecse.rpi.edu/Research/Short_Notes/pnpoly.html
#         x_bnd = np.concatenate(
#             (
#                 allx[0, :-1],
#                 allx[:-1, -1],
#                 allx[-1, nx - 1 : 0 : -1],
#                 allx[ny - 1 : 0 : -1, 0],
#             )
#         )
#         y_bnd = np.concatenate(
#             (
#                 ally[0, :-1],
#                 ally[:-1, -1],
#                 ally[-1, nx - 1 : 0 : -1],
#                 ally[ny - 1 : 0 : -1, 0],
#             )
#         )
#         assert not np.isnan(x_bnd).any(), f"Invalid x boundary: {x_bnd}."
#         assert not np.isnan(y_bnd).any(), f"Invalid y boundary: {y_bnd}."
#         assert x_bnd.size == 2 * ny + 2 * nx - 4
#         inside = False
#         for i, (vertxi, vertyi) in enumerate(zip(x_bnd, y_bnd)):
#             vertxj, vertyj = x_bnd[i - 1], y_bnd[i - 1]
#             if (vertyi > y) != (vertyj > y) and x < (vertxj - vertxi) * (y - vertyi) / (
#                 vertyj - vertyi
#             ) + vertxi:
#                 inside = not inside
#         return inside
