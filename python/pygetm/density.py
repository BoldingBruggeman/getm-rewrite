from typing import Optional, Union, Tuple
import numbers

import numpy as np
import xarray

from . import core
from . import _pygetm
from .constants import INTERFACES
from .input import Operator, LazyArray
import pygsw

# Below we package up all equation-of-state/TEOS10 methods in a class.
# This allows users to substitute other methods (e..g, alternative equation of state
# methods) by inheriting from this class and overriding methods.


class Density:
    CP = pygsw.CP  #: specific heat (J kg-1 K-1)

    @staticmethod
    def get_buoyancy_frequency(
        SA: core.Array,
        ct: core.Array,
        p: Optional[core.Array] = None,
        out: Optional[core.Array] = None,
    ) -> core.Array:
        """Calculate the square of the buoyancy frequency at layer interfaces from
        absolute salinity, conservative temperature and pressure at the layer centers.

        Args:
            SA: absolute salinity (g kg-1)
            ct: conservative temperature (degrees Celsius)
            p: pressure (dbar). If not provided, the water depth in m will be used as
                approximate pressure.
            out: array to store buoyancy frequency result in. If not provided,
                a new array will be created.

        Returns:
            array with values for the square of the buoyancy frequency (s-2)
        """
        assert SA.grid is ct.grid
        assert SA.ndim == 3
        if out is None:
            out = SA.grid.array(z=INTERFACES)
        assert out.grid is SA.grid and out.z == INTERFACES
        if p is None:
            p = _pygetm.thickness2center_depth(SA.grid.mask, SA.grid.hn)
        assert p.grid is SA.grid
        pygsw.nsquared(
            SA.grid.mask.all_values,
            SA.grid.hn.all_values,
            SA.all_values,
            ct.all_values,
            p.all_values,
            SA.grid.lat.all_values,
            out.all_values[1:-1, :, :],
        )
        return out

    @staticmethod
    def get_density(
        SA: core.Array,
        ct: core.Array,
        p: Optional[core.Array] = None,
        out: Optional[core.Array] = None,
    ) -> core.Array:
        """Calculate in-situ density from absolute salinity and conservative
        temperature. Inputs can be 2D or 3D.

        Args:
            SA: absolute salinity (g kg-1)
            ct: conservative temperature (degrees Celsius)
            p: pressure (dbar). If not provided, the water depth in m will be used as
                approximate pressure.
            out: array to store density result in. If not provided,
                a new array will be created.

        Returns:
            array with density values (kg m-3)
        """
        assert SA.grid is ct.grid and SA.z == ct.z
        assert out is None or (out.grid is SA.grid and out.z == SA.z)
        assert p is None or (p.grid is SA.grid and p.z == SA.z)
        if out is None:
            out = SA.grid.array(z=SA.z)
        if p is None:
            p = _pygetm.thickness2center_depth(SA.grid.mask, SA.grid.hn)
        pygsw.rho(
            SA.all_values.ravel(),
            ct.all_values.ravel(),
            p.all_values.ravel(),
            out.all_values.ravel(),
        )
        return out

    @staticmethod
    def get_potential_temperature(
        SA: core.Array, ct: core.Array, out: Optional[core.Array] = None
    ) -> core.Array:
        """Calculate potential temperature from absolute salinity and conservative
        temperature. Inputs can be 2D or 3D.

        Args:
            SA: absolute salinity (g kg-1)
            ct: conservative temperature (degrees Celsius)
            out: array to store potential temperature result in. If not provided, a new
                array will be created.

        Returns:
            array with potential temperature values (degrees Celsius)
        """
        assert SA.grid is ct.grid
        if out is None:
            out = SA.grid.array(z=SA.z)
        assert out.grid is SA.grid and out.shape == SA.shape
        pygsw.pt_from_ct(
            SA.all_values.ravel(), ct.all_values.ravel(), out.all_values.ravel()
        )
        return out

    @staticmethod
    def convert_ts(
        salt: core.Array,
        temp: core.Array,
        p: Optional[core.Array] = None,
        in_situ: bool = False,
    ) -> None:
        """Convert practical salinity and potential temperature to absolute salinity
        and conservative temperature. The conversion happens in-place: absolute
        salinity (g kg-1) will replace practical salinity, conservative temperature
        (degrees Celsius) will replace potential temperature.

        Args:
            salt: practical salinity (PSU)
            temp: potential temperature or, if ``in_situ=True``, in-situ temperature
                (degrees Celsius)
            p: pressure (dbar). If not provided, the water depth in m will be used as
                approximate pressure.
            in_situ: ``temp`` is in-situ temperature rather than potential temperature
        """
        assert salt.grid is temp.grid
        grid = salt.grid

        # Prepare additional inputs: pressure, longitude, latitude
        if p is None:
            p = _pygetm.thickness2center_depth(grid.mask, grid.hn)
        assert p.grid is grid
        out = np.empty_like(salt.all_values)
        lon = np.broadcast_to(grid.lon.all_values, out.shape)
        lat = np.broadcast_to(grid.lat.all_values, out.shape)

        # Convert from practical to absolute salinity
        pygsw.sa_from_sp(
            lon.ravel(),
            lat.ravel(),
            p.all_values.ravel(),
            salt.all_values.ravel(),
            out.ravel(),
        )
        salt.fill(out)

        if in_situ:
            # Convert from in-situ to potential temperature
            pygsw.pt0_from_t(
                salt.all_values.ravel(),
                temp.all_values.ravel(),
                p.all_values.ravel(),
                out.ravel(),
            )
            temp.fill(out)

        # Convert from potential to conservative temperature
        pygsw.ct_from_pt(salt.all_values.ravel(), temp.all_values.ravel(), out.ravel())
        temp.fill(out)

    @staticmethod
    def lazy_convert_ts(
        salt: xarray.DataArray,
        temp: xarray.DataArray,
        p: Optional[core.Array] = None,
        lon: Optional[core.Array] = None,
        lat: Optional[core.Array] = None,
        in_situ: bool = False,
    ) -> Tuple[xarray.DataArray, xarray.DataArray]:
        """Lazily convert practical salinity and potential temperature to absolute
        salinity and conservative temperature. The conversion is done only when the
        returned objects are indexed or cast to a :class:`numpy.ndarray`.

        Args:
            salt: practical salinity (PSU)
            temp: potential temperature or, if ``in_situ=True``, in-situ temperature
                (degrees Celsius)
            p: pressure (dbar). If not provided, it will be inferred from the depth
                coordinate of ``salt`` or ``temp``.
            lon: longitude (degrees East). If not provided, it will be inferred from
                the coordinates of ``salt`` or ``temp``.
            lat: latitude (degrees North). If not provided, it will be inferred from
                the coordinates of ``salt`` or ``temp``.
            in_situ: input is in-situ temperature rather than potential temperature

        Returns:
            a tuple with lazy arrays for absolute salinity (g kg-1) and conservative
            temperature (degrees Celsius)
        """
        gridsrc = salt if isinstance(salt, xarray.DataArray) else temp

        def _broadcast_coordinate(
            c: xarray.DataArray,
        ) -> Union[xarray.DataArray, np.ndarray]:
            if c.ndim == gridsrc.ndim:
                return c
            s = [np.newaxis] * gridsrc.ndim
            for i, d in enumerate(gridsrc.dims):
                if d in c.dims:
                    s[i] = slice(None)
            return np.asarray(c)[tuple(s)]

        if lon is None:
            lon = _broadcast_coordinate(gridsrc.getm.longitude)
        if lat is None:
            lat = _broadcast_coordinate(gridsrc.getm.latitude)
        if p is None:
            p = _broadcast_coordinate(gridsrc.getm.z)
        lon = np.asarray(lon)
        lat = np.asarray(lat)
        p = np.asarray(p)
        salt = salt if not isinstance(salt, xarray.DataArray) else salt.variable
        temp = temp if not isinstance(temp, xarray.DataArray) else temp.variable
        csalt = LazyConvert(salt, temp, lon, lat, p, in_situ, True)
        ctemp = LazyConvert(salt, temp, lon, lat, p, in_situ, False)
        return (
            xarray.DataArray(csalt, coords=gridsrc.coords, dims=gridsrc.dims),
            xarray.DataArray(ctemp, coords=gridsrc.coords, dims=gridsrc.dims),
        )


class LazyConvert(Operator):
    def __init__(
        self,
        salt: Union[np.ndarray, numbers.Number, LazyArray, xarray.Variable],
        temp: Union[np.ndarray, numbers.Number, LazyArray, xarray.Variable],
        lon: Union[np.ndarray, numbers.Number, LazyArray, xarray.Variable],
        lat: Union[np.ndarray, numbers.Number, LazyArray, xarray.Variable],
        p: Union[np.ndarray, numbers.Number, LazyArray, xarray.Variable],
        in_situ: bool = False,
        return_salt: bool = True,
    ):
        if np.all(p < 0):
            p = -p
        lon = np.asarray(lon, dtype=float)
        lat = np.asarray(lat, dtype=float)
        p = np.asarray(p, dtype=float)
        super().__init__(salt, temp, lon, lat, p, passthrough=True)
        self.in_situ = in_situ
        self.return_salt = return_salt

    def apply(self, salt, temp, lon, lat, p, dtype=None) -> np.ndarray:
        salt = np.asarray(salt, dtype=float)
        temp = np.asarray(temp, dtype=float)
        salt, temp, lon, lat, p = np.broadcast_arrays(salt, temp, lon, lat, p)
        sa = np.empty_like(salt)
        pygsw.sa_from_sp(lon.ravel(), lat.ravel(), p.ravel(), salt.ravel(), sa.ravel())
        if self.return_salt:
            return sa
        pt0 = temp
        if self.in_situ:
            # Convert from in-situ to potential temperature
            pt0 = np.empty_like(temp)
            pygsw.pt0_from_t(sa.ravel(), temp.ravel(), p.ravel(), pt0.ravel())
        ct = np.empty_like(pt0)
        pygsw.ct_from_pt(sa.ravel(), pt0.ravel(), ct.ravel())
        return ct
