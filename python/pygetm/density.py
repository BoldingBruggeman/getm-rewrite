from typing import Optional

import numpy

from . import core
from . import _pygetm
from .constants import INTERFACES
import pygsw

# Below we package up all equation-of-state/TEOS10 methods in a class.
# This allows users to substitute other methods (e..g, alternative equation of state
# methods) by inheriting from this class and overriding methods.


class Density:
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
            SA: absolute salinity
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
            SA: absolute salinity
            ct: conservative temperature (degrees Celsius)
            p: pressure (dbar). If not provided, the water depth in m will be used as
                approximate pressure.
            out: array to store density result in. If not provided,
                a new array will be created.
        
        Returns:
            array with density values (kg m-3)
        """
        assert SA.grid is ct.grid
        if out is None:
            out = SA.grid.array(z=SA.z)
        assert out.grid is SA.grid
        if p is None:
            p = _pygetm.thickness2center_depth(SA.grid.mask, SA.grid.hn)
        assert p.grid is SA.grid
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
            SA: absolute salinity
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
    def convert_ts(sp: core.Array, pt: core.Array, p: core.Array = None) -> None:
        """Convert practical salinity and potential temperature to absolute salinity
        and conservative temperature. The conversion happens in-place: absolute salinity
        will replace practical salinity, conservative temperature (degrees Celsius)
        will replace potential temperature.
        
        Args:
            sp: practical salinity (PSU)
            pt: potential temperature (degrees Celsius)
            p: pressure (dbar). If not provided, the water depth in m will be used as
                approximate pressure.
        """
        assert sp.grid is pt.grid
        if p is None:
            p = _pygetm.thickness2center_depth(sp.grid.mask, sp.grid.hn)
        assert p.grid is sp.grid
        out = numpy.empty_like(sp.all_values)
        lon = numpy.broadcast_to(sp.grid.lon.all_values, out.shape)
        lat = numpy.broadcast_to(sp.grid.lat.all_values, out.shape)
        pygsw.sa_from_sp(
            lon.ravel(),
            lat.ravel(),
            p.all_values.ravel(),
            sp.all_values.ravel(),
            out.ravel(),
        )
        sp.all_values[...] = out
        pygsw.ct_from_pt(sp.all_values.ravel(), pt.all_values.ravel(), out.ravel())
        pt.all_values[...] = out
