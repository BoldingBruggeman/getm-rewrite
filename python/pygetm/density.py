import numpy

from . import core
from . import _pygetm
from .constants import *
import pygsw

# Below we package up all equation-of-state/TEOS10 methods in a class.
# This allows users to substitute other methods (e..g, alternative equation of state methods)
# by inheriting from this class and overriding methos.

class Density:
    def get_buoyancy_frequency(self, SA: core.Array, ct: core.Array, p: core.Array=None, out: core.Array=None) -> core.Array:
        """Compute the square of the buoyancy frequency at layer interface from absolute salinity, conservative temperature and pressure at the layer centers."""
        assert SA.grid is ct.grid
        assert SA.ndim == 3
        if out is None:
            out = SA.grid.array(z=INTERFACES)
        assert out.grid is SA.grid and out.z == INTERFACES
        if p is None:
            p = _pygetm.thickness2center_depth(SA.grid.mask, SA.grid.hn)
        assert p.grid is SA.grid
        pygsw.nsquared(SA.grid.mask.all_values, SA.grid.hn.all_values, SA.all_values, ct.all_values, p.all_values, SA.grid.lat.all_values, out.all_values[1:-1, :, :])
        return out

    def get_density(self, SA: core.Array, ct: core.Array, p: core.Array=None, out: core.Array=None) -> core.Array:
        """Compute in-situ density from absolute salinity and conservative temperature. Inputs can be 2D or 3D."""
        assert SA.grid is ct.grid
        if out is None:
            out = SA.grid.array(z=SA.z)
        assert out.grid is SA.grid
        if p is None:
            p = _pygetm.thickness2center_depth(SA.grid.mask, SA.grid.hn)
        assert p.grid is SA.grid
        pygsw.rho(SA.all_values.ravel(), ct.all_values.ravel(), p.all_values.ravel(), out.all_values.ravel())
        return out

    def get_potential_temperature(self, SA: core.Array, ct: core.Array, out: core.Array=None) -> core.Array:
        """Compute potential temperature from absolute salinity and conservative temperature. Inputs can be 2D or 3D."""
        assert SA.grid is ct.grid
        if out is None:
            out = SA.grid.array(z=SA.z)
        assert out.grid is SA.grid and out.shape == SA.shape
        pygsw.pt_from_ct(SA.all_values.ravel(), ct.all_values.ravel(), out.all_values.ravel())
        return out

    def convert_ts(self, sp: core.Array, pt: core.Array, p: core.Array=None):
        """Convert practical salinity and potential temperature to absolute salinity and conservative temperature.
        Pressure must be provided as well as practical salinity and potential temperature.
        The conversion happens in-place."""
        assert sp.grid is pt.grid
        if p is None:
            p = _pygetm.thickness2center_depth(sp.grid.mask, sp.grid.hn)
        assert p.grid is sp.grid
        out = numpy.empty_like(sp.all_values)
        lon = numpy.broadcast_to(sp.grid.lon.all_values, out.shape)
        lat = numpy.broadcast_to(sp.grid.lat.all_values, out.shape)
        pygsw.sa_from_sp(lon.ravel(), lat.ravel(), p.all_values.ravel(), sp.all_values.ravel(), out.ravel())
        sp.all_values[...] = out
        pygsw.ct_from_pt(sp.all_values.ravel(), pt.all_values.ravel(), out.ravel())
        pt.all_values[...] = out