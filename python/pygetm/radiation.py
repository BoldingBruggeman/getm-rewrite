from typing import Optional

from . import core
from . import domain
from . import _pygetm
from .constants import FILL_VALUE, INTERFACES, CENTERS

import numpy as np

JERLOV_I = 1
JERLOV_1 = 2
JERLOV_IA = 3
JERLOV_IB = 4
JERLOV_II = 5
JERLOV_III = 6


class Radiation:
    """Base class that provides heating due to shortwave absorption throughout the
    water column. When using this class directly, heating is effectively prescribed,
    not calculated. In this case, the heating per layer defaults to zero; assign to
    swr_abs or call swr_abs.set to change this."""

    def __init__(self, grid: domain.Grid):
        self.grid = grid
        self.logger = grid.domain.root_logger.getChild("radiation")
        self.swr_abs = grid.array(
            name="swr_abs",
            units="W m-2",
            long_name="heating due to shortwave absorption (layer-integrated)",
            z=CENTERS,
            fill_value=FILL_VALUE,
        )
        self.swr_abs.fill(0.0)

    def __call__(self, swr: core.Array):
        """Compute heating due to shortwave radiation throughout the water column"""
        pass


class TwoBand(Radiation):
    """Two-band (visible and non-visible) model for shortwave radiation and absorption
    throughout the water column. It is driven by downwelling shortwave radiation just
    below the water surface, in combination with the fraction of shortwave radiation
    that is non-visible, and the attenuation coefficients for the visible and
    non-visible fractions. All of these can vary spatially, but the non-visible fraction
    and attenuation coefficients can vary only horizontally, not vertically.
    """

    def __init__(self, grid: domain.Grid, jerlov_type: Optional[int] = None):
        super().__init__(grid)

        # Inputs
        self.A = grid.array(
            name="A",
            units="1",
            long_name="non-visible fraction of shortwave radiation",
            fill_value=FILL_VALUE,
        )
        self.kc1 = grid.array(
            name="kc1",
            units="m-1",
            long_name="attenuation of non-visible fraction of shortwave radiation",
            fill_value=FILL_VALUE,
        )
        self.kc2 = grid.array(
            name="kc2",
            units="m-1",
            long_name="attenuation of visible fraction of shortwave radiation",
            fill_value=FILL_VALUE,
        )

        self._first = True
        if jerlov_type:
            self.set_jerlov_type(jerlov_type)

        # Outputs
        self.rad = grid.array(
            name="rad",
            units="W m-2",
            long_name="shortwave radiation",
            fabm_standard_name="downwelling_shortwave_flux",
            z=INTERFACES,
            fill_value=FILL_VALUE,
            attrs=dict(standard_name="downwelling_shortwave_flux_in_sea_water")
        )
        self.par = grid.array(
            name="par",
            units="W m-2",
            long_name="photosynthetically active radiation",
            fabm_standard_name="downwelling_photosynthetic_radiative_flux",
            z=CENTERS,
            fill_value=FILL_VALUE,
            attrs=dict(standard_name="downwelling_photosynthetic_radiative_flux_in_sea_water")
        )
        self.par0 = grid.array(
            name="par0",
            units="W m-2",
            long_name="surface photosynthetically active radiation",
            fabm_standard_name="surface_downwelling_photosynthetic_radiative_flux",
            fill_value=FILL_VALUE,
        )

    def set_jerlov_type(self, jerlov_type: int):
        """Derive non-visible fraction of shortwave radiation (A), attenuation
        coefficients of non-visible shortwave radiation (kc1), and attenuation
        length of visible shortwave radiation (kc2) from the Jerlov water type
        """
        # Note that attentuation in the dictionary below is described by
        # a length scale (g1, g2 in GOTM/old GETM)
        # Its reciprocal is then calculated to set kc1, kc2.
        A, g1, g2 = {
            JERLOV_I: (0.58, 0.35, 23.0),
            JERLOV_1: (0.68, 1.2, 28.0),
            JERLOV_IA: (0.62, 0.6, 20.0),
            JERLOV_IB: (0.67, 1.0, 17.0),
            JERLOV_II: (0.77, 1.5, 14.0),
            JERLOV_III: (0.78, 1.4, 7.9),
        }[jerlov_type]
        self.A.fill(A)
        self.kc1.fill(1.0 / g1)
        self.kc2.fill(1.0 / g2)

    def __call__(self, swr: core.Array):
        """Compute heating due to shortwave radiation throughout the water column"""
        if self._first:
            assert (
                self.A.require_set(self.logger)
                * self.kc1.require_set(self.logger)
                * self.kc2.require_set(self.logger)
            )
            self._first = False

        assert swr.grid is self.grid and not swr.z
        _pygetm.exponential_profile_2band_interfaces(
            self.grid.mask,
            self.grid.hn,
            self.A,
            self.kc1,
            self.kc2,
            top=swr,
            out=self.rad,
        )

        if self.par0.saved or self.par.saved:
            # Visible part of shortwave radiation just below sea surface
            # (i.e., reflection/albedo already accounted for)
            self.par0.all_values[...] = (1.0 - self.A.all_values) * swr.all_values
        if self.par.saved:
            # Visible part of shortwave radiation at layer centers,
            # often used by biogeochemistry
            _pygetm.exponential_profile_1band_centers(
                self.grid.mask, self.grid.hn, self.kc2, top=self.par0, out=self.par
            )

        self.swr_abs.all_values[...] = np.diff(self.rad.all_values, axis=0)

        # all remaining radiation is absorbed at the bottom and
        # injected in the water layer above it
        self.swr_abs.all_values[0, ...] += self.rad.all_values[0, ...]
