from typing import Optional
import enum
import logging

from . import core
from . import _pygetm
from .constants import FILL_VALUE, INTERFACES, CENTERS

import numpy as np


class Jerlov(enum.Enum):
    """Jerlov water types, expressed as fraction of total shortwave in band 1
    (non-visible; fast attenuating), the attenuation length (m) of band 1,
    and the attenuation length (m) of band 2 (visible; slowly attenuating).
    Taken from Table 2 of Paulson & Simpson (1977)
    https://doi.org/10.1175/1520-0485(1977)007%3C0952:IMITUO%3E2.0.CO;2,
    """

    Type_I = (0.58, 0.35, 23.0)
    Type_1 = (0.68, 1.2, 28.0)
    Type_IA = (0.62, 0.6, 20.0)
    Type_IB = (0.67, 1.0, 17.0)
    Type_II = (0.77, 1.5, 14.0)
    Type_III = (0.78, 1.4, 7.9)


# For backward compatibility:
JERLOV_I = Jerlov.Type_I
JERLOV_1 = Jerlov.Type_1
JERLOV_IA = Jerlov.Type_IA
JERLOV_IB = Jerlov.Type_IB
JERLOV_II = Jerlov.Type_II
JERLOV_III = Jerlov.Type_III


class Radiation:
    """Base class that provides heating due to shortwave absorption throughout the
    water column. When using this class directly, heating is effectively prescribed,
    not calculated. In this case, the heating per layer defaults to zero; assign to
    :attr:`swr_abs` or call ``swr_abs.set`` to change this."""

    def initialize(self, grid: core.Grid, logger: logging.Logger):
        self.grid = grid
        self.logger = logger
        self.swr_abs = grid.array(
            name="swr_abs",
            units="W m-2",
            long_name="heating due to shortwave absorption (layer-integrated)",
            z=CENTERS,
            fill_value=FILL_VALUE,
        )
        self.swr_abs.fill(0.0)

    def __call__(self, swr: core.Array, kc2_add: Optional[core.Array] = None):
        """Compute heating due to shortwave radiation throughout the water column"""
        pass


class TwoBand(Radiation):
    def __init__(
        self, jerlov_type: Optional[Jerlov] = None, reflect_at_bottom: bool = False
    ):
        """Two-band (visible and non-visible) model for shortwave radiation and
        absorption throughout the water column. It is driven by downwelling shortwave
        radiation just below the water surface, in combination with the fraction of
        shortwave radiation that is non-visible, and the attenuation coefficients for
        the visible and non-visible fractions. All of these can vary spatially, but
        the non-visible fraction and attenuation coefficients can vary only
        horizontally, not vertically.

        Args:
            jerlov_type: Jerlov water type to infer attenuation coefficients
                :attr:`kc1` and :attr:`kc2` and the non-visible fraction of shortwave
                radiation :attr:`A` from. If not provided, these quantities
                are potentially horizontally and temporally variable; they can be set
                by calling :meth:`pygetm.core.Array.set` on :attr:`kc1`, :attr:`kc2` and
                :attr:`A`.
            reflect_at_bottom: reflect part of the radiation arriving at the bottom,
                and track this stream upward, attenuating it along the way.
                The spatially explicit :attr:`bottom_albedo` can subsequently be set
                by calling :meth:`pygetm.core.Array.set` on it. Radiation that is
                absorbed at the bottom (the fraction that is not reflected) will heat
                the bottom water layer, as the heat budget of sediments is not modelled.
        """
        self.reflect_at_bottom = reflect_at_bottom
        self.initial_jerlov_type = jerlov_type
        self._first = True

    def initialize(self, grid: core.Grid, logger: logging.Logger):
        super().initialize(grid, logger)
        self.swr_abs.attrs["_mask_output"] = True

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

        self._kc = grid.array(z=CENTERS)
        self._rad = grid.array(z=INTERFACES, fill=0.0)
        self._swr = grid.array()

        if self.initial_jerlov_type:
            self.jerlov_type = self.initial_jerlov_type

        # Outputs
        self.rad = grid.array(
            name="rad",
            units="W m-2",
            long_name="shortwave radiation",
            fabm_standard_name="downwelling_shortwave_flux",
            z=INTERFACES,
            fill_value=FILL_VALUE,
            attrs=dict(standard_name="downwelling_shortwave_flux_in_sea_water"),
        )
        self.par = grid.array(
            name="par",
            units="W m-2",
            long_name="photosynthetically active radiation",
            fabm_standard_name="downwelling_photosynthetic_radiative_flux",
            z=CENTERS,
            fill_value=FILL_VALUE,
            attrs=dict(
                standard_name="downwelling_photosynthetic_radiative_flux_in_sea_water"
            ),
        )
        self.par0 = grid.array(
            name="par0",
            units="W m-2",
            long_name="surface photosynthetically active radiation",
            fabm_standard_name="surface_downwelling_photosynthetic_radiative_flux",
            fill_value=FILL_VALUE,
        )
        if self.reflect_at_bottom:
            self.bottom_albedo = grid.array(
                name="bottom_albedo",
                units="1",
                long_name="bottom albedo",
                fill_value=FILL_VALUE,
            )
            self.bottom_albedo.fill(0.0)
            self.rad_up = grid.array(
                name="rad_up",
                units="W m-2",
                long_name="upwelling shortwave radiation",
                z=INTERFACES,
                fill_value=FILL_VALUE,
            )
            self._swr_up = grid.array()

    def set_jerlov_type(self, jerlov_type: Jerlov):
        """Derive the non-visible fraction of shortwave radiation (:attr:`A`), the
        attenuation coefficient of non-visible shortwave radiation (:attr:`kc1`),
        and the attenuation coefficient of visible shortwave radiation (:attr:`kc2`)
        from the Jerlov water type. These three coefficients will be constant in
        time and space.
        """
        # Attentuation in the Jerlov values is described by a length scale (g1, g2
        # in GOTM/legacy GETM). Its reciprocal is calculated to set kc1, kc2.
        A, g1, g2 = jerlov_type.value
        kc1 = 1.0 / g1
        kc2 = 1.0 / g2
        self.A.fill(A)
        self.kc1.fill(kc1)
        self.kc2.fill(kc2)
        self.logger.info(
            f"Using Jerlov {jerlov_type.name.replace('_', ' ')} attenuation:"
            f" A={A}, kc1={kc1:.4f} m-1, kc2={kc2:.4f} m-1"
        )

    jerlov_type = property(fset=set_jerlov_type)

    def __call__(self, swr: core.Array, kc2_add: Optional[core.Array] = None):
        """Compute heating due to shortwave radiation throughout the water column

        Args:
            swr: net downwelling shortwave radiation just below the water surface
                (i.e., what is left after reflection).
            kc2_add: additional depth-varying attenuation (m-1) of second waveband,
                typically associated with chlorophyll or suspended matter
        """
        if self._first:
            assert (
                self.A.require_set(self.logger)
                * self.kc1.require_set(self.logger)
                * self.kc2.require_set(self.logger)
            )
            self._first = False

        assert swr.grid is self.grid and not swr.z
        assert kc2_add is None or (kc2_add.grid is self.grid and kc2_add.z == CENTERS)

        # Non-visible band - downward
        self._kc.all_values[...] = self.kc1.all_values
        self._swr.all_values[...] = self.A.all_values * swr.all_values
        _pygetm.exponential_profile_1band_interfaces(
            self.grid.mask, self.grid.hn, self._kc, self._swr, up=False, out=self.rad
        )

        # Non-visible band - upward
        if self.reflect_at_bottom:
            np.multiply(
                self.rad.all_values[0, ...],
                self.bottom_albedo.all_values,
                out=self._swr_up.all_values,
            )
            _pygetm.exponential_profile_1band_interfaces(
                self.grid.mask,
                self.grid.hn,
                self._kc,
                self._swr_up,
                up=True,
                out=self.rad_up,
            )

        # Visible band - downward
        self._kc.all_values[...] = self.kc2.all_values
        if kc2_add is not None:
            self._kc.all_values += kc2_add.all_values
        self._swr.all_values[...] = swr.all_values - self._swr.all_values
        _pygetm.exponential_profile_1band_interfaces(
            self.grid.mask, self.grid.hn, self._kc, self._swr, up=False, out=self._rad
        )

        # Total downward
        self.rad.all_values += self._rad.all_values
        self.swr_abs.all_values[...] = np.diff(self.rad.all_values, axis=0)

        # all remaining radiation is absorbed at the bottom and
        # injected in the water layer above it
        self.swr_abs.all_values[0, ...] += self.rad.all_values[0, ...]

        # Visible band and total - upward
        if self.reflect_at_bottom:
            np.multiply(
                self._rad.all_values[0, ...],
                self.bottom_albedo.all_values,
                out=self._swr_up.all_values,
            )
            _pygetm.exponential_profile_1band_interfaces(
                self.grid.mask,
                self.grid.hn,
                self._kc,
                self._swr_up,
                up=True,
                out=self._rad,
            )
            self.rad_up.all_values += self._rad.all_values
            self.swr_abs.all_values[0, ...] -= self.rad_up.all_values[0, ...]
            self.swr_abs.all_values -= np.diff(self.rad_up.all_values, axis=0)

        if self.par0.saved:
            # Visible part of shortwave radiation just below sea surface
            # (i.e., reflection/albedo already accounted for)
            self.par0.all_values[...] = self._swr.all_values
        if self.par.saved:
            # Visible part of shortwave radiation at layer centers,
            # often used by biogeochemistry
            _pygetm.exponential_profile_1band_centers(
                self.grid.mask, self.grid.hn, self._kc, top=self._swr, out=self.par
            )
