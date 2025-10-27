import logging

import numpy as np

import pygetm.core
import pygetm.airsea
from pygetm.constants import FILL_VALUE, CellType


class Base:
    """Simple ice model that assumes a cell is completely ice covered when its
    surface temperature drops below freezing. At that point, surface heat
    fluxes in that cell are clipped to positive values (= no further cooling).
    Surface momentum fluxes for that same cell are switched off altogether."""

    def initialize(self, grid: pygetm.core.Grid, logger: logging.Logger):
        self.logger = logger
        self.grid = grid
        self.ice = grid.array(
            name="ice",
            long_name="ice cover",
            units="1",
            fill_value=FILL_VALUE,
            attrs=dict(
                standard_name="sea_ice_area_fraction", _valid_at=(CellType.BOUNDARY,)
            ),
            fabm_standard_name="ice_area_fraction",
        )
        self.ice.fill(0.0)


class Prescribed(Base):
    def initialize(self, grid: pygetm.core.Grid, logger: logging.Logger):
        super().initialize(grid, logger)
        self.has_ice = False
        self.ice_free = np.full(grid.H.all_values.shape, 1.0)

    def __call__(
        self,
        macro: bool,
        ct_sf: pygetm.core.Array,
        sa_sf: pygetm.core.Array,
        airsea: pygetm.airsea.Fluxes,
    ):
        if macro:
            f = np.where(self.ice.all_mask, 0.0, self.ice.all_values)
            self.has_ice = f.any()

            # 1st order freezing point approximation based on
            # gsw_mod_freezing_poly_coefficients
            ct_freezing = 0.017947064327968736 - 0.06076099099929818 * sa_sf.all_values

            if self.has_ice:
                self.ice_free[...] = 1.0 - f

                # Set temperature in ice-covered fraction to freezing temperature
                ct_sf.all_values += f * (ct_freezing - ct_sf.all_values)

            np.maximum(
                ct_sf.all_values,
                ct_freezing,
                out=ct_sf.all_values,
                where=ct_sf.all_mask,
            )

            # Allow outward surface heat flux [cooling] only in ice-free area
            cooling = airsea.shf.all_values < 0
            airsea.shf.all_values[cooling] *= self.ice_free[cooling]

        if self.has_ice:
            airsea.taux.all_values *= self.ice_free
            airsea.tauy.all_values *= self.ice_free


class Ice(Base):
    """Simple ice model that assumes a cell is completely ice covered when its
    surface temperature drops below freezing. At that point, surface heat
    fluxes in that cell are clipped to positive values (= no further cooling).
    Surface momentum fluxes for that same cell are switched off altogether."""

    def initialize(self, grid: pygetm.core.Grid, logger: logging.Logger):
        super().initialize(grid, logger)
        self.has_ice = False
        self.covered = np.full(grid.H.all_values.shape, False)

    def __call__(
        self,
        macro: bool,
        ct_sf: pygetm.core.Array,
        sa_sf: pygetm.core.Array,
        airsea: pygetm.airsea.Fluxes,
    ):
        if macro:
            # 1st order freezing point approximation based on
            # gsw_mod_freezing_poly_coefficients
            ct_freezing = 0.017947064327968736 - 0.06076099099929818 * sa_sf.all_values
            unmasked = ~ct_sf.all_mask
            np.less_equal(
                ct_sf.all_values, ct_freezing, out=self.covered, where=unmasked
            )
            self.has_ice = self.covered.any()
            np.putmask(self.ice.all_values, unmasked, np.where(self.covered, 1.0, 0.0))
            if self.has_ice:
                np.putmask(ct_sf.all_values, self.covered, ct_freezing)
                airsea.shf.all_values[self.covered & (airsea.shf.all_values < 0)] = 0.0

        if self.has_ice:
            airsea.taux.all_values[self.covered] = 0.0
            airsea.tauy.all_values[self.covered] = 0.0
