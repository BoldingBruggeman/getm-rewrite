from . import core
from . import domain
from . import _pygetm
from .constants import FILL_VALUE, INTERFACES, CENTERS

import numpy

class Radiation:
    def __init__(self, grid: domain.Grid):
        self.grid = grid

        self.A = grid.array(name='A', units='1', long_name='non-visible fraction of shortwave radiation', fill_value=FILL_VALUE)
        self.kc1 = grid.array(name='kc1', units='m-1', long_name='attenuation of non-visible fraction of shortwave radiation', fill_value=FILL_VALUE)
        self.kc2 = grid.array(name='kc2', units='1', long_name='attenuation of visible fraction of shortwave radiation', fill_value=FILL_VALUE)

        self.A.fill(0.7)
        self.kc1.fill(1.)
        self.kc2.fill(1./15.)

        self.rad = grid.array(name='rad', units='W m-2', long_name='shortwave radiation', fabm_standard_name='downwelling_shortwave_flux', z=INTERFACES, fill_value=FILL_VALUE)
        self.par = grid.array(name='par', units='W m-2', long_name='photosynthetically active radiation', fabm_standard_name='downwelling_photosynthetic_radiative_flux', z=CENTERS, fill_value=FILL_VALUE)
        self.par0 = grid.array(name='par0', units='W m-2', long_name='surface photosynthetically active radiation', fabm_standard_name='surface_downwelling_photosynthetic_radiative_flux', fill_value=FILL_VALUE)

    def __call__(self, swr: core.Array):
        """Compute downwelling shortwave radiation throughout the water column"""
        assert swr.grid is self.grid and not swr.z
        _pygetm.exponential_profile_2band_interfaces(self.grid.mask, self.grid.hn, self.A, self.kc1, self.kc2, top=swr, out=self.rad)

        if self.par0.saved or self.par.saved:
            # Visible part of shortwave radiation just below sea surface (i.e., reflection/albedo already accounted for)
            self.par0.all_values[...] = (1. - self.A.all_values) * swr.all_values
        if self.par.saved:
            # Visible part of shortwave radiation at layer centers, often used by biogeochemistry
            _pygetm.exponential_profile_1band_centers(self.grid.mask, self.grid.hn, self.kc2, top=self.par0, out=self.par)
