from typing import Optional, Sequence, Tuple
import logging

import numpy
import cftime

import pyairsea
from . import domain
from . import core
from .constants import FILL_VALUE

CPA = 1008.

class Fluxes:
    def __init__(self, domain, logger: Optional[logging.Logger]=None):
        self.logger = logger or logging.getLogger()

        self.taux = domain.T.array(name='tausx', long_name='wind stress in Eastward direction', units='Pa', fill_value=FILL_VALUE)
        self.tauy = domain.T.array(name='tausy', long_name='wind stress in Northward direction', units='Pa', fill_value=FILL_VALUE)
        self.qe = domain.T.array(name='qe', long_name='sensible heat flux', units='W m-2', fill_value=FILL_VALUE)
        self.qh = domain.T.array(name='qh', long_name='latent heat flux', units='W m-2', fill_value=FILL_VALUE)
        self.ql = domain.T.array(name='ql', long_name='net downwelling longwave radiation', units='W m-2', fill_value=FILL_VALUE)
        self.sp = domain.T.array(name='sp', long_name='surface air pressure', units='Pa', fill_value=FILL_VALUE)
        self.swr = domain.T.array(name='swr', long_name='surface net downwelling shortwave radiation', units='W m-2', fill_value=FILL_VALUE, fabm_standard_name='surface_downwelling_shortwave_flux')

        self.taux_U = domain.U.array(name='tausxu', fill_value=FILL_VALUE)
        self.tauy_V = domain.V.array(name='tausyv', fill_value=FILL_VALUE)

        self._ready = False

    def __call__(self, time: cftime.datetime, sst: core.Array) -> None:
        if not self._ready:
            assert self.taux.require_set(self.logger) * self.tauy.require_set(self.logger) * self.sp.require_set(self.logger)
            self._ready = True
        self.taux.update_halos()
        self.taux.interp(self.taux_U)
        self.tauy.update_halos()
        self.tauy.interp(self.tauy_V)

class FluxesFromMeteo(Fluxes):
    def __init__(self, domain, logger: Optional[logging.Logger]=None, longwave_method: int=1, albedo_method: int=1, compute_swr: bool=True):
        super().__init__(domain, logger)
        self.es = domain.T.array(name='es', long_name='vapor pressure at saturation', units='Pa', fill_value=FILL_VALUE)
        self.ea = domain.T.array(name='ea', long_name='vapor pressure', units='Pa', fill_value=FILL_VALUE)
        self.qs = domain.T.array(name='qs', long_name='specific humidity at saturation', units='kg kg-1', fill_value=FILL_VALUE)
        self.qa = domain.T.array(name='qa', long_name='specific humidity', units='kg kg-1', fill_value=FILL_VALUE)
        self.rhoa = domain.T.array(name='rhoa', long_name='air density', units='kg m-3', fill_value=FILL_VALUE)

        self.zen = domain.T.array(name='zen', long_name='zenith angle', units='degrees', fill_value=FILL_VALUE)
        self.albedo = domain.T.array(name='albedo', long_name='albedo', units='1', fill_value=FILL_VALUE)

        self.t2m = domain.T.array(name='t2m', long_name='air temperature @ 2 m', units='degrees_Celsius', fill_value=FILL_VALUE)
        self.d2m = domain.T.array(name='d2m', long_name='dewpoint temperature @ 2 m', units='degrees_Celsius', fill_value=FILL_VALUE)
        self.u10 = domain.T.array(name='u10', long_name='wind speed in Eastward direction @ 10 m', units='m s-1', fill_value=FILL_VALUE)
        self.v10 = domain.T.array(name='v10', long_name='wind speed in Northward direction @ 10 m', units='m s-1', fill_value=FILL_VALUE)
        self.tcc = domain.T.array(name='tcc', long_name='total cloud cover', units='1', fill_value=FILL_VALUE)

        self.w = domain.T.array(name='w', long_name='wind speed', units='m s-1', fabm_standard_name='wind_speed', fill_value=FILL_VALUE)

        self.lon = domain.T.lon
        self.lat = domain.T.lat

        self.cd_mom = domain.T.array()
        self.cd_latent = domain.T.array()
        self.cd_sensible = domain.T.array()

        self.longwave_method = longwave_method
        self.albedo_method = albedo_method
        self.compute_swr = compute_swr

    def humidity(self, sst: core.Array):
        pyairsea.humidity(3, self.d2m.all_values, self.sp.all_values, sst.all_values, self.t2m.all_values, self.es.all_values, self.ea.all_values, self.qs.all_values, self.qa.all_values, self.rhoa.all_values)

    def longwave_radiation(self, sst: core.Array):
        sst_K = sst.all_values + 273.15
        t2m_K = self.t2m.all_values + 273.15
        pyairsea.longwave_radiation(self.longwave_method, self.lat.all_values, sst_K, t2m_K, self.tcc.all_values, self.ea.all_values, self.qa.all_values, self.ql.all_values)

    def shortwave_radiation(self, time: cftime.datetime):
        hh = time.hour + time.minute / 60. + time.second / 3600.
        yday = time.timetuple()[-2]
        pyairsea.solar_zenith_angle(yday, hh, self.lon.all_values, self.lat.all_values, self.zen.all_values)
        pyairsea.shortwave_radiation(yday, self.zen.all_values, self.lon.all_values, self.lat.all_values, self.tcc.all_values, self.swr.all_values)
        pyairsea.albedo_water(self.albedo_method, self.zen.all_values, yday, self.albedo.all_values)
        self.swr.all_values[...] *= 1 - self.albedo.all_values

    def fluxes(self, sst: core.Array):
        sst_K = sst.all_values + 273.15
        t2m_K = self.t2m.all_values + 273.15
        self.w.all_values[...] = numpy.sqrt(self.u10.all_values**2 + self.v10.all_values**2)
        L = 2.5e6 - 0.00234e6 * sst.all_values   # latent heat of vaporization (J/kg) at sea surface, note SST must be in degrees Celsius
        pyairsea.transfer_coefficients(1, sst_K, t2m_K, self.w.all_values, self.cd_mom.all_values, self.cd_latent.all_values, self.cd_sensible.all_values)
        tmp = self.cd_mom.all_values * self.rhoa.all_values * self.w.all_values
        self.taux.all_values[...] = tmp * self.u10.all_values
        self.tauy.all_values[...] = tmp * self.v10.all_values
        self.qe.all_values[...] = -self.cd_sensible.all_values * CPA * self.rhoa.all_values * self.w.all_values * (sst.all_values - self.t2m.all_values)
        self.qh.all_values[...] = -self.cd_latent.all_values * L * self.rhoa.all_values * self.w.all_values * (self.qs.all_values - self.qa.all_values)

    def __call__(self, time: cftime.datetime, sst: core.Array) -> None:
        if not self._ready:
            assert (self.t2m.require_set(self.logger) * self.d2m.require_set(self.logger)
                  * self.u10.require_set(self.logger) * self.v10.require_set(self.logger)
                  * self.tcc.require_set(self.logger) * sst.require_set(self.logger))

        self.humidity(sst)
        self.fluxes(sst)
        self.longwave_radiation(sst)
        if self.compute_swr:
            self.shortwave_radiation(time)
        super().__call__(time, sst)

