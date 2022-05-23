import numpy
import cftime
import enum

import pyairsea
import pygetm.domain
from . import core
from . import parallel
from .constants import FILL_VALUE

CPA = 1008.  #: specific heat capacity of air (J kg-1 K-1)

class HumidityMeasure(enum.IntEnum):
    """Measure used to specify air humidity"""
    RELATIVE_HUMIDITY = 1     #: relative humidity in %
    WET_BULB_TEMPERATURE = 2  #: wet bulb temperature in in degrees Celsius
    DEW_POINT_TEMPERATURE = 3 #: dewpoint temperature in in degrees Celsius
    SPECIFIC_HUMIDITY = 4     #: specific humidity in kg kg-1

class Fluxes:
    """Base class that provides air-water fluxes of heat and momentum, as well as surface air pressure.
    When using this class directly, these fluxes (stresses ``taux`` and ``tauy``, surface heat flux ``shf``, air pressure ``sp``,
    and net downwelling shortwave radiation ``swr``) are prescribed, not calculated."""
    def initialize(self, domain: pygetm.domain.Domain):
        self.logger = domain.root_logger.getChild('airsea')

        self.taux = domain.T.array(name='tausx', long_name='wind stress in Eastward direction', units='Pa', fill_value=FILL_VALUE)
        self.tauy = domain.T.array(name='tausy', long_name='wind stress in Northward direction', units='Pa', fill_value=FILL_VALUE)
        self.shf = domain.T.array(name='shf', units='W m-2', long_name='surface heat flux', fill_value=FILL_VALUE)
        self.sp = domain.T.array(name='sp', long_name='surface air pressure', units='Pa', fill_value=FILL_VALUE, attrs={'_require_halos': True})
        self.swr = domain.T.array(name='swr', long_name='surface net downwelling shortwave radiation', units='W m-2', fill_value=FILL_VALUE, fabm_standard_name='surface_downwelling_shortwave_flux', attrs={'_3d_only': True})

        self.taux_U = domain.U.array(name='tausxu', fill_value=FILL_VALUE, attrs={'_mask_output': True})
        self.tauy_V = domain.V.array(name='tausyv', fill_value=FILL_VALUE, attrs={'_mask_output': True})

        # Forcing variables for processes operating over the 3D/macro timestep; these lag behind the 2D forcing variables defined above.
        self.spo = domain.T.array()
        self.taux_Uo = domain.U.array()
        self.tauy_Vo = domain.V.array()

        self._ready = False

    def __call__(self, time: cftime.datetime, sst: core.Array, sss: core.Array, calculate_heat_flux: bool) -> None:
        """Update surface flux of momentum (``taux_U`` at U grid, ``tauy_V`` at V grid, ``taux`` and ``tauy`` at T grid)
        and optionally the surface heat flux (``shf``) and net downwelling shortwave flux (``swr``)
        """
        if not self._ready:
            assert self.taux.require_set(self.logger) * self.tauy.require_set(self.logger) * self.sp.require_set(self.logger) * (not calculate_heat_flux or self.shf.require_set(self.logger)) * (not calculate_heat_flux or self.swr.require_set(self.logger))
            self._ready = True
        self.taux.update_halos(parallel.Neighbor.RIGHT)
        self.taux.interp(self.taux_U)
        self.tauy.update_halos(parallel.Neighbor.TOP)
        self.tauy.interp(self.tauy_V)

class FluxesFromMeteo(Fluxes):
    """Calculate air-water fluxes of heat and momentum, as well as surface air pressure, using the pyairsea library.
    The heat flux is the sum of the sensible heat flux, the latent heat flux, and net downwelling longwave radiation
    """
    def __init__(self, longwave_method: int=1, albedo_method: int=1, humidity_measure: HumidityMeasure=HumidityMeasure.DEW_POINT_TEMPERATURE, calculate_swr: bool=True):
        self.longwave_method = longwave_method
        self.albedo_method = albedo_method
        self.humidity_measure = humidity_measure
        self.calculate_swr = calculate_swr

    def initialize(self, domain: pygetm.domain.Domain):
        super().initialize(domain)

        self.es = domain.T.array(name='es', long_name='vapor pressure at saturation', units='Pa', fill_value=FILL_VALUE)
        self.ea = domain.T.array(name='ea', long_name='vapor pressure', units='Pa', fill_value=FILL_VALUE)
        self.qs = domain.T.array(name='qs', long_name='specific humidity at saturation', units='kg kg-1', fill_value=FILL_VALUE)
        self.qa = domain.T.array(name='qa', long_name='specific humidity', units='kg kg-1', fill_value=FILL_VALUE)
        if self.humidity_measure == HumidityMeasure.DEW_POINT_TEMPERATURE:
            self.hum = self.d2m = domain.T.array(name='d2m', long_name='dew point temperature @ 2 m', units='degrees_Celsius', fill_value=FILL_VALUE)
        elif self.humidity_measure == HumidityMeasure.RELATIVE_HUMIDITY:
            self.hum = self.rh = domain.T.array(name='rh', long_name='relative humidity @ 2 m', units='%', fill_value=FILL_VALUE)
        elif self.humidity_measure == HumidityMeasure.WET_BULB_TEMPERATURE:
            self.hum = self.wbt = domain.T.array(name='wbt', long_name='wet bulb temperature @ 2 m', units='degrees_Celsius', fill_value=FILL_VALUE)
        else:
            self.hum = self.qa
        self.rhoa = domain.T.array(name='rhoa', long_name='air density', units='kg m-3', fill_value=FILL_VALUE)

        self.zen = domain.T.array(name='zen', long_name='zenith angle', units='degrees', fill_value=FILL_VALUE, attrs={'_3d_only': True})
        self.albedo = domain.T.array(name='albedo', long_name='albedo', units='1', fill_value=FILL_VALUE, attrs={'_3d_only': True})

        self.t2m = domain.T.array(name='t2m', long_name='air temperature @ 2 m', units='degrees_Celsius', fill_value=FILL_VALUE)

        self.u10 = domain.T.array(name='u10', long_name='wind speed in Eastward direction @ 10 m', units='m s-1', fill_value=FILL_VALUE)
        self.v10 = domain.T.array(name='v10', long_name='wind speed in Northward direction @ 10 m', units='m s-1', fill_value=FILL_VALUE)
        self.tcc = domain.T.array(name='tcc', long_name='total cloud cover', units='1', fill_value=FILL_VALUE, attrs={'_3d_only': True})

        self.w = domain.T.array(name='w', long_name='wind speed', units='m s-1', fabm_standard_name='wind_speed', fill_value=FILL_VALUE)

        self.lon = domain.T.lon
        self.lat = domain.T.lat

        self.qe = domain.T.array(name='qe', long_name='latent heat flux', units='W m-2', fill_value=FILL_VALUE, attrs={'_3d_only': True})
        self.qh = domain.T.array(name='qh', long_name='sensible heat flux', units='W m-2', fill_value=FILL_VALUE, attrs={'_3d_only': True})
        self.ql = domain.T.array(name='ql', long_name='net downwelling longwave radiation', units='W m-2', fill_value=FILL_VALUE, attrs={'_3d_only': True})

        self.cd_mom = domain.T.array()
        self.cd_latent = domain.T.array()
        self.cd_sensible = domain.T.array()

    def update_humidity(self, sst: core.Array):
        """Update humidity metrics: saturation vapor pressure ``es``, actual vapor pressure ``ea``,
        saturation specific humidity ``qs``, actual specific humidity ``qa`` and air density ``rhoa``
        """
        pyairsea.humidity(self.humidity_measure, self.hum.all_values, self.sp.all_values, sst.all_values, self.t2m.all_values, self.es.all_values, self.ea.all_values, self.qs.all_values, self.qa.all_values, self.rhoa.all_values)

    def update_longwave_radiation(self, sst: core.Array):
        """Update net downwelling longwave radiation ``ql``
        """
        sst_K = sst.all_values + 273.15
        t2m_K = self.t2m.all_values + 273.15
        pyairsea.longwave_radiation(self.longwave_method, self.lat.all_values, sst_K, t2m_K, self.tcc.all_values, self.ea.all_values, self.qa.all_values, self.ql.all_values)

    def update_shortwave_radiation(self, time: cftime.datetime):
        """Update net downwelling shortwave radiation ``swr``.
        This represents the value just below the water surface (i.e., what is left after reflection).
        """
        hh = time.hour + time.minute / 60. + time.second / 3600.
        yday = time.timetuple()[-2]
        pyairsea.solar_zenith_angle(yday, hh, self.lon.all_values, self.lat.all_values, self.zen.all_values)
        pyairsea.shortwave_radiation(yday, self.zen.all_values, self.lon.all_values, self.lat.all_values, self.tcc.all_values, self.swr.all_values)
        pyairsea.albedo_water(self.albedo_method, self.zen.all_values, yday, self.albedo.all_values)
        self.swr.all_values *= 1. - self.albedo.all_values

    def update_transfer_coefficients(self, sst: core.Array):
        """Update transfer coefficients for momentum (``cd_mom``), latent heat (``cd_latent``) and sensible heat (``cd_sensible``)
        """
        sst_K = sst.all_values + 273.15
        t2m_K = self.t2m.all_values + 273.15
        numpy.sqrt(self.u10.all_values**2 + self.v10.all_values**2, out=self.w.all_values)
        pyairsea.transfer_coefficients(1, sst_K, t2m_K, self.w.all_values, self.cd_mom.all_values, self.cd_sensible.all_values, self.cd_latent.all_values)

    def __call__(self, time: cftime.datetime, sst: core.Array, sss: core.Array, calculate_heat_flux: bool) -> None:
        if not self._ready:
            assert (self.t2m.require_set(self.logger) * self.hum.require_set(self.logger)
                  * self.u10.require_set(self.logger) * self.v10.require_set(self.logger)
                  * self.tcc.require_set(self.logger) * sst.require_set(self.logger))

        self.update_humidity(sst)
        self.update_transfer_coefficients(sst)

        tmp = self.cd_mom.all_values * self.rhoa.all_values * self.w.all_values
        self.taux.all_values[...] = tmp * self.u10.all_values
        self.tauy.all_values[...] = tmp * self.v10.all_values

        if calculate_heat_flux:
            L = 2.5e6 - 0.00234e6 * sst.all_values   # latent heat of vaporization (J/kg) at sea surface, note SST must be in degrees Celsius
            self.qh.all_values[...] = -self.cd_sensible.all_values * CPA * self.rhoa.all_values * self.w.all_values * (sst.all_values - self.t2m.all_values)
            self.qe.all_values[...] = -self.cd_latent.all_values * L * self.rhoa.all_values * self.w.all_values * (self.qs.all_values - self.qa.all_values)
            self.update_longwave_radiation(sst)
            self.shf.all_values[...] = self.qh.all_values + self.qe.all_values + self.ql.all_values
            self.shf.all_values[numpy.logical_and(sst.all_values < -0.0575 * sss.all_values, self.shf.all_values < 0)] = 0.
            if self.calculate_swr:
                self.update_shortwave_radiation(time)

        super().__call__(time, sst, sss, calculate_heat_flux)

