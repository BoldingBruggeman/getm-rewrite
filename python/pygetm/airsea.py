from typing import Optional, Sequence, Tuple

import numpy
import cftime

import pyairsea
from . import domain
from . import core

cpa=1008.

class Fluxes:
    def __init__(self, domain):
        self.taux = domain.T.array(name='tausx', long_name='wind stress in Eastward direction', units='Pa', fill=0.)
        self.tauy = domain.T.array(name='tausy', long_name='wind stress in Northward direction', units='Pa', fill=0.)
        self.qe = domain.T.array(name='qe', long_name='latent heat flux', units='W m-2', fill=0.)
        self.qh = domain.T.array(name='qh', long_name='sensible heat flux', units='W m-2', fill=0.)
        self.sp = domain.T.array(name='sp', long_name='surface air pressure', units='Pa', fill=0.)

        self.taux_U = domain.U.array(name='tausxu')
        self.tauy_V = domain.V.array(name='tausyv')

    def __call__(self, sst: core.Array) -> None:
        pass

class FluxesFromMeteo(Fluxes):
    def __init__(self, domain):
        super().__init__(domain)
        self.es = domain.T.array(name='es', long_name='vapor pressure at saturation', units='Pa')
        self.ea = domain.T.array(name='ea', long_name='vapor pressure', units='Pa')
        self.qs = domain.T.array(name='qs', long_name='specific humidity at saturation', units='kg kg-1')
        self.qa = domain.T.array(name='qa', long_name='specific humidity', units='kg kg-1')
        self.rhoa = domain.T.array(name='rhoa', long_name='air density', units='kg m-3')

        self.t2m = domain.T.array(name='t2m', long_name='air temperature @ 2 m', units='degrees_Celsius')
        self.d2m = domain.T.array(name='d2m', long_name='dewpoint temperature @ 2 m', units='degrees_Celsius')
        self.u10 = domain.T.array(name='u10', long_name='wind speed in Eastward direction @ 10 m', units='m s-1')
        self.v10 = domain.T.array(name='v10', long_name='wind speed in Northward direction @ 10 m', units='m s-1')
        self.tcc = domain.T.array(name='tcc', long_name='total cloud cover', units='1')

        self.w = domain.T.array(name='w', long_name='wind speed', units='m s-1', fabm_standard_name='wind_speed')

        self.cd_mom = domain.T.array()
        self.cd_latent = domain.T.array()
        self.cd_sensible = domain.T.array()

    def __call__(self, sst: core.Array) -> None:
        pyairsea.humidity(3, self.d2m.all_values, self.sp.all_values, sst.all_values, self.t2m.all_values, self.es.all_values, self.ea.all_values, self.qs.all_values, self.qa.all_values, self.rhoa.all_values)

        self.w.all_values[...] = numpy.sqrt(self.u10.all_values**2 + self.v10.all_values**2)
        L = 2.5e6 - 0.00234e6 * sst.all_values   # latent heat of vaporization (J/kg) at sea surface, note SST must be in degrees Celsius
        pyairsea.transfer_coefficients(1, sst.all_values + 273.15, self.t2m.all_values + 273.15, self.w.all_values, self.cd_mom.all_values, self.cd_latent.all_values, self.cd_sensible.all_values)
        tmp = self.cd_mom.all_values * self.rhoa.all_values * self.w.all_values
        self.taux.all_values[...] = tmp * self.u10.all_values
        self.tauy.all_values[...] = tmp * self.v10.all_values
        self.taux.update_halos()
        self.taux.interp(self.taux_U)
        self.tauy.update_halos()
        self.tauy.interp(self.tauy_V)
        self.qe.all_values[...] = -self.cd_sensible.all_values * cpa * self.rhoa.all_values * self.w.all_values * (sst.all_values - self.t2m.all_values)
        self.qh.all_values[...] = -self.cd_latent.all_values * L * self.rhoa.all_values * self.w.all_values * (self.qs.all_values - self.qa.all_values)

def solar_zenith_angle(time: cftime.datetime, grid: domain.Grid, out: Optional[core.Array]=None):
    if out is None:
        out = grid.array(long_name='zenith angle', units='degrees')
    hh = time.hour + time.minute / 60. + time.second / 3600.
    yday = time.timetuple()[-2]
    pyairsea.solar_zenith_angle(yday, hh, grid.lon.all_values, grid.lat.all_values, out.all_values)
    return out

def shortwave_radiation(time: cftime.datetime, grid: domain.Grid, zenith_angle: core.Array, cloud: core.Array, out: Optional[core.Array]=None):
    if out is None:
        out = grid.array(long_name='shortwave radiation', units='W m-2')
    yday = time.timetuple()[-2]
    pyairsea.shortwave_radiation(yday, zenith_angle.all_values, grid.lon.all_values, grid.lat.all_values, cloud.all_values, out.all_values)
    return out

def albedo_water(time: cftime.datetime, grid: domain.Grid, method: int, zenith_angle: core.Array, out: Optional[core.Array]=None):
    if out is None:
        out = grid.array(long_name='albedo', units='1')
    yday = time.timetuple()[-2]
    pyairsea.albedo_water(method, zenith_angle.all_values, yday, out.all_values)
    return out

def longwave_radiation(grid: domain.Grid, method: int, tw: core.Array, ta: core.Array, cloud: core.Array, ea: core.Array, qa: core.Array, out: Optional[core.Array]=None):
    if out is None:
        out = grid.array(long_name='longwave radiation', units='W m-2')
    pyairsea.longwave_radiation(method, grid.lat.all_values, tw.all_values, ta.all_values, cloud.all_values, ea.all_values, qa.all_values, out.all_values)
    return out
