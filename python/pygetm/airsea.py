from typing import Optional, Sequence, Tuple

import cftime

import pyairsea
from . import domain
from . import core

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

def humidity(grid: domain.Grid, method: int, hum: core.Array, airp: core.Array, tw: core.Array, ta: core.Array, out: Optional[Sequence[core.Array]]=None) -> Tuple[core.Array, core.Array, core.Array, core.Array, core.Array]:
    if out is None:
        es = grid.array(long_name='vapor pressure at saturation', units='Pa')
        ea = grid.array(long_name='vapor pressure', units='Pa')
        qs = grid.array(long_name='specific humidity at saturation', units='kg kg-1')
        qa = grid.array(long_name='specific humidity', units='kg kg-1')
        rhoa = grid.array(long_name='air density', units='kg m-3')
    else:
        assert len(out) == 5
        es, ea, qs, qa, rhoa = out
    pyairsea.humidity(method, hum.all_values, airp.all_values, tw.all_values, ta.all_values, es.all_values, ea.all_values, qs.all_values, qa.all_values, rhoa.all_values)
    return es, ea, qs, qa, rhoa

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
    