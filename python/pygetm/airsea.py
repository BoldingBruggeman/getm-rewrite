import numpy as np
import cftime

import pyairsea
import pygetm.domain
from . import core
from .constants import FILL_VALUE, RHO0, TimeVarying
from pyairsea import HumidityMeasure, LongwaveMethod, AlbedoMethod

CPA = 1008.0  #: specific heat capacity of air (J kg-1 K-1)


class Fluxes:
    """Base class that provides air-water fluxes of heat and momentum, as well as
    surface air pressure. When using this class directly, these fluxes (stresses
    :attr:`taux` and :attr:`tauy`, surface heat flux :attr:`shf`, air pressure
    :attr:`sp`, net downwelling shortwave radiation :attr:`swr` and the net
    freshwater flux :attr:`pe`) are prescribed, not calculated.
    """

    def initialize(self, grid: pygetm.domain.Grid):
        self.logger = grid.domain.root_logger.getChild("airsea")

        self.taux = grid.array(
            name="tausx",
            long_name="wind stress in x direction",
            units="Pa",
            fill_value=FILL_VALUE,
            attrs=dict(standard_name="downward_x_stress_at_sea_water_surface"),
        )
        self.tauy = grid.array(
            name="tausy",
            long_name="wind stress in y direction",
            units="Pa",
            fill_value=FILL_VALUE,
            attrs=dict(standard_name="downward_y_stress_at_sea_water_surface"),
        )
        self.shf = grid.array(
            name="shf",
            units="W m-2",
            long_name="surface heat flux",
            fill_value=FILL_VALUE,
        )
        self.sp = grid.array(
            name="sp",
            long_name="surface air pressure",
            units="Pa",
            fill_value=FILL_VALUE,
            attrs=dict(_require_halos=True, standard_name="surface_air_pressure"),
        )
        self.swr = grid.array(
            name="swr",
            long_name="surface net downwelling shortwave radiation",
            units="W m-2",
            fill_value=FILL_VALUE,
            fabm_standard_name="surface_downwelling_shortwave_flux",
            attrs=dict(
                _time_varying=TimeVarying.MACRO,
                standard_name="net_downward_shortwave_flux_at_sea_water_surface",
            ),
        )
        self.pe = grid.array(
            name="pe",
            long_name=(
                "net freshwater flux due to precipitation, condensation,"
                " evaporation",
            ),
            units="m s-1",
            fill_value=FILL_VALUE,
            attrs=dict(_time_varying=TimeVarying.MACRO),
        )
        self.pe.fill(0.0)

        self._ready = False

    def __call__(
        self, time: cftime.datetime, sst: core.Array, calculate_heat_flux: bool,
    ) -> None:
        """Update surface fluxes. Stresses and air pressure are always updated.
        The surface heat flux (:attr:`shf`), net downwelling shortwave flux
        (:attr:`swr`) and net freshwater flux (:attr:`pe`) are updated only if
        ``calculate_heat_flux`` is ``True``.
        """
        if not self._ready:
            assert (
                self.taux.require_set(self.logger)
                * self.tauy.require_set(self.logger)
                * self.sp.require_set(self.logger)
                * (not calculate_heat_flux or self.shf.require_set(self.logger))
                * (not calculate_heat_flux or self.swr.require_set(self.logger))
            )
            self._ready = True


class FluxesFromMeteo(Fluxes):
    """Calculate air-water fluxes of heat and momentum, as well as surface air
    pressure, using the :mod:`pyairsea` library. The heat flux is the sum of the
    sensible heat flux, the latent heat flux, and net downwelling longwave radiation.
    """

    def __init__(
        self,
        longwave_method: LongwaveMethod = LongwaveMethod.CLARK,
        albedo_method: AlbedoMethod = AlbedoMethod.PAYNE,
        humidity_measure: HumidityMeasure = HumidityMeasure.DEW_POINT_TEMPERATURE,
        calculate_swr: bool = True,
        calculate_evaporation: bool = False,
    ):
        self.longwave_method = longwave_method
        self.albedo_method = albedo_method
        self.humidity_measure = humidity_measure
        self.calculate_swr = calculate_swr
        self.calculate_evaporation = calculate_evaporation

    def initialize(self, grid: pygetm.domain.Grid):
        super().initialize(grid)

        self.logger.info("Longwave method: %s" % self.longwave_method.name)
        self.logger.info("Albedo method: %s" % self.albedo_method.name)
        self.logger.info("Humidity measure method: %s" % self.humidity_measure.name)
        if self.calculate_swr:
            self.logger.info(
                "Shortwave radiation calculated from time, location and cloud cover"
            )
        if self.calculate_evaporation:
            self.logger.info("Evaporation calculated from latent heat flux")

        self.es = grid.array(
            name="es",
            long_name="vapor pressure at saturation",
            units="Pa",
            fill_value=FILL_VALUE,
        )
        self.ea = grid.array(
            name="ea", long_name="vapor pressure", units="Pa", fill_value=FILL_VALUE
        )
        self.qs = grid.array(
            name="qs",
            long_name="specific humidity at saturation",
            units="kg kg-1",
            fill_value=FILL_VALUE,
        )
        self.qa = grid.array(
            name="qa",
            long_name="specific humidity",
            units="kg kg-1",
            fill_value=FILL_VALUE,
            attrs=dict(standard_name="specific_humidity"),
        )
        if self.humidity_measure == HumidityMeasure.DEW_POINT_TEMPERATURE:
            self.hum = self.d2m = grid.array(
                name="d2m",
                long_name="dew point temperature @ 2 m",
                units="degrees_Celsius",
                fill_value=FILL_VALUE,
                attrs=dict(standard_name="dew_point_temperature"),
            )
        elif self.humidity_measure == HumidityMeasure.RELATIVE_HUMIDITY:
            self.hum = self.rh = grid.array(
                name="rh",
                long_name="relative humidity @ 2 m",
                units="%",
                fill_value=FILL_VALUE,
                attrs=dict(standard_name="relative_humidity"),
            )
        elif self.humidity_measure == HumidityMeasure.WET_BULB_TEMPERATURE:
            self.hum = self.wbt = grid.array(
                name="wbt",
                long_name="wet bulb temperature @ 2 m",
                units="degrees_Celsius",
                fill_value=FILL_VALUE,
                attrs=dict(standard_name="wet_bulb_temperature"),
            )
        else:
            self.hum = self.qa
        self.rhoa = grid.array(
            name="rhoa",
            long_name="air density",
            units="kg m-3",
            fill_value=FILL_VALUE,
            attrs=dict(standard_name="air_density"),
        )

        self.zen = grid.array(
            name="zen",
            long_name="solar zenith angle",
            units="degree",
            fill_value=FILL_VALUE,
            attrs=dict(
                _time_varying=TimeVarying.MACRO, standard_name="solar_zenith_angle"
            ),
        )
        self.albedo = grid.array(
            name="albedo",
            long_name="albedo",
            units="1",
            fill_value=FILL_VALUE,
            attrs=dict(_time_varying=TimeVarying.MACRO, standard_name="surface_albedo"),
        )

        self.t2m = grid.array(
            name="t2m",
            long_name="air temperature @ 2 m",
            units="degrees_Celsius",
            fill_value=FILL_VALUE,
            attrs=dict(standard_name="air_temperature"),
        )

        self.u10 = grid.array(
            name="u10",
            long_name="wind speed in Eastward direction @ 10 m",
            units="m s-1",
            fill_value=FILL_VALUE,
            attrs=dict(standard_name="eastward_wind"),
        )
        self.v10 = grid.array(
            name="v10",
            long_name="wind speed in Northward direction @ 10 m",
            units="m s-1",
            fill_value=FILL_VALUE,
            attrs=dict(standard_name="northward_wind"),
        )
        self.tcc = grid.array(
            name="tcc",
            long_name="total cloud cover",
            units="1",
            fill_value=FILL_VALUE,
            attrs=dict(
                _time_varying=TimeVarying.MACRO, standard_name="cloud_area_fraction"
            ),
        )

        self.w = grid.array(
            name="w",
            long_name="wind speed",
            units="m s-1",
            fabm_standard_name="wind_speed",
            fill_value=FILL_VALUE,
            attrs=dict(standard_name="wind_speed"),
        )

        self.lon = grid.lon
        self.lat = grid.lat

        self.qe = grid.array(
            name="qe",
            long_name="latent heat flux",
            units="W m-2",
            fill_value=FILL_VALUE,
            attrs=dict(
                _time_varying=TimeVarying.MACRO,
                standard_name="surface_downward_latent_heat_flux",
            ),
        )
        self.qh = grid.array(
            name="qh",
            long_name="sensible heat flux",
            units="W m-2",
            fill_value=FILL_VALUE,
            attrs=dict(
                _time_varying=TimeVarying.MACRO,
                standard_name="surface_downward_sensible_heat_flux",
            ),
        )
        self.ql = grid.array(
            name="ql",
            long_name="net downwelling longwave radiation",
            units="W m-2",
            fill_value=FILL_VALUE,
            attrs=dict(
                _time_varying=TimeVarying.MACRO,
                standard_name="surface_net_downward_longwave_flux",
            ),
        )

        self.tp = grid.array(
            name="tp",
            long_name="total precipitation",
            units="m s-1",
            fill_value=FILL_VALUE,
            attrs=dict(
                _time_varying=TimeVarying.MACRO, standard_name="lwe_precipitation_rate"
            ),
        )
        self.e = grid.array(
            name="e",
            long_name="evaporation minus condensation",
            units="m s-1",
            fill_value=FILL_VALUE,
            attrs=dict(_time_varying=TimeVarying.MACRO),
        )

        self.cd_mom = grid.array()
        self.cd_latent = grid.array()
        self.cd_sensible = grid.array()
        self.grid = grid

    def update_humidity(self, sst: core.Array):
        """Update humidity metrics: saturation vapor pressure :attr:`es`,
        actual vapor pressure :attr:`ea`, saturation specific humidity :attr:`qs`,
        actual specific humidity :attr:`qa` and air density :attr:`rhoa`
        """
        pyairsea.humidity(
            self.humidity_measure,
            self.hum.all_values,
            self.sp.all_values,
            sst.all_values,
            self.t2m.all_values,
            self.es.all_values,
            self.ea.all_values,
            self.qs.all_values,
            self.qa.all_values,
            self.rhoa.all_values,
        )

    def update_longwave_radiation(self, sst: core.Array):
        """Update net downwelling longwave radiation :attr:`ql`
        """
        sst_K = sst.all_values + 273.15
        t2m_K = self.t2m.all_values + 273.15
        pyairsea.longwave_radiation(
            self.longwave_method,
            self.lat.all_values,
            sst_K,
            t2m_K,
            self.tcc.all_values,
            self.ea.all_values,
            self.qa.all_values,
            self.ql.all_values,
        )

    def update_shortwave_radiation(self, time: cftime.datetime):
        """Update net downwelling shortwave radiation :attr:`swr`.
        This represents the value just below the water surface
        (i.e., what is left after reflection).
        """
        hh = time.hour + time.minute / 60.0 + time.second / 3600.0
        yday = time.timetuple().tm_yday
        pyairsea.solar_zenith_angle(
            yday, hh, self.lon.all_values, self.lat.all_values, self.zen.all_values
        )
        pyairsea.shortwave_radiation(
            yday,
            self.zen.all_values,
            self.lon.all_values,
            self.lat.all_values,
            self.tcc.all_values,
            self.swr.all_values,
        )
        pyairsea.albedo_water(
            self.albedo_method, self.zen.all_values, yday, self.albedo.all_values
        )
        self.swr.all_values *= 1.0 - self.albedo.all_values

    def update_transfer_coefficients(self, sst: core.Array):
        """Update transfer coefficients for momentum (:attr:`cd_mom`), latent heat
        (:attr:`cd_latent`) and sensible heat (:attr:`cd_sensible`)
        """

        # Air and water temperature in Kelvin
        sst_K = sst.all_values + 273.15
        t2m_K = self.t2m.all_values + 273.15

        # Wind speed
        np.hypot(self.u10.all_values, self.v10.all_values, out=self.w.all_values)

        pyairsea.transfer_coefficients(
            1,
            sst_K,
            t2m_K,
            self.w.all_values,
            self.cd_mom.all_values,
            self.cd_sensible.all_values,
            self.cd_latent.all_values,
        )

    def __call__(
        self, time: cftime.datetime, sst: core.Array, calculate_heat_flux: bool,
    ) -> None:

        if not self._ready:
            assert (
                self.t2m.require_set(self.logger)
                * self.hum.require_set(self.logger)
                * self.u10.require_set(self.logger)
                * self.v10.require_set(self.logger)
                * (not calculate_heat_flux or self.tcc.require_set(self.logger))
                * sst.require_set(self.logger)
                * (
                    not self.calculate_evaporation
                    or not calculate_heat_flux
                    or self.tp.require_set(self.logger)
                )
            )

        # Air humidity
        self.update_humidity(sst)

        # Transfer coefficients of heat and momentum
        self.update_transfer_coefficients(sst)

        # Momentum flux in x and y direction
        tmp = self.cd_mom.all_values * self.rhoa.all_values * self.w.all_values
        u10, v10 = self.grid.rotate(self.u10, self.v10)
        self.taux.all_values[...] = tmp * u10.all_values
        self.tauy.all_values[...] = tmp * v10.all_values

        if calculate_heat_flux:
            # Latent heat of vaporization (J/kg) at sea surface
            # Note SST must be in degrees Celsius
            L = 2.5e6 - 0.00234e6 * sst.all_values

            # Sensible heat flux
            self.qh.all_values[...] = (
                self.cd_sensible.all_values
                * CPA
                * self.rhoa.all_values
                * self.w.all_values
                * (self.t2m.all_values - sst.all_values)
            )

            # Latent heat flux
            self.qe.all_values[...] = (
                self.cd_latent.all_values
                * L
                * self.rhoa.all_values
                * self.w.all_values
                * (self.qa.all_values - self.qs.all_values)
            )

            # Longwave radiation
            self.update_longwave_radiation(sst)

            # Net heat flux is the sum of sensible, latent, longwave fluxes
            self.shf.all_values[...] = (
                self.qh.all_values + self.qe.all_values + self.ql.all_values
            )

            # Shortwave radiation just below water surface
            if self.calculate_swr:
                self.update_shortwave_radiation(time)

            # Evaporation (derived from latent heat flux)
            # and precipitation minus evaporation
            if self.calculate_evaporation:
                np.multiply(
                    self.qe.all_values, -1.0 / (RHO0 * L), out=self.e.all_values
                )
                np.subtract(
                    self.tp.all_values, self.e.all_values, out=self.pe.all_values
                )

        super().__call__(time, sst, calculate_heat_flux)
