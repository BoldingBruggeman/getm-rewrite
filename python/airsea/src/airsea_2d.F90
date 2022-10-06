
module mod_airsea_2d

   use mod_airsea_variables, only: rk
   use mod_humidity
   use mod_solar_zenith_angle
   use mod_shortwave_radiation
   use mod_longwave_radiation
   use mod_albedo_water
   use mod_kondo

   implicit none

   private

contains

   subroutine solar_zenith_angle_2d(nx, ny, istart, istop, jstart, jstop, yday, hh, dlon, dlat, zenith_angle) bind(c)
      integer,  intent(in), value                :: nx, ny, istart, istop, jstart, jstop, yday
      real(rk), intent(in), value                :: hh
      real(rk), intent(in),    dimension(nx, ny) :: dlon, dlat
      real(rk), intent(inout), dimension(nx, ny) :: zenith_angle

      real(rk) :: sundec

      sundec = solar_declination_angle(yday)
      zenith_angle(istart:istop, jstart:jstop) = solar_zenith_angle(sundec, hh, dlon(istart:istop, jstart:jstop), &
         dlat(istart:istop, jstart:jstop))
   end subroutine

   subroutine shortwave_radiation_2d(nx, ny, istart, istop, jstart, jstop, yday, zenith_angle, dlat, cloud, swr) bind(c)
      integer,  intent(in), value                :: nx, ny, istart, istop, jstart, jstop, yday
      real(rk), intent(in),    dimension(nx, ny) :: zenith_angle, dlat, cloud
      real(rk), intent(inout), dimension(nx, ny) :: swr
      swr(istart:istop, jstart:jstop) = shortwave_radiation(zenith_angle(istart:istop, jstart:jstop), yday, &
         dlat(istart:istop, jstart:jstop), cloud(istart:istop, jstart:jstop))
   end subroutine

   subroutine humidity_2d(nx, ny, istart, istop, jstart, jstop, method, hum, airp, tw, ta, es, ea, qs, qa, rhoa) bind(c)
      integer,  intent(in), value                :: nx, ny, istart, istop, jstart, jstop, method
      real(rk), intent(in),    dimension(nx, ny) :: hum, airp, tw, ta
      real(rk), intent(inout), dimension(nx, ny) :: es, ea, qs, qa, rhoa

      !  saturation vapor pressure - using SST
      !  correction for seawater, following Kraus 1972
      !  correcting for salt water assuming 98% RH
      es(istart:istop, jstart:jstop) = 0.98_rk * saturation_vapor_pressure(tw(istart:istop, jstart:jstop))

     !  saturation specific humidity
      qs(istart:istop, jstart:jstop) = specific_humidity(es(istart:istop, jstart:jstop), airp(istart:istop, jstart:jstop))

      !  must be also calcuated for airtemperature, depending on humidity input
      !  see ../ncdf/ncdf_meteo.F90 for defined constants
      select case (method)
         case (1)  ! relative humidity in % given
   !        get actual vapor pressure from saturation vapor pressure at that air temperature
            ea(istart:istop, jstart:jstop) = 0.01_rk * hum(istart:istop, jstart:jstop) &
               * saturation_vapor_pressure(ta(istart:istop, jstart:jstop))
   !        convert to specific humidity
            qa(istart:istop, jstart:jstop) = specific_humidity(ea(istart:istop, jstart:jstop), airp(istart:istop, jstart:jstop))
         case (2)  ! Specific humidity from wet bulb temperature
   !        calculate the SVP at wet bulb temp then
   !        use the psychrometer formula to get vapour pressure
   !        See Smithsonian Met tables 6th Edition pg 366 eqn 3
   !        Make sure this is in degC
   !        saturation vapor pressure at wet bulb temperature
   !        actual vapor pressure
            ea(istart:istop, jstart:jstop) = saturation_vapor_pressure(hum(istart:istop, jstart:jstop)) &
               - 6.6e-4_rk * (1._rk + 1.15e-3_rk * hum(istart:istop, jstart:jstop)) * airp &
               * (ta(istart:istop, jstart:jstop) - hum(istart:istop, jstart:jstop))
   !        specific humidity in kg/kg
            qa(istart:istop, jstart:jstop) = specific_humidity(ea(istart:istop, jstart:jstop), airp(istart:istop, jstart:jstop))
         case (3)  ! Specific humidity from dew point temperature
            ea(istart:istop, jstart:jstop) = saturation_vapor_pressure(hum(istart:istop, jstart:jstop))
            qa(istart:istop, jstart:jstop) = specific_humidity(ea(istart:istop, jstart:jstop), airp(istart:istop, jstart:jstop))
         case (4)  ! specific humidity in kg/kg is given
            qa(istart:istop, jstart:jstop) = hum(istart:istop, jstart:jstop)
   !        actual water vapor pressure in Pascal
            ea(istart:istop, jstart:jstop) = vapor_pressure(qa(istart:istop, jstart:jstop), airp(istart:istop, jstart:jstop))
         case default
            stop 'humidity_2d()'
      end select

      rhoa(istart:istop, jstart:jstop) = air_density(airp(istart:istop, jstart:jstop), ta(istart:istop, jstart:jstop), &
         qa(istart:istop, jstart:jstop))

   end subroutine

   subroutine longwave_radiation_2d(nx, ny, istart, istop, jstart, jstop, method, dlat, tw, ta, cloud, ea, qa, ql) bind(c)
      integer,  intent(in), value                :: nx, ny, istart, istop, jstart, jstop, method
      real(rk), intent(in),    dimension(nx, ny) :: dlat, tw, ta, cloud, ea, qa
      real(rk), intent(inout), dimension(nx, ny) :: ql

      select case (method)
         case (1)
            call clark(dlat(istart:istop, jstart:jstop), tw(istart:istop, jstart:jstop), ta(istart:istop, jstart:jstop), &
               cloud(istart:istop, jstart:jstop), ea(istart:istop, jstart:jstop), ql(istart:istop, jstart:jstop))
         case (2)
            call hastenrath(dlat(istart:istop, jstart:jstop), tw(istart:istop, jstart:jstop), ta(istart:istop, jstart:jstop), &
               cloud(istart:istop, jstart:jstop), qa(istart:istop, jstart:jstop), ql(istart:istop, jstart:jstop))
         case (3)
            call bignami(dlat(istart:istop, jstart:jstop), tw(istart:istop, jstart:jstop), ta(istart:istop, jstart:jstop), &
               cloud(istart:istop, jstart:jstop), ea(istart:istop, jstart:jstop), ql(istart:istop, jstart:jstop))
         case (4)
            call berliand(tw(istart:istop, jstart:jstop), ta(istart:istop, jstart:jstop), cloud(istart:istop, jstart:jstop), &
               ea(istart:istop, jstart:jstop), ql(istart:istop, jstart:jstop))
         case (5)
            call josey1(tw(istart:istop, jstart:jstop), ta(istart:istop, jstart:jstop), cloud(istart:istop, jstart:jstop), &
               ql(istart:istop, jstart:jstop))
         case (6)
            call josey2(tw(istart:istop, jstart:jstop), ta(istart:istop, jstart:jstop), cloud(istart:istop, jstart:jstop), &
               ea(istart:istop, jstart:jstop), ql(istart:istop, jstart:jstop))
         case default
            stop 'longwave_radiation_2d()'
      end select
   end subroutine

   subroutine albedo_water_2d(nx, ny, istart, istop, jstart, jstop, method, zen, yday, albedo) bind(c)
      integer, intent(in), value                 :: nx, ny, istart, istop, jstart, jstop, method, yday
      real(rk), intent(in),    dimension(nx, ny) :: zen
      real(rk), intent(inout), dimension(nx, ny) :: albedo

      select case (method)
         case (1)
            albedo(istart:istop, jstart:jstop) = albedo_payne(zen(istart:istop, jstart:jstop))
         case (2)
            albedo(istart:istop, jstart:jstop) = albedo_cogley(zen(istart:istop, jstart:jstop), yday)
         case default
            stop 'albedo_2d()'
      end select
   end subroutine

   subroutine transfer_coefficients_2d(nx, ny, istart, istop, jstart, jstop, method, tw, ta, w, cd_mom, cd_sensible, cd_latent) &
      bind(c)
      integer,  intent(in), value                :: nx, ny, istart, istop, jstart, jstop, method
      real(rk), intent(in),    dimension(nx, ny) :: tw, ta, w
      real(rk), intent(inout), dimension(nx, ny) :: cd_mom, cd_latent, cd_sensible
      call kondo(tw(istart:istop, jstart:jstop), ta(istart:istop, jstart:jstop), w(istart:istop, jstart:jstop), &
         cd_mom(istart:istop, jstart:jstop), cd_sensible(istart:istop, jstart:jstop), cd_latent(istart:istop, jstart:jstop))
   end subroutine

end module
!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
