module mod_solar_zenith_angle

   use mod_airsea_variables

   implicit none

   public

contains
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Calculate the solar zenith angle \label{sec:swr}
!
! !INTERFACE:
   elemental real(rk) function solar_zenith_angle(yday, hh, dlon, dlat)
!
! !DESCRIPTION:
!  This subroutine calculates the solar zenith angle as being used both
!  in albedo_water() and shortwave_radiation(). The result is in degrees.
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: yday
   real(rk), intent(in)                :: hh
   real(rk), intent(in)                :: dlon, dlat
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
! !LOCAL VARIABLES:
   real(rk), parameter       :: pi = 3.14159265358979323846_rk
   real(rk), parameter       :: deg2rad = pi / 180._rk
   real(rk), parameter       :: rad2deg = 180._rk / pi
   real(rk), parameter       :: yrdays = 365.25_rk

   real(rk)                  :: rlon, rlat
   real(rk)                  :: th0, th02, th03, sundec
   real(rk)                  :: thsun, coszen
!EOP
!-----------------------------------------------------------------------
!BOC
!  from now on everything in radians
   rlon = deg2rad * dlon
   rlat = deg2rad * dlat

   th0 = 2.0_rk * pi * yday / yrdays
   th02 = 2.0_rk * th0
   th03 = 3.0_rk * th0
!  sun declination :
   sundec = 0.006918_rk - 0.399912_rk * cos(th0) + 0.070257_rk * sin(th0)   &
           - 0.006758_rk * cos(th02) + 0.000907_rk * sin(th02)              &
           - 0.002697_rk * cos(th03) + 0.001480_rk * sin(th03)
!  sun hour angle :
   thsun = (hh - 12.0_rk) * 15.0_rk * deg2rad + rlon

!  cosine of the solar zenith angle :
   coszen = max(0.0_rk, sin(rlat) * sin(sundec) + cos(rlat) * cos(sundec) * cos(thsun))

   solar_zenith_angle = rad2deg * acos(coszen)

   end function solar_zenith_angle
!EOC
end module
!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
