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
   elemental real(rk) function solar_zenith_angle(yday,hh,dlon,dlat)
!
! !DESCRIPTION:
!  This subroutine calculates the solar zenith angle as being used both
!  in albedo_water() and shortwave_radiation(). The result is in degrees.
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: yday
   real(rk), intent(in)                :: hh
   real(rk), intent(in)                :: dlon,dlat
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
! !LOCAL VARIABLES:
   real(rk), parameter       :: pi=3.14159265358979323846_rk
   real(rk), parameter       :: deg2rad=pi/180._rk
   real(rk), parameter       :: rad2deg=180._rk/pi

   real(rk)                  :: rlon,rlat
   real(rk)                  :: yrdays
   real(rk)                  :: th0,th02,th03,sundec
   real(rk)                  :: thsun,coszen
!EOP
!-----------------------------------------------------------------------
!BOC
!  from now on everything in radians
   rlon = deg2rad*dlon
   rlat = deg2rad*dlat

   yrdays=365.25

   th0 = 2.*pi*yday/yrdays
   th02 = 2.*th0
   th03 = 3.*th0
!  sun declination :
   sundec = 0.006918 - 0.399912*cos(th0) + 0.070257*sin(th0)         &
           - 0.006758*cos(th02) + 0.000907*sin(th02)                 &
           - 0.002697*cos(th03) + 0.001480*sin(th03)
!  sun hour angle :
   thsun = (hh-12.)*15.*deg2rad + rlon

!  cosine of the solar zenith angle :
   coszen =sin(rlat)*sin(sundec)+cos(rlat)*cos(sundec)*cos(thsun)
   if (coszen .lt. 0._rk) coszen = 0._rk

   solar_zenith_angle = rad2deg*acos(coszen)

   end function solar_zenith_angle
!EOC
end module
!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
