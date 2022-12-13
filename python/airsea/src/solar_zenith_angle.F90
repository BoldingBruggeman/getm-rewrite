module mod_solar_zenith_angle

   use mod_airsea_variables

   implicit none

   private

   real(rk), parameter :: pi = 3.14159265358979323846_rk
   real(rk), parameter :: deg2rad = pi / 180._rk
   real(rk), parameter :: rad2deg = 180._rk / pi
   real(rk), parameter :: yrdays = 365.25_rk

contains

   subroutine solar_zenith_angle_2d(nx, ny, imin, imax, jmin, jmax, yday, hh, dlon, dlat, zenith_angle) bind(c)
      integer,  intent(in), value                :: nx, ny, imin, imax, jmin, jmax, yday
      real(rk), intent(in), value                :: hh
      real(rk), intent(in),    dimension(nx, ny) :: dlon, dlat
      real(rk), intent(inout), dimension(nx, ny) :: zenith_angle

      real(rk) :: sundec

      sundec = solar_declination_angle(yday)
      zenith_angle(imin:imax, jmin:jmax) = solar_zenith_angle(sundec, hh, dlon(imin:imax, jmin:jmax), dlat(imin:imax, jmin:jmax))
   end subroutine

   elemental real(rk) function solar_declination_angle(yday)
      integer, intent(in) :: yday

      real(rk) :: th0, th02, th03

      th0 = 2.0_rk * pi * yday / yrdays
      th02 = 2.0_rk * th0
      th03 = 3.0_rk * th0

      solar_declination_angle = 0.006918_rk - 0.399912_rk * cos(th0) + 0.070257_rk * sin(th0)   &
              - 0.006758_rk * cos(th02) + 0.000907_rk * sin(th02)                               &
              - 0.002697_rk * cos(th03) + 0.001480_rk * sin(th03)
   end function

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Calculate the solar zenith angle \label{sec:swr}
!
! !INTERFACE:
   elemental real(rk) function solar_zenith_angle(sundec, hh, dlon, dlat)
!
! !DESCRIPTION:
!  This subroutine calculates the solar zenith angle as being used both
!  in albedo_water() and shortwave_radiation(). The result is in degrees.
!
! !INPUT PARAMETERS:
   real(rk), intent(in) :: sundec, hh
   real(rk), intent(in) :: dlon, dlat
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
! !LOCAL VARIABLES:
   real(rk) :: rlon, rlat
   real(rk) :: thsun, coszen
!EOP
!-----------------------------------------------------------------------
!BOC
!  from now on everything in radians
   rlon = deg2rad * dlon
   rlat = deg2rad * dlat

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
