module mod_shortwave_radiation

   use mod_airsea_variables, only: rk

   implicit none

   private

contains

   subroutine shortwave_radiation_2d(nx, ny, imin, imax, jmin, jmax, yday, zenith_angle, dlat, cloud, swr) bind(c)
      integer,  intent(in), value                :: nx, ny, imin, imax, jmin, jmax, yday
      real(rk), intent(in),    dimension(nx, ny) :: zenith_angle, dlat, cloud
      real(rk), intent(inout), dimension(nx, ny) :: swr
      swr(imin:imax, jmin:jmax) = shortwave_radiation(zenith_angle(imin:imax, jmin:jmax), yday, &
         dlat(imin:imax, jmin:jmax), cloud(imin:imax, jmin:jmax))
   end subroutine

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Calculate the short--wave radiation \label{sec:swr}
!
! !INTERFACE:
   elemental real(rk) function shortwave_radiation(zenith_angle, yday, dlat, cloud)
!
! !DESCRIPTION:
!  This subroutine calculates the short--wave net radiation based on
!  solar zenith angle, year day, longitude, latitude, and fractional cloud cover.
!  No corrections for albedo - must be done by calls to {\tt albedo\_water()} and
!  if ice is included {\tt albedo\_ice()}.
!  The basic formula for the short-wave radiation at the surface, $Q_s$,
!  has been taken from \cite{RosatiMiyacoda88}, who adapted the work
!  of \cite{Reed77} and \cite{SimpsonPaulson99}:
!
!  \begin{equation}
!  Q_s=Q_{tot} (1-0.62 C + 0.0019 \beta) (1-\alpha),
!  \end{equation}
!
!  with the total radiation reaching the surface under clear skies,
!  $Q_{tot}$, the fractional cloud cover, $C$, the solar noon altitude,
!  $\beta$, and the albedo, $\alpha$.
!  This piece of code has been taken the MOM-I (Modular Ocean Model)
!  version at the INGV (Istituto Nazionale di Geofisica e Vulcanologia,
!  see {\tt http://www.bo.ingv.it/}).
!
! !INPUT PARAMETERS:
   real(rk), intent(in)                :: zenith_angle
   integer, intent(in)                 :: yday
   real(rk), intent(in)                :: dlat
   real(rk), intent(in)                :: cloud
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   real(rk), parameter       :: pi = 3.14159265358979323846_rk
   real(rk), parameter       :: deg2rad = pi / 180._rk
   real(rk), parameter       :: rad2deg = 180._rk / pi

   real(rk), parameter       :: solar = 1350._rk
   real(rk), parameter       :: eclips = 23.439_rk * deg2rad
   real(rk), parameter       :: tau = 0.7_rk
   real(rk), parameter       :: aozone = 0.09_rk
   real(rk), parameter       :: yrdays = 365._rk

   real(rk)                  :: coszen, sunbet
   real(rk)                  :: qatten, qzer, qdir, qdiff, qtot
   real(rk)                  :: rlon, rlat, eqnx
!EOP
!-----------------------------------------------------------------------
!BOC

!  Calculate downwelling shortwave at sea surface, qtot (W m-2), before reflection [albedo]
!  Rosati & Miyakoda (1988, https://doi.org/10.1175/1520-0485(1988)018<1601:AGCMFU>2.0.CO;2), Eqs 3.3 - 3.7

   coszen = cos(deg2rad * zenith_angle)
   if (coszen <= 0.0) then
      coszen = 0.0
      qatten = 0.0
   else
      qatten = tau**(1._rk / coszen)
   end if

   qzer  = coszen * solar                             ! radiation at top of the atmosphere
   qdir  = qzer * qatten                              ! direct radiation at water surface
   qdiff = ((1._rk - aozone) * qzer - qdir) * 0.5_rk  ! diffuse radiation at water surface
   qtot  =  qdir + qdiff                              ! total downwelling shortwave radiation at water surface

!  from now on everything in radians
   rlat = deg2rad * dlat

   eqnx = (yday - 81._rk) / yrdays * 2.0_rk * pi
!  sin of the solar noon altitude in radians (Rosati & Miyakoda (1988), Eqs 3.9):
   sunbet = sin(rlat) * sin(eclips * sin(eqnx)) + cos(rlat) * cos(eclips * sin(eqnx))
!  solar noon altitude in degrees :
   sunbet = asin(sunbet) * rad2deg

!  Cloud cover correction from:
!  Reed (1977), https://doi.org/10.1175/1520-0485(1977)007%3C0482:OEIOTO%3E2.0.CO;2
!  Simpson and Paulson (1979), https://doi.org/10.1002/qj.49710544412
!  Rosati & Miyakoda (1988), Eq 3.8
#if 1
   shortwave_radiation  = qtot * min(1.0_rk, 1.0_rk - 0.62_rk * cloud + 0.0019_rk * sunbet)
#else
!  original implementation
   if(cloud < 0.3_rk) then
      shortwave_radiation  = qtot
   else
      shortwave_radiation  = qtot * (1.0_rk - 0.62_rk * cloud + 0.0019_rk * sunbet)
   endif
#endif

   end function shortwave_radiation
!EOC

end module
!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
