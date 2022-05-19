module mod_shortwave_radiation

   use mod_airsea_variables, only: rk

   implicit none

   public

contains

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Calculate the short--wave radiation \label{sec:swr}
!
! !INTERFACE:
   elemental real(rk) function shortwave_radiation(zenith_angle,yday,dlon,dlat,cloud)
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
   real(rk), intent(in)                :: dlon,dlat
   real(rk), intent(in)                :: cloud
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   real(rk), parameter       :: pi=3.14159265358979323846_rk
   real(rk), parameter       :: deg2rad=pi/180._rk
   real(rk), parameter       :: rad2deg=180._rk/pi

   real(rk), parameter       :: solar=1350._rk
   real(rk), parameter       :: eclips=23.439_rk*deg2rad
   real(rk), parameter       :: tau=0.7_rk
   real(rk), parameter       :: aozone=0.09_rk

   real(rk)                  :: coszen,sunbet
   real(rk)                  :: qatten,qzer,qdir,qdiff,qtot,qshort
   real(rk)                  :: rlon,rlat,eqnx
   real(rk)                  :: yrdays
!EOP
!-----------------------------------------------------------------------
!BOC
   coszen = cos(deg2rad*zenith_angle)
   if (coszen .le. 0.0) then
      coszen = 0.0
      qatten = 0.0
   else
      qatten = tau**(1._rk/coszen)
   end if

   qzer  = coszen * solar
   qdir  = qzer * qatten
   qdiff = ((1._rk-aozone)*qzer - qdir) * 0.5
   qtot  =  qdir + qdiff

!  from now on everything in radians
   rlon = deg2rad*dlon
   rlat = deg2rad*dlat

   yrdays=365.
   eqnx = (yday-81.)/yrdays*2.*pi
!  sin of the solar noon altitude in radians :
   sunbet=sin(rlat)*sin(eclips*sin(eqnx))+cos(rlat)*cos(eclips*sin(eqnx))
!  solar noon altitude in degrees :
   sunbet = asin(sunbet)*rad2deg

!  radiation as from Reed(1977), Simpson and Paulson(1979)
!  calculates SHORT WAVE FLUX ( watt/m*m )
!  Rosati,Miyakoda 1988 ; eq. 3.8
!  clouds from COADS perpetual data set
#if 1
   qshort  = qtot*(1-0.62*cloud + .0019*sunbet)
   if(qshort .gt. qtot ) then
      qshort  = qtot
   end if
#else
!  original implementation
   if(cloud .lt. 0.3) then
      qshort  = qtot
   else
      qshort  = qtot*(1-0.62*cloud + 0.0019*sunbet)
   endif
#endif
   shortwave_radiation = qshort

   end function shortwave_radiation
!EOC

end module
!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
