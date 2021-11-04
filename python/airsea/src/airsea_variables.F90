!-----------------------------------------------------------------------
!BOP
!
! !MODULE: airsea_variables\label{sec:airsea-variables}
!
! !INTERFACE:
   module mod_airsea_variables
!
! !DESCRIPTION:
!
!  Here, number of public variables in the airsea module is declared.
!
! !USES:
   IMPLICIT NONE

!  default: all is private.
   private
   
   integer, public, parameter :: rk = selected_real_kind(13)
!
! !PUBLIC DATA MEMBERS:
   real(rk), public, parameter         :: cpa=1008.
   real(rk), public, parameter         :: cpw=3985.
   real(rk), public, parameter         :: emiss=0.97
   real(rk), public, parameter         :: bolz=5.67e-8
   real(rk), public, parameter         :: kelvin=273.15
   real(rk), public, parameter         :: const06=0.62198
   real(rk), public, parameter         :: rgas = 287.1    !
   real(rk), public, parameter         :: g = 9.81        ! [m/s2]
   real(rk), public, parameter         :: rho_0 = 1025.   ! [kg/m3]
   real(rk), public, parameter         :: kappa = 0.41    ! von Karman
!   real(rk), public                    :: es
!   real(rk), public                    :: ea
!   real(rk), public                    :: qs
!   real(rk), public, target            :: qa              ! specific humidity (kg/kg)
!   real(rk), public                    :: L
!   real(rk), public                    :: rhoa
!   real(rk), public, target            :: ta              ! 2m air temperature (degree_Celsius)
!   logical, public                     :: rain_impact
!   logical, public                     :: calc_evaporation
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding, Hans Burchard
!
!EOP
!-----------------------------------------------------------------------

   end module

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
