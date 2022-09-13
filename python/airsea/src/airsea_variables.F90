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

   ! From Seawater-Ice-Air (SIA), http://www.teos-10.org/ (https://github.com/TEOS-10/SIA-Fortran, src/Constants_0.F90)
   real(rk), parameter :: gas_constant_molar_si = 8.314472_rk      !MOLAR GAS CONSTANT J MOL-1 K-1,IAPWS 2005
   real(rk), parameter :: molar_mass_h2o_si = 0.018015268_rk       !MOLAR MASS OF H2O IN KG/MOL, IAPWS 2009
   real(rk), parameter :: molar_mass_air_si = 0.02896546_rk        !MOLAR MASS OF DRY AIR IN KG MOL-1,  PICARD ET AL. 2008
   real(rk), parameter :: gas_constant_h2o_si = gas_constant_molar_si / molar_mass_h2o_si            !SPECIFIC GAS CONSTANT OF H2O IN J KG-1 K-1, IAPWS 2005
   real(rk), parameter :: gas_constant_air_si = gas_constant_molar_si / molar_mass_air_si           !SPECIFIC GAS CONSTANT OF AIR IN J KG-1 K-1, PICARD ET AL. 2008

   real(rk), public, parameter         :: cpa=1008._rk
   real(rk), public, parameter         :: cpw=3985._rk
   real(rk), public, parameter         :: emiss=0.97_rk
   real(rk), public, parameter         :: bolz=5.67e-8_rk
   real(rk), public, parameter         :: kelvin=273.15_rk
   real(rk), public, parameter         :: Ra = gas_constant_air_si ! [J kg-1 K-1]
   real(rk), public, parameter         :: Rw = gas_constant_h2o_si ! [J kg-1 K-1]
   real(rk), public, parameter         :: mw_per_ma = molar_mass_h2o_si / molar_mass_air_si
   real(rk), public, parameter         :: g = 9.81_rk        ! [m/s2]
   real(rk), public, parameter         :: rho_0 = 1025._rk   ! [kg/m3]
   real(rk), public, parameter         :: kappa = 0.41_rk    ! von Karman
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
