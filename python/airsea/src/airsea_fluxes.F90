module mod_airsea_fluxes

   use mod_airsea_variables, only: rk

   implicit none

   private

contains
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Wrapper for air-sea fluxes calculations \label{sec:airsea-fluxes}
!
! !INTERFACE:
   subroutine airsea_fluxes(method,sst,airt,u10,v10,precip, &
                            evap,taux,tauy,qe,qh)
!
! !DESCRIPTION:
!  A wrapper around the different methods for calculating momentum
!  fluxes and sensible and latent heat fluxes at the air-sea interface.
!  To have a complete air-sea exchange also the short wave radiation
!  and back-wave radiation must be calculated.
!
! !USES:
   use mod_airsea_variables, only: kelvin
   use mod_airsea_variables, only: cpa,cpw
   use mod_kara,             only: temp_diff
   use mod_kara,             only: kara2000_cd, kara2000_clcs
   use mod_kara,             only: kara2005_cd, kara2005_clcs
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: method
   real(rk), intent(in)                :: sst,airt,u10,v10,precip
!
! !INPUT/OUTPUT PARAMETERS:
   real(rk), intent(inout)             :: evap
!
! !OUTPUT PARAMETERS:
   real(rk), intent(out)               :: taux,tauy,qe,qh
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
! !LOCAL VARIABLES:
   real(rk)                  :: ta,ta_k,tw,tw_k
   real(rk)                  :: cd_mom,cd_latent,cd_sensible
   real(rk)                  :: w,L,tmp
   real(rk)                  :: Wstar,Qstar,Tstar
   real(rk)                  :: td
   integer                   :: egon=0
   real(rk)                  :: qs,qa,rhoa
!EOP
!-----------------------------------------------------------------------
!BOC
   w = sqrt(u10*u10+v10*v10)
   L = (2.5-0.00234*sst)*1.e6
   if (sst .lt. 100.) then
      tw  = sst
      tw_k= sst+kelvin
   else
      tw  = sst-kelvin
      tw_k= sst
   end if

   if (airt .lt. 100.) then
      ta_k  = airt + kelvin
      ta = airt
   else
      ta  = airt - kelvin
      ta_k = airt
   end if



   egon = egon+1

   select case (method)
      case (1) ! Kondo
         call kondo(sst,airt,w,cd_mom,cd_latent,cd_sensible)
         tmp = cd_mom*rhoa*w
         taux = tmp*u10
         tauy = tmp*v10
         qe=-cd_sensible*cpa*rhoa*w*(sst-airt) ! sensible
         qh=-cd_latent*L*rhoa*w*(qs-qa)        ! latent
if (egon .eq. 1) then
write(100,*) 'X',' cdd ',' chd ',' ced '
end if
write(100,*) egon,cd_mom,cd_latent,cd_sensible
!          =           L*rhoa*Wstar*Qstar 
      case (2) ! Fairall et. all
!STDERR sst,airt
!STDERR rhoa,qa,qs
         call fairall(tw_k,ta_k,w,qa,qs,Wstar,Qstar,Tstar)
!STDERR Wstar,Qstar,Tstar
!         tmp=Wstar*Wstar/(Wspeed*Wspeed)
!STDERR 'Cd ',Wstar*Wstar/(w*w)
         tmp=rhoa*Wstar*Wstar/(w+0.1)
         taux = tmp*u10
         tauy = tmp*v10
         qh=L*rhoa*Wstar*Qstar
         qe=cpa*rhoa*Wstar*Tstar
!STDERR 'AAA ',taux,tauy
!STDERR 'AAB ',qh,qe
!stop
      case (3) ! Kara et. all 2000
         call temp_diff(airt,sst,qs,qa,td)
         call kara2000_cd(td, w, cd_mom)
         call kara2000_clcs(td, w,cd_latent,cd_sensible)
!STDERR 'CCC ',egon,cd_mom,cd_latent,cd_sensible
         tmp = cd_mom*rhoa*w
         taux = tmp*u10
         tauy = tmp*v10
         qe=-cd_sensible*cpa*rhoa*w*(sst-airt) ! sensible
         qh=-cd_latent*L*rhoa*w*(qs-qa)        ! latent
if (egon .eq. 1) then
write(100,*) 'X',' cdd ',' chd ',' ced '
end if
write(100,*) egon,cd_mom,cd_latent,cd_sensible
      case (4) ! Kara et. all 2005
         call temp_diff(airt,sst,qs,qa,td)
         call kara2005_cd(td, w, cd_mom)
         call kara2005_clcs(td, w,cd_latent,cd_sensible)
!STDERR 'CCC ',egon,cd_mom,cd_latent,cd_sensible
         tmp = cd_mom*rhoa*w
         taux = tmp*u10
         tauy = tmp*v10
         qe=-cd_sensible*cpa*rhoa*w*(sst-airt) ! sensible
         qh=-cd_latent*L*rhoa*w*(qs-qa)        ! latent
if (egon .eq. 1) then
write(100,*) 'X',' cdd ',' chd ',' ced '
end if
write(100,*) egon,cd_mom,cd_latent,cd_sensible
#if 0
      case (5) ! Kara et. all 2005 - using 10m T and q
         call kara_2005(sst,airt,u10,v10,precip,evap,taux,tauy,qe,qh)
#endif
      case default
   end select

#if 0
!  compute sensible heatflux due to rain fall
   if (rain_impact) then
!     units of qs and qa - should be kg/kg
      rainfall=precip * 1000. ! (convert from m/s to kg/m2/s)
      x1 = 2.11e-5*(ta_k/kelvin)**1.94
      x2 = 0.02411*(1.0+ta*(3.309e-3-1.44e-6*ta))/(rhoa*cpa)
      x3 = qa * L /(rgas * ta_K * ta_K)
      cd_rain = 1.0/(1.0+const06*(x3*L*x1)/(cpa*x2))
      cd_rain = cd_rain*cpw*((tw-ta) + (qs-qa)*L/cpa)
      qe = qe - rainfall * cd_rain
   end if

!  Compute Webb correction (Webb effect) to latent heat flux
   upvel=-1.61*Wstar*Qstar-(1.0+1.61*qa)*Wstar*Tstar/ta_k
   qe=qe-rhoa*L*upvel*qa

!  calculation of evaporation/condensation in m/s
   if (rain_impact .and. calc_evaporation) then
      evap = rhoa/rho_0*Wstar*Qstar
   else
      evap = _ZERO_
   end if

!  Compute momentum flux (N/m2) due to rainfall (kg/m2/s).
!  according to Caldwell and Elliott (1971, JPO)
   if ( rain_impact ) then
      tmp  = 0.85d0 * rainfall
      taux  = taux + tmp * u10
      tauy  = tauy + tmp * v10
   end if
#endif

   end subroutine airsea_fluxes
!EOC
end module
!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
