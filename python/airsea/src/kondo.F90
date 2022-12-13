module mod_kondo

   use mod_airsea_variables, only: rk

   implicit none

   private

contains

   subroutine kondo_2d(nx, ny, istart, istop, jstart, jstop, tw, ta, w, cd_mom, cd_sensible, cd_latent) bind(c)
      integer,  intent(in), value                :: nx, ny, istart, istop, jstart, jstop
      real(rk), intent(in),    dimension(nx, ny) :: tw, ta, w
      real(rk), intent(inout), dimension(nx, ny) :: cd_mom, cd_latent, cd_sensible
      call kondo(tw(istart:istop, jstart:jstop), ta(istart:istop, jstart:jstop), w(istart:istop, jstart:jstop), &
         cd_mom(istart:istop, jstart:jstop), cd_sensible(istart:istop, jstart:jstop), cd_latent(istart:istop, jstart:jstop))
   end subroutine

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Heat and momemtum fluxes according to Kondo \label{sec:kondo}
!
! !INTERFACE:
   elemental subroutine kondo(sst,airt,w,cdd,chd,ced)
!
! !DESCRIPTION:
!  Based on the model sea surface temperature, the wind vector
!  at 10 m height, the air pressure at 2 m, the dry air
!  temperature and the air pressure at 2 m, and the relative
!  humidity (either directly given or recalculated from the
!  wet bulb or the dew point temperature),
!  this routine first computes the transfer coefficients for the surface
!  momentum flux vector, $(\tau_x^s,\tau_y^s)$ ($c_{dd}$),
!  the latent heat flux, $Q_e$, ($c_{ed}$)
!  and the sensible heat flux, $Q_h$,
!  ($c_{hd}$) heat flux according to the \cite{Kondo75}
!  bulk formulae. Afterwards, these fluxes are calculated according
!  to the following formulae:
!
!  \begin{equation}
!  \begin{array}{rcl}
!  \tau_x^s &=& c_{dd} \rho_a W_x W \\ \\
!  \tau_y^s &=& c_{dd} \rho_a W_y W \\ \\
!  Q_e &=& c_{ed} L \rho_a W (q_s-q_a) \\ \\
!  Q_h &=& c_{hd} C_{pa} \rho_a W (T_w-T_a)
!  \end{array}
!  \end{equation}
!
!  with the air density $\rho_a$, the wind speed at 10 m, $W$,
!  the $x$- and the $y$-component of the wind velocity vector,
!  $W_x$ and $W_y$, respectively, the specific evaporation heat of sea water,
!  $L$, the specific saturation humidity, $q_s$, the actual
!  specific humidity $q_a$, the specific heat capacity of air at constant
!  pressure, $C_{pa}$, the sea surface temperature, $T_w$ and the
!  dry air temperature, $T_a$.
!
! !INPUT PARAMETERS:
   real(rk), intent(in)                :: sst,airt,w
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   real(rk), intent(out)               :: cdd,chd,ced
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard and Karsten Bolding
!
! !LOCAL VARIABLES:
   real(rk)                  :: s,s0,x

   ! Numbers from Table AI
   real(rk), parameter, dimension(5) :: ae_d=(/ 0._rk, 0.771_rk, 0.867_rk, 1.2_rk  , 0._rk    /)
   real(rk), parameter, dimension(5) :: ae_h=(/ 0._rk, 0.927_rk, 1.15_rk , 1.17_rk , 1.652_rk /)
   real(rk), parameter, dimension(5) :: ae_e=(/ 0._rk, 0.969_rk, 1.18_rk , 1.196_rk, 1.68_rk  /)

   real(rk), parameter, dimension(5) :: be_d=(/ 1.08_rk , 0.0858_rk, 0.0667_rk, 0.025_rk ,  0.073_rk /)
   real(rk), parameter, dimension(5) :: be_h=(/ 1.185_rk, 0.0546_rk, 0.01_rk  , 0.0075_rk, -0.017_rk /)
   real(rk), parameter, dimension(5) :: be_e=(/ 1.23_rk , 0.0521_rk, 0.01_rk  , 0.008_rk , -0.016_rk /)

   real(rk), parameter, dimension(5) :: ce_h=(/ 0._rk, 0._rk, 0._rk, -0.00045_rk, 0._rk /)
   real(rk), parameter, dimension(5) :: ce_e=(/ 0._rk, 0._rk, 0._rk, -0.0004_rk , 0._rk /)

   real(rk), parameter, dimension(5) :: pe_d=(/ -0.15_rk , 1._rk, 1._rk, 1._rk, 1._rk /)
   real(rk), parameter, dimension(5) :: pe_h=(/ -0.157_rk, 1._rk, 1._rk, 1._rk, 1._rk /)
   real(rk), parameter, dimension(5) :: pe_e=(/ -0.16_rk , 1._rk, 1._rk, 1._rk, 1._rk /)
!EOP
!-----------------------------------------------------------------------
!BOC
!  Stability
   s0 = 0.25_rk * (sst - airt) / (w + 1.0e-10_rk)**2
   s = s0 * abs(s0) / (abs(s0) + 0.01_rk)

! Transfer coefficient for heat and momentum
!   if (w < 0.3) then
   if (w < 2.2_rk) then
      cdd = (be_d(1) * w**pe_d(1)) * 1.0e-3_rk
      chd = (be_h(1) * w**pe_h(1)) * 1.0e-3_rk
      ced = (be_e(1) * w**pe_e(1)) * 1.0e-3_rk
   else if (w < 5.0_rk) then
      cdd = (ae_d(2) + be_d(2) * w) * 1.0e-3_rk
      chd = (ae_h(2) + be_h(2) * w) * 1.0e-3_rk
      ced = (ae_e(2) + be_e(2) * w) * 1.0e-3_rk
   else if (w < 8.0_rk) then
      cdd = (ae_d(3) + be_d(3) * w) * 1.0e-3_rk
      chd = (ae_h(3) + be_h(3) * w) * 1.0e-3_rk
      ced = (ae_e(3) + be_e(3) * w) * 1.0e-3_rk
   else if (w < 25.0_rk) then
      cdd = (ae_d(4) + be_d(4) * w) * 1.0e-3_rk
      chd = (ae_h(4) + be_h(4) * w + ce_h(4) * (w - 8.0_rk)**2) * 1.0e-3_rk
      ced = (ae_e(4) + be_e(4) * w + ce_e(4) * (w - 8.0_rk)**2) * 1.0e-3_rk
   else
      cdd = (ae_d(5) + be_d(5) * w) * 1.0e-3_rk
      chd = (ae_h(5) + be_h(5) * w) * 1.0e-3_rk
      ced = (ae_e(5) + be_e(5) * w) * 1.0e-3_rk
   end if

   if(s < 0._rk) then
      if (s > -3.3_rk) then
         x = 0.1_rk + 0.03_rk * s + 0.9_rk * exp(4.8_rk * s)
      else
         x = 0.0_rk
      end if
      cdd = x * cdd
      chd = x * chd
      ced = x * ced
   else
      cdd = cdd * (1.0_rk + 0.47_rk * sqrt(s))
      chd = chd * (1.0_rk + 0.63_rk * sqrt(s))
      ced = ced * (1.0_rk + 0.63_rk * sqrt(s))
   end if

   end subroutine kondo
!EOC

end module

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
