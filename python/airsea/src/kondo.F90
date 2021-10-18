module mod_kondo

   use mod_airsea_variables, only: rk

   implicit none

   private

contains
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
   real(rk), parameter, dimension(5) :: ae_d=(/ 0., 0.771, 0.867, 1.2  , 0.    /)
   real(rk), parameter, dimension(5) :: ae_h=(/ 0., 0.927, 1.15 , 1.17 , 1.652 /)
   real(rk), parameter, dimension(5) :: ae_e=(/ 0., 0.969, 1.18 , 1.196, 1.68  /)

   real(rk), parameter, dimension(5) :: be_d=(/ 1.08 , 0.0858, 0.0667, 0.025 ,  0.073 /)
   real(rk), parameter, dimension(5) :: be_h=(/ 1.185, 0.0546, 0.01  , 0.0075, -0.017 /)
   real(rk), parameter, dimension(5) :: be_e=(/ 1.23 , 0.0521, 0.01  , 0.008 , -0.016 /)

   real(rk), parameter, dimension(5) :: ce_h=(/ 0., 0., 0., -0.00045, 0. /)
   real(rk), parameter, dimension(5) :: ce_e=(/ 0., 0., 0., -0.0004 , 0. /)

   real(rk), parameter, dimension(5) :: pe_d=(/ -0.15 , 1., 1., 1., 1. /)
   real(rk), parameter, dimension(5) :: pe_h=(/ -0.157, 1., 1., 1., 1. /)
   real(rk), parameter, dimension(5) :: pe_e=(/ -0.16 , 1., 1., 1., 1. /)
!EOP
!-----------------------------------------------------------------------
!BOC
!  Stability
   s0=0.25*(sst-airt)/(w+1.0e-10)**2
   s=s0*abs(s0)/(abs(s0)+0.01)

! Transfer coefficient for heat and momentum
!   if (w .lt. 0.3) then
   if (w .lt. 2.2) then
      x = log(w+0.01)
      cdd=(be_d(1)*exp(pe_d(1)*x))*1.0e-3
      chd=(be_h(1)*exp(pe_h(1)*x))*1.0e-3
      ced=(be_e(1)*exp(pe_e(1)*x))*1.0e-3
   else if (w .lt. 5.0) then
      x = exp(log(w))
      cdd=(ae_d(2)+be_d(2)*x)*1.0e-3
      chd=(ae_h(2)+be_h(2)*x)*1.0e-3
      ced=(ae_e(2)+be_e(2)*x)*1.0e-3
   else if (w .lt. 8.0) then
      x = exp(log(w))
      cdd=(ae_d(3)+be_d(3)*x)*1.0e-3
      chd=(ae_h(3)+be_h(3)*x)*1.0e-3
      ced=(ae_e(3)+be_e(3)*x)*1.0e-3
   else if (w .lt. 25.0) then
      x = exp(log(w))
      cdd=(ae_d(4)+be_d(4)*x)*1.0e-3
      chd=(ae_h(4)+be_h(4)*x+ce_h(4)*(w-8.0)**2)*1.0e-3
      ced=(ae_e(4)+be_e(4)*x+ce_e(4)*(w-8.0)**2)*1.0e-3
   else
      x = exp(log(w))
      cdd=(ae_d(5)+be_d(5)*x)*1.0e-3
      chd=(ae_h(5)+be_h(5)*x)*1.0e-3
      ced=(ae_e(5)+be_e(5)*x)*1.0e-3
   end if

   if(s .lt. 0.) then
      if (s .gt. -3.3) then
         x = 0.1+0.03*s+0.9*exp(4.8*s)
      else
         x = 0.0
      end if
      cdd=x*cdd
      chd=x*chd
      ced=x*ced
   else
      cdd=cdd*(1.0+0.47*sqrt(s))
      chd=chd*(1.0+0.63*sqrt(s))
      ced=ced*(1.0+0.63*sqrt(s))
   end if

   end subroutine kondo
!EOC

end module

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
