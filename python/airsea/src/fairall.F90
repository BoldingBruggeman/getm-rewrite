module mod_fairall

   use mod_airsea_variables, only: kelvin,g,kappa, rk

   implicit none

   private

contains

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Heat and momentum fluxes according to Fairall et al.
!
! !INTERFACE:
   subroutine fairall(sst,airt,w,qa,qs,Wstar,Qstar,Tstar)
!
! !DESCRIPTION:
!  The surface momentum flux vector, $(\tau_x^s,\tau_y^s)$,
!  in [N\,m$^{-2}$],
!  the latent heat flux, $Q_e$,
!  and the sensible heat flux, $Q_h$, both in [W\,m$^{-2}$]
!  are calculated here according to the \cite{Fairalletal96a} bulk
!  formulae, which are build on the Liu-Katsaros-Businger
!  (\cite{Liuetal79}) method.
!  Cool skin and warm layer effects are considered according to the
!  suggestions of \cite{Fairalletal96b}.
!
!  The air temperature {\tt airt} and the sea surface temperature
!  {\tt sst} may be given in Kelvin or Celsius:
!  if they are $>$ 100 - Kelvin is assumed.
!
!  This piece of code has been adapted from the COARE code originally
!  written by David Rutgers and Frank Bradley - see
!  http://www.coaps.fsu.edu/COARE/flux\_algor/flux.html.
!
! !INPUT PARAMETERS:
   real(rk), intent(in)                :: sst,airt,w,qa,qs
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   real(rk), intent(out)               :: Wstar,Qstar,Tstar
!
! !REVISION HISTORY:
!  Original author(s): Adolf Stips
!
! !DEFINED PARAMETERS:
!  Fairall LKB roughness Reynolds number to Von Karman
   real(rk),parameter        :: fdg = 1.0          ! non-dimensional

!  Beta parameter evaluated from Fairall low windspeed turbulence data.
   real(rk),parameter        :: beta = 1.2         ! non-dimensional

!  Zabl      Height (m) of atmospheric boundary layer.
   real(rk),parameter        :: Zabl = 600.0       ! in [m]

   real(rk), parameter       :: r3 = 1.0/3.0
!
!  Liu et al. (1979) look-up table coefficients to compute roughness
!  Reynolds number for temperature (rt) and moisture (rq) as a
!  function of wind Reynolds number (rr):
!
!       rt = Liu_a(:,1) * Rr   ** Liu_b(:,1)    temperature
!       rq = Liu_a(:,2) * Rr   ** Liu_b(:,2)    moisture
!
   real(rk),parameter, dimension(8,2) :: Liu_a = reshape ( &
                 (/ 0.177,  1.376,    1.026,      1.625,   &
                    4.661, 34.904, 1667.190, 588000.0,     &
                    0.292,  1.808,    1.393,      1.956,   &
                    4.994, 30.709, 1448.680, 298000.0 /),  &
                 (/ 8, 2 /) )

   real(rk),parameter, dimension(8,2) :: Liu_b = reshape ( &
                 (/  0.0,    0.929, -0.599, -1.018,        &
                    -1.475, -2.067, -2.907, -3.935,        &
                     0.0,    0.826, -0.528, -0.870,        &
                    -1.297, -1.845, -2.682, -3.616 /),     &
                 (/ 8, 2 /) )

   real(rk),parameter, dimension(9) :: Liu_Rr =            &
                 (/    0.0,  0.11,   0.825,   3.0,         &
                      10.0, 30.0,  100.0,   300.0,         &
                    1000.0 /)
!
!  Height (m) of surface air temperature measurement.
   real(rk), parameter       ::  zt= 2.0
!  Height (m) of surface air humidity measurement
   real(rk), parameter       ::  zq= 2.0
!  Height (m) of surface winds measurement
   real(rk), parameter       ::  zw= 10.0
   integer,  parameter       :: itermax = 20
#ifdef GUSTINESS
   real(rk), parameter       :: wgust=0.2
#else
   real(rk), parameter       :: wgust=0.0
#endif
   real(rk),external         :: psi

! !LOCAL VARIABLES:
   real(rk)                  :: Bf,cff
   real(rk)                  :: delw,delq,delt
   real(rk)                  :: rr,rt,rq
   real(rk)                  :: Ri
   real(rk)                  :: vis_air,ta
   real(rk)                  :: L
   real(rk)                  :: TVStar,wgus
   real(rk)                  :: tpsi,qpsi,wpsi,ZWoL,oL,ZToL,ZQoL,ZoW,ZoT, ZoQ
   integer                   :: ier,iter,k
!EOP
!-----------------------------------------------------------------------
!BOC
!
!  Initialize.
!
   delw=sqrt(w*w+wgust*wgust)
   if (delw .ne. 0.0) then
!-----------------------------------------------------------------------
!     Compute Monin-Obukhov similarity parameters for wind (Wstar),
!     heat (Tstar), and moisture (Qstar), Liu et al. (1979).
!-----------------------------------------------------------------------

!     Kinematic viscosity of dry air (m2/s), Andreas (1989).
      ta = airt-kelvin
      vis_air=1.326e-5*(1.0+ta*(6.542e-3+ta*(8.301e-6-4.84e-9*ta)))

!     Compute latent heat of vaporization (J/kg) at sea surface
      L = (2.501-0.00237*sst)*1.e6
!
!     Assume that wind is measured relative to sea surface and include
!     gustiness.
!     Initialize.
      ier = 0
      delq=qa-qs
      delt=airt-sst

!     Initial guesses for Monin-Obukhov similarity scales.
      ZWoL=0.0
      ZoW=0.0005
      Wstar=0.04*delw
      Tstar=0.04*delt
      Qstar=0.04*delq
      TVstar=Tstar*(1.0+0.61*qa)+0.61*airt*Qstar

!     Compute Richardson number.
      Ri=g*zw*(delt+0.61*airt*delq)/(airt*delw*delw)

!     Fairall computes turbulent fluxes only if Ri< 0.25
      if ( ri .le. 0.25) then
!        Iterate until convergence or when IER is negative.  It usually
!        converges within four iterations.
         do iter=1,itermax
            if ( ier .ge. 0 ) then
!              Compute Monin-Obukhov stability parameter, Z/L.
               oL=g*kappa*TVstar/(airt*(1.0+0.61*qa)*Wstar*Wstar)
               ZWoL=zw*oL
               ZToL=zt*oL
               ZQoL=zq*oL

!              Evaluate stability functions at Z/L.
               wpsi=psi(1,ZWoL)
               tpsi=psi(2,ZToL)
               qpsi=psi(2,ZQoL)

!              Compute wind scaling parameters, Wstar.
               ZoW=0.011*Wstar*Wstar/g+0.11*vis_air/Wstar
               Wstar=delw*kappa/(log(zw/ZoW)-wpsi)

!              Computes roughness Reynolds number for wind (Rr), heat (Rt),
!              and moisture (Rq). Use Liu et al. (1976) look-up table to
!              compute "Rt" and "Rq" as function of "Rr".
               rr=ZoW*Wstar/vis_air
               if ((rr .ge. 0.0).and.(rr .lt. 1000.0)) then
                  do k=1,8
                     if ((liu_rr(k).le.rr).and.(rr .lt. liu_rr(k+1))) then
                        rt=liu_a(k,1)*rr**liu_b(k,1)
                        rq=liu_a(k,2)*rr**liu_b(k,2)
                     end if
                  end do

!                Compute heat and moisture scaling parameters,
!                Tstar and Qstar.
                  cff=vis_air/Wstar
                  ZoT=rt*cff
                  ZoQ=rq*cff
                  cff=kappa*fdg
                  Tstar=(delt)*cff/(log(zt/ZoT)-tpsi)
                  Qstar=(delq)*cff/(log(zq/ZoQ)-qpsi)

!                 Compute gustiness in wind speed.
                  TVstar=Tstar*(1.0+0.61*qa)+0.61*airt*Qstar
                  bf=-g/airt*Wstar*TVstar
                  if (bf .gt. 0) then
                     wgus=beta*(bf*Zabl)**r3
                  else
                     wgus=0._rk
                  end if
                  delw=sqrt(w*w+wgus*wgus)
               else
                  ier = -2
               end if
            end if
         end do
      end if ! Ri < 0.25
   end if  !delw != 0.0

   end subroutine fairall

   elemental function psi(iflag, ZoL)
!=======================================================================
!                                                                      !
!  This function evaluates the stability function, PSI, for wind       !
!  speed (iflag=1) or for air temperature and moisture (iflag=2)       !
!  profiles as function of the stability parameter, ZoL (z/L).         !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Liu, W.T., K.B. Katsaros, and J.A. Businger, 1979:  Bulk          !
!        parameterization of the air-sea exchange of heat and          !
!        water vapor including the molecular constraints at            !
!        the interface, J. Atmos. Sci, 36, 1722-1735.                  !
!                                                                      !
!=======================================================================
!
      real(rk)  :: psi
!
!  Imported variable declarations.
!
   integer,  intent(in) :: iflag
   real(rk), intent(in) :: ZoL
!
!  Local variable declarations.
!
   real(rk), parameter :: r3 = 1.0/3.0
   real(rk), parameter :: sqr3 = 1.7320508
   real(rk), parameter :: pi=3.141592653589
   real(rk)            :: Fw, chic, chik, psic, psik

!  Initialize for the zero "ZoL" case.
!
   psi=0.0
!
!  Unstable conditions.
!
   if (ZoL .lt. 0.0) then
      chik=(1.0-16.0*ZoL)**0.25
      if (iflag .eq. 1) then
         psik=2.0*LOG(0.5*(1.0+chik))+LOG(0.5*(1.0+chik*chik))-   &
              2.0*ATAN(chik)+ 0.5*pi
      else if (iflag .eq. 2) then
            psik=2.0*LOG(0.5*(1.0+chik*chik))
      end if
!
!  For very unstable conditions, use free-convection (Fairall).
!
      chic=(1.0-12.87*ZoL)**r3
      psic=1.5*LOG(r3*(1.0+chic+chic*chic))-                    &
         sqr3*ATAN((1.0+2.0*chic)/sqr3)+ pi/sqr3
!
!  Match Kansas and free-convection forms with weighting Fw.
!
      Fw=1.0/(1.0+ZoL*ZoL)
      psi=Fw*psik+(1.0-Fw)*psic
!
!  Stable conditions.
!
   else if (ZoL .gt. 0.0) then
      psi=-4.7*ZoL
   end if

   end function psi

end module

!EOC
!-----------------------------------------------------------------------
!Copyright (C) 2007 - Adolf Stips
!-----------------------------------------------------------------------
