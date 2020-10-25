! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!! Calculate the time varying sealevel using ...

MODULE getm_sealevel

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use memory_manager
   use logging
   use field_manager
   use getm_domain, only: type_getm_domain

   IMPLICIT NONE

   PRIVATE  ! Private scope by default

!  Module constants

!  Module types and variables
   type, public :: type_getm_sealevel

      class(type_logging), pointer :: logs
      class(type_field_manager), pointer :: fm
      class(type_getm_domain), pointer :: domain

      contains

      procedure :: configuration => sealevel_configuration
      procedure :: initialize => sealevel_initialize
      procedure :: update => sealevel_calculate

   end type type_getm_sealevel

!---------------------------------------------------------------------------

CONTAINS

!---------------------------------------------------------------------------

SUBROUTINE sealevel_configuration(self,logs,fm)

   !! Configure the components belonging to the dynamics

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_sealevel), intent(inout) :: self
   class(type_logging), intent(in), target :: logs
   class(type_field_manager), intent(inout), target :: fm

!  Local constants

!  Local variables
!---------------------------------------------------------------------------
   self%logs => logs
   call self%logs%info('sealevel_configuration()',level=2)
   self%fm => fm
   return
END SUBROUTINE sealevel_configuration

!---------------------------------------------------------------------------

SUBROUTINE sealevel_initialize(self,domain)

   !! Initialize all dynamical components

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_sealevel), intent(inout) :: self
   class(type_getm_domain), intent(in), target :: domain

!  Local constants

!  Local variables
   integer :: stat
   type (type_field), pointer :: f
!---------------------------------------------------------------------------
   call self%logs%info('sealevel_initialize()',level=2)
   self%domain => domain
   return
END SUBROUTINE sealevel_initialize

!---------------------------------------------------------------------------

SUBROUTINE sealevel_calculate(self,dt,U,V,fwf)

   !! Sealevel calculation based on equation
   !! Here, the sea surface elevation is iterated according to the vertically
   !! integrated continuity equation given in (\ref{Elevation}) on page
   !! \pageref{Elevation}.
   !!
   !! When working with the option {\tt SLICE\_MODEL}, the elevations
   !! at $j=2$ are copied to $j=3$.
   !!
   !! Now with consideration of fresh water fluxes (precipitation and
   !! evaporation). Positive for flux into the water.

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_sealevel), intent(inout) :: self
   real(real64), intent(in) :: dt
      !! timestep [s]
#define _U_ self%domain%U%l(1):,self%domain%U%l(2):
   real(real64), dimension(:,:), intent(in) :: U(_U_)
      !! X transports
#undef _U_
#define _V_ self%domain%V%l(1):,self%domain%V%l(2):
   real(real64), dimension(:,:), intent(in) :: V(_V_)
      !! Y transports
#undef _V_
#define _T_ self%domain%T%l(1):,self%domain%T%l(2):
   real(real64), dimension(:,:), intent(in), optional :: fwf(_T_)
      !! surface fresh water sources
#undef _T_

!  Local constants

!  Local variables
   integer :: i,j
!---------------------------------------------------------------------------
   call self%logs%info('sealevel_calculate()',level=2)
   TGrid: associate( TG => self%domain%T )
   UGrid: associate( UG => self%domain%U )
   VGrid: associate( VG => self%domain%V )
   TG%zo = TG%z
   do j=TG%l(2)+1,TG%u(2)
      do i=TG%l(1)+1,TG%u(1)
         if (TG%mask(i,j) > 0) then
            TG%z(i,j)=TG%z(i,j) &
                     -dt*((U(i,j)*UG%dy(i,j)-U(i-1,j  )*UG%dy(i-1,j)) &
                         +(V(i,j)*VG%dx(i,j)-V(i  ,j-1)*VG%dx(i,j-1))) &
                         *TG%inv_area(i,j)
!                         *TG%inv_area(i,j) &
!                         +dt*fwf(i,j)
         end if
      end do
   end do

! U-points
   UG%zo = UG%z
   do j=UG%l(2),UG%u(2)
      do i=UG%l(1),UG%u(1)-1
         if (UG%mask(i,j) > 0) then
            UG%z(i,j)=max(0.25_real64*(TG%zo(i,j)+TG%zo(i+1,j) &
                                      +TG%z(i,j)+TG%z(i+1,j)), &
                                      -UG%H(i,j)+self%domain%Dmin)
         end if
      end do
   end do

! V-points
   VG%zo = VG%z
   do j=VG%l(2),VG%u(2)
      do i=VG%l(1),VG%u(1)-1
         if (VG%mask(i,j) > 0) then
            VG%z(i,j)=max(0.25_real64*(TG%zo(i,j)+TG%zo(i,j+1) &
                                      +TG%z(i,j)+TG%z(i,j+1)), &
                                      -VG%H(i,j)+self%domain%Dmin)
         end if
      end do
   end do
   end associate VGrid
   end associate UGrid
   end associate TGrid
   return
END SUBROUTINE sealevel_calculate

!---------------------------------------------------------------------------

END MODULE getm_sealevel
