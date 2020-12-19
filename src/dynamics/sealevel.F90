! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!! Calculate the time varying sealevel using ...
!! [kaj](|url|/page//science/momentum.html)


MODULE getm_sealevel

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use memory_manager
   use logging
   use field_manager
   use getm_domain, only: type_getm_domain,g

   IMPLICIT NONE

   PRIVATE  ! Private scope by default

!  Module constants

!  Module types and variables
   type, public :: type_getm_sealevel

      class(type_logging), pointer :: logs => null()
      class(type_field_manager), pointer :: fm => null()
      class(type_getm_domain), pointer :: domain

      real(real64), allocatable :: bdyz(:)

      contains

      procedure :: configure => sealevel_configure
      procedure :: initialize => sealevel_initialize
      procedure :: update => sealevel_calculate
      procedure :: boundaries => sealevel_boundaries

   end type type_getm_sealevel

ENUM, BIND(C)
   ENUMERATOR :: zero_gradient=1
   ENUMERATOR :: sommerfeld=2
   ENUMERATOR :: clamped=3
   ENUMERATOR :: flather_elev=4
END ENUM


!---------------------------------------------------------------------------

CONTAINS

!---------------------------------------------------------------------------

SUBROUTINE sealevel_configure(self,logs,fm)

   !! Configure the components belonging to the dynamics

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_sealevel), intent(inout) :: self
   class(type_logging), intent(in), target, optional :: logs
   class(type_field_manager), intent(inout), target, optional :: fm

!  Local constants

!  Local variables
!---------------------------------------------------------------------------
   if (present(logs)) then
      self%logs => logs
      call self%logs%info('sealevel_configuration()',level=2)
   end if
   if (present(fm)) then
      self%fm => fm
   end if
   return
END SUBROUTINE sealevel_configure

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
   if (associated(self%logs)) call self%logs%info('sealevel_initialize()',level=2)
   self%domain => domain
   if (self%domain%nbdyp > 0) then
      allocate(self%bdyz(self%domain%nbdyp),stat=stat)
      self%bdyz=0._real64
   end if
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
#define _U2_ self%domain%U%l(1):,self%domain%U%l(2):
   real(real64), intent(in) :: U(_U2_)
      !! X transports
#undef _U2_
#define _V2_ self%domain%V%l(1):,self%domain%V%l(2):
   real(real64), intent(in) :: V(_V2_)
      !! Y transports
#undef _V2_
#define _T2_ self%domain%T%l(1):,self%domain%T%l(2):
   real(real64), intent(in), optional :: fwf(_T2_)
      !! surface fresh water sources
#undef _T2_

!  Local constants

!  Local variables
   integer :: i,j
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('sealevel_calculate()',level=2)
   TGrid: associate( TG => self%domain%T )
   VGrid: associate( VG => self%domain%V )
   UGrid: associate( UG => self%domain%U )
   TG%zo = TG%z
   do j=TG%l(2)+1,TG%u(2)
      do i=TG%l(1)+1,TG%u(1)
         if (TG%mask(i,j) == 1) then
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
            UG%z(i,j)=max(0.25_real64*(TG%zo(i,j)+TG%zo(i+1,j)+TG%z(i,j)+TG%z(i+1,j)),-UG%H(i,j)+self%domain%Dmin)
         end if
      end do
   end do
   end associate UGrid

   ! V-points
   VG%zo = VG%z
   do j=VG%l(2),VG%u(2)-1
      do i=VG%l(1),VG%u(1)
         if (VG%mask(i,j) > 0) then
            VG%z(i,j)=max(0.25_real64*(TG%zo(i,j)+TG%zo(i,j+1)+TG%z(i,j)+TG%z(i,j+1)),-VG%H(i,j)+self%domain%Dmin)
         end if
      end do
   end do
   end associate VGrid

   ! X-points
   XGrid: associate( XG => self%domain%X )
   XG%zo = XG%z
   do j=XG%l(2),XG%u(2)-1
      do i=XG%l(1),XG%u(1)-1
         if (XG%mask(i,j) > 0) then
            XG%z(i,j)=max(0.125_real64*(TG%zo(i,j)+TG%zo(i+1,j)+TG%zo(i+1,j+1)+TG%zo(i,j+1) &
                                       +TG%z (i,j)+TG%z (i+1,j)+TG%z (i+1,j+1)+TG%z (i,j+1)), &
                          -XG%H(i,j)+self%domain%Dmin)
         end if
      end do
   end do
   end associate XGrid
   end associate TGrid
END SUBROUTINE sealevel_calculate

!---------------------------------------------------------------------------

SUBROUTINE sealevel_boundaries(self,dt,U,V,bdyu,bdyv)

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_sealevel), intent(inout) :: self
   real(real64), intent(in) :: dt
      !! timestep [s]
#define _U2_ self%domain%U%l(1):,self%domain%U%l(2):
   real(real64), intent(in) :: U(_U2_)
      !! X transports
#undef _U2_
#define _V2_ self%domain%V%l(1):,self%domain%V%l(2):
   real(real64), intent(in) :: V(_V2_)
      !! Y transports
#undef _V2_
   real(real64), intent(in) :: bdyu(:),bdyv(:)

!  Local constants

!  Local variables
   real(real64) :: a,fac=1._real64
   integer :: i,j,k,l,n
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('sealevel_boundaries()',level=2)
   Domain: associate( domain => self%domain )
   TGrid: associate( TG => self%domain%T )
   UGrid: associate( UG => self%domain%U )
   VGrid: associate( VG => self%domain%V )
   l=0
   do n=1,domain%nwb
      l=l+1
      k=domain%bdy_index(l)
      i=domain%wi(n)
      do j=domain%wfj(n),domain%wlj(n)
         select case (domain%bdy_2d_type(l))
            case (zero_gradient)
               TG%z(i,j)=TG%z(i+1,j)
            case (sommerfeld)
!              KK-TODO: change DXC to DXU ?!
!                       change D(i,j) to _HALF_*(D(i,j)+D(i+1,j)) ?
               TG%z(i,j)=TG%z(i,j)+dt*sqrt(g*TG%D(i,j))*(TG%z(i+1,j)-TG%z(i,j))/TG%dx(i,j)
            case (clamped)
               TG%z(i,j)=max(fac*self%bdyz(k),-TG%H(i,j)+domain%Dmin)
            case (flather_elev)
               a= sqrt(UG%D(i,j)/g)*(U(i,j)/UG%D(i,j)-bdyu(k))
               TG%z(i,j)=max(fac*(self%bdyz(k)-a),-TG%H(i,j)+domain%Dmin)
         end select
         k= k+1
      end do
   end do

   do n=1,domain%nnb
      l=l+1
      k=domain%bdy_index(l)
      j=domain%nj(n)
      do i=domain%nfi(n),domain%nli(n)
         select case (domain%bdy_2d_type(l))
            case (zero_gradient)
               TG%z(i,j)=TG%z(i,j-1)
            case (sommerfeld)
!              KK-TODO: change DYC to DYVJM1 ?! (not yet in cppdefs.h!)
!                       change D(i,j) to _HALF_*(D(i,j-1)+D(i,j)) ?
               TG%z(i,j)=TG%z(i,j)-dt*sqrt(g*TG%D(i,j))*(TG%z(i,j)-TG%z(i,j-1))/TG%dy(i,j)
            case (clamped)
               TG%z(i,j)=max(fac*self%bdyz(k),-TG%H(i,j)+domain%Dmin)
            case (flather_elev)
               a=sqrt(VG%D(i,j)/g)*(V(i,j-1)/VG%D(i,j-1)-bdyv(k))
               TG%z(i,j)=max(fac*(self%bdyz(k)+a),-TG%H(i,j)+domain%Dmin)
         end select
         k=k+1
      end do
   end do

   do n=1,domain%neb
      l=l+1
      k=domain%bdy_index(l)
      i=domain%ei(n)
      do j=domain%efj(n),domain%elj(n)
         select case (domain%bdy_2d_type(l))
            case (zero_gradient)
               TG%z(i,j)=TG%z(i-1,j)
            case (sommerfeld)
!              KK-TODO: change DXC to DXUIM1 ?! (not yet in cppdefs.h!)
!                       change D(i,j) to _HALF_*(D(i-1,j)+D(i,j)) ?
               TG%z(i,j)=TG%z(i,j)-dt*sqrt(g*TG%D(i,j))*(TG%z(i,j)-TG%z(i-1,j))/TG%dx(i,j)
            case (clamped)
               TG%z(i,j)=max(fac*self%bdyz(k),-TG%H(i,j)+domain%Dmin)
            case (flather_elev)
               a=sqrt(UG%D(i,j)/g)*(U(i-1,j)/UG%D(i-1,j)-bdyu(k))
               TG%z(i,j)=max(fac*(self%bdyz(k)+a),-TG%H(i,j)+domain%Dmin)
         end select
         k=k+1
      end do
   end do

   do n=1,domain%nsb
      l=l+1
      k=domain%bdy_index(l)
      j=domain%sj(n)
      do i=domain%sfi(n),domain%sli(n)
         select case (domain%bdy_2d_type(l))
            case (zero_gradient)
               TG%z(i,j)=TG%z(i,j+1)
            case (sommerfeld)
!              KK-TODO: change DYC to DYV ?!
!                       change D(i,j) to _HALF_*(D(i,j)+D(i,j+1)) ?
               TG%z(i,j)=TG%z(i,j)+dt*sqrt(g*TG%D(i,j))*(TG%z(i,j+1)-TG%z(i,j))/TG%dy(i,j)
            case (clamped)
               TG%z(i,j)=max(fac*self%bdyz(k),-TG%H(i,j)+domain%Dmin)
            case (flather_elev)
               a=sqrt(VG%D(i,j)/g)*(V(i,j)/VG%D(i,j)-bdyv(k))
               TG%z(i,j)=max(fac*(self%bdyz(k)-a),-TG%H(i,j)+domain%Dmin)
         end select
         k=k+1
      end do
   end do
   end associate VGrid
   end associate UGrid
   end associate TGrid
   end associate Domain
END SUBROUTINE sealevel_boundaries

!---------------------------------------------------------------------------

END MODULE getm_sealevel
