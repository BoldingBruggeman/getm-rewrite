! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

SUBMODULE (getm_operators) advection_smod

!KB   use getm_domain, only: type_getm_domain

!-----------------------------------------------------------------------------

INTERFACE
#if 0
   module subroutine upstream_initialize()
!      class(type_advection), intent(inout) :: self
   end subroutine upstream_initialize
   module subroutine upstream_calculate(self,domain,dt,var)
      class(type_advection), intent(inout) :: self
      class(type_getm_domain), intent(in) :: domain
      real(real64), intent(in) :: dt
      real(real64), dimension(:,:,:), intent(inout) :: var
   end subroutine upstream_calculate
   module subroutine init_gvc(self)
      class(type_getm_domain), intent(inout) :: self
   end subroutine init_gvc
   module subroutine do_gvc(self,dt)
      class(type_getm_domain), intent(inout) :: self
      real(real64), intent(in) :: dt
   end subroutine do_gvc
#endif
END INTERFACE

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

module SUBROUTINE advection_initialize(self,scheme)

   !! Initialize the salinity field

   IMPLICIT NONE

   ! Subroutine arguments
   class(type_advection), intent(inout) :: self
   integer, intent(in) :: scheme
!   real(real64), dimension(:,:,:), intent(in) :: var

!  Local constants

!  Local variables
!---------------------------------------------------------------------------
#if 0
   select case (scheme)
      case (1)
         call upstream_initialize()
   end select
#endif
   return
END SUBROUTINE advection_initialize

!---------------------------------------------------------------------------

MODULE PROCEDURE advection_calculate_1d

   IMPLICIT NONE

!  Local constants

!  Local variables
   real(real64), allocatable, save :: h_1d(:)
   integer :: imin, imax
!---------------------------------------------------------------------------

   imin = lbound(f,dim=1)
   imax = ubound(f,dim=1)

   if (allocated(h_1d)) then
      if (imin.LT.lbound(h_1d,dim=1) .OR. ubound(h_1d,dim=1).LT.imax) then
         deallocate(h_1d)
      end if
   end if
   if (.NOT.allocated(h_1d)) allocate(h_1d(imin:imax))

   h_1d(imin:imax) = h(:)
   call advection_directional_split(dt, f, h_1d(imin:imax), arcd1, u, hu, dyu, dxu, az, au)

   return
END PROCEDURE advection_calculate_1d

!---------------------------------------------------------------------------

MODULE PROCEDURE advection_calculate_2d

   IMPLICIT NONE

!  Local constants

!  Local variables
   real(real64), allocatable, save :: h_2d(:,:)
   integer :: imin, imax, jmin, jmax, i, j
!---------------------------------------------------------------------------

   imin = lbound(f,dim=1)
   imax = ubound(f,dim=1)
   jmin = lbound(f,dim=2)
   jmax = ubound(f,dim=2)

   if (allocated(h_2d)) then
      if (imin.LT.lbound(h_2d,dim=1) .OR. ubound(h_2d,dim=1).LT.imax .OR. jmin.LT.lbound(h_2d,dim=2) .OR. ubound(h_2d,dim=2).LT.jmax) then
         deallocate(h_2d)
      end if
   end if
   if (.NOT.allocated(h_2d)) allocate(h_2d(imin:imax,jmin:jmax))

   h_2d(imin:imax,jmin:jmax) = h(:,:)
   do j=jmin,jmax
      call advection_directional_split(dt/2, f(:,j), h_2d(imin:imax,j), arcd1(:,j), u(:,j), hu(:,j), dyu(:,j), dxu(:,j), az(:,j), au(:,j))
   end do
!  TODO: halo-update
   do i=imin,imax
      call advection_directional_split(dt  , f(i,:), h_2d(i,jmin:jmax), arcd1(i,:), v(i,:), hv(i,:), dxv(i,:), dyv(i,:), az(i,:), av(:,j))
   end do
!  TODO: halo-update
   do j=jmin,jmax
      call advection_directional_split(dt/2, f(:,j), h_2d(imin:imax,j), arcd1(:,j), u(:,j), hu(:,j), dyu(:,j), dxu(:,j), az(:,j), au(:,j))
   end do

   return
END PROCEDURE advection_calculate_2d

!---------------------------------------------------------------------------

MODULE PROCEDURE advection_calculate_3d
   IMPLICIT NONE

!  Local constants

!  Local variables
   real(real64), allocatable, save :: h_3d(:,:,:), dzw(:), ones(:)
   integer, allocatable, save      :: ones_int(:)
   integer :: imin, imax, jmin, jmax, kmin, kmax, i, j, k
!---------------------------------------------------------------------------

   imin = lbound(f,dim=1)
   imax = ubound(f,dim=1)
   jmin = lbound(f,dim=2)
   jmax = ubound(f,dim=2)
   kmin = lbound(f,dim=3)
   kmax = ubound(f,dim=3)

   if (allocated(h_3d)) then
      if (imin.LT.lbound(h_3d,dim=1) .OR. ubound(h_3d,dim=1).LT.imax .OR. jmin.LT.lbound(h_3d,dim=2) .OR. ubound(h_3d,dim=2).LT.jmax .OR. kmin.LT.lbound(h_3d,dim=3) .OR. ubound(h_3d,dim=3).LT.kmax) then
         deallocate(h_3d)
      end if
   end if
   if (.NOT.allocated(h_3d)) allocate(h_3d(imin:imax,jmin:jmax,kmin:kmax))

   if (allocated(dzw)) then
      if (kmin.LT.lbound(dzw,dim=1) .OR. ubound(dzw,dim=1).LT.kmax) then
         deallocate(dzw, ones, ones_int)
      end if
   end if
   if (.NOT.allocated(dzw)) then
      allocate(dzw(kmin:kmax), ones(kmin:kmax), ones_int(kmin:kmax))
      ones(:) = 1
      ones_int(:) = 1
   end if

   h_3d(imin:imax,jmin:jmax,kmin:kmax) = h(:,:,:)
   do k=kmin,kmax
      do j=jmin,jmax
         call advection_directional_split(dt/2, f(:,j,k), h_3d(imin:imax,j,k), arcd1(:,j), u(:,j,k), hu(:,j,k), dyu(:,j), dxu(:,j), az(:,j), au(:,j))
      end do
   end do
!  TODO: halo-update
   do k=kmin,kmax
      do i=imin,imax
         call advection_directional_split(dt/2, f(i,:,k), h_3d(i,jmin:jmax,k), arcd1(i,:), v(i,:,k), hv(i,:,k), dxv(i,:), dyv(i,:), az(i,:), av(:,j))
      end do
   end do
   do j=jmin,jmax
      do i=imin,imax
         if (az(i,j).eq.1) then
            dzw(kmin:kmax-1) = 0.5 * ( h_3d(i,j,kmin:kmax-1) + h_3d(i,j,kmin+1:kmax) )
            call advection_directional_split(dt, f(i,j,:), h_3d(i,j,kmin:kmax), ones(kmin:kmax), w(i,j,:), ones(kmin:kmax), ones(kmin:kmax), dzw(kmin:kmax), ones_int(kmin:kmax), ones_int(kmin:kmax))
         end if
      end do
   end do
!  TODO: halo-update
   do k=kmin,kmax
      do i=imin,imax
         call advection_directional_split(dt/2, f(i,:,k), h_3d(i,jmin:jmax,k), arcd1(i,:), v(i,:,k), hv(i,:,k), dxv(i,:), dyv(i,:), az(i,:), av(:,j))
      end do
   end do
!  TODO: halo-update
   do k=kmin,kmax
      do j=jmin,jmax
         call advection_directional_split(dt/2, f(:,j,k), h_3d(imin:imax,j,k), arcd1(:,j), u(:,j,k), hu(:,j,k), dyu(:,j), dxu(:,j), az(:,j), au(:,j))
      end do
   end do

   return
END PROCEDURE advection_calculate_3d

!---------------------------------------------------------------------------

SUBROUTINE advection_directional_split(dt, f, h, arcd1, u, hu, dyu, dxu, az, au)

   IMPLICIT NONE

   ! Subroutine arguments
   real(real64), intent(in) :: dt
   real(real64), intent(in) :: arcd1(:), u(:), hu(:), dyu(:), dxu(:)
   integer, intent(in)      :: az(:), au(:)
   real(real64), intent(inout) :: f(:), h(:)

!  Local constants
   integer, parameter              :: scheme=12

!  Local variables
   real(real64), allocatable, save :: uflux(:), QU(:)
   real(real64)                    :: cfl, fu, fuu, fd, hfo, advn
   integer :: i, imin, imax
!---------------------------------------------------------------------------

! u is velocity!!!

   imin = lbound(f,dim=1)
   imax = ubound(f,dim=1)

   if (allocated(uflux)) then
      if (imin-1.LT.lbound(uflux,dim=1) .OR. ubound(uflux,dim=1).LT.imax) then
         deallocate(uflux, QU)
      end if
   end if
   if (.NOT.allocated(uflux)) allocate(uflux(imin-1:imax), QU(imin-1:imax))

!  the provided velocity MUST be ZERO at land!
   QU   (imin-1) = 0 ! assume staggered u-fields [imin:imax]
   QU   (imin:imax) = u(:) * hu(:) * dyu(:)
   uflux(imin-1:imax) = 0

!  TODO: (vertical) interation!!!
   do i=imin,imax-1
      if (au(i).eq.1 .or. au(i).eq.2) then
         if (u(i) .gt. 0) then
            cfl = u(i)*dt/dxu(i)
            fu  = f(i)
            fuu = fu
            if (i.GT.imin) then
               if (au(i-1).eq.1 .or. au(i-1).eq.2) fuu = f(i-1)
            end if
            fd = f(i+1)
         else
            cfl = -u(i)*dt/dxu(i)
            fu  = f(i+1)
            fuu = fu
            if (i.LT.imax-1) then
               if (au(i+1).eq.1 .or. au(i+1).eq.2) fuu = f(i+2)
            end if
            fd = f(i)
         end if
         fu = adv_reconstruct(scheme,cfl,fuu,fu,fd)
         uflux(i) = QU(i)*fu
      end if
   end do

   do i=imin,imax
      if (az(i).eq.1) then
         hfo = h(i)*f(i)
         h(i) = h(i) - dt*( QU(i)-QU(i-1) )*arcd1(i)
         advn = ( uflux(i)-uflux(i-1) )*arcd1(i)
         f(i) = ( hfo - dt*advn ) / h(i)
      end if
   end do

   return
END SUBROUTINE advection_directional_split

!---------------------------------------------------------------------------

PURE real(real64) FUNCTION adv_reconstruct(scheme,cfl,fuu,fu,fd)

   IMPLICIT NONE

   integer,intent(in)  :: scheme
   real(real64),intent(in) :: cfl,fuu,fu,fd
!
! !DEFINED PARAMETERS:
   integer, parameter :: SUPERBEE=4, P2_PDM=6, SPLMAX13=13, HSIMT=12, MUSCL=5
   real(real64),parameter :: one3rd = 1.0_real64 / 3
   real(real64),parameter :: one6th = 1.0_real64 / 6
! !LOCAL VARIABLES:
   real(real64)           :: ratio,limiter,kappa,x,deltaf,deltafu
!
!---------------------------------------------------------------------------

   adv_reconstruct = fu

   deltaf  = fd - fu
   deltafu = fu - fuu

   if (deltaf*deltafu .gt. 0) then

      ratio = deltafu / deltaf   ! slope ratio

      select case (scheme)
         case (SUPERBEE)
            limiter = MAX( MIN( 2*ratio , 1.0_real64 ) , MIN( ratio , 2.0_real64 ) )
         case (P2_PDM)
            x = one6th*(1-2*cfl)
            limiter = (0.5+x) + (0.5-x)*ratio
            limiter = MIN( 2*ratio/(cfl+1.d-10) , limiter , 2/(1-cfl) )
         case (SPLMAX13)
            limiter = MIN( 2*ratio , one3rd*MAX( 1+2*ratio , 2+ratio ) , 2.0_real64 )
         case (HSIMT) ! Wu and Zhu (2010, OCEMOD)
            kappa = 1 - cfl
            x = 0.25*( kappa - one3rd/kappa )
            limiter = (0.5+x) + (0.5-x)*ratio ! can become negative!!!
            limiter = MAX( 0.0_real64 , limiter )
            limiter = MIN( 2*ratio , limiter , 2.0_real64 )
         case (MUSCL)
            limiter = MIN( 2*ratio , 0.5*(1+ratio) , 2.0_real64 )
         case default
!           UPSTREAM
            limiter = 0
      end select

      adv_reconstruct = fu + 0.5*limiter*(1-cfl)*deltaf

   end if

   return
END FUNCTION adv_reconstruct

!---------------------------------------------------------------------------

END SUBMODULE advection_smod
