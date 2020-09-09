! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

SUBMODULE (getm_operators) advection_smod

!KB   use getm_domain, only: type_getm_domain
!KB   private u_advection, v_advection

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

MODULE subroutine advection_calculate_2d(self,ugrid,u,vgrid,v,dt,tgrid,f)

   IMPLICIT NONE

   class(type_advection), intent(inout) :: self
   type(type_getm_grid), intent(in) :: ugrid, vgrid
   real(real64), intent(in) :: u(:,:), v(:,:)
   real(real64), intent(in) :: dt
   type(type_getm_grid), intent(inout) :: tgrid
   real(real64), intent(inout) :: f(:,:)

!  Local constants

!  Local variables
!---------------------------------------------------------------------------
   call u_advection(tgrid%ill,tgrid%ihl,tgrid%jll,tgrid%jhl, &
                    ugrid%dx,ugrid%dy,ugrid%D,u, &
                    tgrid%inv_area,dt/2,tgrid%D,f)

!KB   call v_advection(dt,vgrid%dx,vgrid%dy,vgrid%D,v,tgrid%A,mask_flux,mask_update,tgrid%D,f)

   call u_advection(tgrid%ill,tgrid%ihl,tgrid%jll,tgrid%jhl, &
                    ugrid%dx,ugrid%dy,ugrid%D,u, &
                    tgrid%inv_area,dt/2,tgrid%D,f)
   return
END SUBROUTINE advection_calculate_2d

!---------------------------------------------------------------------------

MODULE subroutine advection_calculate_3d(self,ugrid,u,vgrid,v,dt,tgrid,f)

   IMPLICIT NONE

   class(type_advection), intent(inout) :: self
   type(type_getm_grid), intent(in) :: ugrid, vgrid
   real(real64), intent(in) :: u(:,:,:), v(:,:,:)
   real(real64), intent(in) :: dt
   type(type_getm_grid), intent(inout) :: tgrid
   real(real64), intent(inout) :: f(:,:,:)

!  Local constants

!  Local variables
   integer :: k
!---------------------------------------------------------------------------
   do k=tgrid%kmin,tgrid%kmax
      call u_advection(tgrid%ill,tgrid%ihl,tgrid%jll,tgrid%jhl, &
                       ugrid%dx,ugrid%dy,ugrid%hn(:,:,k),u(:,:,k), &
                       tgrid%inv_area,dt/2,tgrid%hn(:,:,k),f(:,:,k))
!KB      call v_advection(dt/2,vgrid%dx,vgrid%dy,vgrid%h(:,:,k),v(:,:,k),tgrid%A,mask_flux,mask_update,tgrid%h(:,:,k),f)
!KB      call w_advection(dt,)
!KB      call v_advection(dt/2,vgrid%dx,vgrid%dy,vgrid%h(:,:,k),v(:,:,k),tgrid%A,mask_flux,mask_update,tgrid%h(:,:,k),f)
      call u_advection(tgrid%ill,tgrid%ihl,tgrid%jll,tgrid%jhl, &
                       ugrid%dx,ugrid%dy,ugrid%hn(:,:,k),u(:,:,k), &
                       tgrid%inv_area,dt/2,tgrid%hn(:,:,k),f(:,:,k))
   end do
   return
END SUBROUTINE advection_calculate_3d

!---------------------------------------------------------------------------

SUBROUTINE u_advection(imin,imax,jmin,jmax,dxu,dyu,hu,u,A,dt,h,f)

   IMPLICIT NONE

   ! Subroutine arguments
   integer, intent(in) :: imin,imax,jmin,jmax
   real(real64), intent(in) :: dxu(:,:), dyu(:,:), hu(:,:), u(:,:), A(:,:)
   real(real64), intent(in) :: dt
   real(real64), intent(inout) :: h(:,:), f(:,:)

!  Local constants
   integer, parameter :: scheme=12

!  Local variables
   real(real64), allocatable, save :: uflux(:,:), QU(:,:)
   real(real64)                    :: cfl, fu, fuu, fd, hfo, advn
   integer :: i, j
!---------------------------------------------------------------------------

   ! the provided velocity MUST be ZERO at land!
!KB   QU   (imin-1,:) = 0 ! assume staggered u-fields [imin:imax]
!KB   QU   (imin:imax) = u(:) * hu(:) * dyu(:)
   uflux(imin-1:imax,:) = 0

!  TODO: (vertical) interation!!!
   do j=jmin,jmax
      do i=imin,imax-1
!KB         if (au(i).eq.1 .or. au(i).eq.2) then
            cfl = abs(u(i,j)*dt/dxu(i,j))
            if (u(i,j) .gt. 0) then
               fu  = f(i,j)
               fuu = fu
               if (i.GT.imin) then
!KB                  if (au(i-1).eq.1 .or. au(i-1).eq.2) fuu = f(i-1,j)
               end if
               fd = f(i+1,j)
            else
               fu  = f(i+1,j)
               fuu = fu
               if (i.LT.imax-1) then
!KB                  if (au(i+1).eq.1 .or. au(i+1).eq.2) fuu = f(i+2,j)
               end if
               fd = f(i,j)
            end if
            fu = adv_reconstruct(scheme,cfl,fuu,fu,fd)
!KB            uflux(i,j) = QU(i)*fu
!KB         end if
      end do
   end do

   do j=jmin,jmax
      do i=imin,imax
!KB         if (az(i).eq.1) then
            hfo = h(i,j)*f(i,j)
            h(i,j) = h(i,j) - dt*( QU(i,j)-QU(i-1,j) )*A(i,j)
            advn = ( uflux(i,j)-uflux(i-1,j) )*A(i,j) ! KB
            f(i,j) = ( hfo - dt*advn ) / h(i,j)
!KB         end if
      end do
   end do

   return
END SUBROUTINE u_advection

#if 0
!---------------------------------------------------------------------------

SUBROUTINE v_advection(dt, f, h, arcd1, u, hu, dyu, dxu, mask_flux, mask_update)
!KBSUBROUTINE v_advection(dt,dxv,dyv,hv,v,A,mask_flux,mask_update,h,f)

   IMPLICIT NONE

   ! Subroutine arguments
   real(real64), intent(in) :: dt
   real(real64), intent(in) :: dxv(:,:), dyv(:,:), hv(:,:), v(:,:), A(:,:)
   logical, intent(in)      :: mask_flux(:,:), mask_update(:,:)
   real(real64), intent(inout) :: h(:,:), f(:,:)

!  Local constants
   integer, parameter              :: scheme=12

!  Local variables
   real(real64), allocatable, save :: uflux(:), QU(:)
      real(real64)                    :: cfl, fu, fuu, fd, hfo, advn
   integer :: i,j
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
   do j=jmin,jmax
      do i=imin,imax-1
         if (mask_flux(i)) then
            cfl = abs(u(i)*dt/dxu(i))
            if (u(i) .gt. 0) then
               fu  = f(i)
               fuu = fu
               if (i.GT.imin) then
                  if (mask_flux(i-1)) fuu = f(i-1)
               end if
               fd = f(i+1)
            else
               fu  = f(i+1)
               fuu = fu
               if (i.LT.imax-1) then
                  if (mask_flux(i+1)) fuu = f(i+2)
               end if
               fd = f(i)
            end if
            fu = adv_reconstruct(scheme,cfl,fuu,fu,fd)
            uflux(i) = QU(i)*fu
         end if
      end do
   end do

   do j=jmin,jmax
      do i=imin,imax
         if (mask_update(i)) then
            hfo = h(i,j)*f(i,j)
            h(i,j) = h(i,j) - dt*( QU(i,j)-QU(i,j-1) )*arcd1(i,j)
            advn = ( uflux(i,j)-uflux(i,j-1) )*arcd1(i,j)
            f(i,j) = ( hfo - dt*advn ) / h(i,j)
         end if
      end do
   end do

   return
END SUBROUTINE v_advection
#endif

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
