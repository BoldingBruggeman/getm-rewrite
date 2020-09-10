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

ENUM, BIND(C)
  ENUMERATOR :: SPLIT_ADVECTION_SCHEMES=0
  ENUMERATOR :: SUPERBEE=1
  ENUMERATOR :: P2_PDM=2
  ENUMERATOR :: SPLMAX13=3
  ENUMERATOR :: HSIMT=4
  ENUMERATOR :: MUSCL=5
  ENUMERATOR :: UPSTREAM=6
END ENUM

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

MODULE subroutine advection_calculate_2d(self,scheme,ugrid,u,vgrid,v,dt,tgrid,f)

   IMPLICIT NONE

   class(type_advection), intent(inout) :: self
   integer, intent(in) :: scheme
   type(type_getm_grid), intent(in) :: ugrid, vgrid
   real(real64), intent(in) :: u(:,:), v(:,:)
   real(real64), intent(in) :: dt
   type(type_getm_grid), intent(inout) :: tgrid
   real(real64), intent(inout) :: f(:,:)

!  Local constants

!  Local variables
!---------------------------------------------------------------------------
   select case (scheme)
      case (SUPERBEE)
         call u_advection_superbee(tgrid%ill,tgrid%ihl,tgrid%jll,tgrid%jhl, &
                          ugrid%mask,ugrid%dx,ugrid%dy,ugrid%D,u, &
                          tgrid%mask,tgrid%inv_area,dt/2,tgrid%D,f)
         call v_advection_superbee()
      case (UPSTREAM)
         call u_advection_upstream(tgrid%ill,tgrid%ihl,tgrid%jll,tgrid%jhl, &
                          ugrid%mask,ugrid%dx,ugrid%dy,ugrid%D,u, &
                          tgrid%mask,tgrid%inv_area,dt/2,tgrid%D,f)
         call v_advection_upstream()
   end select
   return
END SUBROUTINE advection_calculate_2d

!---------------------------------------------------------------------------

MODULE subroutine advection_calculate_3d(self,scheme,ugrid,u,vgrid,v,dt,tgrid,f)

   IMPLICIT NONE

   class(type_advection), intent(inout) :: self
   integer, intent(in) :: scheme
   type(type_getm_grid), intent(in) :: ugrid, vgrid
   real(real64), intent(in) :: u(:,:,:), v(:,:,:)
   real(real64), intent(in) :: dt
   type(type_getm_grid), intent(inout) :: tgrid
   real(real64), intent(inout) :: f(:,:,:)

!  Local constants

!  Local variables
   integer :: k
!---------------------------------------------------------------------------
   select case (scheme)
      case (SUPERBEE)
         do k=tgrid%kmin,tgrid%kmax
            call u_advection_superbee(tgrid%ill,tgrid%ihl,tgrid%jll,tgrid%jhl, &
                             ugrid%mask,ugrid%dx,ugrid%dy,ugrid%hn(:,:,k),u(:,:,k), &
                             tgrid%mask,tgrid%inv_area,dt/2,tgrid%hn(:,:,k),f(:,:,k))
            call v_advection_superbee()
!KB      call v_advection_superbee(dt/2,vgrid%dx,vgrid%dy,vgrid%h(:,:,k),v(:,:,k),tgrid%A,mask_flux,mask_update,tgrid%h(:,:,k),f)
            call w_advection_superbee()
!KB      call v_advection_superbee(dt/2,vgrid%dx,vgrid%dy,vgrid%h(:,:,k),v(:,:,k),tgrid%A,mask_flux,mask_update,tgrid%h(:,:,k),f)
            call v_advection_superbee()
            call u_advection_superbee(tgrid%ill,tgrid%ihl,tgrid%jll,tgrid%jhl, &
                             ugrid%mask,ugrid%dx,ugrid%dy,ugrid%hn(:,:,k),u(:,:,k), &
                             tgrid%mask,tgrid%inv_area,dt/2,tgrid%hn(:,:,k),f(:,:,k))
         end do
      case (UPSTREAM)
         do k=tgrid%kmin,tgrid%kmax
            call u_advection_upstream(tgrid%ill,tgrid%ihl,tgrid%jll,tgrid%jhl, &
                             ugrid%mask,ugrid%dx,ugrid%dy,ugrid%hn(:,:,k),u(:,:,k), &
                             tgrid%mask,tgrid%inv_area,dt/2,tgrid%hn(:,:,k),f(:,:,k))
            call v_advection_upstream()
!KB      call v_advection_upstream(dt/2,vgrid%dx,vgrid%dy,vgrid%h(:,:,k),v(:,:,k),tgrid%A,mask_flux,mask_update,tgrid%h(:,:,k),f)
            call w_advection_upstream()
!KB      call v_advection_upstream(dt/2,vgrid%dx,vgrid%dy,vgrid%h(:,:,k),v(:,:,k),tgrid%A,mask_flux,mask_update,tgrid%h(:,:,k),f)
            call v_advection_upstream()
            call u_advection_upstream(tgrid%ill,tgrid%ihl,tgrid%jll,tgrid%jhl, &
                             ugrid%mask,ugrid%dx,ugrid%dy,ugrid%hn(:,:,k),u(:,:,k), &
                             tgrid%mask,tgrid%inv_area,dt/2,tgrid%hn(:,:,k),f(:,:,k))
         end do
   end select
   return
END SUBROUTINE advection_calculate_3d

!---------------------------------------------------------------------------

END SUBMODULE advection_smod
