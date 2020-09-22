! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

SUBMODULE (getm_operators) advection_smod

!-----------------------------------------------------------------------------

!   USE advection_superbee
!   USE advection_upstream

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
   module subroutine u_advection_superbee(imin,imax,jmin,jmax,umask,dxu,dyu,hu,u,tmask,A,dt,h,f)
      integer, intent(in) :: imin,imax,jmin,jmax
      integer, intent(in) :: umask(:,:)
      real(real64), intent(in) :: dxu(:,:), dyu(:,:), hu(:,:), u(:,:)
      integer, intent(in) :: tmask(:,:)
      real(real64), intent(in) :: A(:,:)
      real(real64), intent(in) :: dt
      real(real64), intent(inout) :: h(:,:), f(:,:)
   end subroutine
#endif
END INTERFACE

ENUM, BIND(C)
  ENUMERATOR :: SPLIT_ADVECTION_SCHEMES=0
  ENUMERATOR :: HSIMT=1
  ENUMERATOR :: MUSCL=2
  ENUMERATOR :: P2_PDM=3
  ENUMERATOR :: SPLMAX13=4
  ENUMERATOR :: SUPERBEE=5
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

#if 0
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
#else
MODULE PROCEDURE advection_calculate_2d
#endif
!---------------------------------------------------------------------------
   select case (scheme)
      case (HSIMT)
         call u_advection_hsimt(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax, &
                          ugrid%mask,ugrid%dx,ugrid%dy,ugrid%D,u, &
                          tgrid%mask,tgrid%inv_area,dt/2,tgrid%D,f)
         call v_advection_hsimt(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax, &
                          vgrid%mask,vgrid%dx,vgrid%dy,vgrid%D,v, &
                          tgrid%mask,tgrid%inv_area,dt,tgrid%D,f)
         call u_advection_hsimt(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax, &
                          ugrid%mask,ugrid%dx,ugrid%dy,ugrid%D,u, &
                          tgrid%mask,tgrid%inv_area,dt/2,tgrid%D,f)
      case (MUSCL)
         call u_advection_muscl(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax, &
                          ugrid%mask,ugrid%dx,ugrid%dy,ugrid%D,u, &
                          tgrid%mask,tgrid%inv_area,dt/2,tgrid%D,f)
         call v_advection_muscl(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax, &
                          vgrid%mask,vgrid%dx,vgrid%dy,vgrid%D,v, &
                          tgrid%mask,tgrid%inv_area,dt,tgrid%D,f)
         call u_advection_muscl(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax, &
                          ugrid%mask,ugrid%dx,ugrid%dy,ugrid%D,u, &
                          tgrid%mask,tgrid%inv_area,dt/2,tgrid%D,f)
      case (P2_PDM)
         call u_advection_p2_pdm(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax, &
                          ugrid%mask,ugrid%dx,ugrid%dy,ugrid%D,u, &
                          tgrid%mask,tgrid%inv_area,dt/2,tgrid%D,f)
         call v_advection_p2_pdm(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax, &
                          vgrid%mask,vgrid%dx,vgrid%dy,vgrid%D,v, &
                          tgrid%mask,tgrid%inv_area,dt,tgrid%D,f)
         call u_advection_p2_pdm(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax, &
                          ugrid%mask,ugrid%dx,ugrid%dy,ugrid%D,u, &
                          tgrid%mask,tgrid%inv_area,dt/2,tgrid%D,f)
      case (SPLMAX13)
         call u_advection_splmax13(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax, &
                          ugrid%mask,ugrid%dx,ugrid%dy,ugrid%D,u, &
                          tgrid%mask,tgrid%inv_area,dt/2,tgrid%D,f)
         call v_advection_splmax13(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax, &
                          vgrid%mask,vgrid%dx,vgrid%dy,vgrid%D,v, &
                          tgrid%mask,tgrid%inv_area,dt,tgrid%D,f)
         call u_advection_splmax13(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax, &
                          ugrid%mask,ugrid%dx,ugrid%dy,ugrid%D,u, &
                          tgrid%mask,tgrid%inv_area,dt/2,tgrid%D,f)
      case (SUPERBEE)
         call u_advection_superbee(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax, &
                          ugrid%mask,ugrid%dx,ugrid%dy,ugrid%D,u, &
                          tgrid%mask,tgrid%inv_area,dt/2,tgrid%D,f)
         call v_advection_superbee(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax, &
                          vgrid%mask,vgrid%dx,vgrid%dy,vgrid%D,v, &
                          tgrid%mask,tgrid%inv_area,dt,tgrid%D,f)
         call u_advection_superbee(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax, &
                          ugrid%mask,ugrid%dx,ugrid%dy,ugrid%D,u, &
                          tgrid%mask,tgrid%inv_area,dt/2,tgrid%D,f)
      case (UPSTREAM)
         call u_advection_upstream(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax, &
                          ugrid%mask,ugrid%dx,ugrid%dy,ugrid%D,u, &
                          tgrid%mask,tgrid%inv_area,dt/2,tgrid%D,f)
         call v_advection_upstream(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax, &
                          vgrid%mask,vgrid%dx,vgrid%dy,vgrid%D,v, &
                          tgrid%mask,tgrid%inv_area,dt/2,tgrid%D,f)
   end select
   return
#if 0
END SUBROUTINE advection_calculate_2d
#else
END PROCEDURE advection_calculate_2d
#endif

!---------------------------------------------------------------------------

#if 0
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
#else
MODULE PROCEDURE advection_calculate_3d
#endif
   integer :: k
!---------------------------------------------------------------------------
   select case (scheme)
      case (SUPERBEE)
         do k=tgrid%kmin,tgrid%kmax
            call u_advection_superbee(tgrid%ill,tgrid%ihl,tgrid%jll,tgrid%jhl, &
                             ugrid%mask,ugrid%dx,ugrid%dy,ugrid%hn(:,:,k),u(:,:,k), &
                             tgrid%mask,tgrid%inv_area,dt/2,tgrid%hn(:,:,k),f(:,:,k))
            call v_advection_superbee(tgrid%ill,tgrid%ihl,tgrid%jll,tgrid%jhl, &
                             vgrid%mask,vgrid%dx,vgrid%dy,vgrid%hn(:,:,k),v(:,:,k), &
                             tgrid%mask,tgrid%inv_area,dt/2,tgrid%hn(:,:,k),f(:,:,k))
!            call w_advection_superbee()
            call v_advection_superbee(tgrid%ill,tgrid%ihl,tgrid%jll,tgrid%jhl, &
                             vgrid%mask,vgrid%dx,vgrid%dy,vgrid%hn(:,:,k),v(:,:,k), &
                             tgrid%mask,tgrid%inv_area,dt/2,tgrid%hn(:,:,k),f(:,:,k))
            call u_advection_superbee(tgrid%ill,tgrid%ihl,tgrid%jll,tgrid%jhl, &
                             ugrid%mask,ugrid%dx,ugrid%dy,ugrid%hn(:,:,k),u(:,:,k), &
                             tgrid%mask,tgrid%inv_area,dt/2,tgrid%hn(:,:,k),f(:,:,k))
         end do
      case (UPSTREAM)
         do k=tgrid%kmin,tgrid%kmax
            call u_advection_upstream(tgrid%ill,tgrid%ihl,tgrid%jll,tgrid%jhl, &
                             ugrid%mask,ugrid%dx,ugrid%dy,ugrid%hn(:,:,k),u(:,:,k), &
                             tgrid%mask,tgrid%inv_area,dt/2,tgrid%hn(:,:,k),f(:,:,k))
            call v_advection_upstream(tgrid%ill,tgrid%ihl,tgrid%jll,tgrid%jhl, &
                             vgrid%mask,vgrid%dx,vgrid%dy,vgrid%hn(:,:,k),v(:,:,k), &
                             tgrid%mask,tgrid%inv_area,dt/2,tgrid%hn(:,:,k),f(:,:,k))
!            call w_advection_upstream()
            call v_advection_upstream(tgrid%ill,tgrid%ihl,tgrid%jll,tgrid%jhl, &
                             vgrid%mask,vgrid%dx,vgrid%dy,vgrid%hn(:,:,k),v(:,:,k), &
                             tgrid%mask,tgrid%inv_area,dt/2,tgrid%hn(:,:,k),f(:,:,k))
            call u_advection_upstream(tgrid%ill,tgrid%ihl,tgrid%jll,tgrid%jhl, &
                             ugrid%mask,ugrid%dx,ugrid%dy,ugrid%hn(:,:,k),u(:,:,k), &
                             tgrid%mask,tgrid%inv_area,dt/2,tgrid%hn(:,:,k),f(:,:,k))
         end do
   end select
   return
#if 0
END SUBROUTINE advection_calculate_3d
#else
END PROCEDURE advection_calculate_3d
#endif

!---------------------------------------------------------------------------

END SUBMODULE advection_smod
