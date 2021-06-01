! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

SUBMODULE (getm_operators) advection_smod
!>  @bug
!>  3D advection is not enabled yet
!>  @endbug

   IMPLICIT NONE

!-----------------------------------------------------------------------------

!INTERFACE
!END INTERFACE

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
END SUBROUTINE advection_initialize

!---------------------------------------------------------------------------

MODULE PROCEDURE advection_calculate_2d

   IMPLICIT NONE

!  Local constants

!  Local variables
   real(real64), allocatable :: D(:,:)
!---------------------------------------------------------------------------
   if (scheme == 0) return
   allocate(D,mold=tgrid%D)
   D = tgrid%D

   select case (scheme)
      case (HSIMT)
         call u_advection_hsimt(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax,tgrid%halo, &
                          ugrid%mask,ugrid%dx,ugrid%dy,ugrid%D,u, &
                          tgrid%mask,tgrid%iarea,dt/2,D,f)
         call v_advection_hsimt(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax,tgrid%halo, &
                          vgrid%mask,vgrid%dx,vgrid%dy,vgrid%D,v, &
                          tgrid%mask,tgrid%iarea,dt,D,f)
         call u_advection_hsimt(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax,tgrid%halo, &
                          ugrid%mask,ugrid%dx,ugrid%dy,ugrid%D,u, &
                          tgrid%mask,tgrid%iarea,dt/2,D,f)
      case (MUSCL)
         call u_advection_muscl(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax,tgrid%halo, &
                          ugrid%mask,ugrid%dx,ugrid%dy,ugrid%D,u, &
                          tgrid%mask,tgrid%iarea,dt/2,D,f)
         call v_advection_muscl(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax,tgrid%halo, &
                          vgrid%mask,vgrid%dx,vgrid%dy,vgrid%D,v, &
                          tgrid%mask,tgrid%iarea,dt,D,f)
         call u_advection_muscl(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax,tgrid%halo, &
                          ugrid%mask,ugrid%dx,ugrid%dy,ugrid%D,u, &
                          tgrid%mask,tgrid%iarea,dt/2,D,f)
      case (P2_PDM)
         call u_advection_p2_pdm(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax,tgrid%halo, &
                          ugrid%mask,ugrid%dx,ugrid%dy,ugrid%D,u, &
                          tgrid%mask,tgrid%iarea,dt/2,D,f)
         call v_advection_p2_pdm(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax,tgrid%halo, &
                          vgrid%mask,vgrid%dx,vgrid%dy,vgrid%D,v, &
                          tgrid%mask,tgrid%iarea,dt,D,f)
         call u_advection_p2_pdm(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax,tgrid%halo, &
                          ugrid%mask,ugrid%dx,ugrid%dy,ugrid%D,u, &
                          tgrid%mask,tgrid%iarea,dt/2,D,f)
      case (SPLMAX13)
         call u_advection_splmax13(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax,tgrid%halo, &
                          ugrid%mask,ugrid%dx,ugrid%dy,ugrid%D,u, &
                          tgrid%mask,tgrid%iarea,dt/2,D,f)
         call v_advection_splmax13(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax,tgrid%halo, &
                          vgrid%mask,vgrid%dx,vgrid%dy,vgrid%D,v, &
                          tgrid%mask,tgrid%iarea,dt,D,f)
         call u_advection_splmax13(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax,tgrid%halo, &
                          ugrid%mask,ugrid%dx,ugrid%dy,ugrid%D,u, &
                          tgrid%mask,tgrid%iarea,dt/2,D,f)
      case (SUPERBEE)
         call u_advection_superbee(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax,tgrid%halo, &
                          ugrid%mask,ugrid%dx,ugrid%dy,ugrid%D,u, &
                          tgrid%mask,tgrid%iarea,dt/2,D,f)
         call v_advection_superbee(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax,tgrid%halo, &
                          vgrid%mask,vgrid%dx,vgrid%dy,vgrid%D,v, &
                          tgrid%mask,tgrid%iarea,dt,D,f)
         call u_advection_superbee(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax,tgrid%halo, &
                          ugrid%mask,ugrid%dx,ugrid%dy,ugrid%D,u, &
                          tgrid%mask,tgrid%iarea,dt/2,D,f)
      case (UPSTREAM)
         call u_advection_upstream(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax,tgrid%halo, &
                          ugrid%mask,ugrid%dx,ugrid%dy,ugrid%D,u, &
                          tgrid%mask,tgrid%iarea,dt/2,D,f)
         call v_advection_upstream(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax,tgrid%halo, &
                          vgrid%mask,vgrid%dx,vgrid%dy,vgrid%D,v, &
                          tgrid%mask,tgrid%iarea,dt,D,f)
         call u_advection_upstream(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax,tgrid%halo, &
                          ugrid%mask,ugrid%dx,ugrid%dy,ugrid%D,u, &
                          tgrid%mask,tgrid%iarea,dt/2,D,f)
   end select
END PROCEDURE advection_calculate_2d

!---------------------------------------------------------------------------

MODULE PROCEDURE advection_calculate_3d

   IMPLICIT NONE

!  Local constants

!  Local variables
   integer :: k
   real(real64), allocatable :: hn(:,:,:)
!---------------------------------------------------------------------------
   if (scheme == 0) return
   allocate(hn,mold=tgrid%hn)
   hn = tgrid%hn
#if 0
   select case (scheme)
      case (SUPERBEE)
         do k=tgrid%kmin,tgrid%kmax
            call u_advection_superbee(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax, &
                             ugrid%mask,ugrid%dx,ugrid%dy,ugrid%hn(:,:,k),u(:,:,k), &
                             tgrid%mask,tgrid%iarea,dt/2,hn(:,:,k),f(:,:,k))
            call v_advection_superbee(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax, &
                             vgrid%mask,vgrid%dx,vgrid%dy,vgrid%hn(:,:,k),v(:,:,k), &
                             tgrid%mask,tgrid%iarea,dt/2,hn(:,:,k),f(:,:,k))
         end do
         call w_advection_superbee(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax, &
                                   tgrid%kmax, w, tgrid%mask, dt, hn, f)
         do k=tgrid%kmin,tgrid%kmax
            call v_advection_superbee(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax, &
                             vgrid%mask,vgrid%dx,vgrid%dy,vgrid%hn(:,:,k),v(:,:,k), &
                             tgrid%mask,tgrid%iarea,dt/2,hn(:,:,k),f(:,:,k))
            call u_advection_superbee(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax, &
                             ugrid%mask,ugrid%dx,ugrid%dy,ugrid%hn(:,:,k),u(:,:,k), &
                             tgrid%mask,tgrid%iarea,dt/2,hn(:,:,k),f(:,:,k))
         end do
      case (UPSTREAM)
         do k=tgrid%kmin,tgrid%kmax
            call u_advection_upstream(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax, &
                             ugrid%mask,ugrid%dx,ugrid%dy,ugrid%hn(:,:,k),u(:,:,k), &
                             tgrid%mask,tgrid%iarea,dt/2,hn(:,:,k),f(:,:,k))
            call v_advection_upstream(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax, &
                             vgrid%mask,vgrid%dx,vgrid%dy,vgrid%hn(:,:,k),v(:,:,k), &
                             tgrid%mask,tgrid%iarea,dt/2,hn(:,:,k),f(:,:,k))
         end do
         call w_advection_upstream(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax, &
                                   tgrid%kmax, w, tgrid%mask, dt, hn, f)
         do k=tgrid%kmin,tgrid%kmax
            call v_advection_upstream(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax, &
                             vgrid%mask,vgrid%dx,vgrid%dy,vgrid%hn(:,:,k),v(:,:,k), &
                             tgrid%mask,tgrid%iarea,dt/2,hn(:,:,k),f(:,:,k))
            call u_advection_upstream(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax, &
                             ugrid%mask,ugrid%dx,ugrid%dy,ugrid%hn(:,:,k),u(:,:,k), &
                             tgrid%mask,tgrid%iarea,dt/2,hn(:,:,k),f(:,:,k))
         end do
   end select
#endif
END PROCEDURE advection_calculate_3d

!---------------------------------------------------------------------------

END SUBMODULE advection_smod
