! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

SUBMODULE (getm_operators) advection_smod
!>  @bug
!>  3D advection is not enabled yet
!>  @endbug

   use advection_hsimt
   use advection_muscl
   use advection_p2_pdm
   use advection_splmax13
   use advection_superbee
   use advection_upstream

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

module SUBROUTINE advection_initialize(self, scheme, tgrid)

   !! Initialize the salinity field

   IMPLICIT NONE

   ! Subroutine arguments
   class(type_advection), intent(inout) :: self
   integer, intent(in) :: scheme
   type(type_getm_grid), intent(in) :: tgrid
!   real(real64), dimension(:,:,:), intent(in) :: var

!  Local constants

!  Local variables
!---------------------------------------------------------------------------
   select case (scheme)
   case (HSIMT); allocate(type_advection_hsimt::self%op)
   case (MUSCL); allocate(type_advection_muscl::self%op)
   case (P2_PDM); allocate(type_advection_p2_pdm::self%op)
   case (SPLMAX13); allocate(type_advection_splmax13::self%op)
   case (SUPERBEE); allocate(type_advection_superbee::self%op)
   case (UPSTREAM); allocate(type_advection_upstream::self%op)
   end select
#if 0
   select case (scheme)
      case (1)
         call upstream_initialize()
   end select
#endif
   !allocate(self%D,mold=tgrid%D)
   !allocate(self%hn,mold=tgrid%hn)
END SUBROUTINE advection_initialize

!---------------------------------------------------------------------------

MODULE PROCEDURE advection_calculate_2d

   IMPLICIT NONE

!  Local constants

!  Local variables
   real(real64), allocatable :: D(:,:)
!---------------------------------------------------------------------------
   if (.not. allocated(self%op)) return

   allocate(D,mold=tgrid%D)
   D(:,:) = tgrid%D

   call self%op%u2d(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax,tgrid%halo, &
                    ugrid%mask,ugrid%dx,ugrid%dy,ugrid%D,u, &
                    tgrid%mask,tgrid%iarea,Ah,dt/2,D,f)
   call self%op%v2d(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax,tgrid%halo, &
                    vgrid%mask,vgrid%dx,vgrid%dy,vgrid%D,v, &
                    tgrid%mask,tgrid%iarea,Ah,dt,D,f)
   call self%op%u2d(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax,tgrid%halo, &
                    ugrid%mask,ugrid%dx,ugrid%dy,ugrid%D,u, &
                    tgrid%mask,tgrid%iarea,Ah,dt/2,D,f)

END PROCEDURE advection_calculate_2d

MODULE PROCEDURE advection_calculate_u2d

   IMPLICIT NONE

!---------------------------------------------------------------------------

   call self%op%u2d(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax,tgrid%halo, &
                    ugrid%mask,ugrid%dx,ugrid%dy,ugrid%D,u, &
                    tgrid%mask,tgrid%iarea,Ah,dt,D,f)

END PROCEDURE advection_calculate_u2d


MODULE PROCEDURE advection_calculate_v2d

   IMPLICIT NONE

!---------------------------------------------------------------------------

   call self%op%v2d(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax,tgrid%halo, &
                    vgrid%mask,vgrid%dx,vgrid%dy,vgrid%D,v, &
                    tgrid%mask,tgrid%iarea,Ah,dt,D,f)

END PROCEDURE advection_calculate_v2d

MODULE PROCEDURE advection_calculate_w3d

   IMPLICIT NONE

!---------------------------------------------------------------------------
   call self%op%w3d(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax,tgrid%kmax,tgrid%halo, &
                    w,tgrid%mask,dt,h,f)

END PROCEDURE advection_calculate_w3d

!---------------------------------------------------------------------------

MODULE PROCEDURE advection_calculate_3d

   IMPLICIT NONE

!  Local constants

!  Local variables
   integer :: k
   real(real64), allocatable :: hn(:,:,:)
!---------------------------------------------------------------------------
   if (.not. allocated(self%op)) return

   allocate(hn,mold=tgrid%hn)
   hn(:,:,:) = tgrid%hn
#if 1
   do k=tgrid%kmin,tgrid%kmax
!KB - evt. self%D=tgrid%hn(:,:,k) - to avoid 3d aux. variable
      call self%op%u2d(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax,tgrid%halo, &
                       ugrid%mask,ugrid%dx,ugrid%dy,ugrid%hn(:,:,k),u(:,:,k), &
                       tgrid%mask,tgrid%iarea,Ah,dt/2,hn(:,:,k),f(:,:,k))
      call self%op%v2d(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax,tgrid%halo, &
                       vgrid%mask,vgrid%dx,vgrid%dy,vgrid%hn(:,:,k),v(:,:,k), &
                       tgrid%mask,tgrid%iarea,Ah,dt/2,hn(:,:,k),f(:,:,k))
   end do
!KB   call w_advection_superbee(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax, &
!KB                             tgrid%kmax, w, tgrid%mask, dt, self%hn, f)
   do k=tgrid%kmin,tgrid%kmax
!KB - evt. self%D=tgrid%hn(:,:,k) - to avoid 3d aux. variable
      call self%op%u2d(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax,tgrid%halo, &
                       ugrid%mask,ugrid%dx,ugrid%dy,ugrid%hn(:,:,k),u(:,:,k), &
                       tgrid%mask,tgrid%iarea,Ah,dt/2,hn(:,:,k),f(:,:,k))
      call self%op%v2d(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax,tgrid%halo, &
                       vgrid%mask,vgrid%dx,vgrid%dy,vgrid%hn(:,:,k),v(:,:,k), &
                       tgrid%mask,tgrid%iarea,Ah,dt/2,hn(:,:,k),f(:,:,k))
   end do
#else
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
