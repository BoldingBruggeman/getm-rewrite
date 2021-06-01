! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!>  Some variables are mirrored outside the calculation domain in the
!>  vicinity of the open boundaries. This is to avoid if statements
!>  when calculating e.g. the Coriolis terms and advection.
!>  Shall be checked with the corresponding GETM code.

SUBMODULE (getm_domain) mirror_bdys_smod

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

module SUBROUTINE mirror_bdy_2d(self,grid,f)

  IMPLICIT NONE

! Subroutine arguments
   class(type_getm_domain), intent(inout) :: self
   class(type_getm_grid), intent(in) :: grid
   real(real64), dimension(:,:), intent(inout) :: f(grid%l(1):,grid%l(2):)

! Local constants

! Local variables
   integer :: i,j,n
!-----------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('mirror_bdy_2d()',level=2)

!KB
#define HALO 0
!KB
   select case (grid%grid_type)
      case (1) ! TGRID
         do n = 1,self%NWB
            i = self%wi(n)
            do j = self%wfj(n)-HALO,self%wlj(n)+HALO
               if (grid%mask(i,j) > 1) f(i-1,j) = f(i,j)
            end do
         end do
         do n = 1,self%NNB
            j = self%nj(n)
            do i = self%nfi(n)-HALO,self%nli(n)+HALO
               if (grid%mask(i,j) > 1) f(i,j+1) = f(i,j)
            end do
         end do
         do n = 1,self%NEB
            i = self%ei(n)
            do j = self%efj(n)-HALO,self%elj(n)+HALO
               if (grid%mask(i,j) > 1) f(i+1,j) = f(i,j)
            end do
         end do
         do n = 1,self%NSB
            j = self%sj(n)
            do i = self%sfi(n)-HALO,self%sli(n)+HALO
               if (grid%mask(i,j) > 1) f(i,j-1) = f(i,j)
            end do
         end do
      case (2) ! UGRID
         do n = 1,self%NNB
            j = self%nj(n)
            do i = self%nfi(n)-HALO,self%nli(n)+HALO
               if (grid%mask(i,j) == 3) f(i,j) = f(i,j-1)
             end do
         end do
         do n = 1,self%NSB
            j = self%sj(n)
            do i = self%sfi(n)-HALO,self%sli(n)+HALO
               if (grid%mask(i,j) == 3) f(i,j) = f(i,j+1)
            end do
         end do
      case (3) ! VGRID
         do n = 1,self%NWB
            i = self%wi(n)
            do j = self%wfj(n)-HALO,self%wlj(n)+HALO
               if (grid%mask(i,j) == 3) f(i,j) = f(i+1,j)
            end do
         end do
         do n = 1,self%NEB
            i = self%ei(n)
            do j = self%efj(n)-HALO,self%elj(n)+HALO
               if (grid%mask(i,j) == 3) f(i,j) = f(i-1,j)
            end do
         end do
   end select
END SUBROUTINE mirror_bdy_2d

!-----------------------------------------------------------------------------

module SUBROUTINE mirror_bdy_3d(self,grid,f)

  IMPLICIT NONE

! Subroutine arguments
   class(type_getm_domain), intent(inout) :: self
   class(type_getm_grid), intent(in) :: grid
   real(real64), dimension(:,:,:), intent(inout) :: f(grid%l(1):,grid%l(2):,grid%l(3):)

! Local constants

! Local variables
   integer :: i,j,n
!-----------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('mirror_bdy_3d()',level=2)

   select case (grid%grid_type)
      case (1) ! TGRID
         do n = 1,self%NWB
            i = self%wi(n)
            do j = self%wfj(n)-HALO,self%wlj(n)+HALO
               if (grid%mask(i,j) > 1) f(i-1,j,:) = f(i,j,:)
            end do
         end do
         do n = 1,self%NNB
            j = self%nj(n)
            do i = self%nfi(n)-HALO,self%nli(n)+HALO
               if (grid%mask(i,j) > 1) f(i,j+1,:) = f(i,j,:)
            end do
         end do
         do n = 1,self%NEB
            i = self%ei(n)
            do j = self%efj(n)-HALO,self%elj(n)+HALO
               if (grid%mask(i,j) > 1) f(i+1,j,:) = f(i,j,:)
            end do
         end do
         do n = 1,self%NSB
            j = self%sj(n)
            do i = self%sfi(n)-HALO,self%sli(n)+HALO
               if (grid%mask(i,j) > 1) f(i,j-1,:) = f(i,j,:)
            end do
         end do
      case (2) ! UGRID
         do n = 1,self%NNB
            j = self%nj(n)
            do i = self%nfi(n)-HALO,self%nli(n)+HALO
               if (grid%mask(i,j) == 3) f(i,j,:) = f(i,j-1,:)
             end do
         end do
         do n = 1,self%NSB
            j = self%sj(n)
            do i = self%sfi(n)-HALO,self%sli(n)+HALO
               if (grid%mask(i,j) == 3) f(i,j,:) = f(i,j+1,:)
            end do
         end do
      case (3) ! VGRID
         do n = 1,self%NWB
            i = self%wi(n)
            do j = self%wfj(n)-HALO,self%wlj(n)+HALO
               if (grid%mask(i,j) == 3) f(i,j,:) = f(i+1,j,:)
            end do
         end do
         do n = 1,self%NEB
            i = self%ei(n)
            do j = self%efj(n)-HALO,self%elj(n)+HALO
               if (grid%mask(i,j) == 3) f(i,j,:) = f(i-1,j,:)
            end do
         end do
   end select
END SUBROUTINE mirror_bdy_3d

!---------------------------------------------------------------------------

END SUBMODULE mirror_bdys_smod
