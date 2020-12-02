! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

  !! @note
  !! ssen and sseo for U and V points
  !! @note

SUBMODULE (getm_domain) vertical_coordinates_smod

INTERFACE
   module subroutine init_sigma(self)
      class(type_getm_domain), intent(inout) :: self
   end subroutine init_sigma
   module subroutine do_sigma(self)
      class(type_getm_domain), intent(inout) :: self
   end subroutine do_sigma
   module subroutine init_gvc(self)
      class(type_getm_domain), intent(inout) :: self
   end subroutine init_gvc
   module subroutine do_gvc(self,dt)
      class(type_getm_domain), intent(inout) :: self
      real(real64), intent(in) :: dt
   end subroutine do_gvc
END INTERFACE

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

module SUBROUTINE init_vertical(self)
   !! A wrapper for vertical coordinate calculations

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_domain), intent(inout) :: self

!  Local constants

!  Local variables
!-----------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('init_vertical()',level=2)

   where (self%T%mask > 0)
      self%T%ssen=self%T%z
   end where

   call init_sigma(self)
!KB   call init_gvc(self)

   return
END SUBROUTINE init_vertical

!---------------------------------------------------------------------------

module SUBROUTINE do_vertical(self,dt)
   !! A wrapper for vertical coordinate calculations

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_domain), intent(inout) :: self
   real(real64), intent(in) :: dt

!  Local constants

!  Local variables
   integer :: i,j,k
!-----------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('do_vertical()',level=2)

   call do_sigma(self)
!KB   call do_gvc(self,dt)

   TGrid: associate( TG => self%T )
   do j=TG%jmin,TG%jmax+1
      do i=TG%imin,TG%imax+1
         if (TG%mask(i,j) > 1) then
            TG%zf(i,j,0)=-TG%H(i,j)
            TG%zf(i,j,1)=-TG%H(i,j)+TG%hn(i,j,1)
            TG%zc(i,j,1)=-TG%H(i,j)+0.5_real64*TG%hn(i,j,1)
            do k=2,TG%kmax
               TG%zf(i,j,1)=TG%zf(i,j,k-1)+TG%hn(i,j,k)
               TG%zc(i,j,k)=TG%zc(i,j,k-1)+0.5_real64*TG%hn(i,j,k)
            end do
         end if
      end do
   end do
   end associate TGrid

   UGrid: associate( UG => self%U )
   do j=UG%jmin,UG%jmax+1
      do i=UG%imin,UG%imax+1
         if (UG%mask(i,j) > 1) then
            UG%zf(i,j,0)=-UG%H(i,j)
            UG%zf(i,j,1)=-UG%H(i,j)+UG%hn(i,j,1)
            UG%zc(i,j,1)=-UG%H(i,j)+0.5_real64*UG%hn(i,j,1)
            do k=2,UG%kmax
               UG%zf(i,j,1)=UG%zf(i,j,k-1)+UG%hn(i,j,k)
               UG%zc(i,j,k)=UG%zc(i,j,k-1)+0.5_real64*UG%hn(i,j,k)
            end do
         end if
      end do
   end do
   end associate UGrid

   VGrid: associate( VG => self%V )
   do j=VG%jmin,VG%jmax+1
      do i=VG%imin,VG%imax+1
         if (VG%mask(i,j) > 1) then
            VG%zf(i,j,0)=-VG%H(i,j)
            VG%zf(i,j,1)=-VG%H(i,j)+VG%hn(i,j,1)
            VG%zc(i,j,1)=-VG%H(i,j)+0.5_real64*VG%hn(i,j,1)
            do k=2,VG%kmax
               VG%zf(i,j,1)=VG%zf(i,j,k-1)+VG%hn(i,j,k)
               VG%zc(i,j,k)=VG%zc(i,j,k-1)+0.5_real64*VG%hn(i,j,k)
            end do
         end if
      end do
   end do
   end associate VGrid
END SUBROUTINE do_vertical

!---------------------------------------------------------------------------

END SUBMODULE vertical_coordinates_smod
