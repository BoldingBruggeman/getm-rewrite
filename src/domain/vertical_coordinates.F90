! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

  !! @note
  !! zin and zio for U and V points
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
   module subroutine init_adaptive(self)
      class(type_getm_domain), intent(inout) :: self
   end subroutine init_adaptive
   module subroutine do_adaptive(self,dt)
      class(type_getm_domain), intent(inout) :: self
      real(real64), intent(in) :: dt
   end subroutine do_adaptive
END INTERFACE

ENUM, BIND(C)
   ENUMERATOR :: method_sigma=1
   ENUMERATOR :: method_gvc=2
   ENUMERATOR :: method_adaptive=3
END ENUM

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
      self%T%zin=self%T%z
   end where

   select case (self%method_vertical_coordinates)
      case(method_sigma)
         call init_sigma(self)
!KB         call depths(self,self%T); call depths(self,self%U); call depths(self,self%V)
      case(method_gvc)
         call init_gvc(self)
!KB         call depths(self,self%T); call depths(self,self%U); call depths(self,self%V)
   end select
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

   select case (self%method_vertical_coordinates)
      case(method_sigma)
         call do_sigma(self)
!KB         call depths(self,self%T); call depths(self,self%U); call depths(self,self%V)
      case(method_gvc)
         call do_gvc(self,dt)
!KB         call depths(self,self%T); call depths(self,self%U); call depths(self,self%V)
   end select
END SUBROUTINE do_vertical

!---------------------------------------------------------------------------

subroutine depths(domain,grid)
   !! Cacluate the depth to the cell face and center

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_domain), intent(inout) :: domain
   class(type_getm_grid), intent(inout) :: grid

!  Local constants

!  Local variables
   integer :: i,j,k
!-----------------------------------------------------------------------------
   if (associated(grid%logs)) call grid%logs%info('depths()',level=3)
   do j=grid%jmin,grid%jmax+1
      do i=grid%imin,grid%imax+1
         if (grid%mask(i,j) > 0) then
            grid%zf(i,j,0)=-grid%H(i,j)
            grid%zf(i,j,1)=-grid%H(i,j)+grid%hn(i,j,1)
            grid%zc(i,j,1)=-grid%H(i,j)+0.5_real64*grid%hn(i,j,1)
            do k=2,grid%kmax
               grid%zf(i,j,k)=grid%zf(i,j,k-1)+grid%hn(i,j,k)
               grid%zc(i,j,k)=grid%zc(i,j,k-1)+0.5_real64*(grid%hn(i,j,k-1)+grid%hn(i,j,k))
            end do
         end if
      end do
   end do
   if (domain%nbdyp > 0) then
      call domain%mirror_bdys(grid,grid%zf)
      call domain%mirror_bdys(grid,grid%zc)
   end if
end subroutine depths

!---------------------------------------------------------------------------

END SUBMODULE vertical_coordinates_smod
