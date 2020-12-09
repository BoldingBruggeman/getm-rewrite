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

ENUM, BIND(C)
   ENUMERATOR :: method_sigma=1
   ENUMERATOR :: method_gvc=2
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
      self%T%ssen=self%T%z
   end where

   select case (self%method_vertical_coordinates)
      case(method_sigma)
         call init_sigma(self)
         call depths(self%T); call depths(self%U); call depths(self%V)
      case(method_gvc)
         call init_gvc(self)
         call depths(self%T); call depths(self%U); call depths(self%V)
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
         call depths(self%T); call depths(self%U); call depths(self%V)
      case(method_gvc)
         call do_gvc(self,dt)
         call depths(self%T); call depths(self%U); call depths(self%V)
   end select
END SUBROUTINE do_vertical

!---------------------------------------------------------------------------

subroutine depths(self)
   !! Cacluate the depth to the cell face and center

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_grid), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: i,j,k
!-----------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('depths()',level=3)
   do j=self%jmin,self%jmax+1
      do i=self%imin,self%imax+1
         if (self%mask(i,j) > 0) then
            self%zf(i,j,0)=-self%H(i,j)
            self%zf(i,j,1)=-self%H(i,j)+self%hn(i,j,1)
            self%zc(i,j,1)=-self%H(i,j)+0.5_real64*self%hn(i,j,1)
            do k=2,self%kmax
               self%zf(i,j,k)=self%zf(i,j,k-1)+self%hn(i,j,k)
               self%zc(i,j,k)=self%zc(i,j,k-1)+0.5_real64*(self%hn(i,j,k-1)+self%hn(i,j,k))
            end do
         end if
      end do
   end do
end subroutine depths

!---------------------------------------------------------------------------

END SUBMODULE vertical_coordinates_smod
