! Copyright (C) 2020 Bolding & Bruggeman

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
   call self%logs%info('init_vertical()',level=2)

   where (self%T%mask > 0)
      self%T%ssen=self%T%z
   end where

!KB   call init_sigma(self)
   call init_gvc(self)

   return
END SUBROUTINE init_vertical

!---------------------------------------------------------------------------

module SUBROUTINE do_vertical(self)
   !! A wrapper for vertical coordinate calculations

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_domain), intent(inout) :: self

!  Local constants

!  Local variables
real(real64) :: dt
!-----------------------------------------------------------------------------
   call self%logs%info('do_vertical()',level=2)

!KB   call do_sigma(self)
   call do_gvc(self,dt)

   return
END SUBROUTINE do_vertical

!---------------------------------------------------------------------------

END SUBMODULE vertical_coordinates_smod
