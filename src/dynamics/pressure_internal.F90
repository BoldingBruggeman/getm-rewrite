! Copyright (C) 2020 Bolding & Bruggeman

SUBMODULE (getm_pressure) pressure_internal_smod

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

module SUBROUTINE pressure_internal(self)
   !! Internal/baroclinic pressure gradients

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_pressure), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: i,j
!---------------------------------------------------------------------------
   call self%logs%info('pressure_internal()',level=2)

   return
END SUBROUTINE pressure_internal

!---------------------------------------------------------------------------

END SUBMODULE pressure_internal_smod
