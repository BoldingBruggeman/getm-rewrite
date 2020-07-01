! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

SUBMODULE (getm_momentum) shear_smod

CONTAINS

!---------------------------------------------------------------------------

module SUBROUTINE velocity_shear(self)

   !! Velocity shear

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: rc
!---------------------------------------------------------------------------
   call self%logs%info('velocity_shear()',level=2)

   return
END SUBROUTINE velocity_shear

!---------------------------------------------------------------------------

END SUBMODULE shear_smod
