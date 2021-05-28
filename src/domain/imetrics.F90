! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

SUBMODULE (getm_domain) imetrics_smod

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

module SUBROUTINE imetrics(self)
   !! Allocate all domain related variables

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_domain), intent(inout) :: self

!  Local constants

!  Local variables
!-----------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('metrics()',level=2)
   self%T%iarea = 1._real64/self%T%area
   self%T%idx   = 1._real64/self%T%dx
   self%T%idy   = 1._real64/self%T%dy
   self%U%iarea = 1._real64/self%U%area
   self%U%idx   = 1._real64/self%U%dx
   self%U%idy   = 1._real64/self%U%dy
   self%V%iarea = 1._real64/self%V%area
   self%V%idx   = 1._real64/self%V%dx
   self%V%idy   = 1._real64/self%V%dy
   self%X%iarea = 1._real64/self%X%area
   self%X%idx   = 1._real64/self%X%dx
   self%X%idy   = 1._real64/self%X%dy
END SUBROUTINE imetrics

!---------------------------------------------------------------------------

END SUBMODULE imetrics_smod
