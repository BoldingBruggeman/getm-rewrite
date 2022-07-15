! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard
!! @note
!! check with GETM
!!
!! rru, rrv, Ui, Vi could maybe be moved to here
!!
!! @endnote

SUBMODULE (getm_momentum) slow_momentum_terms_smod

CONTAINS

!---------------------------------------------------------------------------

MODULE SUBROUTINE slow_momentum_terms(self,dt)
   !! Slow terms

   IMPLICIT NONE

   class(type_getm_momentum), intent(inout) :: self
   real(real64), intent(in) :: dt

!  Local constants

!  Local variables
!---------------------------------------------------------------------------
   if(associated(self%logs)) call self%logs%info('slow_momentum_terms()',level=2)
   ! [GETM Scientific Report: eqs. 2.21, 2.22]
!KB   call self%slow_advection(dt,.false.)

   ! [GETM Scientific Report: eqs. 2.18, 2.19]
!KB   if (self%apply_diffusion) call self%slow_diffusion()

   ! [GETM Scientific Report: eqs. 2.22, 2.23]
   if (self%apply_bottom_friction) call self%slow_bottom_friction()

   self%Uio=self%Ui; self%Ui=0._real64
   self%Vio=self%Vi; self%Vi=0._real64
END SUBROUTINE slow_momentum_terms

!---------------------------------------------------------------------------

END SUBMODULE slow_momentum_terms_smod
