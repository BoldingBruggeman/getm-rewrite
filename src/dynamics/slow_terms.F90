! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard
!! @note
!! check with GETM
!!
!! rru, rrv, Ui, Vi could maybe be moved to here
!!
!! @endnote

SUBMODULE (getm_momentum) slow_terms_smod

CONTAINS

!---------------------------------------------------------------------------

MODULE SUBROUTINE slow_terms(self,dt,idpdx,idpdy)
   !! Slow terms

   IMPLICIT NONE

   class(type_getm_momentum), intent(inout) :: self
   real(real64), intent(in) :: dt
#define _T3_ self%domain%T%l(1):,self%domain%T%l(2):,self%domain%T%l(3):
   real(real64), intent(in), optional :: idpdx(_T3_)
   real(real64), intent(in), optional :: idpdy(_T3_)
#undef _T3_

!  Local constants

!  Local variables
   integer :: i,j,k
!---------------------------------------------------------------------------
   if(associated(self%logs)) call self%logs%info('slow_terms()',level=2)
   call self%slow_advection(dt)
   call self%slow_diffusion()
   TGrid: associate( TG => self%domain%T )
   UGrid: associate( UG => self%domain%U )
   VGrid: associate( VG => self%domain%V )
   if (present(idpdx) .and. present(idpdy)) then
      if(associated(self%logs)) call self%logs%info('slow_buoyancy()',level=3)
      do j=TG%jmin,TG%jmax
         do i=TG%imin,TG%imax
            if (UG%mask(i,j) > 0) then
               self%SxB(i,j)=-SUM(idpdx(i,j,1:))
            end if
            if (VG%mask(i,j) > 0) then
               self%SyB(i,j)=-SUM(idpdy(i,j,1:))
            end if
         end do
      end do
   end if
   end associate VGrid
   end associate UGrid
   end associate TGrid
   call self%slow_bottom_friction()
   self%Uio=self%Ui; self%Ui=0._real64
   self%Vio=self%Vi; self%Vi=0._real64
END SUBROUTINE slow_terms

!---------------------------------------------------------------------------

END SUBMODULE slow_terms_smod
