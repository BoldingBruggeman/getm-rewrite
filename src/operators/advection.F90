! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

SUBMODULE (getm_operators) advection_smod

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

#define _NORMAL_ORDER_

module SUBROUTINE advection_initialize(self,var)

   !! Initialize the salinity field

   IMPLICIT NONE

   ! Subroutine arguments
   class(type_advection), intent(inout) :: self
   real(real64), dimension(:,:,:), intent(in) :: var
      !! grid dimensions in case of dynamic memory allocation

!  Local constants

!  Local variables
!---------------------------------------------------------------------------
   return
END SUBROUTINE advection_initialize

!---------------------------------------------------------------------------

module SUBROUTINE advection_calculate(self,mask,dz,dt,cnpar,avmol,nuh,var)

   !! Vertical diffusion

   IMPLICIT NONE

   ! Subroutine arguments
   class(type_advection), intent(inout) :: self
   integer, dimension(:,:), intent(in) :: mask
   real(real64), dimension(:,:,:), intent(in) :: dz
   real(real64), intent(in) :: dt
   real(real64), intent(in) :: cnpar
   real(real64), intent(in) :: avmol
   real(real64), dimension(:,:,:), intent(in) :: nuh
   real(real64), dimension(:,:,:), intent(inout) :: var

!  Local constants

!  Local variables
   integer :: i,j,k
!---------------------------------------------------------------------------

   if (self%kmax == 1) return

   return
END SUBROUTINE advection_calculate

!---------------------------------------------------------------------------

END SUBMODULE advection_smod
