! Copyright (C) 2020 Bolding & Bruggeman

SUBMODULE (getm_pressure) pressure_surface_smod

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

module SUBROUTINE pressure_surface(self,z,sp)
   !! External/barotropic pressure gradients

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_pressure), intent(inout) :: self
   real(real64), dimension(:,:), intent(in) :: z,sp

!  Local constants

!  Local variables
   integer :: i,j
   real(real64) :: zp,zm
   integer :: rc
!KB
   real(real64) :: gammai= 10., min_depth=0.5_real64
!---------------------------------------------------------------------------
   call self%logs%info('pressure_external()',level=2)
   do j=self%domain%U%l(2),self%domain%U%u(2)
      do i=self%domain%U%l(1),self%domain%U%u(1)
         if (self%domain%U%mask(i,j) == 1 .or. self%domain%U%mask(i,j) == 2) then
            zp = max( z(i+1,j) , -self%domain%T%H(i  ,j)+min( min_depth , self%domain%T%D(i+1,j) ) )
            zm = max( z(i  ,j) , -self%domain%T%H(i+1,j)+min( min_depth , self%domain%T%D(i  ,j) ) )
            self%dpdx(i,j) = ( zp - zm + (sp(i+1,j)-sp(i,j))*gammai ) / self%domain%U%dx(i,j)
         end if
      end do
   end do
   do j=self%domain%V%l(2),self%domain%V%u(2)
      do i=self%domain%V%l(1),self%domain%V%u(1)
         if (self%domain%V%mask(i,j) == 1 .or. self%domain%V%mask(i,j) == 2) then
            zp = max( z(i,j+1) , -self%domain%T%H(i  ,j)+min( min_depth , self%domain%T%D(i,j+1) ) )
            zm = max( z(i,j  ) , -self%domain%T%H(i,j+1)+min( min_depth , self%domain%T%D(i,j  ) ) )
            self%dpdy(i,j) = ( zp - zm + (sp(i,j+1)-sp(i,j))*gammai ) / self%domain%V%dy(i,j)
         end if
      end do
   end do
   return
END SUBROUTINE pressure_surface

!---------------------------------------------------------------------------

END SUBMODULE pressure_surface_smod
