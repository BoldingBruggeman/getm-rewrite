! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

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
   real(real64) :: gammai= 10., Dmin=0.5_real64
!---------------------------------------------------------------------------
   call self%logs%info('pressure_external()',level=2)
#define TG self%domain%T
#define UG self%domain%U
   do j=UG%l(2),UG%u(2)
      do i=UG%l(1),UG%u(1)
         if (UG%mask(i,j) == 1 .or. UG%mask(i,j) == 2) then
            zp = max( z(i+1,j) , -TG%H(i  ,j)+min( Dmin , TG%D(i+1,j) ) )
            zm = max( z(i  ,j) , -TG%H(i+1,j)+min( Dmin , TG%D(i  ,j) ) )
            self%dpdx(i,j) = ( zp - zm + (sp(i+1,j)-sp(i,j))*gammai ) / UG%dx(i,j)
         end if
      end do
   end do
#undef UG
#define VG self%domain%V
   do j=VG%l(2),VG%u(2)
      do i=VG%l(1),VG%u(1)
         if (VG%mask(i,j) == 1 .or. VG%mask(i,j) == 2) then
            zp = max( z(i,j+1) , -TG%H(i  ,j)+min( Dmin , TG%D(i,j+1) ) )
            zm = max( z(i,j  ) , -TG%H(i,j+1)+min( Dmin , TG%D(i,j  ) ) )
            self%dpdy(i,j) = ( zp - zm + (sp(i,j+1)-sp(i,j))*gammai ) / VG%dy(i,j)
         end if
      end do
   end do
#undef VG
#undef TG
   return
END SUBROUTINE pressure_surface

!---------------------------------------------------------------------------

END SUBMODULE pressure_surface_smod
