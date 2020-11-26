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
#define _T2_ self%domain%T%l(1):,self%domain%T%l(2):
   real(real64), dimension(:,:), intent(in) :: z(_T2_)
   real(real64), dimension(:,:), intent(in) :: sp(_T2_)
#undef _T2_

!  Local constants

!  Local variables
   integer :: i,j
   real(real64) :: zp,zm
!KB
   real(real64) :: gammai= 1._real64/(9.81_real64*1025._real64)
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('pressure_external()',level=2)
   TGrid: associate( TG => self%domain%T )
   UGrid: associate( UG => self%domain%U )
   do j=UG%jmin,UG%jmax
      do i=UG%imin,UG%imax
         if (UG%mask(i,j) == 1 .or. UG%mask(i,j) == 2) then
            zp = max(z(i+1,j),-TG%H(i  ,j)+min(self%domain%Dmin,TG%D(i+1,j)))
            zm = max(z(i  ,j),-TG%H(i+1,j)+min(self%domain%Dmin,TG%D(i  ,j)))
            self%dpdx(i,j) = (zp-zm+(sp(i+1,j)-sp(i,j))*gammai)/UG%dx(i,j)
         end if
      end do
   end do
   end associate UGrid

   VGrid: associate( VG => self%domain%V )
   do j=VG%jmin,VG%jmax
      do i=VG%imin,VG%imax
         if (VG%mask(i,j) == 1 .or. VG%mask(i,j) == 2) then
            zp = max(z(i,j+1),-TG%H(i  ,j)+min(self%domain%Dmin,TG%D(i,j+1)))
            zm = max(z(i,j  ),-TG%H(i,j+1)+min(self%domain%Dmin,TG%D(i,j  )))
            self%dpdy(i,j) = (zp-zm+(sp(i,j+1)-sp(i,j))*gammai)/VG%dy(i,j)
         end if
      end do
   end do
   end associate VGrid
   end associate TGrid
   return
END SUBROUTINE pressure_surface

!---------------------------------------------------------------------------

END SUBMODULE pressure_surface_smod
