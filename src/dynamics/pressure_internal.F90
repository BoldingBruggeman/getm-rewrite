! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

SUBMODULE (getm_pressure) pressure_internal_smod

INTERFACE
   module subroutine blumberg_mellor(self,buoy)
      class(type_getm_pressure), intent(inout) :: self
#define _T3_ self%domain%T%l(1):,self%domain%T%l(2):,self%domain%T%l(3):
      real(real64), intent(in) :: buoy(_T3_)
#undef _T3_
   end subroutine blumberg_mellor
END INTERFACE

ENUM, BIND(C)
   ENUMERATOR :: method_blumberg_mellor=1
!   ENUMERATOR :: use_parabolic=1
!   ENUMERATOR :: use_gotm=2
END ENUM


!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

module SUBROUTINE pressure_internal(self,buoy)
   !! Internal/baroclinic pressure gradients

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_pressure), intent(inout) :: self
#define _T3_ self%domain%T%l(1):,self%domain%T%l(2):,self%domain%T%l(3):
   real(real64), intent(in) :: buoy(_T3_)
#undef _T3_

!  Local constants

!  Local variables
   integer :: i,j
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('pressure_internal()',level=2)

   select case (self%ip_method)
      case(method_blumberg_mellor)
         call blumberg_mellor(self,buoy)
   end select

#if 0
   self%SxB(i,j)=-SUM(self%idpdx(i,j,1:))
   self%SyB(i,j)=-SUM(self%idpdy(i,j,1:))
#endif
END SUBROUTINE pressure_internal

!---------------------------------------------------------------------------

END SUBMODULE pressure_internal_smod
