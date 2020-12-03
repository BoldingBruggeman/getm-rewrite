! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

SUBMODULE (getm_pressure) pressure_internal_smod

INTERFACE
   module subroutine blumberg_mellor(self,buoy)
      class(type_getm_pressure), intent(inout) :: self
#define _T3_ self%domain%T%l(1):,self%domain%T%l(2):,self%domain%T%l(3):
      real(real64), intent(in) :: buoy(_T3_)
#undef _T3_
   end subroutine blumberg_mellor
   module subroutine init_shchepetkin_mcwilliams(self)
      class(type_getm_pressure), intent(inout) :: self
   end subroutine init_shchepetkin_mcwilliams
   module subroutine shchepetkin_mcwilliams(self,buoy)
      class(type_getm_pressure), intent(inout) :: self
#define _T3_ self%domain%T%l(1):,self%domain%T%l(2):,self%domain%T%l(3):
      real(real64), intent(in) :: buoy(_T3_)
#undef _T3_
   end subroutine shchepetkin_mcwilliams
END INTERFACE

ENUM, BIND(C)
   ENUMERATOR :: method_blumberg_mellor=1
   ENUMERATOR :: method_shchepetkin_mcwilliams=2
!   ENUMERATOR :: use_parabolic=1
!   ENUMERATOR :: use_gotm=2
END ENUM


!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

module SUBROUTINE pressure_internal_initialize(self)
   !! Internal/baroclinic pressure gradients

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_pressure), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: stat
   type (type_field), pointer :: f
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('pressure_internal_initialize()',level=3)

   TGrid: associate( TG => self%domain%T )
#ifndef _STATIC_
   call mm_s('idpdx',self%idpdx,TG%l(1:3),TG%u(1:3),def=0._real64,stat=stat)
   call mm_s('idpdy',self%idpdy,TG%l(1:3),TG%u(1:3),def=0._real64,stat=stat)
#endif
   call self%fm%register('idpdx', 'Pa/m', 'internal pressure gradient - x', &
                         standard_name='', &
                         dimensions=(self%domain%T%dim_3d_ids), &
    !KB                        output_level=output_level_debug, &
                         part_of_state=.false., &
                         category='baroclinicity', field=f)
   call self%fm%send_data('idpdx', self%idpdx(TG%imin:TG%imax,TG%jmin:TG%jmax,TG%kmin:TG%kmax))
   call self%fm%register('idpdy', 'Pa/m', 'internal pressure gradient - y', &
                         standard_name='', &
                         dimensions=(self%domain%T%dim_3d_ids), &
    !KB                        output_level=output_level_debug, &
                         part_of_state=.false., &
                         category='baroclinicity', field=f)
   call self%fm%send_data('idpdy', self%idpdy(TG%imin:TG%imax,TG%jmin:TG%jmax,TG%kmin:TG%kmax))
   end associate TGrid

   select case (self%method_internal_pressure)
      case(method_blumberg_mellor)
      case(method_shchepetkin_mcwilliams)
         call init_shchepetkin_mcwilliams(self)
   end select
END SUBROUTINE pressure_internal_initialize

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

   select case (self%method_internal_pressure)
      case(method_blumberg_mellor)
         call blumberg_mellor(self,buoy)
      case(method_shchepetkin_mcwilliams)
         call shchepetkin_mcwilliams(self,buoy)
   end select

#if 0
   self%SxB(i,j)=-SUM(self%idpdx(i,j,1:))
   self%SyB(i,j)=-SUM(self%idpdy(i,j,1:))
#endif
END SUBROUTINE pressure_internal

!---------------------------------------------------------------------------

END SUBMODULE pressure_internal_smod
