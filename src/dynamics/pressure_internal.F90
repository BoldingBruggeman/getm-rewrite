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

module SUBROUTINE pressure_internal_initialize(self,runtype)
   !! Internal/baroclinic pressure gradients

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_pressure), intent(inout) :: self
   integer, intent(in) :: runtype

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
   if (runtype > 2) then
      if (associated(self%fm)) then
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
      end if

      select case (self%method_internal_pressure)
         case(method_blumberg_mellor)
         case(method_shchepetkin_mcwilliams)
            call init_shchepetkin_mcwilliams(self)
      end select
   end if
   end associate TGrid
END SUBROUTINE pressure_internal_initialize

!-----------------------------------------------------------------------------

module SUBROUTINE pressure_internal(self,buoy,SxB,SyB)
   !! Internal/baroclinic pressure gradients

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_pressure), intent(inout) :: self
#define _T3_ self%domain%T%l(1):,self%domain%T%l(2):,self%domain%T%l(3):
   real(real64), intent(in) :: buoy(_T3_)
#undef _T3_
#define _T2_ self%domain%T%l(1):,self%domain%T%l(2):
   real(real64), intent(inout) :: SxB(_T2_)
   real(real64), intent(inout) :: SyB(_T2_)
#undef _T2_

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
   TGrid: associate( TG => self%domain%T )
   UGrid: associate( UG => self%domain%U )
   VGrid: associate( VG => self%domain%V )
   if(associated(self%logs)) call self%logs%info('slow_buoyancy()',level=3)
   do j=TG%jmin,TG%jmax
      do i=TG%imin,TG%imax
         if (UG%mask(i,j) > 0) then
            ! [GETM Scientific Report: eqs. 2.24]
            SxB(i,j)=-SUM(self%idpdx(i,j,1:))
         end if
         if (VG%mask(i,j) > 0) then
            ! [GETM Scientific Report: eqs. 2.25]
            SyB(i,j)=-SUM(self%idpdy(i,j,1:))
         end if
      end do
   end do
   end associate VGrid
   end associate UGrid
   end associate TGrid
END SUBROUTINE pressure_internal

!---------------------------------------------------------------------------

END SUBMODULE pressure_internal_smod
