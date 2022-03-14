! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!> In this module, the temperature equation is processed by
!> reading in the namelist {\tt temp} and initialising the temperature field
!> (this is done in {\tt init\_temperature}),
!> and calculating the advection-diffusion-equation, which includes
!> penetrating short-wave radiation as source term (see {\tt do\_temperature}).

MODULE getm_pressure

   !! Description:
   !!   < Say what this module contains >
   !!
   !! Current Code Owner: < Name of person responsible for this code >
   !!
   !! Code Description:
   !!   Language: Fortran 90.
   !!   This code is written to JULES coding standards v1.

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use memory_manager
   use logging
   use field_manager
   use getm_domain

   IMPLICIT NONE

   PRIVATE  ! Private scope by default

!  Module constants

!  Module types and variables
   type, public :: type_getm_pressure

      class(type_logging), pointer :: logs => null()
      class(type_field_manager), pointer :: fm => null()
      class(type_getm_domain), pointer :: domain

#ifdef _STATIC_
      real(real64), dimension(E2DFIELD) :: dpdx, dpdy
      real(real64), dimension(I3DFIELD) :: idpdx, idpdy
#else
      real(real64), dimension(:,:), allocatable :: dpdx, dpdy
      real(real64), dimension(:,:,:), allocatable :: idpdx, idpdy
      ! shchepetkin_mcwilliams
      real(real64), allocatable :: dR(:,:,:)
      real(real64), allocatable :: dZ(:,:,:)
      real(real64), allocatable :: P(:,:,:)
      real(real64), allocatable :: FC(:,:)
      real(real64), allocatable :: dZx(:,:)
      real(real64), allocatable :: dRx(:,:)
      real(real64) :: P_time,U_time,V_time
#endif
      integer :: method_internal_pressure=1

      contains

      procedure :: configure => pressure_configure
      procedure :: initialize => pressure_initialize
      procedure :: surface => pressure_surface
      procedure :: internal => pressure_internal
      procedure :: internal_initialize => pressure_internal_initialize

   end type type_getm_pressure

   INTERFACE
      module subroutine pressure_surface(self,z,sp)
         class(type_getm_pressure), intent(inout) :: self
#define _T2_ self%domain%T%l(1):,self%domain%T%l(2):
         real(real64), intent(in) :: z(_T2_)
           !! elevation [m]
         real(real64), intent(in) :: sp(_T2_)
           !! pressure [Pa]
#undef _T2_
      end subroutine pressure_surface
      module subroutine pressure_internal_initialize(self,runtype)
         class(type_getm_pressure), intent(inout) :: self
         integer, intent(in) :: runtype
      end subroutine pressure_internal_initialize
      module subroutine pressure_internal(self,buoy,SxB,SyB)
         class(type_getm_pressure), intent(inout) :: self
#define _T3_ self%domain%T%l(1):,self%domain%T%l(2):,self%domain%T%l(3):
         real(real64), intent(in) :: buoy(_T3_)
#undef _T3_
#define _T2_ self%domain%T%l(1):,self%domain%T%l(2):
         real(real64), intent(inout) :: SxB(_T2_)
         real(real64), intent(inout) :: SyB(_T2_)
#undef _T2_
      end subroutine pressure_internal
   END INTERFACE

!---------------------------------------------------------------------------

CONTAINS

!---------------------------------------------------------------------------

SUBROUTINE pressure_configure(self,logs,fm)

   !! Configure the components belonging to the dynamics

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_pressure), intent(inout) :: self
   class(type_logging), intent(in), target, optional :: logs
   TYPE(type_field_manager), intent(inout), target, optional :: fm

!  Local constants

!  Local variables
!-----------------------------------------------------------------------
   if (present(logs)) then
      self%logs => logs
      call self%logs%info('pressure_configuration()',level=2)
   end if
   if (present(fm)) then
      self%fm => fm
   end if
END SUBROUTINE pressure_configure

!---------------------------------------------------------------------------

SUBROUTINE pressure_initialize(self,runtype,domain)
   !! Initialize all pressure components

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_pressure), intent(inout) :: self
   integer, intent(in) :: runtype
   class(type_getm_domain), intent(in), target :: domain

!  Local constants

!  Local variables
   integer :: stat
   type (type_field), pointer :: f
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('pressure_initialize()',level=2)
   self%domain => domain
   UGrid: associate( UG => self%domain%U )
   VGrid: associate( VG => self%domain%V )
#ifndef _STATIC_
   call mm_s('dpdx',self%dpdx,UG%l(1:2),UG%u(1:2),def=0._real64,stat=stat)
   call mm_s('dpdy',self%dpdy,VG%l(1:2),VG%u(1:2),def=0._real64,stat=stat)
#endif
   if (associated(self%fm)) then
      call self%fm%register('dpdx', 'Pa/m', 'surface pressure gradient - x', &
                            standard_name='', &
                            dimensions=(self%domain%T%dim_2d_ids), &
    !KB                        output_level=output_level_debug, &
                            part_of_state=.false., &
                            category='airsea', field=f)
      call self%fm%send_data('dpdx', self%dpdx(UG%imin:UG%imax,UG%jmin:UG%jmax))
      call self%fm%register('dpdy', 'Pa/m', 'surface pressure gradient - y', &
                            standard_name='', &
                            dimensions=(self%domain%T%dim_2d_ids), &
    !KB                        output_level=output_level_debug, &
                            part_of_state=.false., &
                            category='airsea', field=f)
      call self%fm%send_data('dpdy', self%dpdy(VG%imin:VG%imax,VG%jmin:VG%jmax))

   end if
   call self%internal_initialize(runtype)
   end associate VGrid
   end associate UGrid
END SUBROUTINE pressure_initialize

!---------------------------------------------------------------------------

END MODULE getm_pressure
