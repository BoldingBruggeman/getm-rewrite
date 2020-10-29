! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!> In this module, the temperature equation is processed by
!> reading in the namelist {\tt temp} and initialising the temperature field
!> (this is done in {\tt init\_temperature}),
!> and calculating the advection-diffusion-equation, which includes
!> penetrating short-wave radiation as source term (see {\tt do\_temperature}).

#ifdef _STATIC_
#include "dimensions.h"
#endif

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
#endif

      contains

      procedure :: configuration => pressure_configuration
      procedure :: initialize => pressure_initialize
      procedure :: surface => pressure_surface
      procedure :: internal => pressure_internal
!KB      procedure :: calculate => physics_calculate

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
      module subroutine pressure_internal(self)
         class(type_getm_pressure), intent(inout) :: self
      end subroutine pressure_internal
   END INTERFACE

!---------------------------------------------------------------------------

CONTAINS

!---------------------------------------------------------------------------

SUBROUTINE pressure_configuration(self,logs,fm)

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
END SUBROUTINE pressure_configuration

!---------------------------------------------------------------------------

SUBROUTINE pressure_initialize(self,domain)

   !! Initialize all dynamical components

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_pressure), intent(inout) :: self
   class(type_getm_domain), intent(in), target :: domain

!  Local constants

!  Local variables
   integer :: imin,imax,jmin,jmax,kmax
   integer :: stat
   type (type_field), pointer :: f
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('pressure_initialize()',level=2)
   self%domain => domain
   TGrid: associate( TG => self%domain%T )
#ifndef _STATIC_
   call mm_s('dpdx',self%dpdx,TG%l(1:2),TG%u(1:2),def=0._real64,stat=stat)
   call mm_s('dpdy',self%dpdy,TG%l(1:2),TG%u(1:2),def=0._real64,stat=stat)
   call mm_s('idpdx',self%idpdx,TG%l(1:3),TG%u(1:3),def=0._real64,stat=stat)
   call mm_s('idpdy',self%idpdy,TG%l(1:3),TG%u(1:3),def=0._real64,stat=stat)
#endif
   if (associated(self%fm)) then
      call self%fm%register('dpdx', 'Pa/m', 'surface pressure gradient - x', &
                            standard_name='', &
                            dimensions=(self%domain%T%dim_2d_ids), &
    !KB                        output_level=output_level_debug, &
                            part_of_state=.false., &
                            category='airsea', field=f)
      call self%fm%send_data('dpdx', self%dpdx(TG%imin:TG%imax,TG%jmin:TG%jmax))
      call self%fm%register('dpdy', 'Pa/m', 'surface pressure gradient - y', &
                            standard_name='', &
                            dimensions=(self%domain%T%dim_2d_ids), &
    !KB                        output_level=output_level_debug, &
                            part_of_state=.false., &
                            category='airsea', field=f)
      call self%fm%send_data('dpdy', self%dpdy(TG%imin:TG%imax,TG%jmin:TG%jmax))
   end if
   end associate TGrid
END SUBROUTINE pressure_initialize

!---------------------------------------------------------------------------

END MODULE getm_pressure
