! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!> In this module, the temperature equation is processed by
!> reading in the namelist {\tt temp} and initialising the temperature field
!> (this is done in {\tt init\_temperature}),
!> and calculating the advection-diffusion-equation, which includes
!> penetrating short-wave radiation as source term (see {\tt do\_temperature}).

#ifdef _STATIC_
#include "dimensions.h"
#endif

MODULE getm_physics

   !! Description:
   !!   < Say what this module contains >
   !!
   !! Current Code Owner: < Name of person responsible for this code >
   !!
   !! Code Description:
   !!   Language: Fortran 90.
   !!   This code is written to JULES coding standards v1.

   use iso_fortran_env
   use logging
   use field_manager
   use getm_domain, only: type_getm_grid
   use getm_operators
   use getm_radiation, only: type_radiation
   use getm_salinity, only: type_salinity
   use getm_temperature, only: type_temperature
   use getm_density, only: type_density
!KB   use getm_turbulence, only: type_turbulence

   IMPLICIT NONE

   PRIVATE  ! Private scope by default

!  Module constants
   integer, parameter :: rk = kind(1.d0)

!  Module types and variables
   type, public :: type_getm_physics
      class(type_logging), pointer :: logs
      class(type_field_manager), pointer :: fm

      TYPE(type_radiation), public :: radiation
      TYPE(type_salinity), public :: salinity
      TYPE(type_temperature), public :: temperature
      TYPE(type_density), public :: density
!KB      TYPE(type_turbulence), public :: turbulence

      contains

      procedure :: configure => physics_configure
      procedure :: initialize => physics_initialize
!KB      procedure :: calculate => physics_calculate

   end type type_getm_physics

!---------------------------------------------------------------------------

CONTAINS

!---------------------------------------------------------------------------

SUBROUTINE physics_configure(self,logs,fm)

   !! Configure the components belonging to the physics

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_physics), intent(inout) :: self
   class(type_logging), intent(in), target :: logs
   class(type_field_manager), intent(in), target :: fm

!  Local constants

!  Local variables
!KB   integer :: rc
! ---------------------------------------------------------------------------
   self%logs => logs
   call self%logs%info('physics_configure()',level=1)
!KB   write(*,*) 'This file was compiled by ', compiler_version(), ' using the options ', compiler_options()

   call self%salinity%configuration(logs,fm)
   call self%temperature%configuration(logs,fm)
   call self%density%configure(logs,fm)
!KB   call self%turbulence%configuration(logs)

   call logs%info('done',level=1)
   return
END SUBROUTINE physics_configure

!---------------------------------------------------------------------------

SUBROUTINE physics_initialize(self,grid,advection,vertical_diffusion)

   !! Initialize all dynamical components

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_physics), intent(inout) :: self
   class(type_getm_grid), intent(in), target :: grid
   class(type_advection), intent(in), target :: advection
   class(type_vertical_diffusion), intent(in), target :: vertical_diffusion

!  Local constants

!  Local variables
!KB   integer :: imin,imax,jmin,jmax,kmin,kmax
!KB   integer :: rc
!---------------------------------------------------------------------------
   call self%logs%info('physics_initialize()',level=1)
!   call logs%warn(trim(self%temperature%config%f))

!   imin = grid%imin; imax = grid%imax
!   jmin = grid%jmin; jmax = grid%jmax
!   kmin = grid%kmin; kmax = grid%kmax

!KB   call self%radiation%initialize(logs,grid)
   call self%salinity%initialize(grid,advection,vertical_diffusion)
   call self%temperature%initialize(grid)
   call self%density%initialize(grid)
   call self%density%density(self%salinity%S,self%temperature%T)
   call self%density%buoyancy()
   call self%logs%info('done',level=1)
   return
END SUBROUTINE physics_initialize

!---------------------------------------------------------------------------

SUBROUTINE physics_do_2d(self,logs)

   !! Advection/diffusion of the temperature field

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_physics), intent(inout) :: self
   class(type_logging), intent(in) :: logs


!  Local constants

!  Local variables
   integer :: rc

!---------------------------------------------------------------------------
   call logs%info('physics_do_2d()',level=1)

   return
END SUBROUTINE physics_do_2d

!---------------------------------------------------------------------------

SUBROUTINE physics_do_3d(self,logs)

   !! Advection/diffusion of the temperature field

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_physics), intent(out) :: self
   class(type_logging), intent(in) :: logs

!  Local constants

!  Local variables
   integer :: rc
!---------------------------------------------------------------------------
   call logs%info('physics_do_3d()',level=1)

   return
END SUBROUTINE physics_do_3d

!---------------------------------------------------------------------------

END MODULE getm_physics

