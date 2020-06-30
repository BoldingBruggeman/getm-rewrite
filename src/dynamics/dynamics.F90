! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!> In this module, the temperature equation is processed by
!> reading in the namelist {\tt temp} and initialising the temperature field
!> (this is done in {\tt init\_temperature}),
!> and calculating the advection-diffusion-equation, which includes
!> penetrating short-wave radiation as source term (see {\tt do\_temperature}).

#ifdef _STATIC_
#include "dimensions.h"
#endif

MODULE getm_dynamics

   !! Description:
   !!   < Say what this module contains >
   !!
   !! Current Code Owner: < Name of person responsible for this code >
   !!
   !! Code Description:
   !!   Language: Fortran 90.
   !!   This code is written to JULES coding standards v1.

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use logging
   use field_manager
   use getm_domain
   use getm_sealevel, only: type_getm_sealevel
   use getm_pressure, only: type_getm_pressure
   use getm_momentum, only: type_getm_momentum

   IMPLICIT NONE

   PRIVATE  ! Private scope by default

!  Module constants

!  Module types and variables
   type, public :: type_getm_dynamics

      TYPE(type_logging), public :: logs
      TYPE(type_getm_domain), public :: domain
      TYPE(type_getm_sealevel), public :: sealevel
      TYPE(type_getm_pressure), public :: pressure
      TYPE(type_getm_momentum), public :: momentum

      contains

      procedure :: configure => dynamics_configure
      procedure :: initialize => dynamics_initialize
!KB      procedure :: calculate => dynamics_calculate

   end type type_getm_dynamics

!---------------------------------------------------------------------------

CONTAINS

!---------------------------------------------------------------------------

SUBROUTINE dynamics_configure(self,logs,fm)

   !! Configure the components belonging to the dynamics

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_dynamics), intent(inout) :: self
   class(type_logging), intent(in), target :: logs
   class(type_field_manager), intent(inout), target :: fm

!  Local constants

!  Local variables
!---------------------------------------------------------------------------
!   self%logs => logs
!KB   call self%logs%info('dynamics_configure()',level=1)
   call self%sealevel%configuration(logs,fm)
   call self%pressure%configuration(logs,fm)
   call self%momentum%configuration(logs,fm)
   call self%logs%info('done',level=1)
   return
END SUBROUTINE dynamics_configure

!---------------------------------------------------------------------------

SUBROUTINE dynamics_initialize(self,domain)

   !! Initialize all dynamical components

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_dynamics), intent(inout) :: self
   class(type_getm_domain), intent(inout), target :: domain

!  Local constants

!  Local variables
   integer :: rc
!---------------------------------------------------------------------------
   call self%logs%info('dynamics_initialize()',level=1)
   call self%sealevel%initialize(domain)
   call self%pressure%initialize(domain)
   call self%momentum%initialize(domain)
   call self%logs%info('done',level=1)
   return
END SUBROUTINE dynamics_initialize

!---------------------------------------------------------------------------

SUBROUTINE dynamics_do_2d(logs)

   !! Advection/diffusion of the temperature field

   IMPLICIT NONE

!  Subroutine arguments
   class(type_logging), intent(in) :: logs

!  Local constants

!  Local variables
   integer :: rc
!---------------------------------------------------------------------------
   call logs%info('dynamics_do_2d()',level=1)

   return
END SUBROUTINE dynamics_do_2d

!---------------------------------------------------------------------------

SUBROUTINE dynamics_do_3d(logs)

   !! Advection/diffusion of the temperature field

   IMPLICIT NONE

!  Subroutine arguments
   class(type_logging), intent(in) :: logs

!  Local constants

!  Local variables
   integer :: rc
!---------------------------------------------------------------------------
   call logs%info('dynamics_do_3d()',level=1)

   return
END SUBROUTINE dynamics_do_3d

!---------------------------------------------------------------------------

END MODULE getm_dynamics
