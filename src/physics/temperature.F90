! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!> In this module, the temperature equation is processed by
!> reading in the namelist {\tt temp} and initialising the temperature field
!> (this is done in {\tt init\_temperature}),
!> and calculating the advection-diffusion-equation, which includes
!> penetrating short-wave radiation as source term (see {\tt do\_temperature}).

#ifdef _STATIC_
#include "dimensions.h"
#endif

MODULE getm_temperature

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
   integer, parameter :: tunit = 20

!  Module types and variables
   type, public :: type_temperature_configuration
      character(len=256) :: f = "temperature.nc"
   end type type_temperature_configuration

   type, public :: type_temperature
      !! author: Karsten Bolding
      !! version: v0.1
      !!
      !! Temperature type

      class(type_logging), pointer :: logs
      class(type_getm_grid), pointer :: G
      class(type_field_manager), pointer :: fm => null()
      TYPE(type_temperature_configuration) :: config

#ifdef _STATIC_
   real(real64), dimension(I3DFIELD), target :: T = 10._real64
#else
   real(real64), dimension(:,:,:), allocatable :: T
#endif

      contains

      procedure :: configuration => temperature_configuration
      procedure :: initialize => temperature_initialize
      procedure :: calculate => temperature_calculate
      procedure :: advection => temperature_advection
      procedure :: vertical_diffusion => temperature_vertical_diffusion

   end type type_temperature

CONTAINS

!---------------------------------------------------------------------------

SUBROUTINE temperature_configuration(self,logs,fm)

   !! Configure the the temperature module

   IMPLICIT NONE

!  Subroutine arguments
   class(type_temperature), intent(out) :: self
   class(type_logging), intent(in), target :: logs
   type(type_field_manager), optional, target  :: fm

!  Local constants

!  Local variables
   integer :: rc
!---------------------------------------------------------------------------
   self%logs => logs
   call self%logs%info('temperature_configuration()',level=2)
   call self%logs%info('reading initial temperature from: ',level=3,msg2=trim(self%config%f))
   if (present(fm)) then
      self%fm => fm
   end if
   return
END SUBROUTINE temperature_configuration

!---------------------------------------------------------------------------

SUBROUTINE temperature_initialize(self,grid)

   !! Initialize the temperature field

   IMPLICIT NONE

!  Subroutine arguments
   class(type_temperature), intent(inout) :: self
   class(type_getm_grid), intent(in), target :: grid

!  Local constants

!  Local variables
   integer :: i,j
   integer :: stat
!---------------------------------------------------------------------------
   call self%logs%info('temperature_initialize()',level=2)
   self%G => grid
!   self%T = null()
#ifndef _STATIC_
   call mm_s('T',self%T,self%G%l,self%G%u,def=15._real64,stat=stat)
#endif
   if (associated(self%fm)) then
      call self%fm%register('temp', 'Celsius', 'potential temperature', &
                            standard_name='sea_water_temperature', &
                            category='temperature_and_salinity', &
                            dimensions=(self%G%dim_3d_ids), &
                            fill_value=-9999._real64, &
                            part_of_state=.true.)
      call self%fm%send_data('temp', self%T(grid%imin:grid%imax,grid%jmin:grid%jmax,grid%kmin:grid%kmax))
   end if
   do j=self%G%jmin,self%G%jmax
      do i=self%G%imin,self%G%imax
         if (self%G%mask(i,j) ==  0) then
            self%T(i,j,:) = -9999._real64
         end if
      end do
   end do
   return
END SUBROUTINE temperature_initialize

!---------------------------------------------------------------------------

!> Here, one time step for the temperature equation is performed.
!> First, preparations for the call to the advection schemes are
!> made, i.e.\ calculating the necessary metric coefficients.
!> After the call to the advection schemes, which actually perform
!> the advection (and horizontal diffusion) step as an operational
!> split step, the solar radiation at the interfaces ({\tt rad(k)})
!> is calculated
!> from given surface radiation ({\tt swr\_loc})
!> by means of a double exponential
!> approach, see equation (\ref{Light}) on page \pageref{Light}).
!> An option to reflect part of the short wave radiation that reaches the
!> bottom has been implemented. In very shallow waters - or with very clear
!> waters - a significant part of the incoming radiation will reach the
!> bottom. Setting swr\_bot\_refl\_frac to a value between 0 and 1 will
!> reflect this fraction of what ever the value of SWR is at the bottom.
!> The default value of swr\_bot\_refl\_frac is 0.
!> The reflection is only done if the ratio between the surface and bottom
!> values of SWR is greater than swr\_min\_bot\_frac (default 0.01).
!> Furthermore, the surface heat flux {\tt sfl\_loc} is given a
!> value.
!> The sea surface temperature is limited by the freezing point
!> temperature (as a most primitive sea ice model). The next
!> step is to set up the tri-diagonal matrix for calculating the
!> new temperature by means of a semi-implicit central scheme for the
!> vertical diffusion. Source terms which appear on the right hand sides
!> are due to the divergence of the solar radiation at the interfaces.
!> The subroutine is completed by solving the tri-diagonal linear
!> equation by means of a tri-diagonal solver.


SUBROUTINE temperature_calculate(self)

   !! Advection/diffusion of the temperature field

   IMPLICIT NONE

!  Subroutine arguments
   class(type_temperature), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: rc
!---------------------------------------------------------------------------
   call self%logs%info('temperature_calculate()',level=2)

   return
END SUBROUTINE temperature_calculate

!---------------------------------------------------------------------------

!> Here, one time step for the temperature equation is performed.
!> First, preparations for the call to the advection schemes are
!> made, i.e.\ calculating the necessary metric coefficients.
!> After the call to the advection schemes, which actually perform
!> the advection (and horizontal diffusion) step as an operational
!> split step, the solar radiation at the interfaces ({\tt rad(k)})
!> is calculated
!> from given surface radiation ({\tt swr\_loc})
!> by means of a double exponential
!> approach, see equation (\ref{Light}) on page \pageref{Light}).
!> An option to reflect part of the short wave radiation that reaches the
!> bottom has been implemented. In very shallow waters - or with very clear
!> waters - a significant part of the incoming radiation will reach the
!> bottom. Setting swr\_bot\_refl\_frac to a value between 0 and 1 will
!> reflect this fraction of what ever the value of SWR is at the bottom.
!> The default value of swr\_bot\_refl\_frac is 0.
!> The reflection is only done if the ratio between the surface and bottom
!> values of SWR is greater than swr\_min\_bot\_frac (default 0.01).
!> Furthermore, the surface heat flux {\tt sfl\_loc} is given a
!> value.
!> The sea surface temperature is limited by the freezing point
!> temperature (as a most primitive sea ice model). The next
!> step is to set up the tri-diagonal matrix for calculating the
!> new temperature by means of a semi-implicit central scheme for the
!> vertical diffusion. Source terms which appear on the right hand sides
!> are due to the divergence of the solar radiation at the interfaces.
!> The subroutine is completed by solving the tri-diagonal linear
!> equation by means of a tri-diagonal solver.


SUBROUTINE temperature_advection(self)

   !! Advection/diffusion of the temperature field

   IMPLICIT NONE

!  Subroutine arguments
   class(type_temperature), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: rc
!---------------------------------------------------------------------------
   call self%logs%info('temperature_advection() - NOT CALLED',level=2)

   return
END SUBROUTINE temperature_advection

!---------------------------------------------------------------------------

!> Here, one time step for the temperature equation is performed.
!> First, preparations for the call to the advection schemes are
!> made, i.e.\ calculating the necessary metric coefficients.
!> After the call to the advection schemes, which actually perform
!> the advection (and horizontal diffusion) step as an operational
!> split step, the solar radiation at the interfaces ({\tt rad(k)})
!> is calculated
!> from given surface radiation ({\tt swr\_loc})
!> by means of a double exponential
!> approach, see equation (\ref{Light}) on page \pageref{Light}).
!> An option to reflect part of the short wave radiation that reaches the
!> bottom has been implemented. In very shallow waters - or with very clear
!> waters - a significant part of the incoming radiation will reach the
!> bottom. Setting swr\_bot\_refl\_frac to a value between 0 and 1 will
!> reflect this fraction of what ever the value of SWR is at the bottom.
!> The default value of swr\_bot\_refl\_frac is 0.
!> The reflection is only done if the ratio between the surface and bottom
!> values of SWR is greater than swr\_min\_bot\_frac (default 0.01).
!> Furthermore, the surface heat flux {\tt sfl\_loc} is given a
!> value.
!> The sea surface temperature is limited by the freezing point
!> temperature (as a most primitive sea ice model). The next
!> step is to set up the tri-diagonal matrix for calculating the
!> new temperature by means of a semi-implicit central scheme for the
!> vertical diffusion. Source terms which appear on the right hand sides
!> are due to the divergence of the solar radiation at the interfaces.
!> The subroutine is completed by solving the tri-diagonal linear
!> equation by means of a tri-diagonal solver.


SUBROUTINE temperature_vertical_diffusion(self,logs,rad,nuh)

   !! Advection/diffusion of the temperature field

   IMPLICIT NONE

!  Subroutine arguments
   class(type_temperature), intent(inout) :: self
   class(type_logging), intent(in) :: logs
   real(real64), dimension(:,:,:), intent(in) :: rad, nuh
      !! radiation and diffusivity

!  Local constants

!  Local variables
   integer :: rc

!---------------------------------------------------------------------------
   call logs%info('temperature_vertical_diffusion() - NOT CALLED',level=2)

   return
END SUBROUTINE temperature_vertical_diffusion

!---------------------------------------------------------------------------

END MODULE getm_temperature

