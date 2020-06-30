! Copyright (C) 2020 Bolding & Bruggeman

!> In this module, the salinity equation is processed by
!> reading in the namelist {\tt salt} and initialising the salinity field
!> (this is done in {\tt init\_salinity}),
!> and calculating the advection-diffusion-equation, which includes
!> penetrating short-wave radiation as source term (see {\tt do\_salinity}).

!> @note
!> ckeck dimension order of auxo and auxn
!> ckeck dimension order of a1, a2, a3, a4
!> @endnote

#ifdef _STATIC_
#include "dimensions.h"
#endif

MODULE getm_salinity

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
   use getm_operators

   IMPLICIT NONE

   PRIVATE  ! Private scope by default

!  Module constants
   real(real64), parameter :: avmols = 0._real64

!  Module types and variables
   real(real64) :: cnpar
   type, public :: type_salinity_configuration
      character(len=256) :: f = "salinity.nc"
   end type type_salinity_configuration

   type, public :: type_salinity
      !! author: Karsten Bolding
      !! version: v0.1
      !!
      !! Salinity type

      class(type_logging), pointer :: logs
      class(type_field_manager), pointer :: fm => null()
      class(type_getm_grid), pointer :: grid => null()
      class(type_advection), pointer :: advection => null()
      class(type_vertical_diffusion), pointer :: vertical_diffusion => null()
      TYPE(type_salinity_configuration) :: config

#ifdef _STATIC_
      real(real64), dimension(I3DFIELD) :: S = 10._real64
#else
      real(real64), dimension(:,:,:), allocatable :: S
#endif
!> remember to allocate in initialize

      real(real64) :: cnpar
      real(real64) :: avmolt

      contains

      procedure :: configuration => salinity_configuration
      procedure :: initialize => salinity_initialize
      procedure :: calculate => salinity_calculate

   end type type_salinity

CONTAINS

!---------------------------------------------------------------------------

SUBROUTINE salinity_configuration(self,logs,fm)

   !! Configure the salinity type

   IMPLICIT NONE

!  Subroutine arguments
   class(type_salinity), intent(out) :: self
   class(type_logging), intent(in), target :: logs
   type(type_field_manager), optional, target  :: fm

!  Local constants

!  Local variables
   integer :: rc
!---------------------------------------------------------------------------
   self%logs => logs
   call self%logs%info('salinity_configuration()',level=2)
   call self%logs%info('reading initial salinity from: ',level=3,msg2=trim(self%config%f))

   if (present(fm)) then
      self%fm => fm
   end if
   return
END SUBROUTINE salinity_configuration

!---------------------------------------------------------------------------

SUBROUTINE salinity_initialize(self,grid,advection,vertical_diffusion)

   !! Initialize the salinity field

   IMPLICIT NONE

!  Subroutine arguments
   class(type_salinity), intent(inout) :: self
   class(type_getm_grid), intent(in), target :: grid
   class(type_advection), intent(in), optional, target :: advection
   class(type_vertical_diffusion), intent(in), optional, target :: vertical_diffusion
      !! grid dimensions in case of dynamic memory allocation

!  Local constants

!  Local variables
   integer :: i,j
   integer :: stat
!---------------------------------------------------------------------------
   call self%logs%info('salinity_initialize()',level=2)

   self%grid => grid
   if (present(advection)) then
      self%advection => advection
   end if
   if (present(vertical_diffusion)) then
      self%vertical_diffusion => vertical_diffusion
   end if
#ifndef _STATIC_
   call mm_s('S',self%S,grid%l,grid%u,def=25._real64,stat=stat)
#endif
   if (associated(self%fm)) then
      call self%fm%register('salt', 'g/kg', 'absolute salinity', &
                            standard_name='sea_water_salinity', &
                            category='temperature_and_salinity', &
                            dimensions=(self%grid%dim_3d_ids), &
                            fill_value=-9999._real64, &
                            part_of_state=.true.)
      call self%fm%send_data('salt', self%S(grid%imin:grid%imax,grid%jmin:grid%jmax,grid%kmin:grid%kmax))
   end if
#define G self%grid
   do j=G%jmin,G%jmax
      do i=G%imin,G%imax
         if (G%mask(i,j) ==  0) then
            self%S(i,j,:) = -9999._real64
         end if
      end do
   end do
#undef G
   return
END SUBROUTINE salinity_initialize

!---------------------------------------------------------------------------

SUBROUTINE salinity_calculate(self,dt,nuh)

   !! Advection/diffusion of the salinity field

   IMPLICIT NONE

!  Subroutine arguments
   class(type_salinity), intent(inout) :: self
   real(real64), intent(in) :: dt
   real(real64), dimension(:,:,:), intent(in) :: nuh

!  Local constants

!  Local variables
   integer :: rc
!---------------------------------------------------------------------------
   call self%logs%info('salinity_calculate()',level=2)

#define G self%grid
   call self%advection%calculate(G%mask,G%hn,dt,self%cnpar,self%avmolt,nuh,self%S)
   call self%vertical_diffusion%calculate(G%mask,G%hn,dt,self%cnpar,self%avmolt,nuh,self%S)
#undef G

   return
END SUBROUTINE salinity_calculate

!---------------------------------------------------------------------------

END MODULE getm_salinity
