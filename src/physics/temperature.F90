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
   use getm_operators

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

      class(type_logging), pointer :: logs => null()
      class(type_field_manager), pointer :: fm => null()
      class(type_getm_domain), pointer :: domain
      class(type_advection), pointer :: advection => null()
      class(type_vertical_diffusion), pointer :: vertical_diffusion => null()
      TYPE(type_temperature_configuration) :: config

#ifdef _STATIC_
   real(real64), dimension(I3DFIELD), target :: T = 10._real64
#else
   real(real64), dimension(:,:,:), allocatable :: T
#endif
      integer :: advection_scheme=1
      real(real64) :: cnpar
      real(real64) :: avmolt=0._real64

      contains

      procedure :: configuration => temperature_configuration
      procedure :: initialize => temperature_initialize
      procedure :: calculate => temperature_calculate

   end type type_temperature

CONTAINS

!---------------------------------------------------------------------------

SUBROUTINE temperature_configuration(self,logs,fm)

   !! Configure the the temperature module

   IMPLICIT NONE

!  Subroutine arguments
   class(type_temperature), intent(out) :: self
   class(type_logging), intent(in), target, optional :: logs
   type(type_field_manager), target, optional  :: fm

!  Local constants

!  Local variables
   integer :: rc
!---------------------------------------------------------------------------
   if (present(logs)) then
      self%logs => logs
      call self%logs%info('temperature_configuration()',level=2)
      call self%logs%info('reading initial temperature from: ',level=3,msg2=trim(self%config%f))
   end if
   if (present(fm)) then
      self%fm => fm
   end if
END SUBROUTINE temperature_configuration

!---------------------------------------------------------------------------

SUBROUTINE temperature_initialize(self,domain,advection,vertical_diffusion)

   !! Initialize the temperature field

   IMPLICIT NONE

!  Subroutine arguments
   class(type_temperature), intent(inout) :: self
   class(type_getm_domain), intent(in), target :: domain
   class(type_advection), intent(in), optional, target :: advection
   class(type_vertical_diffusion), intent(in), optional, target :: vertical_diffusion

!  Local constants

!  Local variables
   integer :: i,j
   integer :: stat
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('temperature_initialize()',level=2)
   self%domain => domain
   if (present(advection)) then
      self%advection => advection
   end if
   if (present(vertical_diffusion)) then
      self%vertical_diffusion => vertical_diffusion
   end if
   TGrid: associate( TG => self%domain%T )
#ifndef _STATIC_
   call mm_s('T',self%T,TG%l,TG%u,def=15._real64,stat=stat)
#endif
   if (associated(self%fm)) then
      call self%fm%register('temp', 'Celsius', 'potential temperature', &
                            standard_name='sea_water_temperature', &
                            category='temperature_and_salinity', &
                            dimensions=(TG%dim_3d_ids), &
                            fill_value=-9999._real64, &
                            part_of_state=.true.)
      call self%fm%send_data('temp', self%T(TG%imin:TG%imax,TG%jmin:TG%jmax,TG%kmin:TG%kmax))
   end if
   do j=TG%jmin,TG%jmax
      do i=TG%imin,TG%imax
         if (TG%mask(i,j) ==  0) then
            self%T(i,j,:) = -9999._real64
         end if
      end do
   end do
   end associate TGrid
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


SUBROUTINE temperature_calculate(self,dt,uk,vk,nuh,rad,shf)

   !! Advection/diffusion of the temperature field

   IMPLICIT NONE

!  Subroutine arguments
   class(type_temperature), intent(inout) :: self
   real(real64), intent(in) :: dt
#define _T3_ self%domain%T%l(1):,self%domain%T%l(2):,self%domain%T%l(3):
   real(real64), dimension(:,:,:), intent(in) :: uk(_T3_)
   real(real64), dimension(:,:,:), intent(in) :: vk(_T3_)
   real(real64), dimension(:,:,:), intent(in) :: nuh(_T3_)
   real(real64), dimension(:,:,:), intent(in) :: rad(_T3_)
#undef _T3_
#define _T2_ self%domain%T%l(1):,self%domain%T%l(2):
   real(real64), intent(in) :: shf(_T2_)
#undef _T2_

!  Local constants

!  Local variables
   integer :: i,j,k
   real(real64), allocatable :: ea4(:,:,:)
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('temperature_calculate()',level=2)

   allocate(ea4,mold=self%T) !KB use imin,imax ... instead to be consistent with a1, a2 ... in diffusion solver

   TGrid: associate( TG => self%domain%T )
   UGrid: associate( UG => self%domain%U )
   VGrid: associate( VG => self%domain%V )
   call self%advection%calculate(self%advection_scheme,UG,uk,VG,vk,dt,TG,self%T)
   end associate VGrid
   end associate UGrid
   ! add extra processes - heating due to short wave radiation and surface heat flux (sensible, latent and net longwave)
   do k=TG%kmin,TG%kmax
      do j=TG%jmin,TG%jmax
         do i=TG%imin,TG%imax
            if (TG%mask(i,j) > 0) ea4(i,j,k) = dt*(rad(i,j,k)-rad(i,j,k-1))
         end do
      end do
   end do
   k=TG%kmax
   do j=TG%jmin,TG%jmax
      do i=TG%imin,TG%imax
         if (TG%mask(i,j) > 0) ea4(i,j,k) = ea4(i,j,k)+dt*shf(i,j) !KB some scaling
      end do
   end do
   call self%vertical_diffusion%calculate(dt,self%cnpar,TG%mask,TG%hn,TG%hn,self%avmolt,nuh,self%T,ea4=ea4)
   end associate TGrid
END SUBROUTINE temperature_calculate

!---------------------------------------------------------------------------

END MODULE getm_temperature

!---------------------------------------------------------------------------
