! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!> In this module, the salinity equation is processed by
!> reading in the namelist {\tt salt} and initialising the salinity field
!> (this is done in {\tt init\_salinity}),
!> and calculating the advection-diffusion-equation, which includes
!> penetrating short-wave radiation as source term (see {\tt do\_salinity}).

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

      class(type_logging), pointer :: logs => null()
      class(type_field_manager), pointer :: fm => null()
      class(type_getm_domain), pointer :: domain => null()
      class(type_advection), pointer :: advection => null()
      class(type_vertical_diffusion), pointer :: vertical_diffusion => null()
      TYPE(type_salinity_configuration) :: config

#ifdef _STATIC_
      real(real64), dimension(I3DFIELD) :: S = 10._real64
#else
      real(real64), dimension(:,:,:), allocatable :: S
      real(real64), dimension(:,:), allocatable :: Sbdy
#endif
      integer :: advection_scheme=1
      real(real64) :: cnpar
      real(real64) :: avmols=0._real64

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
   class(type_logging), intent(in), target, optional :: logs
   type(type_field_manager), target, optional  :: fm

!  Local constants

!  Local variables
   integer :: rc
!---------------------------------------------------------------------------
   if (present(logs)) then
      self%logs => logs
      call self%logs%info('salinity_configuration()',level=2)
      call self%logs%info('reading initial salinity from: ',level=3,msg2=trim(self%config%f))
   end if
   if (present(fm)) then
      self%fm => fm
   end if
END SUBROUTINE salinity_configuration

!---------------------------------------------------------------------------

SUBROUTINE salinity_initialize(self,domain,advection,vertical_diffusion)

   !! Initialize the salinity field

   IMPLICIT NONE

!  Subroutine arguments
   class(type_salinity), intent(inout) :: self
   class(type_getm_domain), intent(in), target :: domain
   class(type_advection), intent(in), optional, target :: advection
   class(type_vertical_diffusion), intent(in), optional, target :: vertical_diffusion
      !! grid dimensions in case of dynamic memory allocation

!  Local constants

!  Local variables
   integer :: i,j
   integer :: stat
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('salinity_initialize()',level=2)

   self%domain => domain
   if (present(advection)) then
      self%advection => advection
   end if
   if (present(vertical_diffusion)) then
      self%vertical_diffusion => vertical_diffusion
   end if

   TGrid: associate( TG => self%domain%T )
#ifndef _STATIC_
   call mm_s('S',self%S,TG%l,TG%u,def=25._real64,stat=stat)
   if (domain%nbdy > 0) then
      call mm_s('Sbdy',self%Sbdy,(/TG%l(3),1/),(/TG%u(3),domain%nbdyp/),def=15._real64,stat=stat)
   end if
#endif
   if (associated(self%fm)) then
      call self%fm%register('salt', 'g/kg', 'absolute salinity', &
                            standard_name='sea_water_salinity', &
                            category='temperature_and_salinity', &
                            dimensions=(TG%dim_3d_ids), &
                            fill_value=-9999._real64, &
                            part_of_state=.true.)
      call self%fm%send_data('salt', self%S(TG%imin:TG%imax,TG%jmin:TG%jmax,TG%kmin:TG%kmax))
   end if
   do j=TG%jmin,TG%jmax
      do i=TG%imin,TG%imax
         if (TG%mask(i,j) ==  0) then
            self%S(i,j,:) = -9999._real64
         end if
      end do
   end do
   end associate TGrid
END SUBROUTINE salinity_initialize

!---------------------------------------------------------------------------

SUBROUTINE salinity_calculate(self,dt,uk,vk,nuh)

   !! Advection/diffusion of the salinity field

   IMPLICIT NONE

!  Subroutine arguments
   class(type_salinity), intent(inout) :: self
   real(real64), intent(in) :: dt
#define _T3_ self%domain%T%l(1):,self%domain%T%l(2):,self%domain%T%l(3):
   real(real64), intent(in) :: uk(_T3_)
   real(real64), intent(in) :: vk(_T3_)
   real(real64), intent(in) :: nuh(_T3_)
#undef _T3_

!  Local constants

!  Local variables
   integer :: rc
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('salinity_calculate()',level=2)

   TGrid: associate( TG => self%domain%T )
   UGrid: associate( UG => self%domain%U )
   VGrid: associate( VG => self%domain%V )
   call self%advection%calculate(self%advection_scheme,UG,uk,VG,vk,dt,TG,self%S)
   end associate VGrid
   end associate UGrid
   call self%vertical_diffusion%calculate(dt,self%cnpar,TG%mask,TG%hn,TG%hn,self%avmols,nuh,self%S)

   call self%domain%mirror_bdys(TG,self%S)
   end associate TGrid
END SUBROUTINE salinity_calculate

!---------------------------------------------------------------------------

END MODULE getm_salinity
