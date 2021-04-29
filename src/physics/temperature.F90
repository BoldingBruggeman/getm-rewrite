! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!> In this module, the temperature equation is processed by
!> reading in the namelist {\tt temp} and initialising the temperature field
!> (this is done in {\tt init\_temperature}),
!> and calculating the advection-diffusion-equation, which includes
!> penetrating short-wave radiation as source term (see {\tt do\_temperature}).

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
   real(real64), parameter :: cp=3985._real64
   real(real64), parameter :: rho0=1025._real64

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
      real(real64), dimension(:,:), allocatable :: Tbdy
      real(real64), dimension(:,:), allocatable :: sst
      real(real64), dimension(:,:,:), allocatable :: ea4
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
   if (domain%nbdy > 0) then
      call mm_s('Tbdy',self%Tbdy,(/TG%l(3)-1,1/),(/TG%u(3),domain%nbdyp/),def=15._real64,stat=stat)
   end if
   call mm_s('sst',self%sst,TG%l(1:2),TG%u(1:2),def=15._real64,stat=stat)
   call mm_s('ea4',self%ea4,self%T,def=0._real64,stat=stat)
#endif
   if (associated(self%fm)) then
      call self%fm%register('temp', 'Celsius', 'conservative temperature', &
                            standard_name='sea_water_temperature', &
                            category='temperature_and_salinity', &
                            dimensions=(TG%dim_3d_ids), &
                            fill_value=-9999._real64, &
                            part_of_state=.true.)
      call self%fm%send_data('temp', self%T(TG%imin:TG%imax,TG%jmin:TG%jmax,TG%kmin:TG%kmax))
      call self%fm%register('sst', 'Celsius', 'potential temperature', &
                            standard_name='sea_water_temperature', &
                            category='temperature_and_salinity', &
                            dimensions=(TG%dim_2d_ids), &
                            fill_value=-9999._real64, &
                            part_of_state=.false.)
      call self%fm%send_data('sst', self%sst(TG%imin:TG%imax,TG%jmin:TG%jmax))
   end if
   do j=TG%jmin,TG%jmax
      do i=TG%imin,TG%imax
         if (TG%mask(i,j) ==  0) then
            self%T(i,j,:) = -9999._real64
            self%sst(i,j) = -9999._real64
         end if
      end do
   end do
   end associate TGrid
END SUBROUTINE temperature_initialize

!---------------------------------------------------------------------------

SUBROUTINE temperature_calculate(self,dt,uk,vk,nuh,rad,shf)

   !! Short wave radiation, surface heat fluxes and advection/diffusion of the temperature field

   IMPLICIT NONE

!  Subroutine arguments
   class(type_temperature), intent(inout) :: self
   real(real64), intent(in) :: dt
#define _T3_ self%domain%T%l(1):,self%domain%T%l(2):,self%domain%T%l(3):
   real(real64), dimension(:,:,:), intent(in) :: uk(_T3_)
   real(real64), dimension(:,:,:), intent(in) :: vk(_T3_)
   real(real64), dimension(:,:,:), intent(in) :: nuh(_T3_)
#undef _T3_
#define _T3_ self%domain%T%l(1):,self%domain%T%l(2):,self%domain%T%l(3)-1:
   real(real64), dimension(:,:,:), intent(in) :: rad(_T3_)
#undef _T3_
#define _T2_ self%domain%T%l(1):,self%domain%T%l(2):
   real(real64), intent(in) :: shf(_T2_)
#undef _T2_

!  Local constants

!  Local variables
   integer :: i,j,k
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('temperature_calculate()',level=2)

   TGrid: associate( TG => self%domain%T )
   UGrid: associate( UG => self%domain%U )
   VGrid: associate( VG => self%domain%V )
   call self%advection%calculate(self%advection_scheme,UG,uk,VG,vk,dt,TG,self%T)
   end associate VGrid
   end associate UGrid
   ! add extra processes
   ! short wave radiation and surface heat flux (sensible, latent and net longwave)
   do k=TG%kmin,TG%kmax
      do j=TG%jmin,TG%jmax
         do i=TG%imin,TG%imax
            self%ea4(i,j,k)=0._real64
            if (TG%mask(i,j) > 0) self%ea4(i,j,k) = dt*(rad(i,j,k)-rad(i,j,k-1))/(rho0*cp)
         end do
      end do
   end do
   k=TG%kmax
   do j=TG%jmin,TG%jmax
      do i=TG%imin,TG%imax
         if (TG%mask(i,j) > 0) self%ea4(i,j,k) = self%ea4(i,j,k)+dt*shf(i,j)/(rho0*cp)
      end do
   end do
   !KB max dt = 0.5*dz^2/nuh
   call self%vertical_diffusion%calculate(dt,self%cnpar,TG%mask,TG%hn,TG%hn,self%avmolt,nuh,self%T,ea4=self%ea4)
   end associate TGrid
END SUBROUTINE temperature_calculate

!---------------------------------------------------------------------------

END MODULE getm_temperature
