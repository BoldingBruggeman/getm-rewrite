! Copyright (C) 2020 Bolding & Bruggeman

#ifdef _STATIC_
#include "dimensions.h"
#endif

MODULE getm_mixing

   USE, INTRINSIC :: ISO_FORTRAN_ENV

   use memory_manager
   use logging
   use field_manager
   use getm_domain, only: type_getm_domain
!   use gsw_mod_toolbox, only: gsw_rho

   IMPLICIT NONE

   PRIVATE  ! Private scope by default

!  Module constants

!  Module types and variables

   type, public :: type_mixing_config
      !! author: Karsten Bolding
      !! version: v0.1
      !!
      !! Turbulence configuration

      character(len=256) :: bathymetry = "bathymetry.nc"
      integer :: kurt = 1

   end type type_mixing_config

   type, public :: type_getm_mixing
      !! author: Karsten Bolding
      !! version: v0.1
      !!
      !! Turbulence type

      TYPE(type_mixing_config) :: config

      class(type_logging), pointer :: logs
      class(type_getm_domain), pointer :: domain
      class(type_field_manager), pointer :: fm => null()




#ifdef _STATIC_
      real(real64), dimension(I3DFIELD) :: tke = 10._real64
#else
      real(real64), dimension(:,:,:), allocatable :: num
         !! turbulent viscosity
      real(real64), dimension(:,:,:), allocatable :: nuh
         !! turbulent diffusivity of heat (scalars)
      real(real64), dimension(:,:,:), allocatable :: tke
         !! turbulent kinetic energy
      real(real64), dimension(:,:,:), allocatable :: eps
         !! turbulent kinetic energy
      real(real64), dimension(:,:,:), allocatable :: NN
         !! Brunt-Vaisalla frequency
      real(real64), dimension(:,:,:), allocatable :: SS
         !! shear frequency
#endif
      contains

      procedure :: configuration => mixing_configuration
      procedure :: initialize => mixing_initialize
      procedure :: calculate => mixing_calculate
!      procedure :: SS => shear_frequency
!      procedure :: NN => bouyancy_frequency

   end type type_getm_mixing

   INTERFACE
      module subroutine shear_frequency(self,u,v)
         class(type_getm_mixing), intent(inout) :: self
         real(real64), dimension(:,:,:), intent(in) :: u,v
      end subroutine shear_frequency
      module subroutine bouyancy_frequency(self,bouy)
         class(type_getm_mixing), intent(inout) :: self
         real(real64), dimension(:,:,:), intent(in) :: bouy
      end subroutine bouyancy_frequency
   END INTERFACE

!---------------------------------------------------------------------------

CONTAINS

!---------------------------------------------------------------------------

SUBROUTINE mixing_configuration(self,logs,fm)

   !! Feeds your cats and dogs, if enough food is available. If not enough
   !! food is available, some of your pets will get angry.

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_mixing), intent(inout) :: self
   class(type_logging), intent(in), target :: logs
   type(type_field_manager), optional, target  :: fm

!  Local constants

!  Local variables
   character(len=256) :: str
!---------------------------------------------------------------------------
   self%logs => logs
   call self%logs%info('mixing_configuration()',level=2)
   write(str,'(I04)') self%config%kurt
   call self%logs%info('config->kurt:',msg2=trim(str),level=3)
   if (present(fm)) then
      self%fm => fm
   end if
   return
END SUBROUTINE mixing_configuration

!---------------------------------------------------------------------------

SUBROUTINE mixing_initialize(self,domain)

   !! Feeds your cats and dogs, if enough food is available. If not enough
   !! food is available, some of your pets will get angry.

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_mixing), intent(inout) :: self
   class(type_getm_domain), intent(in) :: domain

!  Local constants

!  Local variables
   integer :: imin,imax,jmin,jmax,kmin,kmax
   integer :: rc
!-----------------------------------------------------------------------------
   call self%logs%info('mixing_initialize()',level=2)

#if 0
#ifndef _STATIC_
   call mm_s('tke',self%tke,self%G%l,self%G%u,def=15._real64,stat=stat)

!KB   call mm_s('tke',self%tke)
!KB   call mm_s('num',self%num)
!KB   call mm_s('nuh',self%nuh)
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
#endif

  return
END SUBROUTINE mixing_initialize

!---------------------------------------------------------------------------

!SUBROUTINE mixing_calculate(self,domain,S,T,p)
SUBROUTINE mixing_calculate(self,logs,SS,NN)

   !! Feeds your cats and dogs, if enough food is available. If not enough
   !! food is available, some of your pets will get angry.

!KB   use gsw_mod_toolbox, only: gsw_rho
!> @note
!> use NN caculation from gws_toolbox
!> @endnote

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_mixing), intent(inout) :: self
!   class(type_getm_domain), intent(in) :: domain
   class(type_logging), intent(in) :: logs
   real(real64), dimension(-1:,-1:,0:), intent(in) :: SS,NN
      !! shear stress and Brunt-Vaisala frequency

!  Local constants

!  Local variables

!-----------------------------------------------------------------------------
   call logs%info('mixing_calculate()',level=2)

  return
END SUBROUTINE mixing_calculate

!---------------------------------------------------------------------------

END MODULE getm_mixing

!>
!>

