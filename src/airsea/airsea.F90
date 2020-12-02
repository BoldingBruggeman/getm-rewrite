! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!#define _STATIC_

!>  This module provides meteo forcing for GETM

MODULE getm_airsea

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use memory_manager
   use logging
   use field_manager
   use getm_domain, only: type_getm_domain

   IMPLICIT NONE

!-----------------------------------------------------------------------------

   PRIVATE  ! Private scope by default

! Module constants
#if 0
   real(real64), parameter :: rearth_default = 6378815.
      !! The earth radius
   real(real64), parameter :: pi             = 3.141592654
      !! $$pi$$
   real(real64), parameter :: deg2rad        = pi/180.
      !! degrees to radians
#endif

! Module types and variables
   type, public :: type_getm_meteo
      !! author: Karsten Bolding

      real(real64), dimension(:), allocatable :: c1,c2
        !! cordinate variables from NetCDF file
      real(real64), dimension(:), allocatable :: lon,lat
        !! longitude and latitude of grid
      real(real64), dimension(:,:), allocatable :: d2m, sp, u10, v10, tcc, t2m
         !! meteorological variables from e.g. ERA5 or similar
   end type type_getm_meteo

   type, public :: type_getm_airsea
      !! author: Karsten Bolding

      class(type_logging), pointer :: logs => null()
      class(type_field_manager), pointer :: fm => null()
      class(type_getm_domain), pointer :: domain

      real(real64) :: taux0=0._real64, tauy0=0._real64
      real(real64) :: swr0=0._real64
      real(real64) :: shf0=0._real64

#ifdef _STATIC_
      real(real64), dimension(E2DFIELD) :: p = 10._real64
#else
      real(real64), dimension(:,:), allocatable :: d2m, sp, u10, v10, tcc, t2m
         !! meteorological variables from e.g. ERA5 or similar
         !! interpolated to the GETM model grid
      real(real64), dimension(:,:), allocatable :: swr, albedo
      real(real64), dimension(:,:), allocatable :: taux, tauy
      real(real64), dimension(:,:), allocatable :: sensible, latent
      real(real64), dimension(:,:), allocatable :: net_long_wave
      real(real64), dimension(:,:), allocatable :: shf
         !! short wave radiation and albedo fields
#endif

      contains

      procedure :: configure => airsea_configure
      procedure :: initialize => airsea_initialize
      procedure :: update => airsea_update

   end type type_getm_airsea

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

SUBROUTINE airsea_configure(self,logs,fm)

   !! Initialize the meteo variables

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_airsea), intent(inout) :: self
   class(type_logging), intent(in), target, optional :: logs
   class(type_field_manager), intent(in), target, optional :: fm

!  Local constants

!  Local variables
   integer :: stat
!-----------------------------------------------------------------------
   if (present(logs)) then
      self%logs => logs
      call self%logs%info('airsea_initialize()',level=1)
   end if
   if (present(fm)) then
      self%fm => fm
   end if
END SUBROUTINE airsea_configure

!-----------------------------------------------------------------------

SUBROUTINE airsea_initialize(self,domain)

   !! Initialize the meteo variables

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_airsea), intent(inout) :: self
   class(type_getm_domain), intent(in), target :: domain

!  Local constants

!  Local variables
   integer :: stat
   type (type_field), pointer :: f
!-----------------------------------------------------------------------
   self%domain => domain
   TGrid: associate( TG => self%domain%T )
#ifndef _STATIC_
   call mm_s('swr',self%swr,TG%l(1:2),TG%u(1:2),def=0._real64,stat=stat)
   call mm_s('albedo',self%albedo,TG%l(1:2),TG%u(1:2),def=0._real64,stat=stat)
   call mm_s('taux',self%taux,TG%l(1:2),TG%u(1:2),def=0._real64,stat=stat)
   call mm_s('tauy',self%tauy,TG%l(1:2),TG%u(1:2),def=0._real64,stat=stat)
   call mm_s('sp',self%sp,TG%l(1:2),TG%u(1:2),def=0._real64,stat=stat)
   call mm_s('shf',self%shf,TG%l(1:2),TG%u(1:2),def=0._real64,stat=stat)
#endif
   if (associated(self%fm)) then
      call self%fm%register('sp', 'Pa', 'atmospheric surface pressure', &
                            standard_name='', &
                            dimensions=(self%domain%T%dim_2d_ids), &
   !KB                         output_level=output_level_debug, &
                            part_of_state=.false., &
                            category='airsea', field=f)
      call self%fm%send_data('sp', self%sp(TG%imin:TG%imax,TG%jmin:TG%jmax))
      call self%fm%register('tausx', 'Pa', 'surface stress - x', &
                            standard_name='', &
                            dimensions=(self%domain%T%dim_2d_ids), &
   !KB                         output_level=output_level_debug, &
                            part_of_state=.false., &
                            category='airsea', field=f)
      call self%fm%send_data('tausx', self%taux(TG%imin:TG%imax,TG%jmin:TG%jmax))
      call self%fm%register('tausy', 'Pa', 'surface stress - y', &
                            standard_name='', &
                            dimensions=(self%domain%T%dim_2d_ids), &
   !KB                         output_level=output_level_debug, &
                            part_of_state=.false., &
                            category='airsea', field=f)
      call self%fm%send_data('tausy', self%tauy(TG%imin:TG%imax,TG%jmin:TG%jmax))
   end if
   self%swr=self%swr0
   self%shf=self%shf0
   self%taux=self%taux0
   self%tauy=self%tauy0
   end associate TGrid
END SUBROUTINE airsea_initialize

!-----------------------------------------------------------------------

SUBROUTINE airsea_update(self,n)

   !! Update the meteo variables

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_airsea), intent(inout) :: self
   integer, intent(in) :: n

!  Local constants
   real(real64) :: p0=101325._real64

!  Local variables
   integer :: i,j
   real(real64) :: dpdx, dpdy
!-----------------------------------------------------------------------
   dpdx = n*1._real64/(self%domain%T%u(1)-self%domain%T%l(1))
   dpdy = n*1._real64/(self%domain%T%u(2)-self%domain%T%l(2))
   do j=self%domain%T%l(2),self%domain%T%u(2)
      do i=self%domain%T%l(1),self%domain%T%u(1)
         self%sp(i,j) = p0+i*dpdx+j*dpdy
         self%sp(i,j) = p0
         self%taux(i,j) = (1._real64+n)/10._real64*self%taux0
         self%tauy(i,j) = (1+n)/10*self%tauy0
         self%taux(i,j) = self%taux0
         self%tauy(i,j) = self%tauy0
         self%swr(i,j) = self%swr0
         self%shf(i,j) = self%shf0
      end do
   end do
END SUBROUTINE airsea_update

!-----------------------------------------------------------------------

END MODULE getm_airsea
