! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!#define _STATIC_

!>  This module provides all input from external files

MODULE getm_bathymetry

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use memory_manager
   use logging
   use input_module
   use getm_domain

   IMPLICIT NONE

!-----------------------------------------------------------------------------

   PRIVATE  ! Private scope by default

! Module constants

! Module types and variables

!   type, extends(type_input_data), public :: type_input_bathymetry

!-----------------------------------------------------------------------------
   type, public :: type_bathymetry

      TYPE(type_netcdf_input) :: depth
      TYPE(type_netcdf_input) :: c1,c2
      TYPE(type_netcdf_input) :: lon,lat
      TYPE(type_netcdf_input) :: xx,yx,convx
      TYPE(type_netcdf_input) :: lonx,latx

   CONTAINS

      procedure :: initialize => initialize_bathymetry

   end type type_bathymetry

CONTAINS

!-----------------------------------------------------------------------------
SUBROUTINE initialize_bathymetry(self,logs,grid,domain_type)
   !! Read the bathymetry from an external file and sets the mask

   IMPLICIT NONE

! Subroutine arguments
   class(type_bathymetry), intent(inout) :: self
   class(type_logging), intent(in) :: logs
      !! grid size in case of dynamic memory allocation
   class(type_getm_grid), intent(inout), target :: grid
      !! The grid type to hold bathymetry and grid size
   integer, intent(inout) :: domain_type

! Local constants

! Local variables
  integer :: stat
!-----------------------------------------------------------------------------
   call logs%info('initialize_bathymetry()',level=1)

!  reading bathymetry and get coordinate names
   call self%depth%initialize()
   self%depth%p2dreal64 => grid%H(grid%imin:grid%imax,grid%jmin:grid%jmax)
   call self%depth%get()
   call self%depth%close()

   if( trim(self%depth%cord_names(1)) == 'x' .and. trim(self%depth%cord_names(2)) == 'y') then
      domain_type = 1
   end if
   if( trim(self%depth%cord_names(1)) == 'lon' .and. trim(self%depth%cord_names(2)) == 'lat') then
      domain_type = 2
   end if
   if( trim(self%depth%cord_names(1)) == 'lon_kurt' .and. trim(self%depth%cord_names(2)) == 'lat_kurt') then
      domain_type = 3
   end if
   call self%depth%print_info()

!  reading the first cordinate
   self%c1%f = self%depth%f
   self%c1%v = self%depth%cord_names(1)
   call self%c1%initialize()
   self%c1%p1dreal64 => grid%c1(grid%imin:grid%imax)
   call self%c1%get()
   call self%c1%close()

!  reading the second cordinate
   self%c2%f = self%depth%f
   self%c2%v = self%depth%cord_names(2)
   call self%c2%initialize()
   self%c2%p1dreal64 => grid%c2(grid%jmin:grid%jmax)
   call self%c2%get()
   call self%c2%close()

   select case (domain_type)
      case(1)
!! @note
!! cartesian: Check for lon and lat grids - will be 2D
!! @endnote
#if 0
         self%lon%f = self%depth%f
         self%lon%v = 'lonc'
         call self%lon%initialize()
!KB         call self%lon%initialize(start=stat)
         call mm_s('lonc',self%pgrid%lon,self%pgrid%l(1:2),self%pgrid%u(1:2),def=-9999._real64,stat=stat)
         self%lon%p2dreal64 => self%pgrid%lon
         call self%lon%get(error_handler=stat)
         call self%lon%close()
         self%lat%f = self%depth%f
         self%lat%v = 'latc'
         call self%lat%initialize()
         call mm_s('latc',self%pgrid%lat,self%pgrid%l(1:2),self%pgrid%u(1:2),def=-9999._real64,stat=stat)
         self%lat%p2dreal64 => self%pgrid%lat
         call self%lat%get(error_handler=stat)
         call self%lat%close()
#endif
      case(2)
!! @note
!! sperical: read additional sperical fields
!! @endnote
      case(3)
!! @note
!! curvilinear: read xx, yy .....
!! @endnote
#if 0
         self%xx%f = self%depth%f
         self%xx%v = 'xx'
         call self%xx%initialize(stat=stat)
         call mm_s('xx',self%pgrid%lon,self%pgrid%l(1:2),self%pgrid%u(1:2),def=-9999._real64,stat=stat)
         self%lon%p2dreal64 => self%pgrid%xx
         call self%xx%get(error_handler=stat)
         call self%xx%close()
#endif
   end select

   where(grid%H > -10._real64) grid%mask = 1

   return
END SUBROUTINE initialize_bathymetry

!-----------------------------------------------------------------------------

END MODULE getm_bathymetry
