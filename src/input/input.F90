! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!#define _STATIC_

!>  This module provides all input from external files

MODULE getm_input

!   use time, only: type_time
   USE, INTRINSIC :: ISO_FORTRAN_ENV
!   use netcdf
!   use input_manager
   use memory_manager
   use logging
!   use getm_domain
!   use getm_physics
!   use getm_dynamics

   IMPLICIT NONE

!-----------------------------------------------------------------------------

   PRIVATE  ! Private scope by default

! Module constants

! Module types and variables
#if 0
!   type, extends(type_input_data), public :: type_input_bathymetry
   type, extends(type_2d_data), public :: type_input_2d_data
      !! author: Karsten Bolding
      !! version: v0.1
      !!
      !! Reading bathymetry

      integer :: ncid
      integer :: id
      integer :: dimids(2)
      integer :: dimlens(2)

      contains

      procedure :: initialize => initialize_2d_data
      procedure :: get => get_2d_data
      procedure :: close => close_2d_data

   end type type_input_2d_data
#endif

! Public members

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

END MODULE getm_input

#if 0
!-----------------------------------------------------------------------------
!SUBROUTINE initialize_2d_data(self,logs,varname,grid)
SUBROUTINE initialize_2d_data(self,logs,varname)
   !! Read the bathymetry from an external file and sets the mask

   IMPLICIT NONE

! Subroutine arguments
   class(type_input_2d_data), intent(inout) :: self
   class(type_logging), intent(in) :: logs
   character(len=256), intent(in) :: varname
!KB   class(type_getm_grid), intent(inout) :: grid
!KB      !! The grid type to hold bathymetry and grid size

! Local constants

! Local variables
!   integer :: i,j
!   integer, dimension(2) :: start, count
!-----------------------------------------------------------------------------
   call logs%info('initialize_2d_data()',level=2)

!   self%time_dependent = .false.

   ! open the NetCDF file with the bathymetry
   call check( nf90_open(trim(self%f), NF90_NOWRITE, self%ncid) )
!   write(*,*) self%ncid

   ! check existence of bathymetry variable
   call check(nf90_inq_varid(self%ncid, trim(varname), self%id))
!   write(*,*) self%id

   ! get dimension ids for the batymetry variable
   call check( nf90_inquire_variable(self%ncid, self%id, dimids = self%dimids(:)) )
!   write(*,*) self%dimids

   ! get dimension lengths
   call check(nf90_inquire_dimension(self%ncid, self%dimids(1), len=self%dimlens(1)))
   call check(nf90_inquire_dimension(self%ncid, self%dimids(2), len=self%dimlens(2)))
!   write(*,*) self%dimlens

   return
END SUBROUTINE initialize_2d_data

!-----------------------------------------------------------------------------

SUBROUTINE get_2d_data(self,logs)
   !! Read the bathymetry from an external file and sets the mask

   IMPLICIT NONE

! Subroutine arguments
   class(type_input_2d_data), intent(inout) :: self
   class(type_logging), intent(in) :: logs
      !! grid size in case of dynamic memory allocation
!KB   class(type_getm_grid), intent(inout) :: grid
!KB      !! The grid type to hold bathymetry and grid size

! Local constants

! Local variables
   integer :: i,j
   integer, dimension(2) :: start, count
!-----------------------------------------------------------------------------
   call logs%info('get_2d_data()',level=2)

!   start = 1
!   count = self%dimlens

   ! call check(nf90_get_var(self%ncid, self%id, grid%H(1:iself%dimlens(1),1:dimlens(2)), start=start, count=count))
   !call check(nf90_get_var(self%ncid, self%id, grid%H(1:self%dimlens(1),1:self%dimlens(2))))

   return
END SUBROUTINE get_2d_data

!-----------------------------------------------------------------------------

SUBROUTINE close_2d_data(self,logs)
   !! Read the bathymetry from an external file and sets the mask

   IMPLICIT NONE

! Subroutine arguments
   class(type_input_2d_data), intent(inout) :: self
   class(type_logging), intent(in) :: logs

! Local constants

! Local variables
!-----------------------------------------------------------------------------
   call logs%info('close_2d_data()',level=2)

   call check( nf90_close(self%ncid) )
   return
END SUBROUTINE close_2d_data

!-----------------------------------------------------------------------------

subroutine check(status)
   integer, intent ( in) :: status

   if(status /= nf90_noerr) then
     print *, trim(nf90_strerror(status))
     stop "Stopped"
   end if
end subroutine check

#endif
