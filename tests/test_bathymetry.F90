! Copyright (C) 2020 Bolding & Bruggeman

!> @note
!> maybe put all grid tests in one big program - test_bathymetry
!> @endnote

PROGRAM test_bathymetry
   !! Testing calculation of time varying depths at S, U and V points

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use logging
   use memory_manager
   use getm_domain
   use getm_bathymetry
   IMPLICIT NONE

!  Local constants

!  Local variables
   TYPE(type_logging) :: logs
   TYPE(type_getm_grid), target :: grid
   TYPE(type_bathymetry) :: bathymetry

   character(len=256) :: varname
   integer :: ndims, dimlens
   integer :: domain_type
!   integer, dimension(:), allocatable :: dimlens

!-----------------------------------------------------------------------

   if (command_argument_count() .ne. 1 ) then
      write(*,*) 'ERROR, must be called like:'
      STOP ' test_bathymetry <name of bathymetry file>'
   end if

   call get_command_argument(1,bathymetry%depth%f)
   bathymetry%depth%v = 'bathymetry'
   !call self%domain%T%configure(self%logs,kmin=kmin,kmax=kmax,halo=(/2, 2, 0/))
   call bathymetry%initialize(logs,grid,domain_type)
   call grid%print(OUTPUT_UNIT)
   call grid%print_info(logs)
   call grid%print_mask(OUTPUT_UNIT)

!-----------------------------------------------------------------------

!CONTAINS 

!-----------------------------------------------------------------------

END PROGRAM test_bathymetry
