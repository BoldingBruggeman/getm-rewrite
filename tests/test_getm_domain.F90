! Copyright (C) 2020 Bolding & Bruggeman

!> @note
!> maybe put all grid tests in one big program - test_bathymetry
!> @endnote

PROGRAM test_getm_domain
   !! Testing calculation of time varying depths at S, U and V points

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use logging
   use memory_manager
   use getm_domain
   use getm_bathymetry
!   use getm_input_bathymetry
   IMPLICIT NONE

!  Local constants
   integer :: halo=2

!  Local variables
   TYPE(type_logging) :: logs
   TYPE(type_getm_grid), pointer :: pgrid
   TYPE(type_bathymetry) :: bathymetry
   class(type_getm_domain), allocatable, target :: domain

   character(len=256) :: varname
   integer :: ndims, dimlens
!   integer :: domain_type
!   integer, dimension(:), allocatable :: dimlens

!-----------------------------------------------------------------------

   if (command_argument_count() .ne. 1 ) then
      write(*,*) 'ERROR, must be called like:'
      STOP ' test_bathymetry <name of bathymetry file>'
   end if

   allocate(domain)

   call get_command_argument(1,bathymetry%depth%f)
   bathymetry%depth%v = 'bathymetry'
   pgrid => domain%T
   call pgrid%create(kmin=-1,kmax=-1)
   call pgrid%create(halo=(/halo, halo, 0/))
   call bathymetry%initialize(logs,pgrid,domain%domain_type)
   call domain%configure(logs)
   call domain%report()

!-----------------------------------------------------------------------

!CONTAINS 

!-----------------------------------------------------------------------

END PROGRAM test_getm_domain
