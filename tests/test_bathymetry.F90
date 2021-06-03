! Copyright (C) 2020 Bolding & Bruggeman

PROGRAM test_bathymetry
   !! Test reading the toppo.nc file created by the Python domain generator
   !! 1) get the size of the bathymetry field.
   !! 2) configure the Arakawa C-grid calculation domain - i.e. allocate all grid related variables
   !! 3) initialize fields for all 4 grids by reading them from the file
   !! 4) print info on the domain as well as the mask for T-points

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use logging
   use getm_domain
   use getm_bathymetry
   IMPLICIT NONE

!  Local constants

!  Local variables
   TYPE(type_logging) :: logs
   TYPE(type_getm_domain) :: domain
   TYPE(type_bathymetry) :: bathymetry

   character(len=256) :: varname
   integer :: ndims, dimlens
   integer :: domain_type
   integer :: imin=1,imax,jmin=1,jmax,kmin=-1,kmax=-1

!-----------------------------------------------------------------------

   if (command_argument_count() .ne. 1 ) then
      write(*,*) 'ERROR, must be called like:'
      STOP ' test_bathymetry <name of bathymetry file>'
   end if

   getdimensions: block
   call get_command_argument(1,bathymetry%H%f)
   call bathymetry%H%open()
   call bathymetry%H%dimlen('x',imax)
   call bathymetry%H%dimlen('y',jmax)
   call bathymetry%H%close()
   end block getdimensions

   call domain%configure(imin,imax,jmin,jmax,kmin,kmax,logs=logs)
   call bathymetry%initialize(logs,domain,domain%domain_type)
   call domain%T%print(OUTPUT_UNIT)
   call domain%T%print_info(logs)
   call domain%T%print_mask(OUTPUT_UNIT)

!-----------------------------------------------------------------------

END PROGRAM test_bathymetry
