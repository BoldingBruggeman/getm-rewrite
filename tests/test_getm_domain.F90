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
   TYPE(type_bathymetry) :: bathymetry
   class(type_getm_domain), allocatable, target :: cart1,cart2,cart3
   class(type_getm_domain), allocatable, target :: sph1,sph2

   character(len=256) :: varname
   integer :: ndims, dimlens
   integer :: imin=1,imax=100,jmin=1,jmax=30
   real(real64) :: dx=500._real64,dy=500._real64
   real(real64) :: dlon=0.25_real64,dlat=0.2_real64
   integer :: i,j
!-----------------------------------------------------------------------

   if (command_argument_count() .ne. 1 ) then
      write(*,*) 'ERROR, must be called like:'
      STOP ' test_bathymetry <name of bathymetry file>'
   end if
   call get_command_argument(1,bathymetry%depth%f)
   bathymetry%depth%v = 'bathymetry'

   allocate(cart1)

   call cart1%configure(imin,imax,jmin,jmax,kmin=-1,kmax=-1,dx=dx,dy=dy)
   cart%T
!KB   call cart1%kaj(dx=dx,dy=dy)

   allocate(cart2)
   call cart2%configure(imin,imax,jmin,jmax,kmin=-1,kmax=-1)
!KB   call cart2%kaj(x=cart1%T%x(:,1),y=cart1%T%y(1,:))

   write(*,*)
   write(*,*) "difference between cartesian 1 and 2 grids:"
   write(*,*) 'cart1-cart2: x  ',sum(cart1%T%x(imin:imax,jmin:jmax)-cart2%T%x(imin:imax,jmin:jmax))
   write(*,*) 'cart1-cart2: y  ',sum(cart1%T%y(imin:imax,jmin:jmax)-cart2%T%y(imin:imax,jmin:jmax))
   write(*,*) 'cart1-cart2: xi ',sum(cart1%T%xi(imin:imax,jmin:jmax)-cart2%T%xi(imin:imax,jmin:jmax))
   write(*,*) 'cart1-cart2: yi ',sum(cart1%T%yi(imin:imax,jmin:jmax)-cart2%T%yi(imin:imax,jmin:jmax))
   write(*,*) 'cart1-cart2: dx ',sum(cart1%T%dx(imin:imax,jmin:jmax)-cart2%T%dx(imin:imax,jmin:jmax))
   write(*,*) 'cart1-cart2: dy ',sum(cart1%T%dy(imin:imax,jmin:jmax)-cart2%T%dy(imin:imax,jmin:jmax))

   allocate(cart3)
   call cart3%configure(imin,imax,jmin,jmax,kmin=-1,kmax=-1)
!KB   call cart3%kaj(xi=cart2%T%xi(:,1),yi=cart2%T%yi(1,:))

   write(*,*)
   write(*,*) "difference between cartesian 1 and 3 grids:"
   write(*,*) 'cart1-cart3: x  ',sum(cart1%T%x(imin:imax,jmin:jmax)-cart3%T%x(imin:imax,jmin:jmax))
   write(*,*) 'cart1-cart3: y  ',sum(cart1%T%y(imin:imax,jmin:jmax)-cart3%T%y(imin:imax,jmin:jmax))
   write(*,*) 'cart1-cart3: xi ',sum(cart1%T%xi(imin:imax,jmin:jmax)-cart3%T%xi(imin:imax,jmin:jmax))
   write(*,*) 'cart1-cart3: yi ',sum(cart1%T%yi(imin:imax,jmin:jmax)-cart3%T%yi(imin:imax,jmin:jmax))
   write(*,*) 'cart1-cart3: dx ',sum(cart1%T%dx(imin:imax,jmin:jmax)-cart3%T%dx(imin:imax,jmin:jmax))
   write(*,*) 'cart1-cart3: dy ',sum(cart1%T%dy(imin:imax,jmin:jmax)-cart3%T%dy(imin:imax,jmin:jmax))

   allocate(sph1)
   call sph1%configure(imin,imax,jmin,jmax,kmin=-1,kmax=-1)
!KB   call sph1%kaj(dlon=dlon,dlat=dlat)

!-----------------------------------------------------------------------

!CONTAINS 

!-----------------------------------------------------------------------

END PROGRAM test_getm_domain
