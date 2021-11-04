! Copyright (C) 2020 Bolding & Bruggeman

!> @note
!> maybe put all grid tests in one big program - test_bathymetry
!> @endnote

PROGRAM test_vertical_coordinates
   !! Testing calculation of time varying depths at S, U and V points

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use logging
   use memory_manager
   use getm_domain
   use getm_sealevel
   IMPLICIT NONE

!  Local variables
   TYPE(type_logging) :: logs
   class(type_getm_domain), allocatable, target :: domain
   class(type_getm_sealevel), allocatable, target :: sealevel

!   character(len=256) :: varname
!   integer :: ndims, dimlens
   integer :: imin=1,imax=201,jmin=1,jmax=101,kmax=30
   real(real64) :: dt=10._real64
   integer :: i,j,n
!-----------------------------------------------------------------------

!   if (command_argument_count() .ne. 1 ) then
!      write(*,*) 'ERROR, must be called like:'
!      STOP ' test_bathymetry <name of bathymetry file>'
!   end if
!   call get_command_argument(1,bathymetry%depth%f)
!   bathymetry%depth%v = 'bathymetry'

   allocate(domain)
   call domain%configure(imin,imax,jmin,jmax,kmin=0,kmax=kmax,logs=logs)
   do j=jmin,jmax
      do i=imin,imax
         domain%T%H(i,j) = 5._real64 + ((i-1)+(j-1))*0.1_real64
      end do
   end do
   domain%T%mask=1 !KB
   call domain%uvx_depths()

#if 0
   write(0,*) domain%T%halo
   write(0,*) size(domain%T%H)
   write(0,*) domain%T%H(1,1),domain%T%H(imax,jmax)
   write(0,*) domain%U%H(1,1),domain%U%H(imax,jmax)
   write(0,*) domain%V%H(1,1),domain%V%H(imax,jmax)
#endif

   call domain%init_vertical()
   call domain%do_vertical(dt)
   write(100,*) sum(domain%T%ho(1,1,1:)),sum(domain%T%hn(1,1,1:))
   write(200,*) sum(domain%T%ho(imax,jmax,1:)),sum(domain%T%hn(imax,jmax,1:))

   allocate(sealevel)
   call sealevel%configure(logs=logs)
   call sealevel%initialize(domain)
   do n=1,100
      domain%T%z = n*0.01_real64
      call sealevel%uvx()
      if (mod(n,10) == 0) then
         call domain%start_3d()
         call domain%do_vertical(dt)

         write(100,*) sum(domain%T%ho(1,1,1:)),sum(domain%T%hn(1,1,1:)),domain%T%H(1,1)+domain%T%z(i,j)
         write(200,*) sum(domain%T%ho(imax,jmax,1:)),sum(domain%T%hn(imax,jmax,1:)),domain%T%H(imax,jmax)+domain%T%z(imax,jmax)

      end if
   end do

!-----------------------------------------------------------------------

!CONTAINS 

!-----------------------------------------------------------------------

END PROGRAM test_vertical_coordinates
