! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

PROGRAM test_advection
   !! Testing advection in a divergence free solenoidal flow field

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use getm_domain, only: type_getm_domain
   use getm_operators, only: type_advection
   IMPLICIT NONE

!  Local constants
   integer, parameter :: imin=1, imax=100, jmin=1, jmax=100, kmin=0, kmax=100
   real(real64), parameter :: tmax=10._real64
   integer, parameter :: Nmax=100

!  Local variables
   type(type_advection) :: advection
   real(real64) :: dt
   real(real64) :: u(imin:imax,jmin:jmax,kmin:kmax), v(imin:imax,jmin:jmax,kmin:kmax)
   real(real64) :: var(imin:imax,jmin:jmax,kmin:kmax)

   integer :: i,j,k,n,i0,j0
   real(real64) :: omega=0.01
   real(real64) :: dx=1._real64,dy=1._real64
!-----------------------------------------------------------------------------

   i0 = (imax-imin+1)/2
   j0 = (jmax-jmin+1)/2
   do j=jmin,jmax
      do i=imin,imax
         u(i,j,:) = -omega*(j-j0)*dy
         v(i,j,:) =  omega*(i-i0)*dx
      end do
   end do

   dt=tmax/Nmax

   ! initial condition
   i=1; j=1
   var(:,:,:)= 0._real64
   do j=j0-20,j0+20
      do i=i-20,i0+20
         var(i,j,:) = 5._real64
      end do
   end do

OPEN(1, FILE='u.dat', FORM='unformatted')
WRITE(1) u
close(1)
OPEN(2, FILE='v.dat', FORM='unformatted')
WRITE(2) v
close(2)
!write(*,*) u(:,:,1)
!write(*,*) v(:,:,1)
!KBwrite(*,*) var

!KB   call advection%initialize(var)

!KB   do n=1,Nmax
!KB!   do n=1,1
!KB      call advection%calculate(domain,dt,var)
!KB         write(100,*) n,var(i,j,55)
!KB   end do

END PROGRAM test_advection
