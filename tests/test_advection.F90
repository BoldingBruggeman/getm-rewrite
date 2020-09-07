! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

PROGRAM test_advection
   !! Testing advection in a divergence free solenoidal flow field

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use logging
   use memory_manager
   use getm_domain, only: type_getm_domain
   use getm_operators, only: type_advection
   IMPLICIT NONE

!  Local constants
   real(real64), parameter :: Lx=100._real64, Ly=100._real64
   integer, parameter :: imin=1, imax=101, jmin=1, jmax=101, kmin=0, kmax=25
   real(real64), parameter :: tmax=10._real64
   integer, parameter :: Nmax=100

!  Local variables
   type(type_getm_domain) :: domain
   type(type_advection) :: advection
   integer :: initial_method=1
   real(real64) :: dt
   real(real64) :: omega=0.01
   real(real64) :: u(imin:imax,jmin:jmax,kmin:kmax), v(imin:imax,jmin:jmax,kmin:kmax)
   real(real64) :: var(imin:imax,jmin:jmax,kmin:kmax)

   integer :: i,j,k,n
   real(real64) :: x0=Lx/2, y0=Ly/2
   real(real64) :: dx=Lx/(imax-imin), dy=Ly/(jmax-jmin)
   integer :: stat
!-----------------------------------------------------------------------------
!KB   call domain%T%configure(imin=imin,imax=imax,jmin=jmin,jmax=jmax,kmin=kmin,kmax=kmax,halo=/0,0,0/)
   call domain%T%configure(imin=imin,imax=imax,jmin=jmin,jmax=jmax,kmin=kmin,kmax=kmax)
   call mm_s('c1',domain%T%c1,domain%T%l(1),domain%T%u(1),stat=stat)
   call mm_s('c2',domain%T%c2,domain%T%l(2),domain%T%u(2),stat=stat)
   call mm_s('H',domain%T%H,domain%T%l(1:2),domain%T%u(1:2),stat=stat)
   call mm_s('mask',domain%T%mask,domain%T%l(1:2),domain%T%u(1:2),stat=stat)
   do i=imin,imax
      domain%T%c1(i) = (i-1)*dx-x0
   end do
   do j=jmin,jmax
      domain%T%c2(j) = (j-1)*dy-y0
   end do
   domain%T%H(imin+1:imax-1,jmin+1:jmax-1) = 25._real64
   where(domain%T%H .gt. 0._real64) domain%T%mask = 1
   domain%domain_type=1
   call domain%configure()
   call domain%initialize()
   call domain%report()
#if 0
write(*,*) domain%T%D(imax/2,jmax/2)
write(*,*) domain%U%H(imax/2,jmax/2)
write(*,*) domain%T%c1
write(*,*) domain%T%c2
write(*,*) domain%U%c1
write(*,*) domain%U%c2
write(*,*) domain%V%c1
write(*,*) domain%V%c2
#endif

   ! divergence fre velocity field
   do j=jmin,jmax
      do i=imin,imax
         u(i,j,:) = -omega*domain%U%c2(j)
         v(i,j,:) =  omega*domain%V%c1(i)
      end do
   end do

   dt=tmax/Nmax

   ! initial condition
   select case(initial_method)
      case (1)
         var(:,:,:)= 1._real64
         var(imin+40:imax-40,jmin+40:imax-40,:) = 5._real64
   end select
write(*,*) var(:,50,1)

!KB   call advection%initialize(var)

!KB   do n=1,Nmax
!KB!   do n=1,1
!KB      call advection%calculate(domain%ugrid,u,domain%vgrid,v,domain%tgrid,var,dt)
!KK      ! This would be 2D advection: CAN YOU ADD HOW TO GET THE METRICS?
!KK      ! masks are logicals (in this test case probably au.eq.1, av.eq.1, az.eq.1)
!KK      call advection%calculate(dt, var(:,:,1), h(:,:,1), arcd1, u(:,:,1), hu(:,:,1), dyu, dxu, v(:,:,1), hv(:,:,1), dxv, dyv, mask_uflux, mask_vflux, mask_update)
!KK      ! TODO: Better would be something like ...
!KK      !call var%advection%calculate(dt)
!KK      ! TODO: ... where the scheme and all velocities, metrics and masks specific for var are already set (for velocities and layer heights pointers to updated arrays)
!KB         write(100,*) n,var(i,j,55)
!KB   end do

#if 0
OPEN(1, FILE='u.dat', FORM='unformatted')
WRITE(1) u
close(1)
OPEN(2, FILE='v.dat', FORM='unformatted')
WRITE(2) v
close(2)
!write(*,*) u(:,:,1)
!write(*,*) v(:,:,1)
!KBwrite(*,*) var
#endif

END PROGRAM test_advection
