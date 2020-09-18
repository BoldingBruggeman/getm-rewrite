! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

PROGRAM test_advection
   !! Testing advection in a divergence free solenoidal flow field

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use datetime_module
   use memory_manager
   use field_manager
   use getm_domain, only: type_getm_domain
   use getm_operators, only: type_advection
   use getm_output
   IMPLICIT NONE

!  Local constants
   real(real64), parameter :: Lx=100._real64, Ly=100._real64
   integer, parameter :: imin=1, imax=101, jmin=1, jmax=101, kmin=0, kmax=25
   real(real64), parameter :: tmax=10._real64
   integer, parameter :: Nmax=100

!  Local variables
   type(type_getm_domain) :: domain
   type(type_advection) :: advection
   TYPE(type_field_manager) :: fm
   TYPE(type_getm_output) :: output
   integer :: initial_method=1
   real(real64) :: omega=0.01
#if 0
   real(real64) :: u(imin:imax,jmin:jmax,kmin:kmax), v(imin:imax,jmin:jmax,kmin:kmax)
   real(real64) :: var(imin:imax,jmin:jmax,kmin:kmax)
#else
   real(real64) :: u(imin:imax,jmin:jmax), v(imin:imax,jmin:jmax)
   real(real64) :: var(imin:imax,jmin:jmax)
#endif
   integer :: i,j,k,n
   real(real64) :: x0=Lx/2, y0=Ly/2
   real(real64) :: dx=Lx/(imax-imin), dy=Ly/(jmax-jmin)
   integer :: scheme=1
   integer :: stat
   type(datetime) :: t
   type(timedelta) :: dt
   real(real64) :: timestep
!-----------------------------------------------------------------------------
   call domain_setup()
   call field_manager_setup()
   call velocity_field()
   call initial_conditions()

   call output%do_output(t)
!KB   call advection%initialize(var)

!KB - time needs a fix
   timestep = tmax/Nmax
   n = 1*(nint(timestep-int(tmax/Nmax)))
   dt= timedelta(seconds=int(tmax/Nmax),milliseconds=100)
!KB   do n=1,Nmax
   do n=1,25
      t=t+dt
      write(*,*) n,' of',Nmax
      call advection%calculate(1,domain%U,u,domain%V,v,timestep,domain%T,var)
      call output%do_output(t)
   end do

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   subroutine domain_setup()
      call domain%T%configure(imin=imin,imax=imax,jmin=jmin,jmax=jmax,kmin=kmin,kmax=kmax)
      call mm_s('c1',domain%T%c1,domain%T%l(1),domain%T%u(1),stat=stat)
      call mm_s('c2',domain%T%c2,domain%T%l(2),domain%T%u(2),stat=stat)
      call mm_s('H',domain%T%H,domain%T%l(1:2),domain%T%u(1:2),stat=stat)
      call mm_s('mask',domain%T%mask,domain%T%l(1:2),domain%T%u(1:2),stat=stat)
write(*,*) x0,y0
      do i=imin,imax
         domain%T%c1(i) = (i-1)*dx-x0
      end do
      do j=jmin,jmax
         domain%T%c2(j) = (j-1)*dy-y0
      end do
      domain%T%H(imin+1:imax-1,jmin+1:jmax-1) = 1._real64
      where(domain%T%H .gt. 0._real64) domain%T%mask = 1
      domain%domain_type=1
      call domain%configure()
      call domain%initialize()
      call domain%report()
   end subroutine domain_setup

!-----------------------------------------------------------------------------

   subroutine field_manager_setup()
      call fm%register_dimension('x',imax-imin+1,id=id_dim_lon)
      call fm%register_dimension('y',jmax-jmin+1,id=id_dim_lat)
      call fm%register_dimension('time',id=id_dim_time)
      call fm%initialize(prepend_by_default=(/id_dim_lon,id_dim_lat/),append_by_default=(/id_dim_time/))
      call fm%register('u','m/s','u-velocity',dimensions=(/id_dim_lon,id_dim_lat/),no_default_dimensions=.true.,data2d=u)
      call fm%register('v','m/s','v-velocity',dimensions=(/id_dim_lon,id_dim_lat/),no_default_dimensions=.true.,data2d=v)
      call fm%register('f','','scalar',data2d=var)
!KB      call fm%list()
      call output%initialize(fm)
   end subroutine field_manager_setup

!-----------------------------------------------------------------------------

   subroutine velocity_field()
      ! divergence free velocity field
      do j=jmin,jmax
         do i=imin,imax
            u(i,j) = -omega*domain%U%c2(j)
            v(i,j) =  omega*domain%V%c1(i)
         end do
      end do
   end subroutine velocity_field

!-----------------------------------------------------------------------------

   subroutine initial_conditions()
      ! initial condition
      select case(initial_method)
         case (1)
#if 0
            var(:,:,:)= 1._real64
            var(imin+40:imax-40,jmin+40:imax-40,:) = 5._real64
#else
            var(:,:)= 1._real64
            var(imin+40:imax-40,jmin+40:imax-40) = 5._real64
#endif
      end select
   end subroutine initial_conditions

END PROGRAM test_advection
