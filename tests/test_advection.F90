! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!! [[advection_kb.F90]]
!! [[advection.F90.template]]

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
   integer, parameter :: Nmax=1260
   integer, parameter :: nsave=1

!  Local variables
   type(type_getm_domain) :: domain
   type(type_advection) :: advection
   TYPE(type_field_manager) :: fm
   TYPE(type_getm_output) :: output
   integer :: initial_method=1
   real(real64) :: omega=0.01
   real(real64) :: u(imin-1:imax+1,jmin-1:jmax+1), v(imin-1:imax+1,jmin-1:jmax+1)
   real(real64) :: var(imin-1:imax+1,jmin-1:jmax+1)
   real(real64) :: x0=Lx/2, y0=Ly/2
   real(real64) :: dx=Lx/(imax-imin), dy=Ly/(jmax-jmin)
   integer :: scheme=1
   integer :: stat
   type(datetime) :: t
   type(timedelta) :: dt
   real(real64) :: timestep
   integer :: i,j,k,n
!-----------------------------------------------------------------------------
   call domain_setup()
   call field_manager_setup()
   call velocity_field()
   call initial_conditions(3)

   call output%do_output(t)
!KB   call advection%initialize(var)

!KB - time needs a fix
   timestep = tmax/Nmax
   timestep = 1._real64
   n = 1*(nint(timestep-int(tmax/Nmax)))
   dt= timedelta(seconds=int(tmax/Nmax),milliseconds=100)
   dt= timedelta(seconds=1,milliseconds=0)
   do n=1,Nmax
!KB      if (n == Nmax) call initial_conditions(2)
      t=t+dt
      write(*,*) n,' of',Nmax
      call advection%calculate(6,domain%U,u,domain%V,v,timestep,domain%T,var)
      if (mod(n,nsave) == 0) call output%do_output(t)
   end do

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   subroutine domain_setup()
      call domain%T%configure(imin=imin,imax=imax,jmin=jmin,jmax=jmax,kmin=kmin,kmax=kmax,halo=(/1,1,0/))
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
      call fm%register('dx','m','grid spacing - x',dimensions=(/id_dim_lon,id_dim_lat/),no_default_dimensions=.true.,data2d=domain%T%dx(imin:imax,jmin:jmax))
      call fm%register('dy','m','grid spacing - y',dimensions=(/id_dim_lon,id_dim_lat/),no_default_dimensions=.true.,data2d=domain%T%dy(imin:imax,jmin:jmax))
      call fm%register('H','m','depth',dimensions=(/id_dim_lon,id_dim_lat/),no_default_dimensions=.true.,data2d=domain%T%H(imin:imax,jmin:jmax))
      call fm%register('u','m/s','u-velocity',dimensions=(/id_dim_lon,id_dim_lat/),no_default_dimensions=.true.,data2d=u(imin:imax,jmin:jmax),fill_value=0._real64)
      call fm%register('v','m/s','v-velocity',dimensions=(/id_dim_lon,id_dim_lat/),no_default_dimensions=.true.,data2d=v(imin:imax,jmin:jmax),fill_value=0._real64)
      call fm%register('D','','depth',data2d=domain%T%D(imin:imax,jmin:jmax),fill_value=-99._real64)
      call fm%register('f','','scalar',data2d=var(imin:imax,jmin:jmax),fill_value=-99._real64)
!KB      call fm%list()
      call output%initialize(fm)
   end subroutine field_manager_setup

!-----------------------------------------------------------------------------

   subroutine velocity_field()
      ! divergence free velocity field
      u = 0._real64
      v = 0._real64
      do j=jmin,jmax
         do i=imin,imax
            if (domain%U%mask(i,j) > 0) then
               u(i,j) = -omega*domain%U%c2(j)
            end if
            if (domain%V%mask(i,j) > 0) then
               v(i,j) =  omega*domain%V%c1(i)
            end if
         end do
      end do
   end subroutine velocity_field

!-----------------------------------------------------------------------------

   subroutine initial_conditions(method)
      integer, intent(in) :: method

      real(real64) :: x,y,x0=-25._real64,y0=0._real64
      ! initial condition
      var = -99._real64
      select case(method)
         case (1)
            where (domain%T%mask > 0)
               var(:,:)= 1._real64
            end where
            var(imin+55:imax-25,jmin+40:imax-40) = 5._real64
         case (2)
            do j=jmin,jmax
               do i=imin,imax
                  if (domain%T%mask(i,j) > 0) then
                     x=domain%T%c1(i)
                     y=domain%T%c2(j)
                     var(i,j) = 1._real64+4._real64*exp( -0.025*((x-x0)**2+(y-y0)**2))
                  end if
               end do
            end do
         case (3)
            do j=jmin,jmax
               do i=imin,imax
                  if (domain%T%mask(i,j) > 0) then
                     x=domain%T%c1(i)
                     y=domain%T%c2(j)
                     var(i,j) = 1._real64+4._real64*exp( -0.025*((x-x0)**2+(y-y0)**2))
                  end if
               end do
            end do
            var(imin+65:imax-15,jmin+40:imax-40) = 5._real64
      end select
   end subroutine initial_conditions

END PROGRAM test_advection
