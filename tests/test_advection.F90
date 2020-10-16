! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!! [[advection_kb.F90]]
!! [[advection.F90.template]]

PROGRAM test_advection
   !! Testing advection in a divergence free solenoidal flow field

!KB   USE, INTRINSIC :: ISO_FORTRAN_ENV
   USE :: ISO_FORTRAN_ENV
   use datetime_module
   use memory_manager
   use field_manager
   use getm_domain, only: type_getm_domain
   use getm_operators, only: type_advection
   use getm_output
   IMPLICIT NONE

!  Local constants
   real(real64), parameter :: Lx=100._real64, Ly=100._real64
   real(real64), parameter :: pi=3.1415926535897932384626433832795_real64
   integer, parameter :: imin=1, imax=100, jmin=1, jmax=100, kmin=0, kmax=25
   integer, parameter :: halo=1
   integer, parameter :: nsave=1

!  Local variables
   type(type_getm_domain) :: domain
   type(type_advection) :: advection
   TYPE(type_field_manager) :: fm
   TYPE(type_getm_output) :: output
   integer :: initial_method=1
   real(real64) :: period=600.
   real(real64) :: cfl=1.0_real64
   real(real64) :: omega
   real(real64) :: umax
   real(real64) :: tmax
   real(real64) :: dt_cfl
   real(real64) :: u(imin-halo:imax+halo,jmin-halo:jmax+halo), v(imin-halo:imax+halo,jmin-halo:jmax+halo)
   real(real64) :: var(imin-halo:imax+halo,jmin-halo:jmax+halo)
   real(real64) :: x0=Lx/2, y0=Ly/2
   real(real64) :: dx=Lx/(imax-imin+1), dy=Ly/(jmax-jmin+1)
   integer :: scheme
   integer :: no_of_revolutions=5
   integer :: stat
   integer :: Nmax
   type(datetime) :: t
   type(timedelta) :: dt
   real(real64) :: timestep
   integer :: i,j,k,n
   character(len=4) :: arg
   real(real64) :: advection_time=0._real64
!-----------------------------------------------------------------------------
   omega=2*pi/period
   call parse_argument()
   call domain_setup()
   call field_manager_setup()
   call velocity_field()
   call initial_conditions(3)
   call output%do_output(t)

   umax=omega*Lx/2
   dt_cfl=cfl*min(dx,dy)/umax
   Nmax=no_of_revolutions*nint(2*pi/omega/dt_cfl)
   tmax=no_of_revolutions*2*pi/omega
   timestep = tmax/Nmax
   dt= timedelta(seconds=int(tmax/Nmax),milliseconds=int(1000*(timestep-int(tmax/Nmax))))
   do n=1,Nmax
      t=t+dt
      write(*,*) n,' of',Nmax
      advection_: block
      real(real64) :: advection_start,advection_stop
      call cpu_time(advection_start)
      call advection%calculate(scheme,domain%U,u,domain%V,v,timestep,domain%T,var)
      call cpu_time(advection_stop)
      advection_time=advection_time+advection_stop-advection_start
      end block advection_
      if (mod(n,nsave) == 0) call output%do_output(t)
   end do
   call domain%cleanup()

   write(*,*) 'advection scheme: ',scheme
   write(*,*) 'Lx, ly:           ',Lx,Ly
   write(*,*) 'dx, dy:           ',dx,dy
   write(*,*) 'omega:            ',omega
   write(*,*) 'max velocity:     ',umax
   write(*,*) 'timestep:         ',dt_cfl
   write(*,*) '# of revolutions: ',no_of_revolutions
   write(*,*) '# of timesteps:   ',Nmax
   write(*,*) 'CPU time:         ',advection_time
   write(*,*) 'CPU time/timestep:',advection_time/Nmax
   write(*,*)
   write(*,*) 'Compiler: ',compiler_version()

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   subroutine parse_argument()
      if (command_argument_count() .lt. 1 ) then
         write(*,*)' test_advection [1-6] <N> | -h'
         stop 0
      end if
      call get_command_argument(1,arg)
      if (trim(arg) == '-h') then
         write(*,*)' test_advection [1-6] <N> | -h'
         write(*,*)'   1 -> HSIMT'
         write(*,*)'   2 -> MUSCL'
         write(*,*)'   3 -> P2_PDM'
         write(*,*)'   4 -> SPLMAX13'
         write(*,*)'   5 -> SUPERBEE'
         write(*,*)'   6 -> UPSTREAM'
         write(*,*)'   N -> number (integer) of revolutions - optional - default=5'
         write(*,*) 'Compiler: ',compiler_version()
         stop 0
      else
         read(arg,*,iostat=stat) scheme
         if (scheme < 1 .or. scheme > 6) then
            write(*,*) 'Error: argument must be in [1-6] - provided',scheme
            stop 1
         end if
         if (command_argument_count() .eq. 2 ) then
            call get_command_argument(2,arg)
            read(arg,*,iostat=stat) no_of_revolutions
         end if
      end if
   end subroutine parse_argument

!-----------------------------------------------------------------------------

   subroutine domain_setup()
      call domain%T%configure(imin=imin,imax=imax,jmin=jmin,jmax=jmax,kmin=kmin,kmax=kmax,halo=(/halo,halo,0/))
      call mm_s('c1',domain%T%c1,domain%T%l(1),domain%T%u(1),stat=stat)
      call mm_s('c2',domain%T%c2,domain%T%l(2),domain%T%u(2),stat=stat)
      call mm_s('H',domain%T%H,domain%T%l(1:2),domain%T%u(1:2),def=-99._real64,stat=stat)
      call mm_s('mask',domain%T%mask,domain%T%l(1:2),domain%T%u(1:2),def=0,stat=stat)
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
                     var(i,j) = 1._real64+4._real64*exp( -0.010*((x-x0)**2+(y-y0)**2))
                  end if
               end do
            end do
            var(imin+65:imax-15,jmin+40:jmax-40) = 5._real64
            var(imin+72:imax-15,jmin+47:jmax-47) = 1._real64
      end select
   end subroutine initial_conditions

END PROGRAM test_advection
