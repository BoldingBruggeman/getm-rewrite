! Copyright (C) 2020 Bolding & Bruggeman

PROGRAM test_vertical_diffusion
   !! Testing the vertical diffusion solver

   USE :: ISO_FORTRAN_ENV
   use datetime_module
   use field_manager
   use getm_output
   use getm_operators, only: type_vertical_diffusion
   IMPLICIT NONE

!  Local constants
   integer, parameter :: imin=1, imax=100, jmin=1, jmax=100
   real(real64), parameter :: depth=100._real64
   real(real64), parameter :: Nsecs=1000._real64
   integer, parameter :: nsave=1
   real(real64), parameter :: cmax=5._real64, a=5._real64, Av=0.1_real64

!  Local variables
   type(type_vertical_diffusion) :: vertical_diffusion
   type(type_field_manager) :: fm
   type(type_getm_output) :: output
   integer :: Nmax=100
   integer :: kmin=1, kmax=100
   integer, allocatable :: mask(:,:)
   real(real64), allocatable :: z(:,:,:)
   real(real64), allocatable :: dz(:,:,:)
   real(real64), allocatable :: nuh(:,:,:)
   real(real64), allocatable :: var(:,:,:)
   real(real64), allocatable :: ana(:,:,:)
   real(real64), allocatable :: diff(:,:,:)
   real(real64) :: conservation

   type(datetime) :: t
   type(timedelta) :: dt
   real(real64) :: timestep=5._real64, cnpar=0.5, avmol, z0
   integer :: i,j,k,n
!-----------------------------------------------------------------------------

   call parse_argument()
   allocate(mask(imin:imax,jmin:jmax))
   mask = 1
   allocate(z(imin:imax,jmin:jmax,kmin:kmax))
   allocate(dz(imin:imax,jmin:jmax,kmin:kmax))
   allocate(nuh(imin:imax,jmin:jmax,kmin:kmax))
   allocate(var(imin:imax,jmin:jmax,kmin:kmax))
   allocate(ana(imin:imax,jmin:jmax,kmin:kmax))
   allocate(diff(imin:imax,jmin:jmax,kmin:kmax))
   call field_manager_setup()
   call initial_conditions(1)
   conservation=sum(var)
   call output%do_output(t)

   Nmax=Nsecs/timestep
   avmol=0._real64
   nuh = Av

   dt= timedelta(seconds=int(timestep),milliseconds=int(1000*(timestep-int(timestep))))

   call vertical_diffusion%initialize(var)

   do n=1,Nmax
      t=t+dt
      write(*,*) n,' of ',Nmax
      call vertical_diffusion%calculate(timestep,cnpar,mask,dz,dz,avmol,nuh,var)
      conservation=sum(var)

      analytical: block
      real(real64) :: A1, B1, C1
      C1 = (4*Av*n*timestep+a**2)
      A1 = cmax*sqrt(a**2/C1)
      do k=1,kmax
         B1 = exp(-(z(1,1,k)-z0)**2/C1) 
         ana(:,:,k) = A1*B1
      end do
      diff=var-ana
      end block analytical

      if (mod(n,nsave) == 0) call output%do_output(t)
   end do

   write(*,*)
   write(*,*) 'time: length=',Nsecs,Nmax,timestep
   write(*,*) 'domain: depth=',depth,'dz=',dz(imax/2,jmax/2,kmax/2)
   write(*,*) 'cnpar: ',cnpar
   write(*,*) 'domain size (i,j,k): ',imax-imin+1,jmax-jmin+1,kmax-kmin+1
   write(*,*) 'timing (pr timestep)'
   write(*,*) 'matrix setup:        ',vertical_diffusion%matrix_time/Nmax
   write(*,*) 'tri-diagonal solver: ',vertical_diffusion%tridiag_time/Nmax
   write(*,*)
   write(*,*) 'Compiler: ',compiler_version()

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   subroutine parse_argument()
      integer :: stat
      character(len=16) :: arg
      if (command_argument_count() .eq. 1 ) then
         call get_command_argument(1,arg)
         if (trim(arg) == '-h') then
            write(*,*)' test_vertical_diffusion kmax (100) timestep (5.0) cnpar (0.5)'
            stop 0
         else
            write(*,*)' test_vertical_diffusion -h'
            stop 0
         end if
      end if
      if (command_argument_count() .eq. 3 ) then
         call get_command_argument(1,arg)
         read(arg,*,iostat=stat) kmax
         call get_command_argument(2,arg)
         read(arg,*,iostat=stat) timestep
         call get_command_argument(3,arg)
         read(arg,*,iostat=stat) cnpar
      end if
   end subroutine parse_argument

!-----------------------------------------------------------------------------


   subroutine field_manager_setup()
      integer, parameter :: i=imax/2,j=jmax/2
      call fm%register_dimension('x',imax-imin+1,id=id_dim_lon)
      call fm%register_dimension('y',jmax-jmin+1,id=id_dim_lat)
      call fm%register_dimension('z',kmax-kmin+1,id=id_dim_z)
      call fm%register_dimension('time',id=id_dim_time)
      call fm%initialize(prepend_by_default=(/id_dim_lon,id_dim_lat,id_dim_z/),append_by_default=(/id_dim_time/))
      call fm%register('nuh','m2/s','difusion',dimensions=(/id_dim_lon,id_dim_lat,id_dim_z/),no_default_dimensions=.true.,data3d=nuh)
      call fm%register('f','-','scalar',data3d=var,fill_value=-99._real64)
      call fm%register('a','-','analytical',data3d=ana,fill_value=-99._real64)
      call fm%register('d','-','difference',data3d=diff,fill_value=-99._real64)
      call fm%register('c','-','conservation',no_default_dimensions=.true.,dimensions=(/id_dim_time/),data0d=conservation,fill_value=-99._real64)
!KB      call fm%list()
      call output%configure(fm=fm)
      call output%initialize()
   end subroutine field_manager_setup

!-----------------------------------------------------------------------------

!   subroutine field_manager_setup()
!   end subroutine field_manager_setup

!-----------------------------------------------------------------------------

   subroutine initial_conditions(method)
      integer, intent(in) :: method
      integer :: k
      dz(:,:,:)=depth/(kmax-kmin+1)
      ! initial condition
      select case(method)
         case (1)
            z(:,:,1)= -depth
            var(:,:,1)= 1._real64
            do k=2,kmax
              z(:,:,k)= z(:,:,k-1)+dz(:,:,k)
            end do
            z0 = z(1,1,kmax/2)
            do k=1,kmax
               ana(:,:,k) = cmax*exp(-(z(:,:,k)-z0)**2/a**2)
            end do
            var=ana
         case (2)
            z(:,:,0)= -kmax*1._real64
            dz(:,:,:) = 1._real64
            var(:,:,0)= 1._real64
            do k=1,kmax
              z(:,:,k)= z(:,:,k-1)+dz(:,:,k)
              var(:,:,k)=tanh(-(z(i,j,k)+50._real64))
            end do
!KB            var(:,:,kmax)= -4._real64
            var(:,:,kmax)= -1._real64
      end select
   end subroutine initial_conditions

!-----------------------------------------------------------------------------

END PROGRAM test_vertical_diffusion

