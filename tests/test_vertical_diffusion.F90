! Copyright (C) 2020 Bolding & Bruggeman

PROGRAM test_vertical_diffusion
   !! Testing calculation of depth at U and V points

   USE :: ISO_FORTRAN_ENV
   use datetime_module
   use field_manager
   use getm_output
   use getm_operators, only: type_vertical_diffusion
   IMPLICIT NONE

!  Local constants
   integer, parameter :: imin=1, imax=100, jmin=1, jmax=100, kmin=1, kmax=100
   real(real64), parameter :: tmax=10._real64
   integer, parameter :: Nmax=100
   integer, parameter :: nsave=1
   real(real64), parameter :: cmax=5._real64, a=5._real64, Kd=10._real64

!  Local variables
   type(type_vertical_diffusion) :: vertical_diffusion
   type(type_field_manager) :: fm
   type(type_getm_output) :: output
   integer :: mask(imin:imax,jmin:jmax)
   real(real64) :: z(imin:imax,jmin:jmax,kmin:kmax)
   real(real64) :: dz(imin:imax,jmin:jmax,kmin:kmax)
   real(real64) :: nuh(imin:imax,jmin:jmax,kmin:kmax)
   real(real64) :: var(imin:imax,jmin:jmax,kmin:kmax)
   real(real64) :: ana(imin:imax,jmin:jmax,kmin:kmax)

   type(datetime) :: t
   type(timedelta) :: dt
   real(real64) :: Nsecs, timestep, cnpar, avmol, z0
   integer :: i,j,k,n
!-----------------------------------------------------------------------------

   call parse_argument()
   call field_manager_setup()
   call initial_conditions(1)
   call output%do_output(t)

   mask = 1
   timestep=0.1_real64
   Nsecs=Nmax*timestep
   cnpar=0.5_real64
   avmol=0._real64
   nuh = Kd

   dt= timedelta(seconds=int(timestep),milliseconds=int(1000*(timestep-int(timestep))))

   call vertical_diffusion%initialize(var)

!   write(*,*) lbound(var,1), ubound(var,1)                                                                                                                                          
!   write(*,*) lbound(var,2), ubound(var,2)                                                                                                                                          
!   write(*,*) lbound(var,3), ubound(var,3)                                                                                                                                          

   do n=1,Nmax
      t=t+dt
      write(*,*) n,' of ',Nmax
      call vertical_diffusion%calculate(mask,dz,timestep,cnpar,avmol,nuh,var)

      analytical: block
      real(real64) :: A1, B1, C1
      C1 = (4*Kd*n*timestep+a**2)
      A1 = cmax*sqrt(a**2/C1)
      do k=1,kmax
         B1 = exp(-(z(1,1,k)-z0)**2/C1) 
         ana(:,:,k) = A1*B1
      end do
      end block analytical

      if (mod(n,nsave) == 0) call output%do_output(t)
   end do

   write(*,*)
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
!KB      call fm%list()
      call output%initialize(fm)
   end subroutine field_manager_setup

!-----------------------------------------------------------------------------

!   subroutine field_manager_setup()
!   end subroutine field_manager_setup

!-----------------------------------------------------------------------------

   subroutine initial_conditions(method)
      integer, intent(in) :: method
      integer :: k
      ! initial condition
      select case(method)
         case (1)
            z(:,:,0)= -kmax*1._real64
            dz(:,:,:) = 1._real64
            var(:,:,0)= 1._real64
            do k=1,kmax
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

