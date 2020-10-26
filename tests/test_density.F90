! Copyright (C) 2020 Bolding & Bruggeman

!> @note
!> test performance of where and explicit loops
!> @endnote

PROGRAM test_density
   !! Testing calculation of time varying depths at S, U and V points

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use gsw_mod_toolbox, only: gsw_rho
   use logging
   use memory_manager
   use getm_domain
   use getm_temperature
   use getm_salinity
   use getm_density
   IMPLICIT NONE

!  Local constants
   integer, parameter :: rk = real64
   integer, parameter :: imin=1, imax=50, jmin=1, jmax=50, kmin=0, kmax=20
   integer, parameter :: nmax = 500
   real(real64) :: t1, t2, t3, t4

!  Local variables
   TYPE(type_logging) :: logs
   TYPE(type_getm_domain), target :: domain
   TYPE(type_temperature) :: temperature
   TYPE(type_salinity) :: salinity
   TYPE(type_density) :: density

   integer :: out_unit = 100
   character(len=*), parameter :: file_name = 'test_density.log'
   character(len=256) :: msg

   integer :: n,stat
   integer, dimension(3) :: l,u
   character(len=16) :: arg
!-----------------------------------------------------------------------
   write(*,*) command_argument_count()

!   if (command_argument_count() .ne. 0 .or. command_argument_count() .ne. 6 ) then
   if (command_argument_count() .ne. 6 ) then
      write(*,*) 'ERROR, must be called like:'
      write(*,*) 'test_density'
      STOP ' test_density imin imax jmin jmax kmin kmax'
   end if

   do n=1,command_argument_count()
      call get_command_argument(n,arg)
      select case (n)
         case(1,3,5)
            read(arg,*,iostat=stat)  domain%T%l((n+1)/2)
         case(2,4,6)
            read(arg,*,iostat=stat)  domain%T%u((n+1)/2)
      end select
   end do

   logs%prepend = ''
   call logs%info('testing density():')
   open(file = trim(file_name), unit = out_unit)
   call logs%costum(out_unit,'testing density():')
   write(out_unit,*) 'nmax:   ',nmax
   write(out_unit,*) 'lbound: ',domain%T%l
   write(out_unit,*) 'ubound: ',domain%T%u

   call domain%configure(imin,imax,jmin,jmax,kmin,kmax,logs=logs)
   domain%T%mask(imin:imax,jmin:jmax) = 1

   call temperature%configuration(logs)
   call temperature%initialize(domain)

   call salinity%configuration(logs)
   call salinity%initialize(domain)

   call density%configure(logs)
   call density%initialize(domain)

   call cpu_time(t1)
   do n=1,nmax
      call density_mask()
   end do
   call cpu_time(t2)
   do n=1,nmax
      call density_loops(1)
   end do
   call cpu_time(t3)
   do n=1,nmax
      call density_loops(2)
   end do
   call cpu_time(t4)

   write(out_unit,*)
   write(msg,'(A,F10.5,A,2F10.5)') 'using masks: ',t2-t1,' using loops(1,2): ',t3-t2,t4-t3
   call logs%costum(out_unit,trim(msg))

   call logs%info('results in: ', msg2=trim(file_name))

CONTAINS 

!----------------------------------------------------------------------

SUBROUTINE density_mask()
   !! Feeds your cats and dogs, if enough food is available. If not enough
   !! food is available, some of your pets will get angry.

   IMPLICIT NONE

!  Subroutine arguments
      !! The number of cats to keep track of.

!  Local constants

!  Local variables
   integer :: i,j
!---------------------------------------------------------------------------
   call logs%info('density_mask()',level=2)

#if 0
   where (mask3d > 1) 
      density%rho = gsw_rho(salinity%S,temperature%T,0._rk)
   end where
#else
   do j=domain%T%l(2),domain%T%u(2)
      do i=domain%T%l(1),domain%T%u(1)
         if (domain%T%mask(i,j) > 0) then
            density%rho(i,j,:) = gsw_rho(salinity%S(i,j,:), &
                                         temperature%T(i,j,:),0._rk)
         end  if
      end do
   end do
#endif
END SUBROUTINE density_mask

!---------------------------------------------------------------------------

SUBROUTINE density_loops(method)
   !! Feeds your cats and dogs, if enough food is available. If not enough
   !! food is available, some of your pets will get angry.

   IMPLICIT NONE

!  Subroutine arguments
   integer, intent(in) :: method
      !! The number of cats to keep track of.

!  Local constants

!  Local variables
   integer :: i,j,k
!---------------------------------------------------------------------------
   call logs%info('density_loops()',level=2)

   if (method == 1) then
      do k=domain%T%l(3),domain%T%u(3)
         do j=domain%T%l(2),domain%T%u(2)
            do i=domain%T%l(1),domain%T%u(1)
               if (domain%T%mask(i,j) > 0) then
                  density%rho(i,j,k) = gsw_rho(salinity%S(i,j,k), &
                                               temperature%T(i,j,k),0._rk)
               end  if
            end do
         end do
      end do
   else
      do j=domain%T%l(2),domain%T%u(2)
         do i=domain%T%l(1),domain%T%u(1)
            if (domain%T%mask(i,j) > 0) then
               do k=domain%T%l(3),domain%T%u(3)
                  density%rho(i,j,k) = gsw_rho(salinity%S(i,j,k), &
                                               temperature%T(i,j,k),0._rk)
               end do
            end  if
         end do
      end do
   end if
END SUBROUTINE density_loops

!---------------------------------------------------------------------------

END PROGRAM test_density
