! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!> Short wave radiation into the water column

MODULE getm_radiation

   !! Description:
   !!   < Say what this module contains >
   !!
   !! Current Code Owner: < Name of person responsible for this code >
   !!
   !! Code Description:
   !!   Language: Fortran 90.
   !!   This code is written to JULES coding standards v1.

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use memory_manager
   use logging
   use field_manager
   use getm_domain

   IMPLICIT NONE

   PRIVATE  ! Private scope by default

!  Module constants

!  Module types and variables
   type, public :: type_radiation_configuration
      character(len=256) :: f = "radiation.nc"
   end type type_radiation_configuration

   type, public :: type_radiation
      !! author: Karsten Bolding
      !! version: v0.1
      !!
      !! Temperature type

      TYPE(type_radiation_configuration) :: config
      class(type_logging), pointer :: logs => null()
      class(type_field_manager), pointer :: fm => null()
      class(type_getm_domain), pointer :: domain

#ifdef _STATIC_
      real(real64), dimension(A3DFIELD) :: rad = 10._real64
#else
      real(real64), dimension(:,:), allocatable :: A, g1, g2
      real(real64), dimension(:,:,:), allocatable :: rad
#endif

      contains

      procedure :: configuration => radiation_configuration
      procedure :: initialize => radiation_initialize
      procedure :: calculate => radiation_calculate

   end type type_radiation

CONTAINS

!---------------------------------------------------------------------------

SUBROUTINE radiation_configuration(self,logs,fm)

   !! Configure the the radiation module

   IMPLICIT NONE

!  Subroutine arguments
   class(type_radiation), intent(out) :: self
   class(type_logging), intent(in), target, optional :: logs
   class(type_field_manager), intent(in), target, optional :: fm

!  Local constants

!  Local variables
!---------------------------------------------------------------------------
   if (present(logs)) then
      self%logs => logs
      call self%logs%info('radiation_configuration()',level=2)
      call self%logs%info('reading initial radiation from: ',level=3,msg2=trim(self%config%f))
   end if
   if (present(fm)) then
      self%fm => fm
   end if
END SUBROUTINE radiation_configuration

!---------------------------------------------------------------------------

SUBROUTINE radiation_initialize(self,domain)
   !! Initialize the radiation field

   IMPLICIT NONE

!  Subroutine arguments
   class(type_radiation), intent(inout) :: self
   class(type_getm_domain), intent(in), target :: domain

!  Local constants

!  Local variables
   integer :: stat
   type (type_field), pointer :: f
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('radiation_initialize()',level=2)
   self%domain => domain
   TGrid: associate( TG => self%domain%T )
#ifndef _STATIC_
   call mm_s('A',self%A,TG%l(1:2),TG%u(1:2),def=0.7_real64,stat=stat)
   call mm_s('g1',self%g1,TG%l(1:2),TG%u(1:2),def=0.4_real64,stat=stat)
   call mm_s('g2',self%g2,TG%l(1:2),TG%u(1:2),def=8._real64,stat=stat)
   call mm_s('rad',self%rad,TG%l+(/0,0,-1/),TG%u,def=0._real64,stat=stat)
#endif
   if (associated(self%fm)) then
      call self%fm%register('rad', 'W/m2', 'short wave radiation', &
                            standard_name='', &
                            dimensions=(TG%dim_3d_ids), &
   !KB                         output_level=output_level_debug, &
                            part_of_state=.false., &
                            category='airsea', field=f)
      call self%fm%send_data('rad', self%rad(TG%imin:TG%imax,TG%jmin:TG%jmax,TG%kmin:TG%kmax))
   end if
   end associate TGrid
END SUBROUTINE radiation_initialize

!---------------------------------------------------------------------------

!> Write description of algorithm for calculation radiation based on surface
!> short-wave radiation

SUBROUTINE radiation_calculate(self,swr,albedo)
   !!

   IMPLICIT NONE

!  Subroutine arguments
   class(type_radiation), intent(inout) :: self
#define _T2_ self%domain%T%l(1):,self%domain%T%l(2):
   real(real64), intent(in) :: swr(_T2_)
   real(real64), intent(in) :: albedo(_T2_)
#undef _T2_

!  Local constants

!  Local variables
   integer :: i,j,k
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('radiation_calculate()',level=2)

   TGrid: associate( TG => self%domain%T )
   self%rad(:,:,TG%u(3)) = (1._real64-albedo(:,:))*swr(:,:)
   do k=TG%u(3)-1,TG%l(3)-1,-1
      do j=TG%l(2),TG%u(2)
         do i=TG%l(1),TG%u(1)
            if (TG%mask(i,j) > 0) then
               self%rad(i,j,k) = self%rad(i,j,k+1)*( &
                                 self%A(i,j) *exp(-TG%hn(i,j,k+1)/self%g1(i,j)) &
                     +(1._real64-self%A(i,j))*exp(-TG%hn(i,j,k+1)/self%g2(i,j)))
            end  if
         end do
      end do
   end do
   end associate TGrid
END SUBROUTINE radiation_calculate

!---------------------------------------------------------------------------

END MODULE getm_radiation

#if 0
            do k=0,kmax
               rad1d(k)=rad(i,j,k)
            end do
!           allow for reflection of SWR from bottom.
            if (swr_bot_refl_frac .gt. _ZERO_ .and. &
                rad1d(0)/rad1d(kmax) .gt. swr_min_bot_frac) then
               swr_refl=rad1d(0)*swr_bot_refl_frac
#if 0
            if (D(i,j) .lt. 0.5) then
               STDERR 'SWR ',i,j,swr_loc,swr_refl
            end if
#endif
               rad1d(0)=rad1d(0)-swr_refl
               zz = _ZERO_
               do k=1,kmax
                  zz=zz+hn(i,j,k)
                  rad1d(k)=rad(i,j,k) - swr_refl &
                         *(A(i,j)*exp(-zz/g1(i,j))+(1-A(i,j))*exp(-zz/g2(i,j)))
               end do
            end if
            do k=0,kmax
               rad1d(k)=rad1d(k)*rho0_cpi                ! note this
            end do

#endif
