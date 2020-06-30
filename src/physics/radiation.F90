! Copyright (C) 2020 Bolding & Bruggeman

!> Short wave radiation into the water column

#ifdef _STATIC_
#include "dimensions.h"
#endif

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

SUBROUTINE radiation_configuration(self,logs)

   !! Configure the the radiation module

   IMPLICIT NONE

!  Subroutine arguments
   class(type_radiation), intent(out) :: self
   class(type_logging), intent(in) :: logs

!  Local constants

!  Local variables
!---------------------------------------------------------------------------
   call logs%info('radiation_configuration()',level=2)
   call logs%info('reading initial radiation from: ',level=3,msg2=trim(self%config%f))

   return
END SUBROUTINE radiation_configuration

!---------------------------------------------------------------------------

SUBROUTINE radiation_initialize(self,logs,grid)
   !! Initialize the radiation field

   IMPLICIT NONE

!  Subroutine arguments
   class(type_radiation), intent(inout) :: self
   class(type_logging), intent(in) :: logs
   class(type_getm_grid), intent(in) :: grid

!  Local constants

!  Local variables
   integer :: stat
!---------------------------------------------------------------------------
   call logs%info('radiation_initialize()',level=2)
!KB   call logs%warn('radiation_initialize()',level=1)
!   call logs%costum(tunit,'testing costum message from radiation',level=2)

#ifndef _STATIC_
   call mm_s('A',self%A,grid%l(1:2),grid%u(1:2),def=0.7_real64,stat=stat)
   call mm_s('g1',self%g1,grid%l(1:2),grid%u(1:2),def=0.4_real64,stat=stat)
   call mm_s('g2',self%g2,grid%l(1:2),grid%u(1:2),def=8._real64,stat=stat)
   call mm_s('rad',self%rad,grid%l,grid%u,def=15._real64,stat=stat)
#endif

   return
END SUBROUTINE radiation_initialize

!---------------------------------------------------------------------------

!> Write description of algorithm for calculation radiation based on surface
!> short-wave radiation

SUBROUTINE radiation_calculate(self,logs,grid,swr,albedo)
   !!

   IMPLICIT NONE

!  Subroutine arguments
   class(type_radiation), intent(out) :: self
   class(type_logging), intent(in) :: logs
   class(type_getm_grid), intent(in) :: grid
   real(real64), dimension(:,:), intent(in) :: swr
   real(real64), dimension(:,:), intent(in) :: albedo

!  Local constants

!  Local variables
   real(real64) :: z
   integer :: i,j,k
!---------------------------------------------------------------------------
   call logs%info('radiation_calculate()',level=2)

   self%rad(:,:,grid%u(3)) = (1._real64-albedo(:,:))*swr(:,:)
   do k=grid%l(3),grid%u(3)-1
      do j=grid%l(2),grid%u(2)
         do i=grid%l(1),grid%u(1)
            if (grid%mask(i,j) > 0) then
               z = grid%zf(i,j,k)
               self%rad(i,j,k) = self%rad(i,j,grid%u(3)) &
                       *(self%A(i,j) &
                       *exp(-z/self%g1(i,j)) &
                       +(1._real64-self%A(i,j))*exp(-z/self%g2(i,j)))
            end  if
         end do
      end do
   end do

   return
END SUBROUTINE radiation_calculate

!---------------------------------------------------------------------------

END MODULE getm_radiation

#if 0

  do j=jmin,jmax
      do i=imin,imax
         if (az(i,j) .ge. 1) then
            swr_loc=swr(i,j)
            rad(i,j,kmax)=swr_loc
            zz = _ZERO_
            do k=kmax-1,0,-1
               zz=zz+hn(i,j,k+1)
               rad(i,j,k)=swr_loc &
                      *(A(i,j)*exp(-zz/g1(i,j))+(1-A(i,j))*exp(-zz/g2(i,j)))
            end do
         end if
      end do
   end do

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
