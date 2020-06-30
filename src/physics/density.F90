! Copyright (C) 2020 Bolding & Bruggeman

!> @note
!> how to deal with passing mask to calculation routines
!>
!> @endnote

#ifdef _STATIC_
#include "dimensions.h"
#endif
MODULE getm_density

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use gsw_mod_toolbox, only: gsw_rho
   use memory_manager
   use logging
   use field_manager
   use getm_domain, only: type_getm_grid
!KB   use getm_variables

   IMPLICIT NONE

   PRIVATE  ! Private scope by default

!  Module constants

!  Module types and variables
   type, public :: type_density
      !! author: Karsten Bolding
      !! version: v0.1
      !!
      !! Density type

      class(type_logging), pointer :: logs
      class(type_getm_grid), pointer :: grid
      class(type_field_manager), pointer :: fm => null()

#ifdef _STATIC_
      real(real64), dimension(A3DFIELD) :: rho, buoy, NN
#else
      real(real64), dimension(:,:,:), allocatable :: rho, buoy, NN
#endif
      contains

      procedure :: configure => density_configure
      procedure :: initialize => density_initialize
      procedure :: density => density_calculate
      procedure :: buoyancy => buoyancy_calculate
      procedure :: brunt_vaisala => brunt_vaisala_calculate

   end type type_density

CONTAINS

SUBROUTINE density_configure(self,logs,fm)

   !! Configure the the density module

   IMPLICIT NONE

!  Subroutine arguments
   class(type_density), intent(inout) :: self
   class(type_logging), intent(in), target :: logs
   type(type_field_manager), optional, target  :: fm

!  Local constants

!  Local variables
   integer :: rc
!---------------------------------------------------------------------------
   self%logs => logs
   call self%logs%info('density_configuration()',level=2)

   if (present(fm)) then
      self%fm => fm
   end if
   return
END SUBROUTINE density_configure

!---------------------------------------------------------------------------

SUBROUTINE density_initialize(self,grid)

   !! Feeds your cats and dogs, if enough food is available. If not enough
   !! food is available, some of your pets will get angry.

   IMPLICIT NONE

!  Subroutine arguments
   class(type_density), intent(inout) :: self
   class(type_getm_grid), intent(in), target :: grid

!  Local constants

!  Local variables
   integer :: i,j
   integer :: stat
!-----------------------------------------------------------------------------
   call self%logs%info('density_initialize()',2)
   self%grid => grid
#ifndef _STATIC_
   call mm_s('rho',self%rho,grid%l,grid%u,def=-9999._real64,stat=stat)
   call mm_s('buoy',self%buoy,self%rho,def=-9999._real64,stat=stat)
   call mm_s('NN',self%NN,self%rho,def=-9999._real64,stat=stat)
#endif
   if (associated(self%fm)) then
      call self%fm%register('rho', 'kg/m3', 'density', &
                            standard_name='sea_water_temperature', &
                            dimensions=(self%grid%dim_3d_ids), &
                            fill_value=-9999._real64, &
                            category='temperature_and_salinity')
      call self%fm%send_data('rho', self%rho)
      call self%fm%register('buoy', 'm/s2', 'buoyancy', &
                            standard_name='', &
                            dimensions=(self%grid%dim_3d_ids), &
                            fill_value=-9999._real64, &
                            category='temperature_and_salinity')
      call self%fm%send_data('buoy', self%buoy)
   end if
!   do j=grid%jmin,grid%jmax
!      do i=grid%imin,grid%imax
!         if (grid%mask(i,j) ==  0) then
!            self%rho(i,j,:) = -9999._real64
!            self%buoy(i,j,:) = -9999._real64
!         end if
!      end do
!   end do
   return
END SUBROUTINE density_initialize

!---------------------------------------------------------------------------

SUBROUTINE density_calculate(self,S,T,p)

   !! Feeds your cats and dogs, if enough food is available. If not enough
   !! food is available, some of your pets will get angry.

!KB   use gsw_mod_toolbox, only: gsw_rho

   IMPLICIT NONE

!  Subroutine arguments
   class(type_density), intent(inout) :: self
   real(real64), dimension(:,:,:), intent(in) :: S,T
      !! Salinity and temperature
   real(real64), dimension(:,:,:), optional, intent(in) :: p
      !! Pressure

!  Local constants

!  Local variables
   integer :: i,j,k
!-----------------------------------------------------------------------------
!KB   call plogs%info('density_calculate()',2)

#if 1
   do k=self%grid%l(3),self%grid%u(3)
      do j=self%grid%l(2),self%grid%u(2)
         do i=self%grid%l(1),self%grid%u(1)
            if (self%grid%mask(i,j) > 0) then
               self%rho(i,j,k) = gsw_rho(S(i,j,k), T(i,j,k), 0._real64)
            end  if
         end do
      end do
   end do
#else
   if (present(p)) then
      self%rho = gsw_rho(S,T,p)
   end if
#endif

  return
END SUBROUTINE density_calculate

!-----------------------------------------------------------------------------

SUBROUTINE buoyancy_calculate(self)
   !! Feeds your cats and dogs, if enough food is available. If not enough
   !! food is available, some of your pets will get angry.

   IMPLICIT NONE

!  Subroutine arguments
   class(type_density), intent(inout) :: self

!  Local constants
   real(real64), parameter :: g = 9.81_real64
      !! g = gravitational acceleration
   real(real64), parameter :: rho_0 = 1025._real64
      !! (/ rho_0 /) = reference density
   real(real64), parameter :: x = g/rho_0

!  Local variables
   integer :: i, j, k
!-----------------------------------------------------------------------------
   do k=self%grid%l(3),self%grid%u(3)
      do j=self%grid%l(2),self%grid%u(2)
         do i=self%grid%l(1),self%grid%u(1)
            if (self%grid%mask(i,j) > 0) then
               self%buoy(i,j,k) = x*(self%rho(i,j,k)-rho_0)
            end  if
         end do
      end do
   end do
  return
END SUBROUTINE buoyancy_calculate

!---------------------------------------------------------------------------
!> @note
!> How to optimize this calculation - the loop order
!>
!> use gsw - either direct or via alpha and beta
!> [link](http://www.teos-10.org/pubs/gsw/html/gsw_Nsquared.html)
!> @endnote

SUBROUTINE brunt_vaisala_calculate(self)

   !!{!./code/brunt_vaisala.md!}

   IMPLICIT NONE

!  Subroutine arguments
   class(type_density), intent(inout) :: self
!KB   real(real64), dimension(:,:,:), intent(inout) :: NN

!  Local constants

!  Local variables
   integer :: i, j, k
   real(real64), allocatable, dimension(:,:,:) :: x

   real(real64) :: small_bvf !!!! KB
   real(real64) :: gravity, rho_0 !!!! KB
   real(real64) :: dz, NNc, NNe, NNn, NNw, NNs
!-----------------------------------------------------------------------------
   call self%logs%info('brunt_vaisala()',2)

#ifndef NO_BAROCLINIC
#define G self%grid
   do j=G%jmin-1,G%jmax+1
      do i=G%imin-1,G%imax+1
         if (G%mask(i,j) .ge. 1 ) then
            self%NN(i,j,G%kmax) = small_bvf
            do k=G%kmax-1,1,-1
               dz=0.5_real64*(G%hn(i,j,k+1)+G%hn(i,j,k))
#define _OLD_BVF_
#ifdef _OLD_BVF_
               NNc =(self%buoy(i,j,k+1)-self%buoy(i,j,k))/dz
#else
               NNc = -gravity / rho_0
               NNc = NNc/dz *(alpha(i,j,k) *(T(i,j,k+1)-T(i,j,k)) &
                         + beta(i,j,k) *(S(i,j,k+1)-S(i,j,k)))
#endif
#undef _OLD_BVF_
               if (abs(NNc) .lt. small_bvf ) then
                  NNc = sign(1._real64,NNc) * small_bvf
               end if
               self%NN(i,j,k)= NNc
            end do
#ifdef _SMOOTH_BVF_VERT_
            if ( kmax .ge. 4 ) then
               below  = NN(i,j,1)
               center = NN(i,j,2)
               above  = NN(i,j,3)
               do k= 2,kmax-2
                  center = 0.5_real64 * center + _QUART_ * (below+above)
                  below  = NN(i,j,k)
                  NN(i,j,k) = center
                  center = NN(i,j,k+1)
                  above  = NN(i,j,k+2)
               end do
            end if
#endif
         end if
      end do
   end do

#ifdef SMOOTH_BVF_HORI
   do j=G%jmin,G%jmax
      do i=G%imin,G%imax
         if (G%mask(i,j) .ge. 1 ) then
            do k=G%kmax-1,1,-1
               NNc = NN(i,j,k)
               if (G%mask(i+1,j) .ge. 1) then
                  NNe= NN(i+1,j,k)
               else
                  NNe=NNc
               end if
               if (G%mask(i-1,j) .ge. 1) then
                  NNw= NN(i-1,j,k)
               else
                  NNw=NNc
               end if
               if (G%mask(i,j+1) .ge. 1) then
                  NNn= NN(i,j+1,k)
               else
                  NNn=NNc
               end if
               if (G%mask(i,j-1) .ge. 1) then
                  NNs= NN(i,j-1,k)
               else
                  NNs=NNc
               end if
               NN(i,j,k) = (1._real64/3._real64)*NNc+(1._real64/6._real64)*(NNe+NNw+NNn+NNs)
            end do
         end if
      end do
   end do
#endif
#endif
#undef G
  return
END SUBROUTINE brunt_vaisala_calculate

!---------------------------------------------------------------------------

END MODULE getm_density
