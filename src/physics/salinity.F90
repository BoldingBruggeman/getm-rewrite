! Copyright (C) 2020 Bolding & Bruggeman

!> In this module, the salinity equation is processed by
!> reading in the namelist {\tt salt} and initialising the salinity field
!> (this is done in {\tt init\_salinity}),
!> and calculating the advection-diffusion-equation, which includes
!> penetrating short-wave radiation as source term (see {\tt do\_salinity}).

!> @note
!> ckeck dimension order of auxo and auxn
!> ckeck dimension order of a1, a2, a3, a4
!> @endnote

#ifdef _STATIC_
#include "dimensions.h"
#endif

MODULE getm_salinity

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
   real(real64), parameter :: avmols = 0._real64

!  Module types and variables
   real(real64) :: cnpar
   type, public :: type_salinity_configuration
      character(len=256) :: f = "salinity.nc"
   end type type_salinity_configuration

   type, public :: type_salinity
      !! author: Karsten Bolding
      !! version: v0.1
      !!
      !! Salinity type

      class(type_logging), pointer :: logs
      class(type_getm_grid), pointer :: G
      class(type_field_manager), pointer :: fm => null()
      TYPE(type_salinity_configuration) :: config

#ifdef _STATIC_
    real(real64), dimension(I3DFIELD) :: S = 10._real64
#else
    real(real64), dimension(:,:,:), allocatable :: S
#endif
!> remember to allocate in initialize
      real(real64), private, dimension(:,:,:), allocatable :: auxo, auxn
      real(real64), private, dimension(:,:,:), allocatable :: a1,a2,a3,a4

      contains

      procedure :: configuration => salinity_configuration
      procedure :: initialize => salinity_initialize
      procedure :: diffusion => salinity_diffusion
      procedure :: calculate => salinity_calculate

   end type type_salinity

CONTAINS

!---------------------------------------------------------------------------

SUBROUTINE salinity_configuration(self,logs,fm)

   !! Configure the salinity type

   IMPLICIT NONE

!  Subroutine arguments
   class(type_salinity), intent(out) :: self
   class(type_logging), intent(in), target :: logs
   type(type_field_manager), optional, target  :: fm

!  Local constants

!  Local variables
   integer :: rc
!---------------------------------------------------------------------------
   self%logs => logs
   call self%logs%info('salinity_configuration()',level=2)
   call self%logs%info('reading initial salinity from: ',level=3,msg2=trim(self%config%f))

   if (present(fm)) then
      self%fm => fm
   end if
   return
END SUBROUTINE salinity_configuration

!---------------------------------------------------------------------------

SUBROUTINE salinity_initialize(self,grid)

   !! Initialize the salinity field

   IMPLICIT NONE

!  Subroutine arguments
   class(type_salinity), intent(inout) :: self
   class(type_getm_grid), intent(in), target :: grid
      !! grid dimensions in case of dynamic memory allocation

!  Local constants

!  Local variables
   integer :: i,j
   integer :: stat
!---------------------------------------------------------------------------
   call self%logs%info('salinity_initialize()',level=2)

   self%G => grid
#ifndef _STATIC_
   call mm_s('S',self%S,grid%l,grid%u,def=25._real64,stat=stat)
   call mm_s('auxo',self%auxo,grid%l,grid%u,def=0._real64,stat=stat)
   call mm_s('auxn',self%auxn,grid%l,grid%u,def=0._real64,stat=stat)
!   call mm_s('a1',self%a1, k, i, j, def=0._real64,stat=stat)
!   call mm_s('a2',self%a2, k, i, j, def=0._real64,stat=stat)
!   call mm_s('a3',self%a3, k, i, j, def=0._real64,stat=stat)
!   call mm_s('a4',self%a4, k, i, j, def=0._real64)
#endif
   if (associated(self%fm)) then
      call self%fm%register('salt', 'g/kg', 'absolute salinity', &
                            standard_name='sea_water_salinity', &
                            category='temperature_and_salinity', &
                            dimensions=(self%G%dim_3d_ids), &
                            fill_value=-9999._real64, &
                            part_of_state=.true.)
      call self%fm%send_data('salt', self%S(grid%imin:grid%imax,grid%jmin:grid%jmax,grid%kmin:grid%kmax))
   end if
   do j=self%G%jmin,self%G%jmax
      do i=self%G%imin,self%G%imax
         if (self%G%mask(i,j) ==  0) then
            self%S(i,j,:) = -9999._real64
         end if
      end do
   end do
   return
END SUBROUTINE salinity_initialize

!---------------------------------------------------------------------------

SUBROUTINE salinity_calculate(self)

   !! Advection/diffusion of the salinity field

   IMPLICIT NONE

!  Subroutine arguments
   class(type_salinity), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: rc
!---------------------------------------------------------------------------
   call self%logs%info('salinity_calculate()',level=2)

   return
END SUBROUTINE salinity_calculate

!---------------------------------------------------------------------------

SUBROUTINE salinity_diffusion(self,logs,dt,nuh)

   !! Advection/diffusion of the salinity field

   IMPLICIT NONE

!  Subroutine arguments
   class(type_salinity), intent(out) :: self
   class(type_logging), intent(in) :: logs
   real(real64), intent(in) :: dt
   real(real64), dimension(:,:,:), intent(in) :: nuh

!  Local constants

!  Local variables
   integer :: i,j,k
   integer :: imin,imax,jmin,jmax,kmin,kmax
   real(real64) :: x
!---------------------------------------------------------------------------
   call logs%info('salinity_diffusion()',level=2)

   imin = self%G%imin; imax = self%G%imax
   jmin = self%G%kmin; jmax = self%G%jmax
   kmin = self%G%kmin; kmax = self%G%kmax

   if (kmax > 1) then  !KB if (grid%S%kmax.gt.1) then
      do k=kmin,kmax-1
         do j=jmin,jmax
            do i=imin,imax
               if (self%g%mask(i,j) ==  1) then
!                 !! Auxilury terms, old and new time level,
                  x = dt*(nuh(i,j,k)+avmols)/ &
                          (self%g%hn(i,j,k+1)+self%g%hn(i,j,k))
                  self%auxo(i,j,k)=2._real64*(1-cnpar)*x
                  self%auxn(i,j,k)=2._real64*   cnpar *x
               end if
            end do
         end do
      end do

!     Matrix elements for surface layer
      k=kmax
      do j=jmin,jmax
         do i=imin,imax
            if (self%g%mask(i,j) ==  1) then
               self%a1(k,i,j)=-self%auxn(i,j,k-1)
               self%a2(k,i,j)=self%g%hn(i,j,k)+self%auxn(i,j,k-1)
               self%a4(k,i,j)=self%S(i,j,k)*(self%g%hn(i,j,k)-self%auxo(i,j,k-1))+self%S(i,j,k-1)*self%auxo(i,j,k-1)
            end if
         end do
      end do

!     Matrix elements for inner layers
!      do k=kmin,kmax-1
      do k=2,kmax-1
         do j=jmin,jmax
            do i=imin,imax
               if (self%g%mask(i,j) ==  1) then
                  self%a3(k,i,j)=-self%auxn(i,j,k  )
                  self%a1(k,i,j)=-self%auxn(i,j,k-1)
                  self%a2(k,i,j)=self%g%hn(i,j,k)+self%auxn(i,j,k)+self%auxn(i,j,k-1)
                  self%a4(k,i,j)=self%S(i,j,k+1)*self%auxo(i,j,k) &
                       +self%S(i,j,k  )*(self%g%hn(i,j,k)-self%auxo(i,j,k)-self%auxo(i,j,k-1)) &
                       +self%S(i,j,k-1)*self%auxo(i,j,k-1)
               end if
            end do
         end do
      end do

!     Matrix elements for bottom layer
      k=1
      do j=jmin,jmax
         do i=imin,imax
            if (self%g%mask(i,j) ==  1) then
               self%a3(k,i,j)=-self%auxn(i,j,k  )
               self%a2(k,i,j)=self%g%hn(i,j,k)+self%auxn(i,j,k)
               self%a4(k,i,j)=self%S(i,j,k+1)*self%auxo(i,j,k)                              &
                    +self%S(i,j,k  )*(self%g%hn(i,j,k)-self%auxo(i,j,k))
            end if
         end do
      end do

   end if

#if 0
   call getm_tridiagonal(kmax,1,kmax,a1,a2,a3,a4,Res)

   do k=1,kmax
      S(i,j,k)=Res(k)
   end do
#endif

   return
END SUBROUTINE salinity_diffusion

!---------------------------------------------------------------------------

END MODULE getm_salinity
