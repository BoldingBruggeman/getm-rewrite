! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!> @note
!> how to deal with passing mask to calculation routines
!>
!> the NN calculation is not correct for non-equidistant vertical grids - as pmid does not reprensent the model
!> interface - does it matter? What is the error?
!>
!> what is S and/or T equations are excluded
!>
!> lat in NN calculation !!!
!>
!> @endnote

MODULE getm_density

   USE, INTRINSIC :: ISO_FORTRAN_ENV
!KB   use gsw_mod_toolbox, only: gsw_rho
   use memory_manager
   use logging
   use field_manager
   use getm_domain, only: type_getm_domain

   IMPLICIT NONE

   PRIVATE  ! Private scope by default

!  Module constants

!  Module types and variables
   type, public :: type_density
      !! author: Karsten Bolding
      !! version: v0.1
      !!
      !! Density type

      class(type_logging), pointer :: logs => null()
      class(type_field_manager), pointer :: fm => null()
      class(type_getm_domain), pointer :: domain

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
   class(type_logging), intent(in), target, optional :: logs
   type(type_field_manager), target, optional  :: fm

!  Local constants

!  Local variables
   integer :: rc
!---------------------------------------------------------------------------
   if (present(logs)) then
      self%logs => logs
      call self%logs%info('density_configuration()',level=2)
   end if
   if (present(fm)) then
      self%fm => fm
   end if
END SUBROUTINE density_configure

!---------------------------------------------------------------------------

SUBROUTINE density_initialize(self,domain)

   !! Feeds your cats and dogs, if enough food is available. If not enough
   !! food is available, some of your pets will get angry.

   IMPLICIT NONE

!  Subroutine arguments
   class(type_density), intent(inout) :: self
   class(type_getm_domain), intent(in), target :: domain

!  Local constants

!  Local variables
   integer :: i,j
   integer :: stat
!-----------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('density_initialize()',2)
   self%domain => domain
   TGrid: associate( TG => self%domain%T )
#ifndef _STATIC_
   call mm_s('rho',self%rho,TG%l,TG%u,def=-9999._real64,stat=stat)
   call mm_s('buoy',self%buoy,self%rho,def=-9999._real64,stat=stat)
   call mm_s('NN',self%NN,TG%l+(/0,0,-1/),TG%u,def=-9999._real64,stat=stat)
#endif
   if (associated(self%fm)) then
      call self%fm%register('rho', 'kg/m3', 'density', &
                            standard_name='sea_water_temperature', &
                            dimensions=(TG%dim_3d_ids), &
                            fill_value=-9999._real64, &
                            category='temperature_and_salinity')
      call self%fm%send_data('rho', self%rho(TG%imin:TG%imax,TG%jmin:TG%jmax,TG%kmin:TG%kmax))
      call self%fm%register('buoy', 'm/s2', 'buoyancy', &
                            standard_name='', &
                            dimensions=(TG%dim_3d_ids), &
                            fill_value=-9999._real64, &
                            category='temperature_and_salinity')
      call self%fm%send_data('buoy', self%buoy(TG%imin:TG%imax,TG%jmin:TG%jmax,TG%kmin:TG%kmax))
      call self%fm%register('NN', 's-2', 'Brunt-Vaisala frequency squared', &
                            standard_name='', &
                            dimensions=(TG%dim_3d_ids), &
                            fill_value=-9999._real64, &
                            category='temperature_and_salinity')
      call self%fm%send_data('NN', self%NN(TG%imin:TG%imax,TG%jmin:TG%jmax,TG%kmin:TG%kmax))
   end if
!   do j=grid%jmin,grid%jmax
!      do i=grid%imin,grid%imax
!         if (grid%mask(i,j) ==  0) then
!            self%rho(i,j,:) = -9999._real64
!            self%buoy(i,j,:) = -9999._real64
!         end if
!      end do
!   end do
   end associate TGrid
END SUBROUTINE density_initialize

!---------------------------------------------------------------------------

SUBROUTINE density_calculate(self,S,T,p)

   !! Feeds your cats and dogs, if enough food is available. If not enough
   !! food is available, some of your pets will get angry.

   use gsw_mod_toolbox, only: gsw_rho

   IMPLICIT NONE

!  Subroutine arguments
   class(type_density), intent(inout) :: self
#define _T3_ self%domain%T%l(1):,self%domain%T%l(2):,self%domain%T%l(3):
   real(real64), intent(in) :: S(_T3_)
      !! absolute salinity
   real(real64), intent(in) :: T(_T3_)
      !! conservative temperature [C|K?]
   real(real64), optional, intent(in) :: p(_T3_)
      !! Pressure [m or ?]
#undef _T3_

!  Local constants

!  Local variables
   integer :: i,j,k
!-----------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('density_calculate()',2)

#if 1
   TGrid: associate( TG => self%domain%T )
   do k=TG%l(3),TG%u(3)
      do j=TG%l(2),TG%u(2)
         do i=TG%l(1),TG%u(1)
            if (TG%mask(i,j) > 0) then
               self%rho(i,j,k) = gsw_rho(S(i,j,k), T(i,j,k), 0._real64)
            end  if
         end do
      end do
   end do
   end associate TGrid
#else
   if (present(p)) then
      self%rho = gsw_rho(S,T,p)
   end if
#endif
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
   real(real64), parameter :: rho0 = 1025._real64
      !! (/ rho0 /) = reference density
   real(real64), parameter :: x = g/rho0

!  Local variables
   integer :: i, j, k
!-----------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('buoyancy_calculate()',3)
   TGrid: associate( TG => self%domain%T )
   do k=TG%l(3),TG%u(3)
      do j=TG%l(2),TG%u(2)
         do i=TG%l(1),TG%u(1)
            if (TG%mask(i,j) > 0) then
               self%buoy(i,j,k) = x*(self%rho(i,j,k)-rho0)
            end  if
         end do
      end do
   end do
   end associate TGrid
END SUBROUTINE buoyancy_calculate

!---------------------------------------------------------------------------
!> @note
!> How to optimize this calculation - the loop order
!>
!> use gsw - either direct or via alpha and beta
!> [link](http://www.teos-10.org/pubs/gsw/html/gsw_Nsquared.html)
!> @endnote

SUBROUTINE brunt_vaisala_calculate(self,S,T)

   !!{!./code/brunt_vaisala.md!}

   use gsw_mod_toolbox, only: gsw_nsquared

   IMPLICIT NONE

!  Subroutine arguments
   class(type_density), intent(inout) :: self
#define _T3_ self%domain%T%l(1):,self%domain%T%l(2):,self%domain%T%l(3):
   real(real64), intent(in) :: S(_T3_)
      !! absolute salinity
   real(real64), intent(in) :: T(_T3_)
      !! conservative temperature [C|K?]
#undef _T3_

!  Local constants

!  Local variables
!KB
#if 0
   integer :: i, j, k
   real(real64), allocatable, dimension(:,:,:) :: x

   real(real64) :: small_bvf !!!! KB
   real(real64) :: gravity, rho0 !!!! KB
   real(real64) :: dz, NNc, NNe, NNn, NNw, NNs
#else
   real(real64), allocatable :: lat(:),pmid(:)
   integer :: i,j
#endif
!-----------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('brunt_vaisala_calculate()',3)
   TGrid: associate( TG => self%domain%T )
   allocate(lat(TG%l(3):TG%u(3)))
   allocate(pmid(TG%l(3):TG%u(3)-1))
   lat=0._real64 !KB!!!!!
   do j=TG%jmin,TG%jmax
      do i=TG%imin,TG%imax
         if (TG%mask(i,j) .ge. 1 ) then
            self%NN(i,j,TG%u(3))=0._real64
            call gsw_nsquared(S(i,j,:),T(i,j,:),TG%zc(i,j,:),lat,self%NN(i,j,TG%l(3):TG%u(3)-1),pmid) !KB check index
            self%NN(i,j,TG%l(3)-1)=0._real64; self%NN(i,j,TG%u(3))=0._real64
         end if
      end do
   end do
   end associate TGrid
END SUBROUTINE brunt_vaisala_calculate

!---------------------------------------------------------------------------

END MODULE getm_density

#if 0
#if 1
#else
   do j=TG%jmin-1,TG%jmax+1
      do i=TG%imin-1,TG%imax+1
         if (TG%mask(i,j) .ge. 1 ) then
            self%NN(i,j,TG%kmax) = small_bvf
            do k=TG%kmax-1,1,-1
               dz=0.5_real64*(TG%hn(i,j,k+1)+TG%hn(i,j,k))
#define _OLD_BVF_
#ifdef _OLD_BVF_
               NNc =(self%buoy(i,j,k+1)-self%buoy(i,j,k))/dz
#else
               NNc = -gravity / rho0
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
   do j=TG%jmin,TG%jmax
      do i=TG%imin,TG%imax
         if (TG%mask(i,j) .ge. 1 ) then
            do k=TG%kmax-1,1,-1
               NNc = NN(i,j,k)
               if (TG%mask(i+1,j) .ge. 1) then
                  NNe= NN(i+1,j,k)
               else
                  NNe=NNc
               end if
               if (TG%mask(i-1,j) .ge. 1) then
                  NNw= NN(i-1,j,k)
               else
                  NNw=NNc
               end if
               if (TG%mask(i,j+1) .ge. 1) then
                  NNn= NN(i,j+1,k)
               else
                  NNn=NNc
               end if
               if (TG%mask(i,j-1) .ge. 1) then
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
#endif
