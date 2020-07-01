! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!!{!./pages/momentum_3d.md!}

SUBMODULE (getm_momentum) momentum_3d_smod

CONTAINS

!---------------------------------------------------------------------------

module SUBROUTINE momentum_x_3d(self,dt,dpdx,taus)

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
      !! GETM momentum type
   real(real64), intent(in) :: dt
      !! timestep [s]
   real(real64), dimension(:,:), intent(in) :: dpdx
      !! surface pressure gradient - including air pressure
   real(real64), dimension(:,:), intent(in) :: taus
      !! surface stress in X-direction

!  Local constants

!  Local variables
   integer :: i,j
   real(real64) :: tausu
!KB
   real(real64) :: Slr, gammai
   real(real64) :: Dmin = 0.5_real64
!KB
!---------------------------------------------------------------------------
   call self%logs%info('velocity_3d_x()',level=2)
#if 0

#endif
   return
END SUBROUTINE momentum_x_3d

!---------------------------------------------------------------------------

module SUBROUTINE momentum_y_3d(self,dt,dpdy,taus)

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
      !! GETM momentum type
   real(real64), intent(in) :: dt
      !! timestep [s]
   real(real64), dimension(:,:), intent(in) :: dpdy
      !! surface pressure gradient - including air pressure
   real(real64), dimension(:,:), intent(in) :: taus
      !! surface stress in X-direction

!  Local constants

!  Local variables
   integer :: i,j
   real(real64) :: tausu
!KB
   real(real64) :: Slr, gammai
   real(real64) :: Dmin = 0.5_real64
!KB
!---------------------------------------------------------------------------
   call self%logs%info('momentum_3d_y()',level=2)
#if 0
#endif
   return
END SUBROUTINE momentum_y_3d

!---------------------------------------------------------------------------

module SUBROUTINE momentum_z_3d(self,dt)

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
      !! GETM momentum type
   real(real64), intent(in) :: dt
      !! timestep [s]

!  Local constants

!  Local variables
   integer :: i,j,k
   real(real64) :: dtm1
!---------------------------------------------------------------------------
   call self%logs%info('momentum_z_3d()',level=2)
   dtm1=1._real64/dt
   do k=self%domain%T%l(3),self%domain%T%u(3)-1
      do j=self%domain%T%l(2)+1,self%domain%T%u(2)
         do i=self%domain%T%l(1)+1,self%domain%T%u(1)
            if (self%domain%T%mask(i,j) .eq. 1) then
                  self%ww(i,j,k) = self%ww(i,j,k-1) &
                              -(self%domain%T%hn(i,j,k)-self%domain%T%ho(i,j,k))*dtm1 &
                              -(self%pk(i,j,k)*self%domain%U%dy(i,j) - self%pk(i-1,j  ,k)*self%domain%U%dy(i-1,j) &
                               +self%qk(i,j,k)*self%domain%V%dx(i,j) - self%qk(i  ,j-1,k)*self%domain%V%dx(i,j-1)) &
                               *self%domain%T%inv_area(i,j)
!KB                              -(self%pk(i,j,k)*DYU - self%pk(i-1,j  ,k)*DYUIM1 &
!KB                               +self%qk(i,j,k)*DXV - self%qk(i  ,j-1,k)*DXVJM1)*ARCD1
            end if
         end do
      end do
   end do
   return
END SUBROUTINE momentum_z_3d

!---------------------------------------------------------------------------

END SUBMODULE momentum_3d_smod

#if 0
#ifdef CALC_HALO_WW
#ifndef SLICE_MODEL
      do j=jmin-1,jmax+HALO
#endif
         do i=imin-1,imax+HALO
#else
#ifndef SLICE_MODEL
      do j=jmin,jmax
#endif
         do i=imin,imax
#endif
            if (az(i,j) .eq. 1) then
               if (k .lt. kmin(i,j)) then
                  ww(i,j,k)= _ZERO_
               else
                  ww(i,j,k) =   ww(i,j,k-1)                             &
                              - ( hn(i,j,k) - ho(i,j,k) )*dtm1          &
                              - (                                       &
                                   pk(i,j,k)*DYU - pk(i-1,j  ,k)*DYUIM1 &
#ifndef SLICE_MODEL
                                 + qk(i,j,k)*DXV - qk(i  ,j-1,k)*DXVJM1 &
#endif
                                )*ARCD1
               end if
            end if
         end do
#ifndef SLICE_MODEL
      end do
#endif
!$OMP END DO
   end do

   return
END SUBROUTINE momentum_z_3d
#endif
