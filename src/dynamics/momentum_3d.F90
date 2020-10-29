! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!!{!./pages/momentum_3d.md!}

SUBMODULE (getm_momentum) momentum_3d_smod

CONTAINS

MODULE PROCEDURE uv_momentum_3d
   !! Solve the 3D momemtum equations

   IMPLICIT NONE

!  Local constants

!  Local variables
   logical, save :: ufirst=.false. ! should likely not have save
!---------------------------------------------------------------------------
   call self%logs%info('uv_momentum_3d()',level=2)
   self%Ui=self%Ui/mode_split
   self%Vi=self%Vi/mode_split

   if(ufirst) then
      call u_3d(self,dt,tausx,dpdx,idpdy,viscosity)
      call v_3d(self,dt,tausy,dpdy,idpdy,viscosity)
      ufirst = .false.
   else
      call v_3d(self,dt,tausy,dpdy,idpdy,viscosity)
      call u_3d(self,dt,tausx,dpdx,idpdy,viscosity)
      ufirst = .true.
   end if
END PROCEDURE uv_momentum_3d

!---------------------------------------------------------------------------

MODULE PROCEDURE w_momentum_3d

   IMPLICIT NONE

!  Local constants

!  Local variables
   integer :: i,j,k
   real(real64) :: dtm1
!---------------------------------------------------------------------------
   call self%logs%info('momentum_z_3d()',level=2)
   dtm1=1._real64/dt
   TGrid: associate( TG => self%domain%T )
   UGrid: associate( UG => self%domain%U )
   VGrid: associate( VG => self%domain%V )
   do k=TG%l(3),TG%u(3)-1
      do j=TG%l(2)+1,TG%u(2)
         do i=TG%l(1)+1,TG%u(1)
            if (TG%mask(i,j) == 1) then
                  self%ww(i,j,k) = self%ww(i,j,k-1) - (TG%hn(i,j,k)-TG%ho(i,j,k))*dtm1 &
                              -(self%pk(i,j,k)*UG%dy(i,j) - self%pk(i-1,j  ,k)*UG%dy(i-1,j) &
                               +self%qk(i,j,k)*VG%dx(i,j) - self%qk(i  ,j-1,k)*VG%dx(i,j-1)) &
                               *TG%inv_area(i,j)
            end if
         end do
      end do
   end do
   end associate VGrid
   end associate UGrid
   end associate TGrid
END PROCEDURE w_momentum_3d

!---------------------------------------------------------------------------

MODULE PROCEDURE uv_advection_3d
   !! 3D velocity advection

   IMPLICIT NONE

!  Local constants

!  Local variables
   integer :: i,j,k
!---------------------------------------------------------------------------
   XGrid: associate( XG => self%domain%X )
   TGrid: associate( TG => self%domain%T )
   UGrid: associate( UG => self%domain%U )
   VGrid: associate( VG => self%domain%V )
   do j=UG%jmin,UG%jmax
      do i=UG%imin,UG%imax
         self%uadvgrid%hn(i,j,:) = TG%hn(i+1,j,:)
         self%vadvgrid%hn(i,j,:) = XG%hn(i,j,:)
         self%pkadv(i,j,:) = 0.5_real64*(self%pk(i,j,:) + self%pk(i+1,j,:))
         self%qkadv(i,j,:) = 0.5_real64*(self%qk(i,j,:) + self%qk(i+1,j,:))
      end do
   end do
   call self%advection%calculate(self%advection_scheme,self%uadvgrid,self%pkadv,self%vadvgrid,self%qkadv,dt,UG,self%pk)

   do j=UG%jmin,UG%jmax
      do i=UG%imin,UG%imax
         self%uadvgrid%hn(i,j,:) = XG%hn(i,j,:)
         self%vadvgrid%hn(i,j,:) = TG%hn(i,j+1,:)
         self%pkadv(i,j,:) = 0.5_real64*(self%pk(i,j,:) + self%pk(i,j+1,:))
         self%qkadv(i,j,:) = 0.5_real64*(self%qk(i,j,:) + self%qk(i,j+1,:))
      end do
   end do
   call self%advection%calculate(self%advection_scheme,self%uadvgrid,self%pkadv,self%vadvgrid,self%qkadv,dt,VG,self%qk)
   end associate VGrid
   end associate UGrid
   end associate TGrid
   end associate XGrid
END PROCEDURE uv_advection_3d

!---------------------------------------------------------------------------

module SUBROUTINE u_3d(self,dt,taus,dpdx,idpdx,viscosity)

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
      !! GETM momentum type
   real(real64), intent(in) :: dt
      !! timestep [s]
#define _T2_ self%domain%T%l(1):,self%domain%T%l(2):
   real(real64), intent(in) :: taus(_T2_)
      !! surface stress in X-direction
   real(real64), intent(in) :: dpdx(_T2_)
      !! surface pressure gradient - including air pressure
#undef _T2_
#define _T3_ self%domain%T%l(1):,self%domain%T%l(2):,self%domain%T%l(3):
   real(real64), intent(in) :: idpdx(_T3_)
      !! internal pressure gradient
   real(real64), intent(in) :: viscosity(_T3_)
#undef _T3_
      !! viscosity

!  Local constants

!  Local variables
   real(real64) :: Vloc
!KB
   real(real64) :: Slr, ip_fac=1._real64
!KB
   real(real64), dimension(:,:,:), allocatable :: ea2,ea4,num
   integer :: i,j,k
!---------------------------------------------------------------------------
   call self%logs%info('u_3d()',level=3)

   allocate(ea2,mold=self%pk) !KB se note in temperature.F90
   allocate(ea4,mold=self%pk)
   allocate(num,mold=viscosity)

   UGrid: associate( UG => self%domain%U )
   do k=UG%kmin,UG%kmax
      do j=UG%jmin,UG%jmax
         do i=UG%imin,UG%imax
            ea2(i,j,k)=0._real64
            ea4(i,j,k)=0._real64
            if (UG%mask(i,j) == 1 .or. UG%mask(i,j) == 2) then

               ! Eddy viscosity at U-points
               num(i,j,k)=0.5_real64*(num(i,j,k)+num(i+1,j,k))

               ! surface pressure
               ea4(i,j,k)=-dt*0.5_real64*(UG%ho(i,j,k)+UG%hn(i,j,k))*g*dpdx(i,j)
#if 0
               ! Coriolis
               ! Espelid et al. [2000], IJNME 49, 1521-1545
#ifdef NEW_CORI
               Vloc=0.25_real64*(work3d(i,j,k)+work3d(i+1,j,k)+work3d(i,j-1,k)+work3d(i+1,j-1,k))*sqrt(huo(i,j,k))
#else
               Vloc=0.25_real64*(self%qk(i,j,k)+self%qk(i+1,j,k)+self%qk(i,j-1,k)+self%qk(i+1,j-1,k))
#endif
#if defined(SPHERICAL) || defined(CURVILINEAR)
               cord_curv=(Vloc*(DYCIP1-DYC)-uu(i,j,k)*(DXX-DXXJM1))   &
                     /huo(i,j,k)*ARUD1
               ex(k)=(cord_curv+coru(i,j))*Vloc
#else
               ex(k)=coru(i,j)*Vloc
#endif
#endif
               ! .... and internal pressure gradient
               ea4(i,j,k)=UG%alpha(i,j)*(ea4(i,j,k)-self%uuEx(i,j,k)+ip_fac*idpdx(i,j,k))
            end if
         end do
      end do
   end do

   ! Additional matrix elements for surface and bottom layer
   do j=UG%jmin,UG%jmax
      do i=UG%imin,UG%imax
         if (UG%mask(i,j) == 1 .or. UG%mask(i,j) == 2) then
            ! surface stress
            ea4(i,j,UG%kmax)=ea4(i,j,UG%kmax)+UG%alpha(i,j)*0.5_real64*(taus(i,j)+taus(i+1,j))/rho0
            ! bottom friction
            ea2(i,j,UG%kmin)=dt*self%rru(i,j)/(0.5_real64*(UG%ho(i,j,UG%kmin)+UG%hn(i,j,UG%kmin)))
         end if
      end do
   end do

!KB scale with huo
   call self%vertical_diffusion%calculate(UG%mask,UG%ho,UG%hn,dt,self%cnpar,self%molecular,num,self%pk,ea2,ea4)
!KB un-scale with huo

#if 0
   ! Transport correction: the integral of the new velocities has to
   ! be the same than the transport calculated by the external mode, Uint.

   ResInt= _ZERO_
   do k=kumin(i,j),kmax
      ResInt=ResInt+Res(k)
   end do
#ifdef MUDFLAT
   Diff=(Uint(i,j)-ResInt)/(ssun(i,j)+HU(i,j))
#else
   Diff=(Uint(i,j)-ResInt)/(ssuo(i,j)+HU(i,j))
#endif

   do k=kumin(i,j),kmax
#ifndef NO_BAROTROPIC
      uu(i,j,k)=Res(k) +hun(i,j,k)*Diff
#else
      uu(i,j,k)=Res(k)
#endif
#endif
   end associate UGrid
END SUBROUTINE u_3d

!---------------------------------------------------------------------------

module SUBROUTINE v_3d(self,dt,taus,dpdy,idpdy,viscosity)

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
      !! GETM momentum type
   real(real64), intent(in) :: dt
      !! timestep [s]
#define _T2_ self%domain%T%l(1):,self%domain%T%l(2):
   real(real64), intent(in) :: taus(_T2_)
      !! surface stress in X-direction
   real(real64), intent(in) :: dpdy(_T2_)
      !! surface pressure gradient - including air pressure
#undef _T2_
#define _T3_ self%domain%T%l(1):,self%domain%T%l(2):,self%domain%T%l(3):
   real(real64), intent(in) :: idpdy(_T3_)
      !! internal pressure gradient
   real(real64), intent(in) :: viscosity(_T3_)
      !! viscosity
#undef _T3_

!  Local constants

!  Local variables
   integer :: i,j
   real(real64) :: tausu
!KB
   real(real64) :: Slr
!KB
!---------------------------------------------------------------------------
   call self%logs%info('v_3d()',level=3)
END SUBROUTINE v_3d

!---------------------------------------------------------------------------

END SUBMODULE momentum_3d_smod
