! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!! @note
!! loop boundaries in w_momentum
!! @endnote


!!{!./pages/momentum_3d.md!}

SUBMODULE (getm_momentum) momentum_3d_smod

   real(real64) :: ip_fac=1._real64 !KB

CONTAINS

!---------------------------------------------------------------------------

MODULE SUBROUTINE uv_initialize_3d(self)
   !! initialize the 3D momentum sub-module

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: stat
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('uv_initialize_3d()',level=2)
   TGrid: associate( TG => self%domain%T )
   UGrid: associate( UG => self%domain%U )
   VGrid: associate( VG => self%domain%V )
   call mm_s('pk',self%pk,UG%l(1:3),UG%u(1:3),def=0._real64,stat=stat)
   call mm_s('qk',self%qk,VG%l(1:3),VG%u(1:3),def=0._real64,stat=stat)
   call mm_s('ww',self%ww,self%pk,def=0._real64,stat=stat)
   call mm_s('pka',self%pka,self%pk,def=0._real64,stat=stat)
   call mm_s('qka',self%qka,self%qk,def=0._real64,stat=stat)
   call mm_s('fpk',self%fpk,self%pk,def=0._real64,stat=stat)
   call mm_s('fqk',self%fqk,self%qk,def=0._real64,stat=stat)
   call mm_s('advpk',self%advpk,self%pk,def=0._real64,stat=stat)
   call mm_s('advqk',self%advqk,self%qk,def=0._real64,stat=stat)
   call mm_s('diffuk',self%diffuk,self%pk,def=0._real64,stat=stat)
   call mm_s('diffvk',self%diffvk,self%qk,def=0._real64,stat=stat)
   call mm_s('uk',self%uk,self%pk,def=0._real64,stat=stat)
   call mm_s('vk',self%vk,self%qk,def=0._real64,stat=stat)
   call mm_s('taus',self%taus,TG%l(1:2),TG%u(1:2),def=0._real64,stat=stat)
   call mm_s('taub',self%taub,TG%l(1:2),TG%u(1:2),def=0._real64,stat=stat)
   call mm_s('taubx',self%taubx,self%U,def=0._real64,stat=stat)
   call mm_s('tauby',self%tauby,self%V,def=0._real64,stat=stat)
   call mm_s('SS',self%SS,TG%l+(/0,0,-1/),TG%u,def=0._real64,stat=stat)
   if (self%advection_scheme > 0) then
      call mm_s('uadvhn',self%uadvgrid%hn,TG%hn,def=0._real64,stat=stat)
      call mm_s('vadvhn',self%vadvgrid%hn,TG%hn,def=0._real64,stat=stat)
   end if
   call mm_s('num',self%num,self%pk,def=0._real64,stat=stat)
   call mm_s('ea2',self%ea2,self%pk,def=0._real64,stat=stat)
   call mm_s('ea4',self%ea4,self%pk,def=0._real64,stat=stat)

   end associate VGrid
   end associate UGrid
   end associate TGrid
END SUBROUTINE uv_initialize_3d

!---------------------------------------------------------------------------

MODULE SUBROUTINE uvw_momentum_3d(self,dt,tausx,tausy,dpdx,dpdy,idpdx,idpdy,viscosity)
   !! Solve the horizontal 3D momemtum equations - pk and qk

   IMPLICIT NONE

   class(type_getm_momentum), intent(inout) :: self
   real(real64), intent(in) :: dt
      !! timestep [s]
#define _T2_ self%domain%T%l(1):,self%domain%T%l(2):
   real(real64), intent(in) :: tausx(_T2_)
     !! surface stress - x
   real(real64), intent(in) :: tausy(_T2_)
     !! surface stress - y
   real(real64), intent(in) :: dpdx(_T2_)
     !! surface pressure (including air pressure) - x-gradient
   real(real64), intent(in) :: dpdy(_T2_)
     !! surface pressure (including air pressure) - y-gradient
#undef _T2_
#define _T3_ self%domain%T%l(1):,self%domain%T%l(2):,self%domain%T%l(3):
   real(real64), intent(in) :: idpdx(_T3_)
     !! internal pressure - x-gradient
   real(real64), intent(in) :: idpdy(_T3_)
     !! internal pressure - y-gradient
   real(real64), intent(in) :: viscosity(_T3_)
     !! viscosity
#undef _T3_

!  Local constants

!  Local variables
   logical, save :: ufirst=.false. !KB should likely not have save
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('uv_momentum_3d()',level=2)
   if (self%apply_bottom_friction) call self%bottom_friction_3d()
   if(ufirst) then
      call pk_3d(self,dt,tausx,dpdx,idpdy,viscosity)
      call self%coriolis_fpk()
      call qk_3d(self,dt,tausy,dpdy,idpdy,viscosity)
      call self%coriolis_fqk()
      ufirst = .false.
   else
      call qk_3d(self,dt,tausy,dpdy,idpdy,viscosity)
      call self%coriolis_fqk()
      call pk_3d(self,dt,tausx,dpdx,idpdy,viscosity)
      call self%coriolis_fpk()
      ufirst = .true.
   end if
   call self%w_momentum_3d(dt)
   call self%velocities_3d()
   call self%uv_advection_3d(dt)
#if 0
   call self%uv_diffusion_3d(dt) !KB - makes model go wrong
#else
if (associated(self%logs)) call self%logs%info('*** missing uv_diffusion_3d() ***',level=0)
#endif
   call self%shear_frequency(viscosity)
   call self%stresses(tausx,tausy)
END SUBROUTINE uvw_momentum_3d

!---------------------------------------------------------------------------

SUBROUTINE pk_3d(self,dt,taus,dpdx,idpdx,viscosity)
   !! solve the 3D momentum equation in the local x-direction

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
      !! viscosity
#undef _T3_

!  Local constants

!  Local variables
   integer :: i,j,k
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('pk_3d()',level=3)
   UGrid: associate( UG => self%domain%U )
   do k=UG%kmin,UG%kmax
      do j=UG%jmin,UG%jmax
         do i=UG%imin,UG%imax
            self%num(i,j,k)=0._real64
            self%ea2(i,j,k)=0._real64
            self%ea4(i,j,k)=0._real64
            if (UG%mask(i,j) == 1 .or. UG%mask(i,j) == 2) then
               ! Eddy viscosity at U-points
               self%num(i,j,k)=0.5_real64*(viscosity(i,j,k)+viscosity(i+1,j,k))
               ! Coriolis, advection/diffusion and internal pressure
!KB               self%ea4(i,j,k)=dt*UG%alpha(i,j)*(self%fqk(i,j,k)-self%uuEx(i,j,k)+ip_fac*idpdx(i,j,k))
               ! check signs for components
#ifndef _UPDATE_ADV_DIFF_
               self%ea4(i,j,k)=dt*UG%alpha(i,j)*(self%fqk(i,j,k)-self%advpk(i,j,k)-self%diffuk(i,j,k)+ip_fac*idpdx(i,j,k))
#else
               self%ea4(i,j,k)=dt*UG%alpha(i,j)*(self%fqk(i,j,k)+ip_fac*idpdx(i,j,k))
#endif
            end if
         end do
      end do
   end do

   ! Additional matrix elements for inner, surface and bottom layer
   do j=UG%jmin,UG%jmax
      do i=UG%imin,UG%imax
         if (UG%mask(i,j) == 1 .or. UG%mask(i,j) == 2) then
            ! surface stress
            k=UG%kmax
            self%ea4(i,j,k)=self%ea4(i,j,k)+dt*UG%alpha(i,j)*0.5_real64*(taus(i,j)+taus(i+1,j))/rho0
            ! external pressure
            self%ea4(i,j,1:)=self%ea4(i,j,1:)-dt*0.5_real64*(UG%ho(i,j,1:)+UG%hn(i,j,1:))*g*dpdx(i,j)
            ! bottom friction
            k=UG%kmin
            self%ea2(i,j,k)=dt*self%rru(i,j)
         end if
      end do
   end do

   !KB - could call with uk instead of doing scaling??
   do j=UG%jmin,UG%jmax
      do i=UG%imin,UG%imax
         if (UG%mask(i,j) == 1 .or. UG%mask(i,j) == 2) then
            self%pk(i,j,:)=self%pk(i,j,:)/UG%ho(i,j,:)
         end if
      end do
   end do
   call self%vertical_diffusion%calculate(dt,self%cnpar,UG%mask,UG%ho,UG%hn,self%molecular,self%num,self%pk,ea2=self%ea2,ea4=self%ea4)
   do j=UG%jmin,UG%jmax
      do i=UG%imin,UG%imax
         if (UG%mask(i,j) == 1 .or. UG%mask(i,j) == 2) then
            self%pk(i,j,:)=self%pk(i,j,:)*UG%hn(i,j,:)
         end if
      end do
   end do

   ! Transport correction: the integral of the new velocities has to
   ! be the same as the transport calculated by the external mode, Uint.
   do j=UG%jmin,UG%jmax
      do i=UG%imin,UG%imax
         if (UG%mask(i,j) == 1 .or. UG%mask(i,j) == 2) then
            self%pk(i,j,1:) = self%pk(i,j,1:)+UG%hn(i,j,1:)*(self%Ui(i,j)-sum(self%pk(i,j,1:)))/UG%D(i,j)
         end if
      end do
   end do
   end associate UGrid
END SUBROUTINE pk_3d

!---------------------------------------------------------------------------

SUBROUTINE qk_3d(self,dt,taus,dpdy,idpdy,viscosity)
   !! solve the 3D momentum equation in the local y-direction

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
   integer :: i,j,k
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('qk_3d()',level=3)
   VGrid: associate( VG => self%domain%V )
   do k=VG%kmin,VG%kmax
      do j=VG%jmin,VG%jmax
         do i=VG%imin,VG%imax
            self%num(i,j,k)=0._real64
            self%ea2(i,j,k)=0._real64
            self%ea4(i,j,k)=0._real64
            if (VG%mask(i,j) == 1 .or. VG%mask(i,j) == 2) then
               ! Eddy viscosity at V-points
               self%num(i,j,k)=0.5_real64*(viscosity(i,j,k)+viscosity(i,j+1,k))
               ! Coriolis, advection/diffusion and internal pressure
!KB               self%ea4(i,j,k)=dt*VG%alpha(i,j)*(-self%fpk(i,j,k)-self%vvEx(i,j,k)+ip_fac*idpdy(i,j,k))
               ! check signs for components
#ifndef _UPDATE_ADV_DIFF_
               self%ea4(i,j,k)=dt*VG%alpha(i,j)*(-self%fpk(i,j,k)-self%advqk(i,j,k)-self%diffvk(i,j,k)+ip_fac*idpdy(i,j,k))
#else
               self%ea4(i,j,k)=dt*VG%alpha(i,j)*(-self%fpk(i,j,k)+ip_fac*idpdy(i,j,k))
#endif
            end if
         end do
      end do
   end do

   ! Additional matrix elements for inner, surface and bottom layer
   do j=VG%jmin,VG%jmax
      do i=VG%imin,VG%imax
         if (VG%mask(i,j) == 1 .or. VG%mask(i,j) == 2) then
            ! surface stress
            k=VG%kmax
            self%ea4(i,j,k)=self%ea4(i,j,k)+dt*VG%alpha(i,j)*0.5_real64*(taus(i,j)+taus(i,j+1))/rho0
            ! external pressure
            self%ea4(i,j,1:)=self%ea4(i,j,1:)-dt*0.5_real64*(VG%ho(i,j,1:)+VG%hn(i,j,1:))*g*dpdy(i,j)
            ! bottom friction
            k=VG%kmin
!KB            self%ea2(i,j,k)=dt*self%rrv(i,j)/(0.5_real64*(VG%ho(i,j,k)+VG%hn(i,j,k)))
            self%ea2(i,j,k)=dt*self%rrv(i,j)
         end if
      end do
   end do

   do j=VG%jmin,VG%jmax
      do i=VG%imin,VG%imax
         if (VG%mask(i,j) == 1 .or. VG%mask(i,j) == 2) then
            self%qk(i,j,:)=self%qk(i,j,:)/VG%ho(i,j,:)
         end if
      end do
   end do
   call self%vertical_diffusion%calculate(dt,self%cnpar,VG%mask,VG%ho,VG%hn,self%molecular,self%num,self%qk,ea2=self%ea2,ea4=self%ea4)
   do j=VG%jmin,VG%jmax
      do i=VG%imin,VG%imax
         if (VG%mask(i,j) == 1 .or. VG%mask(i,j) == 2) then
            self%qk(i,j,:)=self%qk(i,j,:)*VG%hn(i,j,:)
         end if
      end do
   end do

   ! Transport correction: the integral of the new velocities has to
   ! be the same as the transport calculated by the external mode, Vi
   do j=VG%jmin,VG%jmax
      do i=VG%imin,VG%imax
         if (VG%mask(i,j) == 1 .or. VG%mask(i,j) == 2) then
            self%qk(i,j,1:) = self%qk(i,j,1:)+VG%hn(i,j,1:)*(self%Vi(i,j)-sum(self%qk(i,j,1:)))/VG%D(i,j)
         end if
      end do
   end do
   end associate VGrid
END SUBROUTINE qk_3d

!---------------------------------------------------------------------------

MODULE SUBROUTINE w_momentum_3d(self,dt)
   !! Solve the vertical momemtum equation - w

   IMPLICIT NONE

   class(type_getm_momentum), intent(inout) :: self
   real(real64), intent(in) :: dt
      !! timestep [s]

!  Local constants

!  Local variables
   real(real64) :: dtm1
   integer :: i,j,k
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('w_momentum_3d()',level=3)
   dtm1=1._real64/dt
   TGrid: associate( TG => self%domain%T )
   UGrid: associate( UG => self%domain%U )
   VGrid: associate( VG => self%domain%V )
   do k=TG%l(3),TG%u(3)-1
      do j=TG%l(2)+1,TG%u(2)
         do i=TG%l(1)+1,TG%u(1)
            if (TG%mask(i,j) == 1) then
                  self%ww(i,j,k) = self%ww(i,j,k-1)-(TG%hn(i,j,k)-TG%ho(i,j,k))*dtm1 &
                              -(self%pk(i,j,k)*UG%dy(i,j)-self%pk(i-1,j  ,k)*UG%dy(i-1,j) &
                               +self%qk(i,j,k)*VG%dx(i,j)-self%qk(i  ,j-1,k)*VG%dx(i,j-1)) &
                               *TG%inv_area(i,j)
            end if
         end do
      end do
   end do
   end associate VGrid
   end associate UGrid
   end associate TGrid
END SUBROUTINE w_momentum_3d

!---------------------------------------------------------------------------

END SUBMODULE momentum_3d_smod
