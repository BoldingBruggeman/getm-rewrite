! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!!{!./code/velocity_diffusion.md!}

SUBMODULE (getm_momentum) uv_diffusion_smod

CONTAINS

!---------------------------------------------------------------------------

MODULE SUBROUTINE uv_diffusion_2d(self,dt)

   !! Velocity shear

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
   real(real64), intent(in) :: dt
      !! timestep [s]

!  Local constants

!  Local variables
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('uv_diffusion_2d()',level=3)
   TGrid: associate( TG => self%domain%T )
   UGrid: associate( UG => self%domain%U )
   VGrid: associate( VG => self%domain%V )
   if (self%Am0 > 0._real64) then !KB - another check is required
      where(UG%mask > 0) self%u1 = self%U/UG%D
      where(VG%mask > 0) self%v1 = self%V/VG%D
      call diffusion_driver(self,self%u1,self%v1,TG%D,UG%D,VG%D,self%diffu1,self%diffv1)
   end if
   if (self%An_method > 0) then
      call numerical_damping(self,self%U,self%V)
   end if
#ifdef _APPLY_ADV_DIFF_
   self%U=dt*(self%diffu1+self%dampU)
   self%V=dt*(self%diffv1+self%dampV)
#endif
   end associate VGrid
   end associate UGrid
   end associate TGrid
END SUBROUTINE uv_diffusion_2d

!---------------------------------------------------------------------------

MODULE SUBROUTINE uv_diffusion_3d(self,dt)

   !! Velocity shear

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
   real(real64), intent(in) :: dt
      !! timestep [s]

!  Local constants

!  Local variables
   integer :: k
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('diffusion_3d()',level=3)
   TGrid: associate( TG => self%domain%T )
   UGrid: associate( UG => self%domain%U )
   VGrid: associate( VG => self%domain%V )
   if (self%Am0 > 0._real64) then
      where(UG%mask > 0) self%uk(:,:,k) = self%pk(:,:,k)/UG%hn(:,:,k)
      where(VG%mask > 0) self%vk(:,:,k) = self%qk(:,:,k)/VG%hn(:,:,k)
      do k=1,TG%kmax
         call diffusion_driver(self,self%vk(:,:,k),self%vk(:,:,k), &
                               TG%hn(:,:,k),UG%hn(:,:,k),VG%hn(:,:,k), &
                               self%diffuk(:,:,k),self%diffvk(:,:,k))
      end do
#ifdef _APPLY_ADV_DIFF_
      self%pk=dt*self%diffuk
      self%qk=dt*self%diffvk
#endif
   end if
   end associate VGrid
   end associate UGrid
   end associate TGrid

END SUBROUTINE uv_diffusion_3d

!---------------------------------------------------------------------------

MODULE SUBROUTINE slow_diffusion(self)

   !! Diffusion of time averaged depth integrated transports and slow diffusion update

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: i,j
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('slow_diffusion()',level=3)
   TGrid: associate( TG => self%domain%T )
   UGrid: associate( UG => self%domain%U )
   VGrid: associate( VG => self%domain%V )
   if (self%Am0 > 0._real64) then
      where(UG%mask > 0) self%u1 = self%U/UG%D
      where(VG%mask > 0) self%v1 = self%V/VG%D
      call diffusion_driver(self,self%u1,self%v1,TG%D,UG%D,VG%D,self%SxD,self%SyD)
   end if
   if (self%An_method > 0) then
      call numerical_damping(self,self%Ui,self%Vi)
   else
      self%dampU=0._real64
      self%dampV=0._real64
   end if

   do j=UG%jmin,UG%jmax
      do i=UG%imin,UG%imax
         if (UG%mask(i,j) .ge. 1) then
            self%SxD(i,j)=SUM(self%diffuk(i,j,1:))-self%SxD(i,j)-self%dampU(i,j)
         end if
      end do
   end do
   do j=VG%jmin,VG%jmax
      do i=VG%imin,VG%imax
         if (VG%mask(i,j) .ge. 1) then
            self%SyD(i,j)=SUM(self%diffvk(i,j,1:))-self%SyD(i,j)-self%dampv(i,j)
         end if
      end do
   end do
   end associate VGrid
   end associate UGrid
   end associate TGrid
END SUBROUTINE slow_diffusion

!---------------------------------------------------------------------------

MODULE SUBROUTINE numerical_damping(self,U,V)

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
#define _U2_ self%domain%U%l(1):,self%domain%U%l(2):
#define _V2_ self%domain%V%l(1):,self%domain%V%l(2):
   real(real64), dimension(:,:), intent(in) :: U(_U2_)
   real(real64), dimension(:,:), intent(in) :: v(_V2_)
#undef _U2_
#undef _V2_

!  Local constants

!  Local variables
   integer :: i,j
!---------------------------------------------------------------------------
   call self%logs%info('numerical_damping()',level=2)
   TGrid: associate( TG => self%domain%T )
   XGrid: associate( XG => self%domain%X )
   UGrid: associate( UG => self%domain%U )
   ! Central for dx(2*Am*dx(U/DU))
   do j=TG%jmin,TG%jmax
      do i=TG%imin,TG%imax+1 ! work2d defined on T-points
         self%work2d(i,j)=0._real64
         if (TG%mask(i,j) == 1) then
            self%work2d(i,j)=self%An(i,j)*TG%dy(i,j)*(U(i,j)-U(i-1,j))/TG%dx(i,j)
         end if
      end do
   end do
   do j=UG%jmin,UG%jmax ! dampU defined on U-points
      do i=UG%imin,UG%imax
         self%dampU(i,j)=0._real64
         if (UG%mask(i,j) == 1 .or. UG%mask(i,j) == 2) then
            self%dampU(i,j)=-(self%work2d(i+1,j)-self%work2d(i,j))*UG%inv_area(i,j)
         end if
      end do
   end do

   ! Central for dy(Am*(dy(U/DU)+dx(V/DV)))
   do j=XG%jmin-1,XG%jmax ! work2d defined on X-points
      do i=XG%imin,XG%imax
         self%work2d(i,j)=0._real64
! XGrids must be fixed
#if 0
         if (XG%mask(i,j) > 0) then
            self%work2d(i,j)=AnX(i,j)*(U(i,j+1)-U(i,j))*XG%dx(i,j)/XG%dy(i,j)
         end if
#endif
      end do
   end do
   do j=UG%jmin,UG%jmax ! dampU defined on U-points
      do i=UG%imin,UG%imax
         if (UG%mask(i,j) == 1 .or. UG%mask(i,j) == 2) then
            self%dampU(i,j)=self%dampU(i,j)-(self%work2d(i,j)-self%work2d(i,j-1))*UG%inv_area(i,j)
         end if
      end do
   end do
   end associate UGrid

   ! Central for dx(Am*(dy(U/DU)+dx(V/DV)))
   VGrid: associate( VG => self%domain%V )
   do j=XG%jmin,XG%jmax
      do i=XG%imin-1,XG%imax ! work2d defined on X-points
         self%work2d(i,j)=0._real64
! XGrids must be fixed
#if 0
         if (XG%mask(i,j) > 0) then
            self%work2d(i,j)=AnX(i,j)*(V(i+1,j)-V(i,j))*XG%dy(i,j)/XG%dx(i,j)
         end if
#endif
      end do
   end do
   do j=VG%jmin,VG%jmax ! dampV defined on V-points
      do i=VG%imin,VG%imax
         self%dampV(i,j)=0._real64
         if (VG%mask(i,j) == 1 .or. VG%mask(i,j) == 2) then
            self%dampV(i,j)=-(self%work2d(i,j)-self%work2d(i-1,j))*VG%inv_area(i,j)
         end if
      end do
   end do

   ! Central for dy(2*Am*dy(V/DV))
   do j=TG%jmin,TG%jmax+1 ! work2d defined on T-points
      do i=TG%imin,TG%imax
         self%work2d(i,j)=0._real64
         if (TG%mask(i,j)  ==  1) then
            self%work2d(i,j)=self%An(i,j)*TG%dx(i,j)*(V(i,j)-V(i,j-1))/TG%dy(i,j)
         end if
      end do
   end do
   do j=VG%jmin,VG%jmax ! dampV defined on V-points
      do i=VG%imin,VG%imax
         if (VG%mask(i,j) == 1 .or. VG%mask(i,j) == 2) then
            self%dampV(i,j)=self%dampV(i,j)-(self%work2d(i,j+1)-self%work2d(i,j))*VG%inv_area(i,j)
         end if
      end do
   end do
   end associate VGrid
   end associate XGrid
   end associate TGrid
END SUBROUTINE numerical_damping

!---------------------------------------------------------------------------

MODULE SUBROUTINE diffusion_driver(self,u,v,D,DU,DV,diffu,diffv)
   !! Driver routine for velocity diffusion
   !! @note
   !! this routine should likely be in a submodule
   !! @endnote

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
#define _T2_ self%domain%T%l(1):,self%domain%T%l(2):
#define _U2_ self%domain%U%l(1):,self%domain%U%l(2):
#define _V2_ self%domain%V%l(1):,self%domain%V%l(2):
   real(real64), dimension(:,:), intent(in) :: u(_U2_)
   real(real64), dimension(:,:), intent(in) :: v(_V2_)
   real(real64), dimension(:,:), intent(in) :: D(_T2_)
   real(real64), dimension(:,:), intent(in) :: DU(_U2_)
   real(real64), dimension(:,:), intent(in) :: DV(_V2_)
   real(real64), dimension(:,:), intent(inout) :: diffu(_U2_)
   real(real64), dimension(:,:), intent(inout) :: diffv(_V2_)
#undef _V2_
#undef _U2_
#undef _T2_

!  Local constants

!  Local variables
   integer :: i,j
!---------------------------------------------------------------------------
   call self%logs%info('diffusion_driver()',level=2)
   TGrid: associate( TG => self%domain%T )
   XGrid: associate( XG => self%domain%X )
   UGrid: associate( UG => self%domain%U )
   ! Central for dx(2*Am*dx(U/DU))
   do j=TG%jmin,TG%jmax
      do i=TG%imin,TG%imax+1 ! work2d defined on T-points
         self%work2d(i,j)=0._real64
         if (TG%mask(i,j) == 1) then
            self%work2d(i,j)=2._real64*self%Am(i,j)*TG%dy(i,j)*D(i,j)*(u(i,j)-u(i-1,j))/TG%dx(i,j)
         end if
      end do
      do i=UG%imin,UG%imax ! diffu defined on U-points
         if (UG%mask(i,j) == 1 .or. UG%mask(i,j) == 2) then
            diffu(i,j)=-(self%work2d(i+1,j)-self%work2d(i,j))*UG%inv_area(i,j)
         end if
      end do
   end do

   ! Central for dy(Am*(dy(U/DU)+dx(V/DV)))
   do j=XG%jmin-1,XG%jmax ! work2d defined on X-points
      do i=XG%imin,XG%imax
         self%work2d(i,j)=0._real64
! XGrids must be fixed
#if 0
         if (XG%mask(i,j) > 0) then
            self%work2d(i,j)=self%Am(i,j)*0.5_real64*(DU(i,j)+DU(i,j+1))*XG%dx(i,j)  &
                            *((u(i,j+1)-u(i,j))/XG%dy(i,j)+(v(i+1,j)-v(i,j))/XG%dx(i,j) )
         end if
#endif
      end do
   end do
   do j=UG%jmin,UG%jmax ! diffu defined on U-points
      do i=UG%imin,UG%imax
         if (UG%mask(i,j) == 1 .or. UG%mask(i,j) == 2) then
            diffu(i,j)=diffu(i,j)-(self%work2d(i,j)-self%work2d(i,j-1))*UG%inv_area(i,j)
         end if
      end do
   end do
   end associate UGrid

   ! Central for dx(Am*(dy(U/DU)+dx(V/DV)))
   VGrid: associate( VG => self%domain%V )
   do j=XG%jmin,XG%jmax
      do i=XG%imin-1,XG%imax ! work2d defined on X-points
         self%work2d(i,j)=0._real64
! XGrids must be fixed
#if 0
         if (XG%mask(i,j) > 0) then
            self%work2d(i,j)=self%Am(i,j)*0.5_real64*(DV(i,j)+DV(i+1,j))*XG%dy(i,j) &
                            *((u(i,j+1)-u(i,j))/XG%dy(i,j)+(v(i+1,j)-v(i,j))/XG%dx(i,j) )
         end if
#endif
      end do
   end do
   do j=VG%jmin,VG%jmax ! diffv defined on V-points
      do i=VG%imin,VG%imax ! diffv defined on V-points
         if (VG%mask(i,j) == 1 .or. VG%mask(i,j) == 2) then
            diffv(i,j)=-(self%work2d(i,j)-self%work2d(i-1,j))*VG%inv_area(i,j)
         end if
      end do
   end do

   ! Central for dy(2*Am*dy(V/DV))
   do j=TG%jmin,TG%jmax+1 ! work2d defined on T-points
      do i=TG%imin,TG%imax
         self%work2d(i,j)=0._real64
         if (TG%mask(i,j) == 1) then
            self%work2d(i,j)=2._real64*self%Am(i,j)*TG%dx(i,j)*XG%D(i,j)*(V(i,j)-V(i,j-1))/TG%dy(i,j)
         end if
      end do
   end do
   do j=VG%jmin,VG%jmax ! diffv defined on V-points
      do i=VG%imin,VG%imax
         if (VG%mask(i,j) == 1 .or. VG%mask(i,j) == 2) then
            diffv(i,j)=diffv(i,j)-(self%work2d(i,j+1)-self%work2d(i,j))*VG%inv_area(i,j)
         end if
      end do
   end do
   end associate VGrid
   end associate XGrid
   end associate TGrid
END SUBROUTINE diffusion_driver

!---------------------------------------------------------------------------

END SUBMODULE uv_diffusion_smod
