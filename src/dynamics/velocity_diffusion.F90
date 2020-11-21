! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!!{!./code/velocity_diffusion.md!}

SUBMODULE (getm_momentum) diffusion_smod

CONTAINS

!---------------------------------------------------------------------------

module SUBROUTINE uv_diffusion_2d(self)

   !! Velocity shear

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self

!  Local constants

!  Local variables
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('diffusion_2d()',level=2)
   TGrid: associate( TG => self%domain%T )
   UGrid: associate( UG => self%domain%U )
   VGrid: associate( VG => self%domain%V )
   if (self%Am > 0._real64 .or. self%An_method > 0) then
      call diffusion_driver(self,self%U,self%V,TG%D,UG%D,VG%D,self%diffu1,self%diffv1)
   end if
   end associate VGrid
   end associate UGrid
   end associate TGrid
END SUBROUTINE uv_diffusion_2d

!---------------------------------------------------------------------------

module SUBROUTINE uivi_diffusion_2d(self)

   !! Velocity shear

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self

!  Local constants

!  Local variables
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('diffusion_2d()',level=2)
   TGrid: associate( TG => self%domain%T )
   UGrid: associate( UG => self%domain%U )
   VGrid: associate( VG => self%domain%V )
   if (self%Am > 0._real64 .or. self%An_method > 0) then
      call diffusion_driver(self,self%Ui,self%Vi,TG%D,UG%D,VG%D,self%SxD,self%SyD)
   end if
   end associate VGrid
   end associate UGrid
   end associate TGrid
END SUBROUTINE uivi_diffusion_2d

!---------------------------------------------------------------------------

module SUBROUTINE uv_diffusion_3d(self)

   !! Velocity shear

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: k
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('diffusion_2d()',level=2)
   TGrid: associate( TG => self%domain%T )
   UGrid: associate( UG => self%domain%U )
   VGrid: associate( VG => self%domain%V )
   if (self%Am > 0._real64 .or. self%An_method > 0) then
      do k=1,TG%kmax
         call diffusion_driver(self,self%pk(:,:,k),self%qk(:,:,k), &
                               TG%hn(:,:,k),UG%hn(:,:,k),VG%hn(:,:,k), &
                               self%diffuk(:,:,k),self%diffvk(:,:,k))
      end do
   end if
   end associate VGrid
   end associate UGrid
   end associate TGrid

END SUBROUTINE uv_diffusion_3d

!---------------------------------------------------------------------------

module SUBROUTINE diffusion_driver(self,U,V,D,DU,DV,diffu,diffv,hsd_u,hsd_v)
   !! Driver routine for velocity diffusion
   !! @note
   !! this routine should likely be in a submodule
   !! @endnote

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
#define _T2_ self%domain%T%l(1):,self%domain%T%l(2):
#define _U2_ self%domain%T%l(1):,self%domain%T%l(2):
#define _V2_ self%domain%T%l(1):,self%domain%T%l(2):
   real(real64), dimension(:,:), intent(in) :: U(_U2_),V(_V2_),D(_T2_),DU(_U2_),DV(_V2_)
   real(real64), dimension(:,:), intent(inout) :: diffu(_U2_),diffv(_V2_)
   real(real64), dimension(:,:), intent(out),optional :: hsd_u(_U2_),hsd_v(_V2_)
#undef _V2_
#undef _U2_
#undef _T2_

!  Local constants

!  Local variables
   real(real64), dimension(:,:), allocatable :: work2d,Anx,Any,An
   real(real64) :: Am !!!!!!!
   logical :: use_Am
   integer :: i,j
!---------------------------------------------------------------------------
   call self%logs%info('diffusion_driver()',level=2)

#if 0
   use_Am = (Am .gt. _ZERO_)

   if (present(hsd_u)) then
      hsd_u = diffu
   end if
   if (present(hsd_v)) then
      hsd_v = diffv
   end if
#endif

   TGrid: associate( TG => self%domain%T )
   XGrid: associate( XG => self%domain%X )
   UGrid: associate( UG => self%domain%U )
   ! Central for dx(2*Am*dx(U/DU))
   do j=TG%jmin,TG%jmax
      do i=TG%imin,TG%imax+1 ! work2d defined on T-points
         work2d(i,j)=0._real64
         if (TG%mask(i,j) .eq. 1) then
            if(use_Am) then
               work2d(i,j)=2._real64*Am*TG%dy(i,j)*D(i,j)               &
                       *(U(i,j)/DU(i,j)-U(i-1,j)/DU(i-1,j))/TG%dx(i,j)
            end if
            if (self%An_method .gt. 0) then
               work2d(i,j)=work2d(i,j)+An(i,j)*TG%dy(i,j)*(U(i,j)-U(i-1,j))/TG%dx(i,j)
            end if
         end if
      end do
      do i=UG%imin,UG%imax ! diffu defined on U-points
         if (UG%mask(i,j).eq.1 .or. UG%mask(i,j).eq.2) then
!KB            diffu(i,j)=diffu(i,j)-(work2d(i+1,j)-work2d(i  ,j))*UG%mask%inv_area(i,j)
         end if
      end do
   end do

   ! Central for dy(Am*(dy(U/DU)+dx(V/DV)))
   do j=XG%jmin-1,XG%jmax ! work2d defined on X-points
      do i=XG%imin,XG%imax
         work2d(i,j)=0._real64
         if (XG%mask(i,j) .ge. 1) then
            if(use_Am) then
               work2d(i,j)=self%Am*0.5_real64*(DU(i,j)+DU(i,j+1))*XG%dx(i,j)  &
                       *((U(i,j+1)/DU(i,j+1)-U(i,j)/DU(i,j))/XG%dy(i,j) &
                        +(V(i+1,j)/DV(i+1,j)-V(i,j)/DV(i,j))/XG%dx(i,j) )
            end if
            if (self%An_method .gt. 0) then
               work2d(i,j)=work2d(i,j)+AnX(i,j)*(U(i,j+1)-U(i,j))*XG%dx(i,j)/XG%dy(i,j)
            end if
         end if
      end do
   end do
   do j=UG%jmin,UG%jmax ! diffu defined on U-points
      do i=UG%imin,UG%imax
         if (UG%mask(i,j).eq.1 .or. UG%mask(i,j).eq.2) then
!KB            diffu(i,j)=diffu(i,j)-(work2d(i,j  )-work2d(i,j-1))*UG%mask%inv_area(i,j)
         end if
      end do
   end do
   end associate UGrid

   ! Central for dx(Am*(dy(U/DU)+dx(V/DV)))
   VGrid: associate( VG => self%domain%V )
   do j=XG%jmin,XG%jmax
      do i=XG%imin-1,XG%imax ! work2d defined on X-points
         work2d(i,j)=0._real64
         if (XG%mask(i,j) .ge. 1) then
            if(use_Am) then
               work2d(i,j)=Am*0.5_real64*(DV(i,j)+DV(i+1,j))*XG%dy(i,j) &
                       *((U(i,j+1)/DU(i,j+1)-U(i,j)/DU(i,j))/XG%dy(i,j) &
                        +(V(i+1,j)/DV(i+1,j)-V(i,j)/DV(i,j))/XG%dx(i,j) )
            end if
            if (self%An_method .gt. 0) then
               work2d(i,j)=work2d(i,j)+AnX(i,j)*(V(i+1,j)-V(i,j))*XG%dy(i,j)/XG%dx(i,j)
            end if
         end if
      end do
      !!!!!!!!!!
      do i=VG%imin,VG%imax ! diffv defined on V-points
         if (VG%mask(i,j).eq.1 .or. VG%mask(i,j).eq.2) then
!KB            diffv(i,j)=diffv(i,j)-(work2d(i  ,j)-work2d(i-1,j))*VG%mask%inv_area(i,j)
         end if
      end do
   end do

   ! Central for dy(2*Am*dy(V/DV))
   do j=TG%jmin,TG%jmax+1 ! work2d defined on T-points
      do i=TG%imin,TG%imax
         work2d(i,j)=0._real64
         if (TG%mask(i,j) .eq. 1) then
            if(use_Am) then
               work2d(i,j)=2._real64*Am*TG%dx(i,j)*D(i,j)               &
                       *(V(i,j)/DV(i,j)-V(i,j-1)/DV(i,j-1))/TG%dy(i,j)
            end if
            if (self%An_method .gt. 0) then
               work2d(i,j)=work2d(i,j)+An(i,j)*TG%dx(i,j)*(V(i,j)-V(i,j-1))/TG%dy(i,j)
            end if
         end if
      end do
   end do
   do j=VG%jmin,VG%jmax ! diffv defined on V-points
      do i=VG%imin,VG%imax
         if (VG%mask(i,j).eq.1 .or. VG%mask(i,j).eq.2) then
!KB            diffv(i,j)=diffv(i,j)-(work2d(i,j+1)-work2d(i,j  ))*VG%inv_area(i,j)
         end if
      end do
   end do
   end associate VGrid
   end associate XGrid
   end associate TGrid

#if 0
   if (present(hsd_u)) then
      hsd_u = diffu - hsd_u
   end if
   if (present(hsd_v)) then
      hsd_v = diffv - hsd_v
   end if
#endif
   return
END SUBROUTINE diffusion_driver

!---------------------------------------------------------------------------

END SUBMODULE diffusion_smod
