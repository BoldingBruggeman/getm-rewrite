! Copyright (C) 2020 Bolding & Bruggeman

!!{!./code/velocity_diffusion.md!}

SUBMODULE (getm_momentum) diffusion_smod

CONTAINS

!---------------------------------------------------------------------------

module SUBROUTINE diffusion_2d(self)

   !! Velocity shear

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: rc
!---------------------------------------------------------------------------
   call self%logs%info('diffusion_2d()',level=2)

   return
END SUBROUTINE diffusion_2d

!---------------------------------------------------------------------------

module SUBROUTINE diffusion_3d(self)

   !! Velocity shear

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: rc
!---------------------------------------------------------------------------
   call self%logs%info('diffusion_2d()',level=2)

   return
END SUBROUTINE diffusion_3d

!---------------------------------------------------------------------------

module SUBROUTINE diffusion_driver(self,An_method,UEx,VEx,U,V,D,DU,DV,hsd_u,hsd_v)
   !! Driver routine for velocity diffusion
   !! @note
   !! this routine should likely be in a submodule
   !! @endnote

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
   integer, intent(in) :: An_method
   real(real64), dimension(:,:), intent(in),optional  :: U,V,D,DU,DV
   real(real64), dimension(:,:), intent(inout)        :: UEx,VEx
   real(real64), dimension(:,:), intent(out),optional :: hsd_u,hsd_v

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
      hsd_u = UEx
   end if
   if (present(hsd_v)) then
      hsd_v = VEx
   end if
#endif

   ! Central for dx(2*Am*dx(U/DU))
#define TG self%domain%T
#define UG self%domain%U
   do j=TG%jmin,TG%jmax
      do i=TG%imin,TG%imax+1 ! work2d defined on T-points
         work2d(i,j)=0._real64
         if (TG%mask(i,j) .eq. 1) then
            if(use_Am) then
               work2d(i,j)=2._real64*Am*TG%dy(i,j)*D(i,j)               &
                       *(U(i,j)/DU(i,j)-U(i-1,j)/DU(i-1,j))/TG%dx(i,j)
            end if
            if (An_method .gt. 0) then
               work2d(i,j)=work2d(i,j)+An(i,j)*TG%dy(i,j)*(U(i,j)-U(i-1,j))/TG%dx(i,j)
            end if
         end if
      end do
      do i=UG%imin,UG%imax      ! UEx defined on U-points
         if (UG%mask(i,j).eq.1 .or. UG%mask(i,j).eq.2) then
!KB            UEx(i,j)=UEx(i,j)-(work2d(i+1,j)-work2d(i  ,j))*UG%mask%inv_area(i,j)
         end if
      end do
   end do
#undef TG

   ! Central for dy(Am*(dy(U/DU)+dx(V/DV)))
#define XG self%domain%X
   do j=XG%jmin-1,XG%jmax ! work2d defined on X-points
      do i=XG%imin,XG%imax
         work2d(i,j)=0._real64
         if (XG%mask(i,j) .ge. 1) then
            if(use_Am) then
               work2d(i,j)=Am*0.5_real64*(DU(i,j)+DU(i,j+1))*XG%dx(i,j)  &
                       *((U(i,j+1)/DU(i,j+1)-U(i,j)/DU(i,j))/XG%dy(i,j) &
                        +(V(i+1,j)/DV(i+1,j)-V(i,j)/DV(i,j))/XG%dx(i,j) )
            end if
            if (An_method .gt. 0) then
               work2d(i,j)=work2d(i,j)+AnX(i,j)*(U(i,j+1)-U(i,j))*XG%dx(i,j)/XG%dy(i,j)
            end if
         end if
      end do
   end do
   do j=UG%jmin,UG%jmax !UEx defined on U-points
      do i=UG%imin,UG%imax
         if (UG%mask(i,j).eq.1 .or. UG%mask(i,j).eq.2) then
!KB            UEx(i,j)=UEx(i,j)-(work2d(i,j  )-work2d(i,j-1))*UG%mask%inv_area(i,j)
         end if
      end do
   end do
#undef XG
#undef UG

   ! Central for dx(Am*(dy(U/DU)+dx(V/DV)))
#define XG self%domain%X
#define VG self%domain%V
   do j=XG%jmin,XG%jmax
      do i=XG%imin-1,XG%imax ! work2d defined on X-points
         work2d(i,j)=0._real64
         if (XG%mask(i,j) .ge. 1) then
            if(use_Am) then
               work2d(i,j)=Am*0.5_real64*(DV(i,j)+DV(i+1,j))*XG%dy(i,j) &
                       *((U(i,j+1)/DU(i,j+1)-U(i,j)/DU(i,j))/XG%dy(i,j) &
                        +(V(i+1,j)/DV(i+1,j)-V(i,j)/DV(i,j))/XG%dx(i,j) )
            end if
            if (An_method .gt. 0) then
               work2d(i,j)=work2d(i,j)+AnX(i,j)*(V(i+1,j)-V(i,j))*XG%dy(i,j)/XG%dx(i,j)
            end if
         end if
      end do
      !!!!!!!!!!
      do i=VG%imin,VG%imax ! VEx defined on V-points
         if (VG%mask(i,j).eq.1 .or. VG%mask(i,j).eq.2) then
!KB            VEx(i,j)=VEx(i,j)-(work2d(i  ,j)-work2d(i-1,j))*VG%mask%inv_area(i,j)
         end if
      end do
   end do
#undef XG

   ! Central for dy(2*Am*dy(V/DV))
#define TG self%domain%T
#define VG self%domain%V
   do j=TG%jmin,TG%jmax+1 ! work2d defined on T-points
      do i=TG%imin,TG%imax
         work2d(i,j)=0._real64
         if (TG%mask(i,j) .eq. 1) then
            if(use_Am) then
               work2d(i,j)=2._real64*Am*TG%dx(i,j)*D(i,j)               &
                       *(V(i,j)/DV(i,j)-V(i,j-1)/DV(i,j-1))/TG%dy(i,j)
            end if
            if (An_method .gt. 0) then
               work2d(i,j)=work2d(i,j)+An(i,j)*TG%dx(i,j)*(V(i,j)-V(i,j-1))/TG%dy(i,j)
            end if
         end if
      end do
   end do
#undef TG
   do j=VG%jmin,VG%jmax ! VEx defined on V-points
      do i=VG%imin,VG%imax
         if (VG%mask(i,j).eq.1 .or. VG%mask(i,j).eq.2) then
!KB            VEx(i,j)=VEx(i,j)-(work2d(i,j+1)-work2d(i,j  ))*VG%inv_area(i,j)
         end if
      end do
   end do
#undef VG

#if 0
   if (present(hsd_u)) then
      hsd_u = UEx - hsd_u
   end if
   if (present(hsd_v)) then
      hsd_v = VEx - hsd_v
   end if
#endif
   return
END SUBROUTINE diffusion_driver


!---------------------------------------------------------------------------

END SUBMODULE diffusion_smod
