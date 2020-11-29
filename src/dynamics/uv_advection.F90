! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

SUBMODULE (getm_momentum) uv_advection_smod

CONTAINS

!KB      real(real64), dimension(:,:), allocatable :: Ua,Va

!---------------------------------------------------------------------------

MODULE SUBROUTINE uv_advection_2d(self,dt)

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
   real(real64), intent(in) :: dt
      !! timestep [s]

!  Local constants

!  Local variables
   integer :: i,j
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('uv_advection_2d()',level=2)
   XGrid: associate( XG => self%domain%X )
   TGrid: associate( TG => self%domain%T )
   UGrid: associate( UG => self%domain%U )
   VGrid: associate( VG => self%domain%V )
   do j=UG%jmin,UG%jmax
      do i=UG%imin,UG%imax
         self%uadvgrid%D(i,j)  = TG%D(i+1,j)
         self%vadvgrid%D(i,j)  = XG%D(i,j)
         self%Ua(i,j) = 0.5_real64*(self%U(i,j) + self%U(i+1,j))
         self%Va(i,j) = 0.5_real64*(self%V(i,j) + self%V(i+1,j))
      end do
   end do
   if (self%store_advection) self%advU=self%U
   call self%advection%calculate(self%advection_scheme,self%uadvgrid,self%Ua,self%vadvgrid,self%Va,dt,UG,self%U)
   if (self%store_advection) self%advU=(self%U-self%advU)/dt
   do j=UG%jmin,UG%jmax
      do i=UG%imin,UG%imax
         self%uadvgrid%D(i,j)  = XG%D(i,j)
         self%vadvgrid%D(i,j)  = TG%D(i,j+1)
         self%Ua(i,j) = 0.5_real64*(self%U(i,j) + self%U(i,j+1))
         self%Va(i,j) = 0.5_real64*(self%V(i,j) + self%V(i,j+1))
      end do
   end do
   if (self%store_advection) self%advV=self%V
   call self%advection%calculate(self%advection_scheme,self%uadvgrid,self%Ua,self%vadvgrid,self%Va,dt,VG,self%V)
   if (self%store_advection) self%advV=(self%V-self%advV)/dt
   end associate VGrid
   end associate UGrid
   end associate TGrid
   end associate XGrid
END SUBROUTINE uv_advection_2d

!---------------------------------------------------------------------------

MODULE SUBROUTINE uivi_advection_2d(self,dt)

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
   real(real64), intent(in) :: dt
      !! timestep [s]

!  Local constants

!  Local variables
   integer :: i,j
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('uivi_advection_2d()',level=2)
   XGrid: associate( XG => self%domain%X )
   TGrid: associate( TG => self%domain%T )
   UGrid: associate( UG => self%domain%U )
   VGrid: associate( VG => self%domain%V )
   do j=UG%jmin,UG%jmax
      do i=UG%imin,UG%imax
         self%uadvgrid%D(i,j)  = TG%H(i+1,j)+TG%ssen(i+1,j) !KB TG%D(i+1,j)
         self%vadvgrid%D(i,j)  = XG%D(i,j) ! Knut
         self%Ua(i,j) = 0.5_real64*(self%Ui(i,j) + self%Ui(i+1,j))
         self%Va(i,j) = 0.5_real64*(self%Vi(i,j) + self%Vi(i+1,j))
      end do
   end do
   self%SxA=self%Ui
   call self%advection%calculate(self%advection_scheme,self%uadvgrid,self%Ua,self%vadvgrid,self%Va,dt,UG,self%Ui)
   self%SxA=(self%Ui-self%SxA)/dt
   do j=UG%jmin,UG%jmax
      do i=UG%imin,UG%imax
         self%uadvgrid%D(i,j)  = XG%D(i,j) ! Knut
         self%vadvgrid%D(i,j)  = TG%H(i,j+1)+TG%ssen(i,j+1) !KB TG%D(i,j+1)
         self%Ua(i,j) = 0.5_real64*(self%Ui(i,j) + self%Ui(i,j+1))
         self%Va(i,j) = 0.5_real64*(self%Vi(i,j) + self%Vi(i,j+1))
      end do
   end do
   self%SyA=self%Vi
   call self%advection%calculate(self%advection_scheme,self%uadvgrid,self%Ua,self%vadvgrid,self%Va,dt,VG,self%Vi)
   self%SyA=(self%Vi-self%SyA)/dt
   end associate VGrid
   end associate UGrid
   end associate TGrid
   end associate XGrid
END SUBROUTINE uivi_advection_2d

!---------------------------------------------------------------------------

MODULE SUBROUTINE uv_advection_3d(self,dt)
   !! 3D velocity advection

   IMPLICIT NONE

   class(type_getm_momentum), intent(inout) :: self
   real(real64), intent(in) :: dt
      !! timestep [s]

!  Local constants

!  Local variables
   integer :: i,j,k
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('uv_advection_3d()',level=2)
   XGrid: associate( XG => self%domain%X )
   TGrid: associate( TG => self%domain%T )
   UGrid: associate( UG => self%domain%U )
   VGrid: associate( VG => self%domain%V )
   do j=UG%jmin,UG%jmax
      do i=UG%imin,UG%imax
         self%uadvgrid%hn(i,j,:) = TG%hn(i+1,j,:)
         self%vadvgrid%hn(i,j,:) = XG%hn(i,j,:)
         self%pka(i,j,:) = 0.5_real64*(self%pk(i,j,:) + self%pk(i+1,j,:))
         self%qka(i,j,:) = 0.5_real64*(self%qk(i,j,:) + self%qk(i+1,j,:))
      end do
   end do
   call self%advection%calculate(self%advection_scheme,self%uadvgrid,self%pka,self%vadvgrid,self%qka,dt,UG,self%pk)
   do j=UG%jmin,UG%jmax
      do i=UG%imin,UG%imax
         self%uadvgrid%hn(i,j,:) = XG%hn(i,j,:)
         self%vadvgrid%hn(i,j,:) = TG%hn(i,j+1,:)
         self%pka(i,j,:) = 0.5_real64*(self%pk(i,j,:) + self%pk(i,j+1,:))
         self%qka(i,j,:) = 0.5_real64*(self%qk(i,j,:) + self%qk(i,j+1,:))
      end do
   end do
   call self%advection%calculate(self%advection_scheme,self%uadvgrid,self%pka,self%vadvgrid,self%qka,dt,VG,self%qk)
   end associate VGrid
   end associate UGrid
   end associate TGrid
   end associate XGrid
END SUBROUTINE uv_advection_3d

!---------------------------------------------------------------------------

END SUBMODULE uv_advection_smod
