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
   if (self%advection_scheme == 0) return
   if (associated(self%logs)) call self%logs%info('uv_advection_2d()',level=3)
   XGrid: associate( XG => self%domain%X )
   TGrid: associate( TG => self%domain%T )
   UGrid: associate( UG => self%domain%U )
   VGrid: associate( VG => self%domain%V )
   do j=UG%jmin-1,UG%jmax
      do i=UG%imin-1,UG%imax
         self%uuadvgrid%D(i,j)  = TG%D(i+1,j)
         self%uvadvgrid%D(i,j)  = XG%D(i,j)
         self%ua(i,j) = 0.5_real64*(self%U(i,j) + self%U(i+1,j)) / self%uuadvgrid%D(i,j)
         self%va(i,j) = 0.5_real64*(self%V(i,j) + self%V(i+1,j)) / self%uvadvgrid%D(i,j)
      end do
   end do
#ifndef _APPLY_ADV_DIFF_
   self%advU=self%U
#endif
   where(UG%mask > 0) self%u1 = self%U/UG%D
   call self%advection%calculate(self%uuadvgrid,self%ua,self%uvadvgrid,self%va,dt,UG,self%u1)
#ifdef _APPLY_ADV_DIFF_
   where(UG%mask > 0) self%U = self%u1*UG%D
#else
   where(UG%mask > 0) self%advU=(self%u1*UG%D-self%advU)/dt
#endif

   do j=VG%jmin-1,VG%jmax
      do i=VG%imin-1,VG%imax
         self%vuadvgrid%D(i,j)  = XG%D(i,j)
         self%vvadvgrid%D(i,j)  = TG%D(i,j+1)
         self%ua(i,j) = 0.5_real64*(self%U(i,j) + self%U(i,j+1)) / self%vuadvgrid%D(i,j)
         self%va(i,j) = 0.5_real64*(self%V(i,j) + self%V(i,j+1)) / self%vvadvgrid%D(i,j)
      end do
   end do
#ifndef _APPLY_ADV_DIFF_
   self%advV=self%V
#endif
   where(VG%mask > 0) self%v1 = self%V/VG%D
   call self%advection%calculate(self%vuadvgrid,self%ua,self%vvadvgrid,self%va,dt,VG,self%v1)
#ifdef _APPLY_ADV_DIFF_
   where(VG%mask > 0) self%V = self%v1*VG%D
#else
   where(VG%mask > 0) self%advV=(self%v1*VG%D-self%advV)/dt
#endif

   end associate VGrid
   end associate UGrid
   end associate TGrid
   end associate XGrid
END SUBROUTINE uv_advection_2d

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
   if (self%advection_scheme == 0) return
   if (associated(self%logs)) call self%logs%info('uv_advection_3d()',level=3)
   XGrid: associate( XG => self%domain%X )
   TGrid: associate( TG => self%domain%T )
   UGrid: associate( UG => self%domain%U )
   VGrid: associate( VG => self%domain%V )
   do j=UG%jmin-1,UG%jmax
      do i=UG%imin-1,UG%imax
         self%uuadvgrid%hn(i,j,:) = TG%hn(i+1,j,:)
         self%uvadvgrid%hn(i,j,:) = XG%hn(i,j,:)
         self%pka(i,j,:) = 0.5_real64*(self%pk(i,j,:) + self%pk(i+1,j,:))
         self%qka(i,j,:) = 0.5_real64*(self%qk(i,j,:) + self%qk(i+1,j,:))
      end do
   end do
   self%advpk=self%pk
   call self%advection%calculate(self%uuadvgrid,self%pka,self%uvadvgrid,self%qka,dt,UG,self%pk)
   self%advpk=(self%pk-self%advpk)/dt
   do j=VG%jmin-1,UG%jmax
      do i=VG%imin-1,UG%imax
         self%vuadvgrid%hn(i,j,:) = XG%hn(i,j,:)
         self%vvadvgrid%hn(i,j,:) = TG%hn(i,j+1,:)
         self%pka(i,j,:) = 0.5_real64*(self%pk(i,j,:) + self%pk(i,j+1,:))
         self%qka(i,j,:) = 0.5_real64*(self%qk(i,j,:) + self%qk(i,j+1,:))
      end do
   end do
   self%advqk=self%qk
   call self%advection%calculate(self%vuadvgrid,self%pka,self%vvadvgrid,self%qka,dt,VG,self%qk)
   self%advqk=(self%qk-self%advqk)/dt
   end associate VGrid
   end associate UGrid
   end associate TGrid
   end associate XGrid
END SUBROUTINE uv_advection_3d

!---------------------------------------------------------------------------

MODULE SUBROUTINE slow_advection(self,dt)

   !! Advection of time averaged depth integrated transports and slow advection update

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
   real(real64), intent(in) :: dt
      !! timestep [s]

!  Local constants

!  Local variables
   integer :: i,j
!---------------------------------------------------------------------------
   if (self%advection_scheme == 0) return
   if (associated(self%logs)) call self%logs%info('slow_advection()',level=3)
   XGrid: associate( XG => self%domain%X )
   TGrid: associate( TG => self%domain%T )
   UGrid: associate( UG => self%domain%U )
   VGrid: associate( VG => self%domain%V )

   ! [GETM Scientific Report: eq. 2.20]
   do j=UG%jmin-1,UG%jmax
      do i=UG%imin-1,UG%imax
         self%uuadvgrid%D(i,j)  = TG%H(i+1,j)+TG%zin(i+1,j) !KB TG%D(i+1,j)
         self%uvadvgrid%D(i,j)  = XG%D(i,j) ! Knut
         self%Ua(i,j) = 0.5_real64*(self%Ui(i,j) + self%Ui(i+1,j))
         self%Va(i,j) = 0.5_real64*(self%Vi(i,j) + self%Vi(i+1,j))
      end do
   end do
   self%SxA=self%Ui
   call self%advection%calculate(self%uuadvgrid,self%Ua,self%uvadvgrid,self%Va,dt,UG,self%Ui)
   do j=UG%jmin,UG%jmax
      do i=UG%imin,UG%imax
         if (UG%mask(i,j) .ge. 1) then
            self%SxA(i,j)=SUM(self%advpk(i,j,1:))-(self%Ui(i,j)-self%SxA(i,j))/dt
         end if
      end do
   end do

   ! [GETM Scientific Report: eq. 2.21]
   do j=VG%jmin-1,VG%jmax
      do i=VG%imin-1,VG%imax
         self%vuadvgrid%D(i,j)  = XG%D(i,j) ! Knut
         self%vvadvgrid%D(i,j)  = TG%H(i,j+1)+TG%zin(i,j+1) !KB TG%D(i,j+1)
         self%Ua(i,j) = 0.5_real64*(self%Ui(i,j) + self%Ui(i,j+1))
         self%Va(i,j) = 0.5_real64*(self%Vi(i,j) + self%Vi(i,j+1))
      end do
   end do
   self%SyA=self%Vi
   call self%advection%calculate(self%vuadvgrid,self%Ua,self%vvadvgrid,self%Va,dt,VG,self%Vi)
   do j=VG%jmin,VG%jmax
      do i=VG%imin,VG%imax
         if (VG%mask(i,j) .ge. 1) then
            self%SyA(i,j)=SUM(self%advqk(i,j,1:))-(self%Vi(i,j)-self%SyA(i,j))/dt
         end if
      end do
   end do
   end associate VGrid
   end associate UGrid
   end associate TGrid
   end associate XGrid
END SUBROUTINE slow_advection

!---------------------------------------------------------------------------

END SUBMODULE uv_advection_smod
