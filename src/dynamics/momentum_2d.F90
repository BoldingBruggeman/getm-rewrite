! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!!{!./code/dynamics/momentum_2d.md!}

SUBMODULE (getm_momentum) momentum_2d_smod

CONTAINS

!---------------------------------------------------------------------------

MODULE SUBROUTINE uv_initialize_2d(self)

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: stat
   integer :: i,j
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('uv_initialize_2d()',level=2)
   UGrid: associate( UG => self%domain%U )
   VGrid: associate( VG => self%domain%V )
   call mm_s('U',self%U,UG%l(1:2),UG%u(1:2),def=0._real64,stat=stat)
   call mm_s('V',self%V,VG%l(1:2),VG%u(1:2),def=0._real64,stat=stat)
   call mm_s('Ui',self%Ui,self%U,def=0._real64,stat=stat)
   call mm_s('Vi',self%Vi,self%V,def=0._real64,stat=stat)
   call mm_s('Uio',self%Uio,self%U,def=0._real64,stat=stat)
   call mm_s('Vio',self%Vio,self%V,def=0._real64,stat=stat)
   call mm_s('Ua',self%Ua,self%U,def=0._real64,stat=stat)
   call mm_s('Va',self%Va,self%V,def=0._real64,stat=stat)
   call mm_s('fU',self%fU,self%U,def=0._real64,stat=stat)
   call mm_s('fV',self%fV,self%V,def=0._real64,stat=stat)
   call mm_s('advU',self%advU,self%U,def=0._real64,stat=stat)
   call mm_s('advV',self%advV,self%V,def=0._real64,stat=stat)
   call mm_s('diffu1',self%diffu1,self%U,def=0._real64,stat=stat)
   call mm_s('diffv1',self%diffv1,self%V,def=0._real64,stat=stat)
   call mm_s('dampU',self%dampU,self%U,def=0._real64,stat=stat)
   call mm_s('dampV',self%dampV,self%V,def=0._real64,stat=stat)
   call mm_s('SxA',self%SxA,self%U,def=0._real64,stat=stat)
   call mm_s('SyA',self%SyA,self%V,def=0._real64,stat=stat)
   call mm_s('SxB',self%SxB,self%U,def=0._real64,stat=stat)
   call mm_s('SyB',self%SyB,self%V,def=0._real64,stat=stat)
   call mm_s('SxD',self%SxD,self%U,def=0._real64,stat=stat)
   call mm_s('SyD',self%SyD,self%V,def=0._real64,stat=stat)
   call mm_s('SxF',self%SxF,self%U,def=0._real64,stat=stat)
   call mm_s('SyF',self%SyF,self%V,def=0._real64,stat=stat)
   call mm_s('Am',self%Am,self%U,def=0._real64,stat=stat)
   call mm_s('An',self%An,self%U,def=0._real64,stat=stat)
   call mm_s('ru',self%ru,self%U,def=0._real64,stat=stat)
   call mm_s('rv',self%rv,self%V,def=0._real64,stat=stat)
   call mm_s('rru',self%rru,self%U,def=0._real64,stat=stat)
   call mm_s('rrv',self%rrv,self%V,def=0._real64,stat=stat)
   call mm_s('u1',self%u1,self%U,def=0._real64,stat=stat)
   call mm_s('v1',self%v1,self%V,def=0._real64,stat=stat)

   if (self%Am0 > 0_real64) self%Am=self%Am0

   if (self%advection_scheme > 0) then
      TGrid: associate( TG => self%domain%T )
      XGrid: associate( XG => self%domain%X )
      ! Grids for U and V advection - updates of time varying fields in advection calling routine
      call mm_s('uuadvmask',self%uuadvgrid%mask,TG%mask,def=0,stat=stat)
      call mm_s('uuadvdx',self%uuadvgrid%dx,TG%dx,def=0._real64,stat=stat)
      call mm_s('uuadvdy',self%uuadvgrid%dy,TG%dy,def=0._real64,stat=stat)
      call mm_s('uuadvD',self%uuadvgrid%D,TG%D,def=0._real64,stat=stat)

      call mm_s('uvadvmask',self%uvadvgrid%mask,TG%mask,def=0,stat=stat)
      call mm_s('uvadvdx',self%uvadvgrid%dx,TG%dx,def=0._real64,stat=stat)
      call mm_s('uvadvdy',self%uvadvgrid%dy,TG%dy,def=0._real64,stat=stat)
      call mm_s('uvadvD',self%uvadvgrid%D,TG%D,def=0._real64,stat=stat)

      call mm_s('vuadvmask',self%vuadvgrid%mask,TG%mask,def=0,stat=stat)
      call mm_s('vuadvdx',self%vuadvgrid%dx,TG%dx,def=0._real64,stat=stat)
      call mm_s('vuadvdy',self%vuadvgrid%dy,TG%dy,def=0._real64,stat=stat)
      call mm_s('vuadvD',self%vuadvgrid%D,TG%D,def=0._real64,stat=stat)

      call mm_s('vvadvmask',self%vvadvgrid%mask,TG%mask,def=0,stat=stat)
      call mm_s('vvadvdx',self%vvadvgrid%dx,TG%dx,def=0._real64,stat=stat)
      call mm_s('vvadvdy',self%vvadvgrid%dy,TG%dy,def=0._real64,stat=stat)
      call mm_s('vvadvD',self%vvadvgrid%D,TG%D,def=0._real64,stat=stat)

      call mm_s('Ua',self%Ua,self%U,def=0._real64,stat=stat)
      call mm_s('Va',self%Va,self%V,def=0._real64,stat=stat)

      do j=UG%l(2),UG%u(2)-1
         do i=UG%l(1),UG%u(1)-1
            if (UG%mask(i,j) /= 0 .and. UG%mask(i+1,j) /= 0) self%uuadvgrid%mask(i,j) = 1
            self%uuadvgrid%dx(i,j) = TG%dx(i+1,j)
            self%uuadvgrid%dy(i,j) = TG%dy(i+1,j)

            if (UG%mask(i,j) /= 0 .and. UG%mask(i,j+1) /= 0) self%uvadvgrid%mask(i,j) = 1
            self%uvadvgrid%dx(i,j) = XG%dx(i,j)
            self%uvadvgrid%dy(i,j) = XG%dy(i,j)

            if (VG%mask(i,j) /= 0 .and. VG%mask(i+1,j) /= 0) self%vuadvgrid%mask(i,j) = 1
            self%vuadvgrid%dx(i,j) = XG%dx(i,j)
            self%vuadvgrid%dy(i,j) = XG%dy(i,j)

            if (VG%mask(i,j) /= 0 .and. VG%mask(i,j+1) /= 0) self%vvadvgrid%mask(i,j) = 1
            self%vvadvgrid%dx(i,j) = TG%dx(i,j+1)
            self%vvadvgrid%dy(i,j) = TG%dy(i,j+1)
         end do
      end do
      end associate XGrid
      end associate TGrid
   end if
   end associate VGrid
   end associate UGrid
END SUBROUTINE uv_initialize_2d

!---------------------------------------------------------------------------

MODULE SUBROUTINE uv_momentum_2d(self,runtype,dt,tausx,tausy,dpdx,dpdy)
   !! Solve the horizontal 2D momemtum equations - U, V
   !! [GETM Scientific Report: eqs. 2.16, 2.17]

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
   integer, intent(in) :: runtype
      !! model runtype
   real(real64), intent(in) :: dt
      !! timestep [s]
#define _U2_ self%domain%U%l(1):,self%domain%U%l(2):
   real(real64), intent(in) :: tausx(_U2_)
      !! surface stress in local x-direction
#undef _U2_
#define _V2_ self%domain%V%l(1):,self%domain%V%l(2):
   real(real64), intent(in) :: tausy(_V2_)
      !! surface stress in local y-direction
#undef _V2_
#define _U2_ self%domain%U%l(1):,self%domain%U%l(2):
   real(real64), intent(in) :: dpdx(_U2_)
      !! surface pressure gradient - including air pressure
#undef _U2_
#define _V2_ self%domain%V%l(1):,self%domain%V%l(2):
   real(real64), intent(in) :: dpdy(_V2_)
      !! surface pressure gradient - including air pressure
#undef _V2_

!  Local constants

!  Local variables
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('uv_momentum_2d()',level=2)
   if (self%apply_bottom_friction) call self%bottom_friction_2d(runtype)
   if (self%advection_scheme > 0) call self%uv_advection_2d(dt)
   if (self%apply_diffusion) call self%uv_diffusion_2d(dt)
   if(self%ufirst) then
      call u_2d(self,dt,tausx,dpdx)
      call self%coriolis_fu()
      call v_2d(self,dt,tausy,dpdy)
      call self%coriolis_fv()
      self%ufirst = .false.
   else
      call v_2d(self,dt,tausy,dpdy)
      call self%coriolis_fv()
      call u_2d(self,dt,tausx,dpdx)
      call self%coriolis_fu()
      self%ufirst = .true.
   end if
END SUBROUTINE uv_momentum_2d

!---------------------------------------------------------------------------

MODULE SUBROUTINE u_2d(self,dt,tausx,dpdx)

!!{!./code/dynamics/u_2d.md!}

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
      !! GETM momentum type
   real(real64), intent(in) :: dt
      !! timestep [s]
#define _U2_ self%domain%U%l(1):,self%domain%U%l(2):
   real(real64), intent(in) :: tausx(_U2_)
      !! surface stress in local x-direction
   real(real64), intent(in) :: dpdx(_U2_)
      !! surface pressure gradient - including air pressure
#undef _U2_

!  Local constants

!  Local variables
   real(real64) :: Slr, rho0i
   integer :: i,j
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('u_2d()',level=3)
   rho0i=1._real64/rho0
   UGrid: associate( UG => self%domain%U )
   do j=UG%jmin,UG%jmax
      do i=UG%imin,UG%imax
         if (UG%mask(i,j) == 1 .or. UG%mask(i,j) == 2) then
            if (self%U(i,j) .gt. 0._real64) then
               Slr = min( self%SxF(i,j) , 0._real64 )
            else
               Slr = max( self%SxF(i,j) , 0._real64 )
            end if
            ! [GETM Scientific Report: eqs. 2.14, 2.16]
            self%U(i,j) = (self%U(i,j) + dt * (-g * UG%D(i,j) * dpdx(i,j) & ! note SxF is multiplied by alpha
                           + UG%alpha(i,j) * (tausx(i,j) * rho0i + self%fV(i,j) &
#ifndef _APPLY_ADV_DIFF_
                           + self%advU(i,j) + self%diffu1(i,j) + self%dampU(i,j) &
#endif
                           + self%SxA(i,j) + self%SxB(i,j) + self%SxD(i,j) + Slr))) &
                          / (1._real64 + dt * self%ru(i,j) / UG%D(i,j))
            !self%Ui(i,j)=self%Ui(i,j)+self%U(i,j)   ! JB now done in Python, needs to include halos for slow advection
         end if
      end do
   end do
   call self%domain%mirror_bdys(UG,self%U)
   end associate UGrid
END SUBROUTINE u_2d

!---------------------------------------------------------------------------

MODULE SUBROUTINE v_2d(self,dt,tausy,dpdy)

!!{!./code/dynamics/v_2d.md!}

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
      !! GETM momentum type
   real(real64), intent(in) :: dt
      !! timestep [s]
#define _V2_ self%domain%V%l(1):,self%domain%V%l(2):
   real(real64), intent(in) :: tausy(_V2_)
      !! surface stress in local y-direction
   real(real64), intent(in) :: dpdy(_V2_)
      !! surface pressure gradient - including air pressure
#undef _U2_

!  Local constants

!  Local variables
   real(real64) :: Slr,rho0i
   integer :: i,j
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('v_2d()',level=3)
   rho0i=1._real64/rho0
   VGrid: associate( VG => self%domain%V )
   do j=VG%jmin,VG%jmax
      do i=VG%imin,VG%imax
         if (VG%mask(i,j) == 1 .or. VG%mask(i,j) == 2) then
            if (self%V(i,j) .gt. 0._real64) then
               Slr = min( self%SyF(i,j) , 0._real64 )
            else
               Slr = max( self%SyF(i,j) , 0._real64 )
            end if
            ! [GETM Scientific Report: eqs. 2.15, 2.17]
            self%V(i,j) = (self%V(i,j) + dt * (-g * VG%D(i,j) * dpdy(i,j) & ! note SyF is multiplied by alpha
                           + VG%alpha(i,j) * (tausy(i,j) * rho0i - self%fU(i,j) &
#ifndef _APPLY_ADV_DIFF_
                           + self%advV(i,j) + self%diffv1(i,j) + self%dampV(i,j) &
#endif
                           + self%SyA(i,j) + self%SyB(i,j) + self%SyD(i,j) + Slr))) &
                          / (1._real64 + dt * self%rv(i,j) / VG%D(i,j))
            !self%Vi(i,j)=self%Vi(i,j)+self%V(i,j)   ! JB now done in Python, needs to include halos for slow advection
         end if
      end do
   end do
   call self%domain%mirror_bdys(VG,self%V)
   end associate VGrid
END SUBROUTINE v_2d

!---------------------------------------------------------------------------

END SUBMODULE momentum_2d_smod
