! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!!{!./pages/momentum_2d.md!}

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
   call mm_s('zub',self%zub,self%U,def=0._real64,stat=stat)
   call mm_s('zub0',self%zub0,self%U,def=0._real64,stat=stat)
   call mm_s('zvb',self%zvb,self%V,def=0._real64,stat=stat)
   call mm_s('zvb0',self%zvb0,self%V,def=0._real64,stat=stat)
   call mm_s('u1',self%u1,self%U,def=0._real64,stat=stat)
   call mm_s('v1',self%v1,self%V,def=0._real64,stat=stat)

   if (self%Am0 > 0_real64) self%Am=self%Am0

!KB
!if (self%advection_scheme > 0) then
   TGrid: associate( TG => self%domain%T )
   XGrid: associate( XG => self%domain%X )
   ! Grids for U and V advection - updates of time varying fields in advection calling routine
   call mm_s('uadvmask',self%uadvgrid%mask,TG%mask,def=0,stat=stat)
   call mm_s('uadvdx',self%uadvgrid%dx,TG%dx,def=0._real64,stat=stat)
   call mm_s('uadvdy',self%uadvgrid%dy,TG%dy,def=0._real64,stat=stat)
   call mm_s('uadvD',self%uadvgrid%D,TG%D,def=0._real64,stat=stat)

   call mm_s('vadvmask',self%vadvgrid%mask,TG%mask,def=0,stat=stat)
   call mm_s('vadvdx',self%vadvgrid%dx,TG%dx,def=0._real64,stat=stat)
   call mm_s('vadvdy',self%vadvgrid%dy,TG%dy,def=0._real64,stat=stat)
   call mm_s('vadvD',self%vadvgrid%D,TG%D,def=0._real64,stat=stat)

   call mm_s('Ua',self%Ua,self%U,def=0._real64,stat=stat)
   call mm_s('Va',self%Va,self%V,def=0._real64,stat=stat)

   do j=UG%jmin,UG%jmax
      do i=UG%imin-1,UG%imax
         self%uadvgrid%mask(i,j) = TG%mask(i+1,j) ! check this
         self%uadvgrid%dx(i,j) = TG%dx(i+1,j)
         self%uadvgrid%dy(i,j) = TG%dy(i+1,j)
         self%vadvgrid%dx(i,j) = XG%dx(i,j)
         self%vadvgrid%dy(i,j) = XG%dy(i,j)
      end do
   end do
   do j=UG%jmin-1,UG%jmax
      do i=UG%imin,UG%imax
         self%uadvgrid%dx(i,j) = XG%dx(i,j)
         self%uadvgrid%dy(i,j) = XG%dy(i,j)
         self%vadvgrid%mask(i,j) = TG%mask(i,j+1) ! check this
         self%vadvgrid%dx(i,j) = TG%dx(i,j+1)
         self%vadvgrid%dy(i,j) = TG%dy(i,j+1)
      end do
   end do
   end associate XGrid
   end associate TGrid
!end if
   end associate VGrid
   end associate UGrid
END SUBROUTINE uv_initialize_2d

!---------------------------------------------------------------------------

MODULE SUBROUTINE uv_momentum_2d(self,runtype,dt,tausx,tausy,dpdx,dpdy)

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
   integer, intent(in) :: runtype
      !! model runtype
   real(real64), intent(in) :: dt
      !! timestep [s]
#define _T2_ self%domain%T%l(1):,self%domain%T%l(2):
   real(real64), intent(in) :: tausx(_T2_),tausy(_T2_)
      !! surface stresses
   real(real64), intent(in) :: dpdx(_T2_), dpdy(_T2_)
      !! surface pressure gradient - including air pressure
#undef _T2_

!  Local constants

!  Local variables
   logical :: ufirst=.false.
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('uv_momentum_2d()',level=2)
   call self%bottom_friction_2d(runtype)
   call self%uv_advection_2d(dt)
   call self%uv_diffusion_2d(dt)
   if(ufirst) then
      call u_2d(self,dt,tausx,dpdx)
      call self%coriolis_fu()
      call v_2d(self,dt,tausy,dpdy)
      call self%coriolis_fv()
      ufirst = .false.
   else
      call v_2d(self,dt,tausy,dpdy)
      call self%coriolis_fv()
      call u_2d(self,dt,tausx,dpdx)
      call self%coriolis_fu()
      ufirst = .true.
   end if
END SUBROUTINE uv_momentum_2d

!---------------------------------------------------------------------------

SUBROUTINE u_2d(self,dt,taus,dpdx)

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

!  Local constants

!  Local variables
   real(real64) :: tausu
   real(real64) :: Slr
   integer :: i,j
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('u_2d()',level=3)
   UGrid: associate( UG => self%domain%U )
   do j=UG%jmin,UG%jmax
      do i=UG%imin,UG%imax
         if (UG%mask(i,j) == 1 .or. UG%mask(i,j) == 2) then
            tausu = 0.5_real64 * ( taus(i,j) + taus(i+1,j) )
            if (self%U(i,j) .gt. 0._real64) then
               Slr = max( self%SxF(i,j) , 0._real64 )
            else
               Slr = min( self%SxF(i,j) , 0._real64 )
            end if
            self%U(i,j)=(self%U(i,j)-dt*(g*UG%D(i,j)*dpdx(i,j) & ! (2.16) - note SxF is multiplied by alpha
                        +UG%alpha(i,j)*(-tausu/rho0-self%fV(i,j) &
#ifndef _APPLY_ADV_DIFF_
                        +self%advU(i,j)-self%diffu1(i,j)-self%dampU(i,j) &
#endif
                        +self%SxA(i,j)-self%SxB(i,j)+self%SxD(i,j)+Slr))) &
                        /(1._real64+dt*self%ru(i,j)/UG%D(i,j))
            self%Ui(i,j)=self%Ui(i,j)+self%U(i,j)
         end if
      end do
   end do
   end associate UGrid
END SUBROUTINE u_2d

!---------------------------------------------------------------------------

SUBROUTINE v_2d(self,dt,taus,dpdy)

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
      !! GETM momentum type
   real(real64), intent(in) :: dt
      !! timestep [s]
#define _T2_ self%domain%T%l(1):,self%domain%T%l(2):
   real(real64), intent(in) :: taus(_T2_)
      !! surface stress in Y-direction
   real(real64), intent(in) :: dpdy(_T2_)
      !! surface pressure gradient - including air pressure
#undef _T2_

!  Local constants

!  Local variables
   real(real64) :: tausv
   real(real64) :: Slr
   integer :: i,j
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('v_2d()',level=3)
   VGrid: associate( VG => self%domain%V )
   do j=VG%jmin,VG%jmax
      do i=VG%imin,VG%imax
         if (VG%mask(i,j) == 1 .or. VG%mask(i,j) == 2) then
            tausv = 0.5_real64 * ( taus(i,j) + taus(i,j+1) )
            if (self%V(i,j) .gt. 0._real64) then
               Slr = max( self%SyF(i,j) , 0._real64 )
            else
               Slr = min( self%SyF(i,j) , 0._real64 )
            end if
            self%V(i,j)=(self%V(i,j)-dt*(g*VG%D(i,j)*dpdy(i,j) & ! (2.17) - note SxF is multiplied by alpha
                        +VG%alpha(i,j)*(-tausv/rho0+self%fU(i,j) &
#ifndef _APPLY_ADV_DIFF_
                        +self%advV(i,j)-self%diffv1(i,j)-self%dampV(i,j) &
#endif
                        +self%SyA(i,j)-self%SyB(i,j)+self%SyD(i,j)+Slr))) &
                        /(1._real64+dt*self%rv(i,j)/VG%D(i,j))
            self%Vi(i,j)=self%Vi(i,j)+self%V(i,j)
         end if
      end do
   end do
   end associate VGrid
END SUBROUTINE v_2d

!---------------------------------------------------------------------------

END SUBMODULE momentum_2d_smod
