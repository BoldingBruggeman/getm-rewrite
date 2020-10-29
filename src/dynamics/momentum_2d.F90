! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!!{!./pages/momentum_2d.md!}

!> @note
!> self%UEx(i,j)+self%SlUx(i,j)+Slr
!> and
!> self%VEx(i,j)+self%SlVx(i,j)+Slr

!> should be replaced with:

!> SxA, SyA, SxB, SyB, SxD, SyD, SxF, SyF

!> to be consistent with the old GETM documentation
!>
!> @endnote

SUBMODULE (getm_momentum) momentum_2d_smod

CONTAINS

!---------------------------------------------------------------------------

MODULE PROCEDURE uv_momentum_2d

   IMPLICIT NONE

!  Local constants

!  Local variables
   logical :: ufirst=.false.
!---------------------------------------------------------------------------
   call self%logs%info('uv_momentum_2d()',level=2)
   if(ufirst) then
      call u_2d(self,dt,tausx,dpdx)
      call v_2d(self,dt,tausy,dpdy)
      ufirst = .false.
   else
      call v_2d(self,dt,tausy,dpdy)
      call u_2d(self,dt,tausx,dpdx)
      ufirst = .true.
   end if
END PROCEDURE uv_momentum_2d

!---------------------------------------------------------------------------

MODULE PROCEDURE uv_advection_2d

   IMPLICIT NONE

!  Local constants

!  Local variables
   integer :: i,j
!---------------------------------------------------------------------------
   XGrid: associate( XG => self%domain%X )
   TGrid: associate( TG => self%domain%T )
   UGrid: associate( UG => self%domain%U )
   VGrid: associate( VG => self%domain%V )
   do j=UG%jmin,UG%jmax
      do i=UG%imin,UG%imax
         self%uadvgrid%D(i,j)  = TG%D(i+1,j)
         self%vadvgrid%D(i,j)  = XG%D(i,j)
         self%Uadv(i,j) = 0.5_real64*(self%U(i,j) + self%U(i+1,j))
         self%Vadv(i,j) = 0.5_real64*(self%V(i,j) + self%V(i+1,j))
      end do
   end do
   call self%advection%calculate(self%advection_scheme,self%uadvgrid,self%Uadv,self%vadvgrid,self%Vadv,dt,UG,self%U)

   do j=UG%jmin,UG%jmax
      do i=UG%imin,UG%imax
         self%uadvgrid%D(i,j)  = XG%D(i,j)
         self%vadvgrid%D(i,j)  = TG%D(i,j+1)
         self%Uadv(i,j) = 0.5_real64*(self%U(i,j) + self%U(i,j+1))
         self%Vadv(i,j) = 0.5_real64*(self%V(i,j) + self%V(i,j+1))
      end do
   end do
   call self%advection%calculate(self%advection_scheme,self%uadvgrid,self%Uadv,self%vadvgrid,self%Vadv,dt,VG,self%V)
   end associate VGrid
   end associate UGrid
   end associate TGrid
   end associate XGrid
END PROCEDURE uv_advection_2d

!---------------------------------------------------------------------------

module SUBROUTINE u_2d(self,dt,taus,dpdx)

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
   integer :: i,j
   real(real64) :: tausu, Uloc, cord_curv
!KB
   real(real64) :: Slr
!---------------------------------------------------------------------------
   call self%logs%info('u_2d()',level=3)
   UGrid: associate( UG => self%domain%U )
   do j=UG%jmin,UG%jmax
      do i=UG%imin,UG%imax
         if (UG%mask(i,j) == 1 .or. UG%mask(i,j) == 2) then
            tausu = 0.5_real64 * ( taus(i,j) + taus(i+1,j) )
            if (self%U(i,j) .gt. 0._real64) then
               Slr = max( self%Slru(i,j) , 0._real64 )
            else
               Slr = min( self%Slru(i,j) , 0._real64 )
            end if
            self%U(i,j)=(self%U(i,j) &
                        -dt*(g*UG%D(i,j)*dpdx(i,j)+UG%alpha(i,j) &
                        *(-tausu/rho0-self%fV(i,j)+self%UEx(i,j)+self%SlUx(i,j)+Slr))) &
                        /(1._real64+dt*self%ru(i,j)/UG%D(i,j))
            self%Ui(i,j)=self%Ui(i,j)+self%U(i,j)
         end if
      end do
   end do

   ! Semi-implicit treatment of Coriolis force for V-momentum eq.
   XGrid: associate( XG => self%domain%X )
   TGrid: associate( TG => self%domain%T )
   VGrid: associate( VG => self%domain%V )
   do j=VG%jmin,VG%jmax
      do i=VG%imin,VG%imax
         if(VG%mask(i,j) .ge. 1) then
#ifdef NEW_CORI
            ! Espelid et al. [2000], IJNME 49, 1521-1545
            Uloc= &
             ( self%U(i,j  )/sqrt(UG%D(i,j  ))+ self%U(i-1,j  )/sqrt(UG%D(i-1,j  ))  &
             + self%U(i,j+1)/sqrt(UG%D(i,j+1))+ self%U(i-1,j+1)/sqrt(UG%D(i-1,j+1))) &
               *0.25_real64*sqrt(VG%D(i,j))
#else
            Uloc=0.25_real64*( self%U(i-1,j)+self%U(i,j)+self%U(i-1,j+1)+self%U(i,j+1))
#endif
!KB check XG%??
            cord_curv=(self%V(i,j)*(XG%dy(i,j)-XG%dy(i-1,j)) &
                      -Uloc*(TG%dx(i,j+1)-TG%dx(i,j))) &
                      /VG%D(i,j)*VG%inv_area(i,j)
!            cord_curv=(self%V(i,j)*(XG%dy(i,j)-XG%dy(i-1,j))-Uloc*(TG%dx(i,j+1)-TG%dx(i,j)))/VG%D(i,j)*VG%inv_area(i,j)
            self%fU(i,j)=(cord_curv+VG%cor(i,j))*Uloc
         else
            self%fU(i,j)= 0._real64
         end if
!KB
            self%fU(i,j)= 0._real64
      end do
   end do
   end associate VGrid
   end associate TGrid
   end associate XGrid
   end associate UGrid
END subroutine u_2d

!---------------------------------------------------------------------------
module SUBROUTINE v_2d(self,dt,taus,dpdy)

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
   integer :: i,j
   real(real64) :: tausv, Vloc, cord_curv
!KB
   real(real64) :: Slr
!---------------------------------------------------------------------------
   call self%logs%info('v_2d()',level=3)

   VGrid: associate( VG => self%domain%V )
   do j=VG%jmin,VG%jmax
      do i=VG%imin,VG%imax
         if (VG%mask(i,j) == 1 .or. VG%mask(i,j) == 2) then
            tausv = 0.5_real64 * ( taus(i,j) + taus(i,j+1) )
            if (self%V(i,j) .gt. 0._real64) then
               Slr = max( self%Slrv(i,j) , 0._real64 )
            else
               Slr = min( self%Slrv(i,j) , 0._real64 )
            end if
            self%V(i,j)=(self%V(i,j) &
                        -dt*(g*VG%D(i,j)*dpdy(i,j)+VG%alpha(i,j) &
                        *(-tausv/rho0-self%fU(i,j)+self%VEx(i,j)+self%SlVx(i,j)+Slr))) &
                        /(1._real64+dt*self%rv(i,j)/VG%D(i,j))
            self%Vi(i,j)=self%Vi(i,j)+self%V(i,j)
         end if
      end do
   end do

   ! Semi-implicit treatment of Coriolis force for U-momentum eq.
   XGrid: associate( XG => self%domain%X )
   TGrid: associate( TG => self%domain%T )
   UGrid: associate( UG => self%domain%U )
   do j=UG%jmin,UG%jmax
      do i=UG%imin,UG%imax
         if(UG%mask(i,j) .ge. 1) then
#ifdef NEW_CORI
            ! Espelid et al. [2000], IJNME 49, 1521-1545
            Vloc= &
             ( self%V(i,j  )/sqrt(VG%D(i,j  ))+ self%V(i-1,j  )/sqrt(VG%D(i-1,j  ))  &
             + self%V(i,j+1)/sqrt(VG%D(i,j+1))+ self%V(i-1,j+1)/sqrt(VG%D(i-1,j+1))) &
               *0.25_real64*sqrt(UG%D(i,j))
#else
            Vloc=0.25_real64*( self%V(i-1,j)+self%V(i,j)+self%V(i-1,j+1)+self%V(i,j+1))
#endif
            cord_curv=(Vloc*(TG%dy(i+1,j)-TG%dy(i,j))) &
                      +self%U(i,j)*(XG%dx(i,j)-XG%dx(i,j-1)) &
                      /UG%D(i,j)*UG%inv_area(i,j)
            self%fV(i,j)=(cord_curv+UG%cor(i,j))*Vloc
         else
            self%fV(i,j)= 0._real64
         end if
!KB
            self%fV(i,j)= 0._real64
      end do
   end do
   end associate UGrid
   end associate TGrid
   end associate XGrid
   end associate VGrid
END subroutine v_2d

!---------------------------------------------------------------------------

END SUBMODULE momentum_2d_smod
