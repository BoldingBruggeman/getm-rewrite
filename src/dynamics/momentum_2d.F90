! Copyright (C) 2020 Bolding & Bruggeman

!!{!./pages/momentum_2d.md!}

SUBMODULE (getm_momentum) momentum_2d_smod

CONTAINS

!---------------------------------------------------------------------------

module SUBROUTINE momentum_x_2d(self,dt,dpdx,taus)

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
      !! GETM momentum type
   real(real64), intent(in) :: dt
      !! timestep [s]
   real(real64), dimension(:,:), intent(in) :: dpdx
      !! surface pressure gradient - including air pressure
   real(real64), dimension(:,:), intent(in) :: taus
      !! surface stress in X-direction

!  Local constants

!  Local variables
   integer :: i,j
   real(real64) :: tausu
!KB
   real(real64) :: Slr, gammai
   real(real64) :: min_depth = 0.5_real64
!KB
!---------------------------------------------------------------------------
   call self%logs%info('velocity_2d_x()',level=2)
#if 1
   do j=self%domain%U%l(2),self%domain%U%u(2)
      do i=self%domain%U%l(1),self%domain%U%u(1)
         if (self%domain%U%mask(i,j) == 1 .or. self%domain%U%mask(i,j) == 2) then
            tausu = 0.5_real64 * ( taus(i,j) + taus(i+1,j) )
            if (self%U(i,j) .gt. 0._real64) then
               Slr = max( self%Slru(i,j) , 0._real64 )
            else
               Slr = min( self%Slru(i,j) , 0._real64 )
            end if
            !KB dry_u
            self%U(i,j)=(self%U(i,j) &
                        -dt*(g*self%domain%U%D(i,j)*dpdx(i,j)+self%dry_u(i,j) &
                        *(-tausu/rho_0-self%fV(i,j)+self%UEx(i,j)+self%SlUx(i,j)+Slr))) &
                        /(1._real64+dt*self%ru(i,j)/self%domain%U%D(i,j))
            self%Uint(i,j)=self%Uint(i,j)+self%U(i,j)
         end if
      end do
   end do
#if 0
   ! Semi-implicit treatment of Coriolis force for V-momentum eq.
   do j=jmin,jmax
      do i=imin,imax
         if(self%V%av(i,j) .ge. 1) then
! Espelid et al. [2000], IJNME 49, 1521-1545
#ifdef NEW_CORI
            Uloc= &
             ( U(i,j  )/sqrt(DU(i,j  ))+ U(i-1,j  )/sqrt(DU(i-1,j  ))  &
             + U(i,j+1)/sqrt(DU(i,j+1))+ U(i-1,j+1)/sqrt(DU(i-1,j+1))) &
               *_QUART_*sqrt(DV(i,j))
#else
            Uloc=_QUART_*( U(i-1,j)+U(i,j)+U(i-1,j+1)+U(i,j+1))
#endif
#if defined(SPHERICAL) || defined(CURVILINEAR)
            cord_curv=(V(i,j)*(DYX-DYXIM1)-Uloc*(DXCJP1-DXC)) &
                      /DV(i,j)*ARVD1
            fU(i,j)=(cord_curv+corv(i,j))*Uloc
#else
            fU(i,j)=corv(i,j)*Uloc
#endif
         else
            fU(i,j)= _ZERO_
         end if
      end do
   end do
#endif
#endif
   return
END SUBROUTINE momentum_x_2d

!---------------------------------------------------------------------------

module SUBROUTINE momentum_y_2d(self,dt,dpdy,taus)

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
      !! GETM momentum type
   real(real64), intent(in) :: dt
      !! timestep [s]
   real(real64), dimension(:,:), intent(in) :: dpdy
      !! surface pressure gradient - including air pressure
   real(real64), dimension(:,:), intent(in) :: taus
      !! surface stress in X-direction

!  Local constants

!  Local variables
   integer :: i,j
   real(real64) :: tausu
!KB
   real(real64) :: Slr, gammai
   real(real64) :: min_depth = 0.5_real64
!KB
!---------------------------------------------------------------------------
#if 0
   call self%logs%info('momentum_2d_y()',level=2)

   do j=domain%U%l(2),domain%U%u(2)
      do i=domain%U%l(1),domain%U%u(1)
         if (domain%U%mask(i,j) == 1 .or. domain%U%mask(i,j) == 2) then
            tausu = 0.5_real64 * ( taus(i,j) + taus(i+1,j) )
            if (self%U(i,j) .gt. 0._real64) then
               Slr = max( self%Slru(i,j) , 0._real64 )
            else
               Slr = min( self%Slru(i,j) , 0._real64 )
            end if
            self%U(i,j)=(self%U(i,j) &
                        -dt*(g*domain%U%D(i,j)*dpdx(i,j)+self%dry_u(i,j) &
                        *(-tausu/rho_0-self%fV(i,j)+self%UEx(i,j)+self%SlUx(i,j)+Slr))) &
                        /(1._real64+dt*self%ru(i,j)/domain%U%D(i,j))
            self%Vint(i,j)=self%Vint(i,j)+self%V(i,j)
         end if
      end do
   end do
#if 0
   ! Semi-implicit treatment of Coriolis force for V-momentum eq.
   do j=jmin,jmax
      do i=imin,imax
         if(av(i,j) .ge. 1) then
! Espelid et al. [2000], IJNME 49, 1521-1545
#ifdef NEW_CORI
            Uloc= &
             ( U(i,j  )/sqrt(DU(i,j  ))+ U(i-1,j  )/sqrt(DU(i-1,j  ))  &
             + U(i,j+1)/sqrt(DU(i,j+1))+ U(i-1,j+1)/sqrt(DU(i-1,j+1))) &
               *_QUART_*sqrt(DV(i,j))
#else
            Uloc=_QUART_*( U(i-1,j)+U(i,j)+U(i-1,j+1)+U(i,j+1))
#endif
#if defined(SPHERICAL) || defined(CURVILINEAR)
            cord_curv=(V(i,j)*(DYX-DYXIM1)-Uloc*(DXCJP1-DXC)) &
                      /DV(i,j)*ARVD1
            fU(i,j)=(cord_curv+corv(i,j))*Uloc
#else
            fU(i,j)=corv(i,j)*Uloc
#endif
         else
            fU(i,j)= _ZERO_
         end if
      end do
   end do
#endif
#endif
   return
END SUBROUTINE momentum_y_2d

!---------------------------------------------------------------------------

END SUBMODULE momentum_2d_smod
