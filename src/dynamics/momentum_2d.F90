! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

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
   real(real64) :: tausu, Uloc, cord_curv
!KB
   real(real64) :: Slr, gammai
   real(real64) :: Dmin = 0.5_real64
!KB
!---------------------------------------------------------------------------
   call self%logs%info('velocity_2d_x()',level=2)
#define UG self%domain%U
   do j=UG%l(2),UG%u(2)
      do i=UG%l(1),UG%u(1)
         if (UG%mask(i,j) == 1 .or. UG%mask(i,j) == 2) then
            tausu = 0.5_real64 * ( taus(i,j) + taus(i+1,j) )
            if (self%U(i,j) .gt. 0._real64) then
               Slr = max( self%Slru(i,j) , 0._real64 )
            else
               Slr = min( self%Slru(i,j) , 0._real64 )
            end if
            !KB dry_u
            self%U(i,j)=(self%U(i,j) &
                        -dt*(g*UG%D(i,j)*dpdx(i,j)+self%dry_u(i,j) &
                        *(-tausu/rho0-self%fV(i,j)+self%UEx(i,j)+self%SlUx(i,j)+Slr))) &
                        /(1._real64+dt*self%ru(i,j)/UG%D(i,j))
            self%Ui(i,j)=self%Ui(i,j)+self%U(i,j)
         end if
      end do
   end do

   ! Semi-implicit treatment of Coriolis force for V-momentum eq.
#define TG self%domain%T
#define VG self%domain%V
#define XG self%domain%X
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
            cord_curv=(self%V(i,j)* &
                       (XG%dy(i,j)-XG%dy(i-1,j)) &
                      -Uloc*(TG%dx(i,j+1)-TG%dx(i,j))) &
                      /VG%D(i,j)*VG%inv_area(i,j)
!            cord_curv=(self%V(i,j)*(XG%dy(i,j)-XG%dy(i-1,j))-Uloc*(TG%dx(i,j+1)-TG%dx(i,j)))/VG%D(i,j)*VG%inv_area(i,j)
            self%fU(i,j)=(cord_curv+VG%cor(i,j))*Uloc
         else
            self%fU(i,j)= 0._real64
         end if
      end do
   end do
#undef TG
#undef UG
#undef VG
#undef XG
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
   real(real64) :: tausv, Vloc, cord_curv
!KB
   real(real64) :: Slr, gammai
   real(real64) :: Dmin = 0.5_real64
!KB
!---------------------------------------------------------------------------
   call self%logs%info('momentum_2d_y()',level=2)

#define VG self%domain%V
   do j=VG%l(2),VG%u(2)
      do i=VG%l(1),VG%u(1)
         if (VG%mask(i,j) == 1 .or. VG%mask(i,j) == 2) then
            tausv = 0.5_real64 * ( taus(i,j) + taus(i,j+1) )
            if (self%V(i,j) .gt. 0._real64) then
               Slr = max( self%Slrv(i,j) , 0._real64 )
            else
               Slr = min( self%Slrv(i,j) , 0._real64 )
            end if
            !KB dry_u
            self%V(i,j)=(self%V(i,j) &
                        -dt*(g*VG%D(i,j)*dpdy(i,j)+self%dry_v(i,j) &
                        *(-tausv/rho0-self%fU(i,j)+self%VEx(i,j)+self%SlVx(i,j)+Slr))) &
                        /(1._real64+dt*self%rv(i,j)/VG%D(i,j))
            self%Vi(i,j)=self%Vi(i,j)+self%V(i,j)
         end if
      end do
   end do

   ! Semi-implicit treatment of Coriolis force for U-momentum eq.
#define TG self%domain%T
#define UG self%domain%U
#define XG self%domain%X
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
                      +self%U(i,j) &
                      *(XG%dx(i,j)-XG%dx(i,j-1)) &
                      /UG%D(i,j)*UG%inv_area(i,j)
            self%fV(i,j)=(cord_curv+UG%cor(i,j))*Vloc
         else
            self%fV(i,j)= 0._real64
         end if
      end do
   end do
#undef TG
#undef UG
#undef VG
#undef XG
   return
END SUBROUTINE momentum_y_2d

!---------------------------------------------------------------------------

END SUBMODULE momentum_2d_smod
