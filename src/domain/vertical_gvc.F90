! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!! @note
!! ddl and ddu == 0
!! @endnote

SUBMODULE (getm_domain : vertical_coordinates_smod) vertical_gvc_smod

   real(real64), dimension(:), allocatable :: ga
   real(real64), dimension(:), allocatable :: beta
   real(real64), dimension(:), allocatable :: sigma

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

MODULE SUBROUTINE init_gvc(self)
   !! Initializing the general vertical coordinates
   !! Following [GETM Scientific Report: pages 25-28]

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_domain), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: stat
   integer :: k,kk
!-----------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('init_vertical_gvc()',level=3)

   allocate(beta(0:self%T%kmax),stat=stat)     ! dimensionless beta-coordinate
   if (stat /= 0) STOP 'coordinates: Error allocating (beta)'
   allocate(ga(0:self%T%kmax),stat=stat)
   if (stat /= 0) stop 'coordinates: Error allocating (ga)'
   allocate(sigma(0:self%T%kmax),stat=stat)    ! dimensionless sigma-coordinate
   if (stat /= 0) STOP 'coordinates: Error allocating (sigma)'
   do k=0,self%T%kmax
      ga(k) = k
   end do

   beta(0)=  -1._real64
   sigma(0)= -1._real64
   if (self%ddu < 0._real64) self%ddu=0._real64
   if (self%ddl < 0._real64) self%ddl=0._real64

   do k=1,self%T%kmax
      beta(k)=tanh((self%ddl+self%ddu)*k/float(self%T%kmax)-self%ddl)+tanh(self%ddl)
      beta(k)=beta(k)/(tanh(self%ddl)+tanh(self%ddu))-1._real64
      sigma(k)=k/float(self%T%kmax)-1._real64
   end do
   if (self%gamma_surf) then
      kk=self%T%kmax
   else
      kk=1
   end if

   call scaling(self%T,kk,self%Dmin,self%Dgamma)
   call scaling(self%U,kk,self%Dmin,self%Dgamma)
   call scaling(self%V,kk,self%Dmin,self%Dgamma)
   call scaling(self%X,kk,self%Dmin,self%Dgamma)

   call do_gvc(self,0._real64)
   self%T%ho=self%T%hn
   self%U%ho=self%U%hn
   self%V%ho=self%V%hn
   self%X%ho=self%X%hn
END SUBROUTINE init_gvc

!---------------------------------------------------------------------------

MODULE SUBROUTINE do_gvc(self,dt)
   !! A wrapper for vertical coordinate calculations

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_domain), intent(inout) :: self
   real(real64), intent(in) :: dt

!  Local constants

!  Local variables
   real(real64) :: cord_relax=0.0_real64
   real(real64) :: HH
   real(real64) :: zz
   real(real64) :: r
   integer :: i,j,k
   real(real64) :: Lstart,Lstop
!-----------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('do_gvc()',level=3)
   call update(self%T,self%Dmin,self%Dgamma,dt)
   call update(self%U,self%Dmin,self%Dgamma,dt)
   call update(self%V,self%Dmin,self%Dgamma,dt)
   call update(self%X,self%Dmin,self%Dgamma,dt)
END SUBROUTINE do_gvc

!-----------------------------------------------------------------------------

SUBROUTINE scaling(grid,kk,Dmin,Dgamma)
   !! A wrapper for vertical coordinate calculations

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_grid), intent(inout) :: grid
   integer, intent(in) :: kk
   real(real64), intent(in) :: Dmin,Dgamma

!  Local constants

!  Local variables
   integer :: stat
   real(real64) :: HH, alpha
   integer :: i,j,k,l(3)
!-----------------------------------------------------------------------------
   l = grid%l+(/0,0,-1/)
   call mm_s('gga',grid%gga,l,grid%u,stat=stat)
   if (stat /= 0) stop 'coordinates(scaling): Error allocating memory (gga)'

   do j=grid%l(2),grid%u(2)
      do i=grid%l(1),grid%u(1)
         HH=max(grid%zio(i,j)+grid%H(i,j),Dmin)
         alpha=min(((beta(kk)-beta(kk-1))-Dgamma/HH*(sigma(kk)-sigma(kk-1))) &
                  /((beta(kk)-beta(kk-1))-(sigma(kk)-sigma(kk-1))),1._real64)
         grid%gga(i,j,0)=-1._real64
         do k=1,grid%kmax
            grid%gga(i,j,k)=alpha*sigma(k)+(1._real64-alpha)*beta(k)
!            if (grid%gga(i,j,k) .lt. grid%gga(i,j,k-1)) then
!               STDERR kk,(beta(kk)-beta(kk-1)),(sigma(kk)-sigma(kk-1))
!               STDERR Dgamma,HH
!               STDERR alpha
!               STDERR k-1,grid%gga(i,j,k-1),beta(k-1),sigma(k-1)
!               STDERR k,grid%gga(i,j,k),beta(k),sigma(k)
!               stop 'coordinates'
!            end if
         end do
         do k=grid%kmax,1,-1
            grid%gga(i,j,k)=grid%gga(i,j,k)-grid%gga(i,j,k-1)
         end do
         grid%gga(i,j,0)=0._real64
      end do
   end do
END SUBROUTINE scaling

!-----------------------------------------------------------------------------

SUBROUTINE update(grid,Dmax,Dgamma,dt)
   !! A wrapper for vertical coordinate calculations

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_grid), intent(inout) :: grid
   real(real64), intent(in) :: Dmax,Dgamma,dt

!  Local constants

!  Local variables
   real(real64) :: HH,zz,r,cord_relax=0._real64
   integer :: i,j,k
!KB   real(real64) :: Lstart,Lstop
!-----------------------------------------------------------------------------
!KB   call cpu_time(Lstart)
!JB   grid%ho=grid%hn
   do j=grid%l(2),grid%u(2)
      do i=grid%l(1),grid%u(1)
         if (grid%mask(i,j) > 0) then
            HH=grid%zin(i,j)+grid%H(i,j)
            if (HH .lt. Dgamma) then
               do k=1,grid%kmax
                  grid%hn(i,j,k)=HH/grid%kmax
               end do
            else
               if (dt > 0) then
                  r=cord_relax/dt*grid%H(i,j)/Dmax
               else
                  r=0._real64
               end if
               zz=-grid%H(i,j)
!KB - check for r
               do k=1,grid%kmax-1
                  grid%hn(i,j,k)=(grid%ho(i,j,k)*r+HH*(grid%gga(i,j,k)-grid%gga(i,j,k-1)))/(r+1._real64)
                  zz=zz+grid%hn(i,j,k)
               end do
               grid%hn(i,j,grid%kmax)=grid%zin(i,j)-zz
            end if
         end if
      end do
   end do
!KB   call cpu_time(Lstop)
!KB   write(57,*) Lstop-Lstart
END SUBROUTINE update

!---------------------------------------------------------------------------

END SUBMODULE vertical_gvc_smod
