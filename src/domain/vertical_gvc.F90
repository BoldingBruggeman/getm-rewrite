! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

SUBMODULE (getm_domain : vertical_coordinates_smod) vertical_gvc_smod

!!{!./code/domain/vertical_gvc.md!}

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
   integer :: k,l(3)
!-----------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('init_vertical_gvc()',level=3)

   allocate(self%beta(0:self%T%kmax),stat=stat)     ! dimensionless beta-coordinate
   if (stat /= 0) STOP 'coordinates: Error allocating (beta)'
   allocate(self%sigma(0:self%T%kmax),stat=stat)    ! dimensionless sigma-coordinate
   if (stat /= 0) STOP 'coordinates: Error allocating (sigma)'

   self%beta(0)=  -1._real64
   self%sigma(0)= -1._real64
   if (self%ddu < 0._real64) self%ddu=0._real64
   if (self%ddl < 0._real64) self%ddl=0._real64

   do k=1,self%T%kmax
      self%beta(k)=tanh((self%ddl+self%ddu)*k/float(self%T%kmax)-self%ddl)+tanh(self%ddl)
      self%beta(k)=self%beta(k)/(tanh(self%ddl)+tanh(self%ddu))-1._real64
      self%sigma(k)=k/float(self%T%kmax)-1._real64
   end do
   if (self%gamma_surf) then
      self%kk=self%T%kmax
   else
      self%kk=1
   end if

   l = self%T%l+(/0,0,-1/)
   call mm_s('T%gga',self%T%gga,l,self%T%u,stat=stat)
   if (stat /= 0) stop 'init_gvc: Error allocating memory (self%T%gga)'
   l = self%U%l+(/0,0,-1/)
   call mm_s('U%gga',self%U%gga,l,self%U%u,stat=stat)
   if (stat /= 0) stop 'init_gvc: Error allocating memory (self%U%gga)'
   l = self%V%l+(/0,0,-1/)
   call mm_s('V%gga',self%V%gga,l,self%V%u,stat=stat)
   if (stat /= 0) stop 'init_gvc: Error allocating memory (self%V%gga)'
   l = self%X%l+(/0,0,-1/)
   call mm_s('X%gga',self%X%gga,l,self%X%u,stat=stat)
   if (stat /= 0) stop 'init_gvc: Error allocating memory (self%X%gga)'

   call do_gvc(self,0._real64)
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
   real(real64) :: Dmax=0._real64 !KB
!-----------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('do_gvc()',level=3)
   call update(self%T,self%beta,self%sigma,self%kk,self%Dmin,Dmax,self%Dgamma,dt)
   call update(self%U,self%beta,self%sigma,self%kk,self%Dmin,Dmax,self%Dgamma,dt)
   call update(self%V,self%beta,self%sigma,self%kk,self%Dmin,Dmax,self%Dgamma,dt)
   call update(self%X,self%beta,self%sigma,self%kk,self%Dmin,Dmax,self%Dgamma,dt)
END SUBROUTINE do_gvc

!-----------------------------------------------------------------------------

SUBROUTINE update(grid,beta,sigma,kk,Dmin,Dmax,Dgamma,dt)
   !! A wrapper for vertical coordinate calculations

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_grid), intent(inout) :: grid
   integer, intent(in) :: kk
   real(real64), intent(in) :: beta(0:grid%kmax),sigma(0:grid%kmax),Dmin,Dmax,Dgamma,dt

!  Local variables
   real(real64) :: HH,alpha,zz
   real(real64) :: r,relax=0._real64
   integer :: i,j,k
!-----------------------------------------------------------------------------
   do j=grid%l(2),grid%u(2)
      do i=grid%l(1),grid%u(1)
         if (grid%mask(i,j) > 0) then
            HH=max(grid%zin(i,j)+grid%H(i,j),Dmin)
            if (HH < Dgamma) then
               grid%gga(i,j,:)=sigma
            else
               alpha=min(((beta(kk)-beta(kk-1))-Dgamma/HH*(sigma(kk)-sigma(kk-1))) &
                        /((beta(kk)-beta(kk-1))-(sigma(kk)-sigma(kk-1))),1._real64)
               do k=1,grid%kmax
                  grid%gga(i,j,k)=alpha*sigma(k)+(1._real64-alpha)*beta(k)
               end do
            end if
            grid%gga(i,j,0)=-1._real64
            do k=grid%kmax,1,-1
               grid%gga(i,j,k)=grid%gga(i,j,k)-grid%gga(i,j,k-1)
            end do
            grid%gga(i,j,0)=0._real64
         end if
      end do
   end do

   if (relax > 0._real64 .and. dt > 0._real64) then
      do j=grid%l(2),grid%u(2)
         do i=grid%l(1),grid%u(1)
            if (grid%mask(i,j) > 0) then
               r=relax/dt*grid%H(i,j)/Dmax
               grid%hn(i,j,:)=(grid%ho(i,j,:)*r &
                              +grid%D(i,j)*grid%gga(i,j,1:grid%kmax))/(r+1._real64)
            end if
         end do
      end do
   else
      do j=grid%l(2),grid%u(2)
         do i=grid%l(1),grid%u(1)
            if (grid%mask(i,j) > 0) then
#if 1
               grid%hn(i,j,:)=grid%D(i,j)*grid%gga(i,j,1:grid%kmax)
#else
               HH=max(grid%zin(i,j)+grid%H(i,j),Dmin)
               grid%hn(i,j,:)=HH*grid%gga(i,j,1:grid%kmax)
#endif
            end if
         end do
      end do
   end if
END SUBROUTINE update

!---------------------------------------------------------------------------

END SUBMODULE vertical_gvc_smod
