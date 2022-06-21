! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

SUBMODULE (getm_momentum) velocities_smod

CONTAINS

!---------------------------------------------------------------------------

MODULE SUBROUTINE velocities_2d(self)
   !! 2D velocities

   IMPLICIT NONE

   class(type_getm_momentum), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: i,j
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('velocities_2d()',level=3)

   UGrid: associate( UG => self%domain%U )
   where (UG%mask > 0)
      self%u1 = self%U/UG%D
   else where
      self%u1 = 0._real64
   end where
   end associate UGrid

   VGrid: associate( VG => self%domain%V )
   where (VG%mask > 0)
      self%v1 = self%V/VG%D
   else where
      self%v1 = 0._real64
   end where
   end associate VGrid
END SUBROUTINE velocities_2d

!---------------------------------------------------------------------------

MODULE SUBROUTINE velocities_3d(self)
   !! 3D velocities

   IMPLICIT NONE

   class(type_getm_momentum), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: i,j,k
!---------------------------------------------------------------------------
   if(associated(self%logs)) call self%logs%info('velocities_3d()',level=3)
   UGrid: associate( UG => self%domain%U )
   do k=UG%l(3),UG%u(3)
      do j=UG%l(2),UG%u(2)
         do i=UG%l(1),UG%u(1)
            if (UG%mask(i,j) > 0) then
               self%uk(i,j,k) = self%pk(i,j,k)/UG%hn(i,j,k)
            else
               self%uk(i,j,k) = 0._real64
            end  if
         end do
      end do
   end do
   end associate UGrid

   VGrid: associate( VG => self%domain%V )
   do k=VG%l(3),VG%u(3)
      do j=VG%l(2),VG%u(2)
         do i=VG%l(1),VG%u(1)
            if (VG%mask(i,j) > 0) then
               self%vk(i,j,k) = self%qk(i,j,k)/VG%hn(i,j,k)
            else
               self%vk(i,j,k) = 0._real64
            end  if
         end do
      end do
   end do
   end associate VGrid
END SUBROUTINE velocities_3d

!---------------------------------------------------------------------------

MODULE SUBROUTINE shear_frequency(self,num)
   !!{!./code/velocity_shear_frequency.md!}

   IMPLICIT NONE

   class(type_getm_momentum), intent(inout) :: self
#define _T3_ self%domain%T%l(1):,self%domain%T%l(2):,self%domain%T%l(3):
   real(real64), intent(in) :: num(_T3_)
#undef _T3_

!  Local constants

!  Local variables
   integer :: i,j,k
!---------------------------------------------------------------------------
   if(associated(self%logs)) call self%logs%info('stresses()',level=2)
   TGrid: associate( TG => self%domain%T )
   UGrid: associate( UG => self%domain%U )
   VGrid: associate( VG => self%domain%V )
   do j=TG%jmin,TG%jmax
      do i=TG%imin,TG%imax
         if (TG%mask(i,j) .ge. 1 ) then
            self%SS(i,j,TG%kmax) = 0._real64
            do k=1,TG%kmax-1
#ifndef NEW_SS
               ! This is an older version which we should keep here.
               self%SS(i,j,k)=0.5_real64* ( &
                   ((self%uk(i  ,j,k+1)-self%uk(i  ,j,k)) &
                   /(0.5_real64*(UG%hn(i  ,j,k+1)+UG%hn(i  ,j,k))))**2 &
                  +((self%uk(i-1,j,k+1)-self%uk(i-1,j,k)) &
                   /(0.5_real64*(UG%hn(i-1,j,k+1)+UG%hn(i-1,j,k))))**2 &
                  +((self%vk(i,j  ,k+1)-self%vk(i,j  ,k)) &
                   /(0.5_real64*(VG%hn(i,j  ,k+1)+VG%hn(i,j  ,k))))**2 &
                  +((self%vk(i,j-1,k+1)-self%vk(i,j-1,k)) &
                   /(0.5_real64*(VG%hn(i,j-1,k+1)+VG%hn(i,j-1,k))))**2 &
                                     )
#else
               ! This version should better conserve energy.
               self%SS(i,j,k)=0.5_real64* &
                  ( &
                   (self%uk(i  ,j,k+1)-self%uk(i  ,j,k))**2 &
                   /(UG%hn(i  ,j,k+1)+UG%hn(i  ,j,k))*(num(i,j,k)+num(i+1,j,k)) &
                  +(self%uk(i-1,j,k+1)-self%uk(i-1,j,k))**2 &
                   /(UG%hn(i-1,j,k+1)+UG%hn(i-1,j,k))*(num(i-1,j,k)+num(i,j,k)) &
                  +(self%vk(i,j  ,k+1)-self%vk(i,j  ,k))**2 &
                   /(VG%hn(i,j  ,k+1)+VG%hn(i,j  ,k))*(num(i,j,k)+num(i,j+1,k)) &
                  +(self%vk(i,j-1,k+1)-self%vk(i,j-1,k))**2 &
                   /(VG%hn(i,j-1,k+1)+VG%hn(i,j-1,k))*(num(i,j-1,k)+num(i,j,k)) &
                  )/(0.5_real64*(TG%hn(i,j,k)+TG%hn(i,j,k+1)))/num(i,j,k)
#endif
             end do
         end if
      end do
   end do
   end associate VGrid
   end associate UGrid
   end associate TGrid
END SUBROUTINE shear_frequency

!---------------------------------------------------------------------------

MODULE SUBROUTINE stresses(self,tausx,tausy)
   !! Bottom stress

   IMPLICIT NONE
   class(type_getm_momentum), intent(inout) :: self
#define _T2_ self%domain%T%l(1):,self%domain%T%l(2):
   real(real64), intent(in) :: tausx(_T2_)
   real(real64), intent(in) :: tausy(_T2_)
#undef _T2_

!  Local constants

!  Local variables
   integer :: i,j,k
!---------------------------------------------------------------------------
   if(associated(self%logs)) call self%logs%info('stresses()',level=3)
   k=1 !note
!  x-component of bottom momentum flux at U-points
!  (square of shear velocity = bottom stress in Pa divided by density rho0)
   UGrid: associate( UG => self%domain%U )
   do j=UG%l(2)+1,UG%u(2) !KB loop boundaries
      do i=UG%l(1)+1,UG%u(1)
         if (UG%mask(i,j) > 0) then
            !k = kumin(i,j) ! bottom index
            self%taubx(i,j) = -self%uk(i,j,k)*self%rru(i,j) ! momentum flux
         end if
      enddo
   enddo
   end associate UGrid

!  y-component of bottom momentum flux at V-points
!  (square of shear velocity = bottom stress in Pa divided by density rho0)
   VGrid: associate( VG => self%domain%V )
   do j=VG%l(2)+1,VG%u(2) !KB loop boundaries
      do i=VG%l(1)+1,VG%u(1)
         if (VG%mask(i,j) > 0) then
            !k          = kvmin(i,j) ! bottom index
            self%tauby(i,j) = -self%vk(i,j,k)*self%rrv(i,j) ! momentum flux
         end if
      enddo
   enddo
   end associate VGrid

!  stress magnitude
   TGrid: associate( TG => self%domain%T )
   do j=TG%l(2)+1,TG%u(2) !KB loop boundaries
      do i=TG%l(1)+1,TG%u(1)
         if (TG%mask(i,j) > 0) then
            ! total surface stress at T-points (square of shear velocity = surface stress in Pa divided by density rho0)
            self%ustar2_s(i,j)=sqrt(tausx(i,j)**2 + tausy(i,j)**2)/rho0

            ! total bottom stress at T-points (square of shear velocity = bottom stress in Pa divided by density rho0)
            self%ustar2_b(i,j)=sqrt(0.5_real64*( &
#if 1
                      (self%taubx(i-1,j  ))**2+(self%taubx(i,j))**2 &
                     +(self%tauby(i  ,j-1))**2+(self%tauby(i,j))**2))
#else
                      (self%uk(i-1,j  ,k)*self%rru(i-1,j  ))**2 &
                     +(self%uk(i  ,j  ,k)*self%rru(i  ,j  ))**2 &
                     +(self%vk(i  ,j-1,k)*self%rrv(i  ,j-1))**2 &
                     +(self%vk(i  ,j  ,k)*self%rrv(i  ,j  ))**2))
#endif
         end if
      end do
   end do
   end associate TGrid
END SUBROUTINE stresses

#if 0

!---------------------------------------------------------------------------

MODULE PROCEDURE destag_velocities_2d

   !! Velocity shear

   IMPLICIT NONE

!  Local constants

!  Local variables
   integer :: i,j
!---------------------------------------------------------------------------
   if(associated(self%logs)) call self%logs%info('destag_velocities_2d()',level=2)
END PROCEDURE destag_velocities_2d

!---------------------------------------------------------------------------

MODULE PROCEDURE destag_velocities_3d

   !! Velocity shear

   IMPLICIT NONE

!  Local constants

!  Local variables
   integer :: i,j
!---------------------------------------------------------------------------
   if(associated(self%logs)) call self%logs%info('destag_velocities_3d()',level=2)
END PROCEDURE destag_velocities_3d
#endif

!---------------------------------------------------------------------------

END SUBMODULE velocities_smod
