! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

SUBMODULE (getm_momentum) velocities_smod

CONTAINS

!---------------------------------------------------------------------------

MODULE PROCEDURE velocities_2d
   !! 2D velocities

   IMPLICIT NONE

!  Local constants

!  Local variables
   integer :: i,j
!---------------------------------------------------------------------------
   call self%logs%info('velocities_2d()',level=2)

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

END PROCEDURE velocities_2d

!---------------------------------------------------------------------------

MODULE PROCEDURE velocities_3d
   !! 3D velocities

   IMPLICIT NONE

!  Local constants

!  Local variables
   integer :: i,j,k
!---------------------------------------------------------------------------
   call self%logs%info('velocities_3d()',level=2)

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
END PROCEDURE velocities_3d

!---------------------------------------------------------------------------

MODULE PROCEDURE velocity_shear_frequency

   !!{!./code/velocity_shear_frequency.md!}

   IMPLICIT NONE

!  Local constants

!  Local variables
   integer :: i,j,k
!---------------------------------------------------------------------------
   call self%logs%info('stresses()',level=2)
! Just use already calculated velocities
#if 1
   TGrid: associate( TG => self%domain%T )
   UGrid: associate( UG => self%domain%U )
   VGrid: associate( VG => self%domain%V )
   do j=TG%jmin,TG%jmax
      do i=TG%imin,TG%imax
         if (TG%mask(i,j) .ge. 1 ) then
            SS(i,j,TG%kmax) = 0._real64
            do k=1,TG%kmax-1
#ifndef NEW_SS
               ! This is an older version which we should keep here.
               SS(i,j,k)=0.5_real64* ( &
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
               SS(i,j,k)=0.5_real64* &
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
#else
! This is an older version which we should keep here.
#ifndef NEW_SS
              SS(i,j,k)=_HALF_* (                                             &
                   ( (pk(i,j,k+1)/hun(i,j,k+1)-pk(i,j,k)/hun(i,j,k))          &
                   /(_HALF_*(hun(i,j,k+1)+hun(i,j,k))) )**2                   &
                +  ( (pk(i-1,j,k+1)/hun(i-1,j,k+1)-pk(i-1,j,k)/hun(i-1,j,k))  &
                   /(_HALF_*(hun(i-1,j,k+1)+hun(i-1,j,k))) )**2               &
                +  ( (qk(i,j,k+1)/hvn(i,j,k+1)-qk(i,j,k)/hvn(i,j,k))          &
                   /(_HALF_*(hvn(i,j,k+1)+hvn(i,j,k))) )**2                   &
                +  ( (qk(i,j-1,k+1)/hvn(i,j-1,k+1)-qk(i,j-1,k)/hvn(i,j-1,k))  &
                   /(_HALF_*(hvn(i,j-1,k+1)+hvn(i,j-1,k))) )**2               &
                            )
#else
! This version should better conserve energy.
              SS(i,j,k)=0.5_real64* ( &
                   (pk(i,j,k+1)/hun(i,j,k+1)-pk(i,j,k)/hun(i,j,k))**2 &
                   /(0.25_real64*(hun(i,j,k+1)+hun(i,j,k))*(num(i,j,k)+num(i+1,j,k)) &
               +  (pk(i-1,j,k+1)/hun(i-1,j,k+1)-pk(i-1,j,k)/hun(i-1,j,k))**2  &
                  /(_HALF_*(hun(i-1,j,k+1)+hun(i-1,j,k)))                     &
                   *_HALF_*(num(i-1,j,k)+num(i,j,k))                          &
                +  (qk(i,j,k+1)/hvn(i,j,k+1)-qk(i,j,k)/hvn(i,j,k))**2         &
                   /(_HALF_*(hvn(i,j,k+1)+hvn(i,j,k)))                        &
                    *_HALF_*(num(i,j,k)+num(i,j+1,k))                         &
                +  (qk(i,j-1,k+1)/hvn(i,j-1,k+1)-qk(i,j-1,k)/hvn(i,j-1,k))**2 &
                   /(_HALF_*(hvn(i,j-1,k+1)+hvn(i,j-1,k)))                    &
                    *_HALF_*(num(i,j-1,k)+num(i,j,k))                         &
                            )/(_HALF_*(hn(i,j,k)+hn(i,j,k+1)))/num(i,j,k)
#endif
            end do
         end if
      end do
   end do
#endif
END PROCEDURE velocity_shear_frequency

!---------------------------------------------------------------------------

MODULE PROCEDURE stresses

   !! Bottom stress

   IMPLICIT NONE

!  Local constants

!  Local variables
   integer :: i,j,k
!---------------------------------------------------------------------------
   call self%logs%info('stresses()',level=2)

   k=1
!  x-component of bottom momentum flux at U-points
   UGrid: associate( UG => self%domain%U )
   do j=UG%l(2)+1,UG%u(2)
      do i=UG%l(1)+1,UG%u(1)
         if (UG%mask(i,j) > 0) then
            !k = kumin(i,j) ! bottom index
            self%taubx(i,j) = -self%uk(i,j,k)*self%rru(i,j) ! momentum flux
         end if
      enddo
   enddo
   end associate UGrid

!  y-component of bottom momentum flux at V-points
   VGrid: associate( VG => self%domain%V )
   do j=VG%l(2)+1,VG%u(2)
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
   do j=TG%l(2)+1,TG%u(2)
      do i=TG%l(1)+1,TG%u(1)
         if (TG%mask(i,j) > 0) then
#if 0
         ! lower indices at U- and V-points
         ku1=kumin(i-1,j  )
         ku2=kumin(i  ,j  )
         kv1=kvmin(i  ,j-1)
         kv2=kvmin(i  ,j  )
#endif
         ! total bottom stress at T-points
         self%taub(i,j)=sqrt(0.5_real64*( &
                   (self%uk(i-1,j  ,k)*self%rru(i-1,j  ))**2 &
                  +(self%uk(i  ,j  ,k)*self%rru(i  ,j  ))**2 &
                  +(self%vk(i  ,j-1,k)*self%rrv(i  ,j-1))**2 &
                  +(self%vk(i  ,j  ,k)*self%rrv(i  ,j  ))**2))
         end if
      end do
   end do
   end associate TGrid
END PROCEDURE stresses

#if 0

!---------------------------------------------------------------------------

MODULE PROCEDURE destag_velocities_2d

   !! Velocity shear

   IMPLICIT NONE

!  Local constants

!  Local variables
   integer :: i,j
!---------------------------------------------------------------------------
   call self%logs%info('destag_velocities_2d()',level=2)

END PROCEDURE destag_velocities_2d

!---------------------------------------------------------------------------

MODULE PROCEDURE destag_velocities_3d

   !! Velocity shear

   IMPLICIT NONE

!  Local constants

!  Local variables
   integer :: i,j
!---------------------------------------------------------------------------
   call self%logs%info('destag_velocities_3d()',level=2)

END PROCEDURE destag_velocities_3d
#endif

!---------------------------------------------------------------------------

END SUBMODULE velocities_smod
