! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

SUBMODULE (getm_momentum) velocities_smod

CONTAINS

!---------------------------------------------------------------------------

module SUBROUTINE velocities_2d(self)
   !! 2D velocities

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: i,j
!---------------------------------------------------------------------------
   call self%logs%info('velocities_2d()',level=2)

#define UG self%domain%U
   where (UG%mask > 0)
      self%u1 = self%U/UG%D
   end where
#undef UG

#define VG self%domain%V
   where (VG%mask > 0)
      self%v1 = self%V/VG%D
   end where
#undef VG

   return
END SUBROUTINE velocities_2d

!---------------------------------------------------------------------------

module SUBROUTINE velocities_3d(self)
   !! 3D velocities

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: i,j,k
!---------------------------------------------------------------------------
   call self%logs%info('velocities_3d()',level=2)
#define UG self%domain%U
   do k=UG%l(3),UG%u(3)
      do j=UG%l(2),UG%u(2)
         do i=UG%l(1),UG%u(1)
            if (UG%mask(i,j) > 0) then
               self%uk(i,j,k) = self%pk(i,j,k)/UG%hn(i,j,k)
            end  if
         end do
      end do
   end do
#undef UG

#define VG self%domain%V
   do k=VG%l(3),VG%u(3)
      do j=VG%l(2),VG%u(2)
         do i=VG%l(1),VG%u(1)
            if (VG%mask(i,j) > 0) then
               self%vk(i,j,k) = self%qk(i,j,k)/VG%hn(i,j,k)
            end  if
         end do
      end do
   end do
#undef VG

   return
END SUBROUTINE velocities_3d

!---------------------------------------------------------------------------

module SUBROUTINE velocity_shear_frequency(self,num,SS)

   !!{!./code/velocity_shear_frequency.md!}

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
   real(real64), dimension(:,:,:), intent(in) :: num
   real(real64), dimension(:,:,:), intent(inout) :: SS

!  Local constants

!  Local variables
   integer :: i,j,k
!---------------------------------------------------------------------------
   call self%logs%info('stresses()',level=2)
! Just use already calculated velocities
#if 1
#define TG self%domain%T
#define UG self%domain%U
#define VG self%domain%V
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
#undef TG
#undef UG
#undef VG
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
   return
END SUBROUTINE velocity_shear_frequency

!---------------------------------------------------------------------------

module SUBROUTINE stresses(self)

   !! Bottom stress

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: i,j,k
!---------------------------------------------------------------------------
   call self%logs%info('stresses()',level=2)

   k=1
!  x-component of bottom momentum flux at U-points
#define UG self%domain%U
   do j=UG%l(2)+1,UG%u(2)
      do i=UG%l(1)+1,UG%u(1)
         if (UG%mask(i,j) > 0) then
            !k = kumin(i,j) ! bottom index
            self%taubx(i,j) = -self%uk(i,j,k)*self%rru(i,j) ! momentum flux
         end if
      enddo
   enddo
#undef UG

!  y-component of bottom momentum flux at V-points
#define VG self%domain%V
   do j=VG%l(2)+1,VG%u(2)
      do i=VG%l(1)+1,VG%u(1)
         if (VG%mask(i,j) > 0) then
            !k          = kvmin(i,j) ! bottom index
            self%tauby(i,j) = -self%vk(i,j,k)*self%rrv(i,j) ! momentum flux
         end if
      enddo
   enddo
#undef VG

!  stress magnitude
#define TG self%domain%T
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
#undef TG
   return
END SUBROUTINE stresses

#if 0

!---------------------------------------------------------------------------

module SUBROUTINE destag_velocities_2d(self)

   !! Velocity shear

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: i,j
!---------------------------------------------------------------------------
   call self%logs%info('destag_velocities_2d()',level=2)

   return
END SUBROUTINE destag_velocities_2d

!---------------------------------------------------------------------------

module SUBROUTINE destag_velocities_3d(self)

   !! Velocity shear

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: i,j
!---------------------------------------------------------------------------
   call self%logs%info('destag_velocities_3d()',level=2)

   return
END SUBROUTINE destag_velocities_3d
#endif

!---------------------------------------------------------------------------

END SUBMODULE velocities_smod
