! Copyright (C) 2020 Bolding & Bruggeman

  !! @note
  !! 2D arrays must be allocated
  !! min_depth must be passed
  !! @endnote

SUBMODULE (getm_momentum) slow_bottom_friction_smod

CONTAINS

!---------------------------------------------------------------------------

module SUBROUTINE slow_bottom_friction(self)

   !! Velocity shear

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: i,j,k
   real(real64) :: uloc,vloc,HH
   real(real64), dimension(:,:), allocatable :: ruu,rvv,Ui,Vi
!---------------------------------------------------------------------------
   call self%logs%info('slow_bottom_friction()',level=2)

#define UG self%domain%U
   do j=UG%l(2),UG%u(2)
      do i=UG%l(1),UG%u(1)
         if (UG%mask(i,j) .ge. 1) then
            Ui(i,j)=self%Uinto(i,j)/(UG%sseo(i,j)+UG%H(i,j))
         else
            Ui(i,j)=0._real64
         end if
      end do
   end do
#undef UG

#define VG self%domain%V
   do j=VG%l(2),VG%u(2)
      do i=VG%l(1),VG%u(1)
         if (VG%mask(i,j) .ge. 1) then
            Vi(i,j)=self%Vinto(i,j)/(VG%sseo(i,j)+VG%H(i,j))
         else
            Vi(i,j)=0._real64
         end if
      end do
   end do
#undef VG

#define UG self%domain%U
   do j=UG%l(2),UG%u(2)
      do i=UG%l(1),UG%u(1)
         if (UG%mask(i,j) .ge. 1) then
            HH=max(UG%ssen(i,j)+UG%H(i,j),self%domain%min_depth)
            ruu(i,j)=(self%zub(i,j)+0.5_real64*HH)/self%zub(i,j)
#if 0
            if (ruu(i,j) .le. 1._real64) then
               LEVEL1 i,j,UG%ssuo(i,j)
               LEVEL1 'Bottom xfriction coefficient infinite.'
               FATAL 'slow_bottom_friction()'
!KB               call getm_error('slow_bottom_friction()','Bottom xfriction coefficient infinite.')
               STOP ! Just if getm_error doesnt actually halt
            end if
#endif
            ruu(i,j)=(kappa/log(ruu(i,j)))**2
            vloc = 0.25_real64*(Vi(i,j)+Vi(i+1,j)+Vi(i,j-1)+Vi(i+1,j-1))
            self%ru(i,j) = ruu(i,j)*sqrt(Ui(i,j)**2+vloc**2)
         else
            self%ru(i,j)=0._real64
         end if
      end do
   end do
#undef UG

#define VG self%domain%V
   do j=VG%l(2),VG%u(2)
      do i=VG%l(1),VG%u(1)
         if (VG%mask(i,j) .ge. 1) then
            HH=max(VG%ssen(i,j)+VG%H(i,j),self%domain%min_depth)
            ruu(i,j)=(self%zvb(i,j)+0.5_real64*HH)/self%zvb(i,j)
#if 0
            if (ruu(i,j) .le. 1._real64) then
               LEVEL1 i,j,UG%ssuo(i,j)
               LEVEL1 'Bottom xfriction coefficient infinite.'
               FATAL 'slow_bottom_friction()'
!KB               call getm_error('slow_bottom_friction()','Bottom xfriction coefficient infinite.')
               STOP ! Just if getm_error doesnt actually halt
            end if
#endif
            ruu(i,j)=(kappa/log(ruu(i,j)))**2
            uloc = 0.25_real64*(Ui(i,j)+Ui(i+1,j)+Ui(i,j-1)+Ui(i+1,j-1))
            self%rv(i,j) = ruu(i,j)*sqrt(uloc**2+Vi(i,j)**2)
         else
            self%rv(i,j)=0._real64
         end if
      end do
   end do
#undef VG
   return
END SUBROUTINE slow_bottom_friction

!---------------------------------------------------------------------------

END SUBMODULE slow_bottom_friction_smod
