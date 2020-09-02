! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard
   !! This routine which is called at every micro time step updates all
   !! necessary depth related information. These are the water depths in the
   !! T-, U- and V-points, {\tt D}, {\tt DU} and {\tt DV}, respectively,
   !! and the drying value $\alpha$ defined in equation (\ref{alpha})
   !! on page \pageref{alpha} in the T-, the U- and the V-points
   !! ({\tt dry\_z}, {\tt dry\_u} and {\tt dry\_v}).
   !!
   !! When working with the option {\tt SLICE\_MODEL}, the water depths in the
   !! V-points are mirrored from $j=2$ to $j=1$ and $j=3$.
   !!
   !! @note
   !! use mask or not
   !! check dry calculation
   !! calculate depth in X-points - used in e.g. diffusion eq.
   !! shall z and zo be part of domain - likely yes
   !! @endnote

SUBMODULE (getm_domain) depth_update_smod

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

module SUBROUTINE depth_update(self)

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_domain), intent(inout) :: self

!  Local constants

!  Local variables
   real(real64) :: x
   integer :: i,j
!---------------------------------------------------------------------------
   call self%logs%info('depth_update()',level=1)

#define USE_MASK

!  T-points
#define TG self%T
   do j=TG%l(2),TG%u(2)
      do i=TG%l(1),TG%u(1)
#ifdef USE_MASK
      if (TG%mask(i,j) > 0) then
#endif
         TG%D(i,j) = TG%z(i,j)+TG%H(i,j)
!         dry_z(i,j)=max(0._real64,min(1._real64,(self%T%D(i,j)-_HALF_*d1)*d2i))
#ifdef USE_MASK
      end if
#endif
      end do
   end do

!  U-points
#define UG self%U
   do j=UG%l(2),UG%u(2)
      do i=UG%l(1),UG%u(1)-1
#ifdef USE_MASK
      if (UG%mask(i,j) > 0) then
#endif
         x=max(0.25_real64*(TG%zo(i,j)+TG%zo(i+1,j)+TG%z(i,j)+TG%z(i+1,j)),-UG%H(i,j)+self%Dmin)
         UG%D(i,j) = x+UG%H(i,j)
!         dry_u(i,j) = max(0._real64,min(1._real64,(self%U%DU(i,j)-d1)*d2i))
#ifdef USE_MASK
      end if
#endif
      end do
   end do
#undef UG

!  V-points
#define VG self%V
   do j=VG%l(2),VG%u(2)-1
      do i=VG%l(1),VG%u(1)
#ifdef USE_MASK
      if (VG%mask(i,j) > 0) then
#endif
         x=max(0.25_real64*(TG%zo(i,j)+TG%zo(i,j+1)+TG%z(i,j)+TG%z(i,j+1)),-VG%H(i,j)+self%Dmin)
         VG%D(i,j) = x+VG%H(i,j)
!         dry_v(i,j) = max(0._real64,min(1._real64,(self%V%DV(i,j)-d1)*d2i))
#ifdef USE_MASK
      end if
#endif
      end do
   end do
#undef VG
#undef TG

#if 0
   d1  = 2*Dmin
   d2i = _ONE_/(Dcrit-2*Dmin)
   do j=jmin-HALO,jmax+HALO
      do i=imin-HALO,imax+HALO
         if (az(i,j) .gt. 0) then
            dry_z(i,j)=max(_ZERO_,min(_ONE_,(D(i,j)-_HALF_*d1)*d2i))
         end if
         if (au(i,j) .gt. 0) then
            dry_u(i,j) = max(_ZERO_,min(_ONE_,(DU(i,j)-d1)*d2i))
         end if
         if (av(i,j) .gt. 0) then
            dry_v(i,j) = max(_ZERO_,min(_ONE_,(DV(i,j)-d1)*d2i))
         end if
     end do
  end do
#endif
   call self%logs%info('done',level=1)
   return
END SUBROUTINE depth_update

!-----------------------------------------------------------------------------

END SUBMODULE depth_update_smod
