! Copyright (C) 2020 Bolding & Bruggeman
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
   !! @endnote

SUBMODULE (getm_domain) depth_update_smod

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

module SUBROUTINE depth_update(self,z,zo)

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_domain), intent(inout) :: self
   real(real64), dimension(:,:), intent(in) :: z,zo

!  Local constants

!  Local variables
   real(real64) :: x
   integer :: i,j
!---------------------------------------------------------------------------
   call self%logs%info('depth_update()',level=1)

#define USE_MASK

!  T-points
   do j=self%T%l(2),self%T%u(2)
      do i=self%T%l(1),self%T%u(1)
#ifdef USE_MASK
      if (self%T%mask(i,j) > 0) then
#endif
         self%T%D(i,j) = z(i,j)+self%T%H(i,j)
!         dry_z(i,j)=max(0._real64,min(1._real64,(self%T%D(i,j)-_HALF_*d1)*d2i))
#ifdef USE_MASK
      end if
#endif
      end do
   end do

!  U-points
   do j=self%U%l(2),self%U%u(2)
      do i=self%U%l(1),self%U%u(1)-1
#ifdef USE_MASK
      if (self%U%mask(i,j) > 0) then
#endif
         x=max(0.25_real64*(zo(i,j)+zo(i+1,j)+z(i,j)+z(i+1,j)),-self%U%H(i,j)+min_depth)
         self%U%D(i,j) = x+self%U%H(i,j)
!         dry_u(i,j) = max(0._real64,min(1._real64,(self%U%DU(i,j)-d1)*d2i))
#ifdef USE_MASK
      end if
#endif
      end do
   end do

!  V-points
   do j=self%V%l(2),self%V%u(2)-1
      do i=self%V%l(1),self%V%u(1)
#ifdef USE_MASK
      if (self%V%mask(i,j) > 0) then
#endif
         x=max(0.25_real64*(zo(i,j)+zo(i,j+1)+z(i,j)+z(i,j+1)),-self%V%H(i,j)+min_depth)
         self%V%D(i,j) = x+self%V%H(i,j)
!         dry_v(i,j) = max(0._real64,min(1._real64,(self%V%DV(i,j)-d1)*d2i))
#ifdef USE_MASK
      end if
#endif
      end do
   end do

#if 0
   d1  = 2*min_depth
   d2i = _ONE_/(crit_depth-2*min_depth)
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
