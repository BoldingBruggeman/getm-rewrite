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
   !! dry_{z,u,v} ---> alpha (consistent with old GETM documentation)
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
   TGrid: associate( TG => self%T )
   do j=TG%l(2),TG%u(2)
      do i=TG%l(1),TG%u(1)
#ifdef USE_MASK
      if (TG%mask(i,j) > 0) then
#endif
         TG%D(i,j) = TG%z(i,j)+TG%H(i,j)
         TG%alpha(i,j)=max(0._real64,min(1._real64,(self%T%D(i,j)-self%Dmin)/(self%Dcrit-self%Dmin)))
#ifdef USE_MASK
      end if
#endif
      end do
   end do

   UGrid: associate( UG => self%U )
   do j=UG%l(2),UG%u(2)
      do i=UG%l(1),UG%u(1)-1
#ifdef USE_MASK
      if (UG%mask(i,j) > 0) then
#endif
         x=max(0.25_real64*(TG%zo(i,j)+TG%zo(i+1,j)+TG%z(i,j)+TG%z(i+1,j)),-UG%H(i,j)+self%Dmin)
         UG%D(i,j) = x+UG%H(i,j)
         UG%alpha(i,j)=max(0._real64,min(1._real64,(self%U%D(i,j)-self%Dmin)/(self%Dcrit-self%Dmin)))
#ifdef USE_MASK
      end if
#endif
      end do
   end do

   VGrid: associate( VG => self%V )
   do j=VG%l(2),VG%u(2)-1
      do i=VG%l(1),VG%u(1)
#ifdef USE_MASK
      if (VG%mask(i,j) > 0) then
#endif
         x=max(0.25_real64*(TG%zo(i,j)+TG%zo(i,j+1)+TG%z(i,j)+TG%z(i,j+1)),-VG%H(i,j)+self%Dmin)
         VG%D(i,j) = x+VG%H(i,j)
         VG%alpha(i,j)=max(0._real64,min(1._real64,(self%V%D(i,j)-self%Dmin)/(self%Dcrit-self%Dmin)))
#ifdef USE_MASK
      end if
#endif
      end do
   end do

   XGrid: associate( XG => self%X )
   do j=XG%l(2)+1,XG%u(2)
      do i=XG%l(1)+1,XG%u(1)
#ifdef USE_MASK
      if (XG%mask(i,j) > 0) then
#endif
         XG%D(i,j)=max(0.25_real64*(UG%D(i,j)+UG%D(i+1,j)+VG%D(i,j)+VG%D(i,j+1)),-XG%H(i,j)+self%Dmin)
!KB
         XG%alpha(i,j)=1._real64
#ifdef USE_MASK
      end if
#endif
      end do
   end do
   end associate XGrid
   end associate VGrid
   end associate UGrid
   end associate TGrid
   call self%logs%info('done',level=1)
END SUBROUTINE depth_update

!-----------------------------------------------------------------------------

END SUBMODULE depth_update_smod
