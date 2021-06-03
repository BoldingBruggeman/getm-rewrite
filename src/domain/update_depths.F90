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
   !! note loop boundaries when using [UVX]G%z directly
   !! @endnote

SUBMODULE (getm_domain) update_depths_smod

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

module SUBROUTINE update_depths(self)

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_domain), intent(inout) :: self

!  Local constants

!  Local variables
   real(real64) :: x
   integer :: i,j
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('update_depths()',level=2)

#define USE_MASK
   TGrid: associate( TG => self%T )
   do j=TG%l(2),TG%u(2)
      do i=TG%l(1),TG%u(1)
#ifdef USE_MASK
      if (TG%mask(i,j) > 0) then
#endif
         TG%D(i,j) = TG%H(i,j)+TG%z(i,j)
         TG%alpha(i,j)=max(0._real64,min(1._real64,(self%T%D(i,j)-self%Dmin)/(self%Dcrit-self%Dmin)))
#ifdef USE_MASK
      end if
#endif
      end do
   end do
   call self%mirror_bdys(TG,TG%D)

   UGrid: associate( UG => self%U )
   do j=UG%l(2),UG%u(2)
      do i=UG%l(1),UG%u(1)
#ifdef USE_MASK
      if (UG%mask(i,j) > 0) then
#endif
         UG%D(i,j) = UG%H(i,j)+UG%z(i,j)
         UG%alpha(i,j)=max(0._real64,min(1._real64,(self%U%D(i,j)-self%Dmin)/(self%Dcrit-self%Dmin)))
#ifdef USE_MASK
      end if
#endif
      end do
   end do
   call self%mirror_bdys(UG,UG%D)
   end associate UGrid

   VGrid: associate( VG => self%V )
   do j=VG%l(2),VG%u(2)
      do i=VG%l(1),VG%u(1)
#ifdef USE_MASK
      if (VG%mask(i,j) > 0) then
#endif
         VG%D(i,j) = VG%H(i,j)+VG%z(i,j)
         VG%alpha(i,j)=max(0._real64,min(1._real64,(self%V%D(i,j)-self%Dmin)/(self%Dcrit-self%Dmin)))
#ifdef USE_MASK
      end if
#endif
      end do
   end do
   call self%mirror_bdys(VG,VG%D)
   end associate VGrid

   XGrid: associate( XG => self%X )
   do j=XG%l(2),XG%u(2)
      do i=XG%l(1),XG%u(1)
#ifdef USE_MASK
      if (XG%mask(i,j) > 0) then
#endif
         XG%D(i,j) = XG%H(i,j)+XG%z(i,j)
         XG%alpha(i,j)=1._real64 !KB
#ifdef USE_MASK
      end if
#endif
      end do
   end do
!KB - not sure if necessary (or correct)   call self%mirror_bdys(XG,XG%D)
   end associate XGrid
   end associate TGrid
END SUBROUTINE update_depths

!-----------------------------------------------------------------------------

END SUBMODULE update_depths_smod
