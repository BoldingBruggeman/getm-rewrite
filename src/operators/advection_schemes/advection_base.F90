module advection_base

   USE, INTRINSIC :: ISO_FORTRAN_ENV

   implicit none

   private

   type, abstract, public :: type_advection_base
   contains
      procedure(u2d), deferred :: u2d
      procedure(v2d), deferred :: v2d
      procedure(w3d), deferred :: w3d
   end type

   interface
SUBROUTINE u2d(self,imin,imax,jmin,jmax,halo,umask,dxu,dyu,hu,u,tmask,iA,dt,h,f)
   import type_advection_base, real64
   ! Subroutine arguments
   class (type_advection_base), intent(inout) :: self
   integer, intent(in) :: imin,imax,jmin,jmax
   integer, intent(in) :: halo(2)
#define _A_  imin-halo(1):imax+halo(1),jmin-halo(2):jmax+halo(2)
   integer, intent(in) :: umask(_A_)
   real(real64), intent(in) :: dxu(_A_)
   real(real64), intent(in) :: dyu(_A_)
   real(real64), intent(in) :: hu(_A_)
   real(real64), intent(in) :: u(_A_)
   integer, intent(in) :: tmask(_A_)
   real(real64), intent(in) :: iA(_A_)
   real(real64), intent(in) :: dt
   real(real64), intent(inout) :: h(_A_)
   real(real64), intent(inout) :: f(_A_)
#undef _A_
END SUBROUTINE

SUBROUTINE v2d(self,imin,imax,jmin,jmax,halo,vmask,dxv,dyv,hv,v,tmask,iA,dt,h,f)
   import type_advection_base, real64
   ! Subroutine arguments
   class (type_advection_base), intent(inout) :: self
   integer, intent(in) :: imin,imax,jmin,jmax
   integer, intent(in) :: halo(2)
#define _A_  imin-halo(1):imax+halo(1),jmin-halo(2):jmax+halo(2)
   integer, intent(in) :: vmask(_A_)
   real(real64), intent(in) :: dxv(_A_)
   real(real64), intent(in) :: dyv(_A_)
   real(real64), intent(in) :: hv(_A_)
   real(real64), intent(in) :: v(_A_)
   integer, intent(in) :: tmask(_A_)
   real(real64), intent(in) :: iA(_A_)
   real(real64), intent(in) :: dt
   real(real64), intent(inout) :: h(_A_)
   real(real64), intent(inout) :: f(_A_)
#undef _A_
END SUBROUTINE

SUBROUTINE w3d(self,imin, imax, jmin, jmax, kmax, w , tmask, dt, h, f)
   import type_advection_base, real64
   ! Subroutine arguments
   class (type_advection_base), intent(inout) :: self
   integer, intent(in) :: imin, imax, jmin, jmax, kmax
   real(real64), intent(in) :: w(imin-1:imax+1,jmin-1:jmax+1,0:kmax)
   integer, intent(in) :: tmask(imin-1:imax+1,jmin-1:jmax+1)
   real(real64), intent(in) :: dt
   real(real64), target, intent(inout) :: h(imin-1:imax+1,jmin-1:jmax+1,0:kmax), f(imin-1:imax+1,jmin-1:jmax+1,0:kmax)
END SUBROUTINE

   end interface

end module
