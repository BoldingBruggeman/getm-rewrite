module advection_base

   USE, INTRINSIC :: ISO_FORTRAN_ENV

   implicit none

   private

   type, abstract, public :: type_advection_base
   contains
      procedure(u2d), deferred, nopass :: u2d
      procedure(v2d), deferred, nopass :: v2d
      procedure(w3d), deferred, nopass :: w3d
   end type

   interface
SUBROUTINE u2d(imin,imax,jmin,jmax,halo,umask,dxu,dyu,hu,u,tmask,iA,dt,h,f)
   import real64
   ! Subroutine arguments
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

SUBROUTINE v2d(imin,imax,jmin,jmax,halo,vmask,dxv,dyv,hv,v,tmask,iA,dt,h,f)
   import real64
   ! Subroutine arguments
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

SUBROUTINE w3d(imin,imax,jmin,jmax,kmax,halo,w,tmask,dt,h,f)
   import real64
   ! Subroutine arguments
   integer, intent(in) :: imin,imax,jmin,jmax,kmax
   integer, intent(in) :: halo(2)
#define _2D_  imin-halo(1):imax+halo(1),jmin-halo(2):jmax+halo(2)
#define _3D_  imin-halo(1):imax+halo(1),jmin-halo(2):jmax+halo(2),0:kmax
   real(real64), intent(in) :: w(_3D_)
   integer, intent(in) :: tmask(_2D_)
   real(real64), intent(in) :: dt
   real(real64), target, intent(inout) :: h(_3D_)
   real(real64), target, intent(inout) :: f(_3D_)
#undef _2D_
#undef _3D_
END SUBROUTINE

   end interface

end module
