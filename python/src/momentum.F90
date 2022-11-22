module c_momentum

   use iso_c_binding, only: c_int, c_double
   use getm_momentum, only: rho0, kappa
   use getm_domain, only: g

   implicit none

   real(c_double), parameter :: rho0i = 1.0_c_double / rho0

contains

   SUBROUTINE c_horizontal_diffusion(imin, imax, jmin, jmax, halox, haloy, umask, vmask, &
   idxu,dyu,idyv,dxv,Ah_u,Ah_v,tmask,iA,dt,f,df) bind(c)
   ! Subroutine arguments
   integer(c_int), value, intent(in) :: imin,imax,jmin,jmax
   integer(c_int), value, intent(in) :: halox
   integer(c_int), value, intent(in) :: haloy
#define _A_  imin-halox:imax+halox,jmin-haloy:jmax+haloy
   integer(c_int), intent(in) :: umask(_A_)
   integer(c_int), intent(in) :: vmask(_A_)
   real(c_double), intent(in) :: idxu(_A_)
   real(c_double), intent(in) :: dyu(_A_)
   real(c_double), intent(in) :: idyv(_A_)
   real(c_double), intent(in) :: dxv(_A_)
   real(c_double), intent(in) :: Ah_u(_A_)
   real(c_double), intent(in) :: Ah_v(_A_)
   integer(c_int), intent(in) :: tmask(_A_)
   real(c_double), intent(in) :: iA(_A_)
   real(c_double), value, intent(in) :: dt
   real(c_double), intent(inout) :: f(_A_)
   real(c_double), intent(inout) :: df(_A_)
#undef _A_

!  Local constants

!  Local variables
#ifdef _AUTOMATIC_
   real(c_double) :: uflux(imin-halox:imax+halox,jmin-haloy:jmax+haloy)
   real(c_double) :: vflux(imin-halox:imax+halox,jmin-haloy:jmax+haloy)
#else
   real(c_double), allocatable :: uflux(:,:), vflux(:,:)
#endif
   integer :: i, j
!---------------------------------------------------------------------------
#ifndef _AUTOMATIC_
   allocate(uflux(imin-halox:imax+halox,jmin-haloy:jmax+haloy))
   allocate(vflux(imin-halox:imax+halox,jmin-haloy:jmax+haloy))
#endif

   do j=jmin,jmax
      do i=imin-1,imax
         if (umask(i,j) == 1) then
            uflux(i,j) = Ah_u(i,j) * (f(i+1,j)-f(i,j)) * idxu(i,j) * dyu(i,j)
         else
            uflux(i,j) = 0.0_c_double
         end if
      end do
   end do

   do j=jmin-1,jmax
      do i=imin,imax
         if (vmask(i,j) == 1) then
            vflux(i,j) = Ah_v(i,j) * (f(i,j+1)-f(i,j)) * idyv(i,j) * dxv(i,j)
         else
            vflux(i,j) = 0.0_c_double
         end if
      end do
   end do

   do j=jmin,jmax
      do i=imin,imax
         if (tmask(i,j) == 1) then
            df(i,j) = dt * ( uflux(i,j)-uflux(i-1,j) + vflux(i,j)-vflux(i,j-1) ) * iA(i,j)
         end if
      end do
   end do
   END SUBROUTINE c_horizontal_diffusion

   SUBROUTINE c_bottom_friction(nx, ny, mask, u, v, D, z0b, z0b_min, avmmol, ru, iterate) bind(c)
      integer(c_int), value, intent(in) :: nx, ny
      integer(c_int), intent(in) :: mask(nx, ny)
      real(c_double), intent(in) :: u(nx, ny)
      real(c_double), intent(in) :: v(nx, ny)
      real(c_double), intent(in) :: D(nx, ny)
      real(c_double), intent(inout) :: z0b(nx, ny)
      real(c_double), intent(in) :: z0b_min(nx, ny)
      real(c_double), intent(in), value :: avmmol
      real(c_double), intent(out) :: ru(nx, ny)
      integer(c_int), value, intent(in) :: iterate

      real(c_double) :: sqrtcd(nx, ny)

      where (mask == 1) sqrtcd = kappa / log(1.0_c_double + 0.5_c_double * D / z0b)
      if (iterate /= 0) then
         where (mask == 1)
            z0b = min(D, z0b_min + 0.1_c_double * avmmol / max(avmmol,sqrtcd*sqrt(u*u + v*v)))
            sqrtcd = kappa / log(1.0_c_double + 0.5_c_double * D / z0b)
         end where
      end if
      where (mask == 1) ru = sqrtcd * sqrtcd * sqrt(u*u + v*v)
   END SUBROUTINE

   SUBROUTINE c_collect_3d_momentum_sources(imin, imax, jmin, jmax, kmax, halox, haloy, mask, alpha, ho, hn, &
      dp, cor, adv, diff, idp, taus, rr, dt, ea2, ea4) bind(c)
      integer(c_int), value, intent(in) :: imin, imax, jmin, jmax, kmax
      integer(c_int), value, intent(in) :: halox
      integer(c_int), value, intent(in) :: haloy
#define _A2D_  imin-halox:imax+halox,jmin-haloy:jmax+haloy
#define _A3D_  imin-halox:imax+halox,jmin-haloy:jmax+haloy,kmax
      integer(c_int), intent(in) :: mask(_A2D_)
      real(c_double), intent(in) :: alpha(_A2D_)
      real(c_double), intent(in) :: ho(_A3D_)
      real(c_double), intent(in) :: hn(_A3D_)
      real(c_double), intent(in) :: dp(_A2D_)
      real(c_double), intent(in) :: cor(_A3D_)
      real(c_double), intent(in) :: adv(_A3D_)
      real(c_double), intent(in) :: diff(_A3D_)
      real(c_double), intent(in) :: idp(_A3D_)
      real(c_double), intent(in) :: taus(_A2D_)
      real(c_double), intent(in) :: rr(_A2D_)
      real(c_double), value, intent(in) :: dt
      real(c_double), intent(inout) :: ea2(_A3D_)
      real(c_double), intent(inout) :: ea4(_A3D_)
#undef _A2D_
#undef _A3D_

      integer :: i, j, k

      do k=1,kmax
         do j=jmin,jmax
            do i=imin-halox,imax+halox
               ea4(i,j,k)=dt*(-0.5_c_double * (ho(i,j,k)+hn(i,j,k)) * g * dp(i,j) &
                              +alpha(i,j)*(cor(i,j,k) + adv(i,j,k) + diff(i,j,k) + idp(i,j,k)) &
                              )
            end do
         end do
      end do

      ! Additional matrix elements for surface and bottom layer
      do j=jmin,jmax
         do i=imin-halox,imax+halox
            ! surface stress
            ea4(i,j,kmax) = ea4(i,j,kmax) + dt * alpha(i,j) * taus(i,j) * rho0i

            ! bottom friction
            ea2(i,j,1) = -dt * rr(i,j)
         end do
      end do
   END SUBROUTINE

   SUBROUTINE c_advance_2d_transport(imin, imax, jmin, jmax, halox, haloy, mask, alpha, D, &
      dp, taus, cor, adv, diff, damp, SA, SB, SD, SF, r, dt, U) bind(c)
      integer(c_int), value, intent(in) :: imin, imax, jmin, jmax
      integer(c_int), value, intent(in) :: halox
      integer(c_int), value, intent(in) :: haloy
#define _A2D_  imin-halox:imax+halox,jmin-haloy:jmax+haloy
      integer(c_int), intent(in) :: mask(_A2D_)
      real(c_double), intent(in) :: alpha(_A2D_)
      real(c_double), intent(in) :: D(_A2D_)
      real(c_double), intent(in) :: dp(_A2D_)
      real(c_double), intent(in) :: taus(_A2D_)
      real(c_double), intent(in) :: cor(_A2D_)
      real(c_double), intent(in) :: adv(_A2D_)
      real(c_double), intent(in) :: diff(_A2D_)
      real(c_double), intent(in) :: damp(_A2D_)
      real(c_double), intent(in) :: SA(_A2D_)
      real(c_double), intent(in) :: SB(_A2D_)
      real(c_double), intent(in) :: SD(_A2D_)
      real(c_double), intent(in) :: SF(_A2D_)
      real(c_double), intent(in) :: r(_A2D_)
      real(c_double), value, intent(in) :: dt
      real(c_double), intent(inout) :: U(_A2D_)
#undef _A2D_

      integer :: i, j
      real(c_double) :: Slr

      do j=jmin,jmax
         do i=imin,imax
            if (mask(i,j) == 1) then
               if (U(i,j) > 0._c_double) then
                  Slr = min(SF(i,j), 0._c_double)
               else
                  Slr = max(SF(i,j), 0._c_double)
               end if
               ! [GETM Scientific Report: eqs. 2.14, 2.16]
               U(i,j) = (U(i,j) + dt * (-g * D(i,j) * dp(i,j) & ! note SxF is multiplied by alpha
                              + alpha(i,j) * (taus(i,j) * rho0i + cor(i,j) &
                              + adv(i,j) + diff(i,j) + damp(i,j) &
                              + SA(i,j) + SB(i,j) + SD(i,j) + Slr))) &
                             / (1.0_c_double + dt * r(i,j) / D(i,j))
            end if
         end do
      end do
   END SUBROUTINE

end module
