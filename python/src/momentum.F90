module c_momentum

   use iso_c_binding, only: c_int, c_double
   use getm_momentum, only: rho0, kappa
   use getm_domain, only: g

   implicit none

   real(c_double), parameter :: rho0i = 1.0_c_double / rho0

contains

   subroutine c_horizontal_diffusion(imin, imax, jmin, jmax, halox, haloy, umask, vmask, &
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
   end subroutine

   subroutine c_bottom_friction(nx, ny, mask, u, v, D, z0b, z0b_min, avmmol, ru, iterate) bind(c)
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
   end subroutine

   subroutine c_collect_3d_momentum_sources(nx, ny, nz, halox, haloy, mask, alpha, ho, hn, &
      dp, cor, adv, diff, idp, taus, rr, dt, ea2, ea4) bind(c)
      integer(c_int), value, intent(in) :: nx, ny, nz
      integer(c_int), value, intent(in) :: halox
      integer(c_int), value, intent(in) :: haloy
      integer(c_int), intent(in) :: mask(nx, ny)
      real(c_double), intent(in) :: alpha(nx, ny)
      real(c_double), intent(in) :: ho(nx, ny, nz)
      real(c_double), intent(in) :: hn(nx, ny, nz)
      real(c_double), intent(in) :: dp(nx, ny)
      real(c_double), intent(in) :: cor(nx, ny, nz)
      real(c_double), intent(in) :: adv(nx, ny, nz)
      real(c_double), intent(in) :: diff(nx, ny, nz)
      real(c_double), intent(in) :: idp(nx, ny, nz)
      real(c_double), intent(in) :: taus(nx, ny)
      real(c_double), intent(in) :: rr(nx, ny)
      real(c_double), value, intent(in) :: dt
      real(c_double), intent(inout) :: ea2(nx, ny, nz)
      real(c_double), intent(inout) :: ea4(nx, ny, nz)

      integer :: i, j, k

      do k = 1, nz
         do j = 1+haloy, ny-haloy
            do i = 1, nx
               ea4(i,j,k)=dt*(-0.5_c_double * (ho(i,j,k)+hn(i,j,k)) * g * dp(i,j) &
                              +alpha(i,j)*(cor(i,j,k) + adv(i,j,k) + diff(i,j,k) + idp(i,j,k)) &
                              )
            end do
         end do
      end do

      ! Additional matrix elements for surface and bottom layer
      do j = 1+haloy, ny-haloy
         do i = 1, nx
            ! surface stress
            ea4(i,j,nz) = ea4(i,j,nz) + dt * alpha(i,j) * taus(i,j) * rho0i

            ! bottom friction
            ea2(i,j,1) = -dt * rr(i,j)
         end do
      end do
   end subroutine

   subroutine c_advance_2d_transport(nx, ny, halox, haloy, mask, alpha, D, &
      dp, taus, cor, adv, diff, damp, SA, SB, SD, SF, r, dt, U) bind(c)
      integer(c_int), value, intent(in) :: nx, ny
      integer(c_int), value, intent(in) :: halox
      integer(c_int), value, intent(in) :: haloy
      integer(c_int), intent(in) :: mask(nx, ny)
      real(c_double), intent(in) :: alpha(nx, ny)
      real(c_double), intent(in) :: D(nx, ny)
      real(c_double), intent(in) :: dp(nx, ny)
      real(c_double), intent(in) :: taus(nx, ny)
      real(c_double), intent(in) :: cor(nx, ny)
      real(c_double), intent(in) :: adv(nx, ny)
      real(c_double), intent(in) :: diff(nx, ny)
      real(c_double), intent(in) :: damp(nx, ny)
      real(c_double), intent(in) :: SA(nx, ny)
      real(c_double), intent(in) :: SB(nx, ny)
      real(c_double), intent(in) :: SD(nx, ny)
      real(c_double), intent(in) :: SF(nx, ny)
      real(c_double), intent(in) :: r(nx, ny)
      real(c_double), value, intent(in) :: dt
      real(c_double), intent(inout) :: U(nx, ny)

      integer :: i, j
      real(c_double) :: Slr

      do j = 1+haloy, ny-haloy
         do i = 1+halox, nx-halox
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
   end subroutine

   subroutine c_w_momentum_3d(nx, ny, nz, imin, imax, jmin, jmax, mask, dyu, dxv, iarea, ho, hn, pk, qk, dt, w) bind(c)
      integer(c_int), value, intent(in) :: nx, ny, nz, imin, imax, jmin, jmax
      integer(c_int), intent(in) :: mask(nx, ny)
      real(c_double), intent(in) :: dyu(nx, ny), dxv(nx, ny), iarea(nx, ny)
      real(c_double), intent(in) :: ho(nx, ny, nz), hn(nx, ny, nz)
      real(c_double), intent(in) :: pk(nx, ny, nz), qk(nx, ny, nz)
      real(c_double), value, intent(in) :: dt
      real(c_double), intent(inout) :: w(nx, ny, 0:nz)

      integer :: i,j,k
      real(c_double) :: idt

      idt = 1.0_c_double / dt
      do k=1,nz-1
         do j=jmin,jmax
            do i=imin,imax
               if (mask(i,j) == 1) then
                     w(i,j,k) = w(i,j,k-1)-(hn(i,j,k)-ho(i,j,k))*idt &
                                 -(pk(i,j,k)*dyu(i,j)-pk(i-1,j  ,k)*dyu(i-1,j) &
                                  +qk(i,j,k)*dxv(i,j)-qk(i  ,j-1,k)*dxv(i,j-1)) &
                                  *iarea(i,j)
               end if
            end do
         end do
      end do
   end subroutine

   subroutine c_shear_frequency(nx, ny, nz, imin, imax, jmin, jmax, mask, h, hu, hv, uk, vk, num, SS) bind(c)
      integer(c_int), value, intent(in) :: nx, ny, nz, imin, imax, jmin, jmax
      integer(c_int), intent(in) :: mask(nx, ny)
      real(c_double), intent(in) :: h(nx, ny, nz), hu(nx, ny, nz), hv(nx, ny, nz)
      real(c_double), intent(in) :: uk(nx, ny, nz), vk(nx, ny, nz), num(nx, ny, 0:nz)
      real(c_double), intent(inout) :: SS(nx, ny, 0:nz)

      integer :: i,j,k

      do k = 1, nz-1
         do j = jmin, jmax
            do i = imin, imax
               if (mask(i,j) > 0) then
#ifndef NEW_SS
                  ! This is an older version which we should keep here.
                  SS(i,j,k) = 0.5_c_double * ( &
                       ((uk(i  ,j,k+1)-uk(i  ,j,k)) / (0.5_c_double*(hu(i  ,j,k+1)+hu(i  ,j,k))))**2 &
                     + ((uk(i-1,j,k+1)-uk(i-1,j,k)) / (0.5_c_double*(hu(i-1,j,k+1)+hu(i-1,j,k))))**2 &
                     + ((vk(i,j  ,k+1)-vk(i,j  ,k)) / (0.5_c_double*(hv(i,j  ,k+1)+hv(i,j  ,k))))**2 &
                     + ((vk(i,j-1,k+1)-vk(i,j-1,k)) / (0.5_c_double*(hv(i,j-1,k+1)+hv(i,j-1,k))))**2 &
                                         )
#else
                  ! This version should better conserve energy.
                  SS(i,j,k) = 0.5_c_double* &
                     ( &
                       (uk(i  ,j,k+1)-uk(i  ,j,k))**2 / (hu(i  ,j,k+1)+hu(i  ,j,k)) * (num(i,  j,  k)+num(i+1,j,  k)) &
                     + (uk(i-1,j,k+1)-uk(i-1,j,k))**2 / (hu(i-1,j,k+1)+hu(i-1,j,k)) * (num(i-1,j,  k)+num(i,  j,  k)) &
                     + (vk(i,j  ,k+1)-vk(i,j  ,k))**2 / (hv(i,j  ,k+1)+hv(i,j  ,k)) * (num(i,  j,  k)+num(i,  j+1,k)) &
                     + (vk(i,j-1,k+1)-vk(i,j-1,k))**2 / (hv(i,j-1,k+1)+hv(i,j-1,k)) * (num(i,  j-1,k)+num(i,  j,  k)) &
                     ) / (0.5_c_double*(h(i,j,k)+h(i,j,k+1))) / num(i,j,k)
#endif
               end if
            end do
         end do
      end do
   end subroutine

end module
