! Copyright (C) 2024 Bolding & Bruggeman and Hans Burchard

module m_tridiagonal

   use iso_c_binding, only: c_int, c_double

   implicit none

!#define _TIMING_
#ifdef _TIMING_
   real(c_double) :: start_, stop_, time_=0._c_double
#endif

contains

#define _USE_3D_
! doing a simple test using 3D local variables is about 4% faster then
! using 1D - KB 2024-12-16

subroutine c_tridiagonal(nx, ny, nz, halox, haloy, dt, mask, nu, var) bind(c)

!  Subroutine arguments
   integer(c_int), intent(in), value :: nx, ny, nz
   integer(c_int), intent(in), value :: halox, haloy
   real(c_double), intent(in), value :: dt
#define _D2_ -halox+1:nx+halox,-haloy+1:ny+haloy
   integer(c_int), intent(in) :: mask(_D2_)
   real(c_double), intent(in) :: nu(_D2_,1:nz)
   real(c_double), intent(inout) :: var(_D2_,0:nz)
#undef _D2_

 !  Local variables
#ifdef _USE_3D_
   real(c_double),allocatable :: a1(:,:,:)
   real(c_double),allocatable :: a2(:,:,:)
   real(c_double),allocatable :: a3(:,:,:)
   real(c_double),allocatable :: a4(:,:,:)
#else
   real(c_double) :: a1(0:nz)
   real(c_double) :: a2(0:nz)
   real(c_double) :: a3(0:nz)
   real(c_double) :: a4(0:nz)
#endif
   real(c_double) :: ru(0:nz)
   real(c_double) :: qu(0:nz)
   integer :: imin=1, jmin=1, imax, jmax, kmax
   integer(c_int) :: i,j,k
   real(c_double) :: x

!-----------------------------------------------------------------------
   imax=nx; jmax=ny; kmax=nz

#ifdef _USE_3D_
#define _D2_ -halox+1:nx+halox,-haloy+1:ny+haloy
   allocate(a1(_D2_,0:nz))
   allocate(a2(_D2_,0:nz))
   allocate(a3(_D2_,0:nz))
   allocate(a4(_D2_,0:nz))
#undef _D2_
#endif

#ifdef _TIMING_
   call cpu_time(start_)
#endif
#ifdef _USE_3D_
   ! picked directly from Bjarnes code
   a2(:,:,kmax) = 1._c_double
   a1(:,:,kmax) = 0._c_double
   a4(:,:,kmax) = 0._c_double
   x = 2._c_double*dt
   do k=1,kmax-1
      do j=jmin,jmax
         do i=imin,imax
            !if (mask(i,j) /= 1) cycle
            if (mask(i,j) < 1) cycle

            a1(i,j,k) = -x*nu(i,j,k)
            a3(i,j,k) = -x*nu(i,j,k+1)
            a2(i,j,k) = 1._c_double + x*nu(i,j,k) + x*nu(i,j,k+1)
            a4(i,j,k) = var(i,j,k)
         end do
      end do
   end do
   a3(:,:,0) = 0._c_double
   a2(:,:,0) = 1._c_double
   a4(:,:,0) = -1._c_double

   ! solve system
   do j=jmin,jmax
      do i=imin,imax
         !if (mask(i,j) /= 1) cycle
         if (mask(i,j) < 1) cycle

         ru(kmax)=a1(i,j,kmax)/a2(i,j,kmax)
         qu(kmax)=a4(i,j,kmax)/a2(i,j,kmax)
         do k=kmax-1,1,-1
            ru(k)=a1(i,j,k)/(a2(i,j,k)-a3(i,j,k)*ru(k+1))
            qu(k)=(a4(i,j,k)-a3(i,j,k)*qu(k+1))/(a2(i,j,k)-a3(i,j,k)*ru(k+1))
         end do
         qu(0)=(a4(i,j,0)-a3(i,j,0)*qu(1)) &
                   /(a2(i,j,0)-a3(i,j,0)*ru(1))
         var(i,j,0)=qu(0)
         var(i,j,0) = -1._c_double
         do k=1,kmax
            var(i,j,k)=qu(k)-ru(k)*var(i,j,k-1)
         end do
         if (var(i,j,kmax) > 0._c_double) then
            write(89,*) i,j,k,var(i,j,kmax)
         end if
         !var(i,j,kmax) = 0._c_double
      end do
   end do

#else

   ! picked directly from Bjarnes code
   a2(kmax) = 1._c_double
   a1(kmax) = 0._c_double
   a4(kmax) = 0._c_double
   a3(0) = 0._c_double
   a2(0) = 1._c_double
   a4(0) = -1._c_double
   x = 2._c_double*dt
   do j=jmin,jmax
      do i=imin,imax
         !if (mask(i,j) /= 1) cycle
         if (mask(i,j) /= 1) cycle

         do k=1,kmax-1
            a1(k) = -x*nu(i,j,k)
            a3(k) = -x*nu(i,j,k+1)
            a2(k) = 1._c_double + x*nu(i,j,k) + x*nu(i,j,k+1)
            a4(k) = var(i,j,k)
         end do

         ! solve system
         ru(kmax)=a1(kmax)/a2(kmax)
         qu(kmax)=a4(kmax)/a2(kmax)
         do k=kmax-1,1,-1
            ru(k)=a1(k)/(a2(k)-a3(k)*ru(k+1))
            qu(k)=(a4(k)-a3(k)*qu(k+1))/(a2(k)-a3(k)*ru(k+1))
         end do
         qu(0)=(a4(0)-a3(0)*qu(1))/(a2(0)-a3(0)*ru(1))
         var(i,j,0)=qu(0)
         do k=1,kmax-1
            var(i,j,k)=qu(k)-ru(k)*var(i,j,k-1)
         end do
         var(i,j,kmax) = 0._c_double
      end do
   end do
#endif
#ifdef _TIMING_
   call cpu_time(stop_)
   time_ = time_+(stop_-start_)
#ifdef _USE_3D_
   write(9,*) '3D ',time_
#else
   write(9,*) '1D ',time_
#endif
#endif

end subroutine c_tridiagonal

!-----------------------------------------------------------------------

end module
