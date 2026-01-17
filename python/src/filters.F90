module m_filters

   use iso_c_binding, only: c_int, c_double

   implicit none

contains

   subroutine c_horizontal_filter(nx, ny, nz, halox, haloy, &
                                  mask, w, var) bind(c)

     !! Simple horizontal filter where the central value is weighted by
     !! w and the average of the - active - neighbors by 1-w.

      integer(c_int), intent(in), value :: nx, ny, nz
      integer(c_int), intent(in), value :: halox, haloy
#define _D2_  -halox+1:nx+halox,-haloy+1:ny+haloy
      integer(c_int), intent(in) :: mask(_D2_)
      real(c_double), intent(in), value :: w
      real(c_double), intent(inout) :: var(_D2_,nz)
#undef _D2_

      real(c_double), allocatable :: x(:,:,:)
      integer :: n1, n2, n3, n4, n5
      integer :: rc
      integer :: imin=1, jmin=1, imax, jmax, kmax
      integer :: i, j, k

      imin = max(1, -halox+2)
      jmin = max(1, -haloy+2)
      imax = min(nx, nx+halox-1)
      jmax = min(ny, ny+haloy-1)
      kmax=nz

      allocate(x, source=var, stat=rc)
      if (rc /= 0) stop 'c_horizontal_filter: Error allocating x'
      do j = jmin, jmax
         do i = imin, imax
            if (mask(i,j) < 1) cycle

            n1=min(1,mask(i-1,j))
            n2=min(1,mask(i+1,j))
            n3=min(1,mask(i,j-1))
            n4=min(1,mask(i,j+1))
            n5=n1+n2+n3+n4
            if (n5 /= 0) then
               var(i,j,1:kmax) = (1._c_double-w)*var(i,j,1:kmax)
               if (n1 /= 0) var(i,j,1:kmax) = var(i,j,1:kmax) + (w/n5)*x(i-1,j,1:kmax)
               if (n2 /= 0) var(i,j,1:kmax) = var(i,j,1:kmax) + (w/n5)*x(i+1,j,1:kmax)
               if (n3 /= 0) var(i,j,1:kmax) = var(i,j,1:kmax) + (w/n5)*x(i,j-1,1:kmax)
               if (n4 /= 0) var(i,j,1:kmax) = var(i,j,1:kmax) + (w/n5)*x(i,j+1,1:kmax)
            end if
         end do
      end do
   end subroutine

   subroutine c_vertical_filter(nx, ny, nz, halox, haloy, &
                                nfilter, mask, w, var) bind(c)

     !! Simple vartical filter where the central value is weighted by
     !! w and the two neighbors by 1-w/2.
     !! The filter will be applied nfilter times.

      integer(c_int), intent(in), value :: nx, ny, nz
      integer(c_int), intent(in), value :: halox, haloy
      integer(c_int), intent(in), value :: nfilter
#define _D2_  -halox+1:nx+halox,-haloy+1:ny+haloy
      integer(c_int), intent(in) :: mask(_D2_)
      real(c_double), intent(in), value :: w
      real(c_double), intent(inout) :: var(_D2_, nz)
#undef _D2_

      real(c_double) :: col(nz)
      real(c_double) :: wc,wn
      integer :: imin=1, jmin=1, imax, jmax, kmax
      integer :: i, j, k, n

!----------------------------------------------------------------------
      imax=nx; jmax=ny; kmax=nz

      wc = 1._c_double-w
      wn = w/2._c_double
      do n = 1, nfilter
         do j = jmin, jmax
            do i = imin, imax
               if (mask(i,j) < 1) cycle

               !col(kmax) = var(i,j,kmax)
               !col(kmin) = var(i,j,kmin)
               col(1) = var(i,j,1)
               do k = 2, kmax-1
                  col(k) = wc*var(i,j,k)+wn*(var(i,j,k-1)+var(i,j,k+1))
               end do
               var(i,j,1:kmax-1) = col(1:kmax-1)
            end do
         end do
      end do
   end subroutine

end module
