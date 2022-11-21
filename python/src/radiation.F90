module c_radiation

   use iso_c_binding, only: c_ptr, c_int, c_double

   implicit none

contains

   subroutine c_exponential_profile_1band_centers(nx, ny, nz, istart, istop, jstart, jstop, mask, h, kc, top, out) bind(c)
      integer(c_int), intent(in), value :: nx, ny, nz, istart, istop, jstart, jstop
      integer(c_int), intent(in)        :: mask(nx, ny)
      real(c_double), intent(in)        :: h(nx, ny, nz), kc(nx, ny, nz), top(nx, ny)
      real(c_double), intent(inout)     :: out(nx, ny, nz)

      integer :: k
      real(c_double) :: hkc(nx, ny), cumkc(nx, ny)

      cumkc = 0._c_double
      do k=nz,1,-1
         hkc = h(:,:,k) * kc(:,:,k)
         where (mask /= 0) out(:,:,k) = top * exp(-(cumkc + 0.5_c_double * hkc))
         cumkc = cumkc + hkc
      end do
   end subroutine

   subroutine c_exponential_profile_1band_interfaces(nx, ny, nz, istart, istop, jstart, jstop, mask, &
      h, kc, initial, up, out) bind(c)
      integer(c_int), intent(in), value :: nx, ny, nz, istart, istop, jstart, jstop, up
      integer(c_int), intent(in)        :: mask(nx, ny)
      real(c_double), intent(in)        :: h(nx, ny, nz), kc(nx, ny, nz), initial(nx, ny)
      real(c_double), intent(inout)     :: out(nx, ny, 0:nz)

      integer :: k, kstart, kstop, kstep, h_offset
      real(c_double) :: cumkc(nx, ny)

      if (up /= 0) then
         kstart = 0
         kstop = nz
         kstep = 1
         h_offset = 0
      else
         kstart = nz
         kstop = 0
         kstep = -1
         h_offset = 1
      end if

      where (mask /= 0) out(:,:,kstart) = initial
      cumkc = 0._c_double
      do k = kstart + kstep, kstop, kstep
         cumkc = cumkc + kc(:,:,k+h_offset) * h(:,:,k+h_offset)
         where (mask /= 0) out(:,:,k) = initial * exp(-cumkc)
      end do
   end subroutine

   subroutine c_exponential_profile_2band_interfaces(nx, ny, nz, istart, istop, jstart, jstop, mask, &
      h, f1, kc1, kc2, initial, up, out) bind(c)
      integer(c_int), intent(in), value :: nx, ny, nz, istart, istop, jstart, jstop, up
      integer(c_int), intent(in)        :: mask(nx, ny)
      real(c_double), intent(in)        :: h(nx, ny, nz), f1(nx, ny), kc1(nx, ny), kc2(nx, ny), initial(nx, ny)
      real(c_double), intent(inout)     :: out(nx, ny, 0:nz)

      integer :: k, kstart, kstop, kstep, h_offset
      real(c_double) :: cumh(nx, ny)

      if (up /= 0) then
         kstart = 0
         kstop = nz
         kstep = 1
         h_offset = 0
      else
         kstart = nz
         kstop = 0
         kstep = -1
         h_offset = 1
      end if

      where (mask /= 0) out(:,:,kstart) = initial
      cumh = 0._c_double
      do k = kstart + kstep, kstop, kstep
         cumh = cumh + h(:,:,k+h_offset)
         where (mask /= 0) out(:,:,k) = initial * (f1  * exp(-kc1 * cumh) + (1._c_double - f1) * exp(-kc2 * cumh))
      end do
   end subroutine

end module
