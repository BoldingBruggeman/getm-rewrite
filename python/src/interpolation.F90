module c_interpolation

   use iso_c_binding, only: c_int, c_double

   implicit none

contains

   subroutine grid_interp_x(nx, ny, nz, source, target, ioffset) bind(c)
      integer(c_int), intent(in), value :: nx, ny, nz, ioffset
      real(c_double), intent(in)        :: source(nx, ny, nz)
      real(c_double), intent(inout)     :: target(nx, ny, nz)

      integer :: i, j, k

      do k = 1, nz
         do j = 1, ny
            do i = 1, nx - 1
               target(i + ioffset, j, k) = 0.5_c_double * (source(i, j, k) + source(i + 1, j, k))
            end do
         end do
      end do
   end subroutine

   subroutine grid_interp_y(nx, ny, nz, source, target, joffset) bind(c)
      integer(c_int), intent(in), value :: nx, ny, nz, joffset
      real(c_double), intent(in)        :: source(nx, ny, nz)
      real(c_double), intent(inout)     :: target(nx, ny, nz)

      integer :: i, j, k

      do k = 1, nz
         do j = 1, ny - 1
            do i = 1, nx
               target(i, j + joffset, k) = 0.5_c_double * (source(i, j, k) + source(i, j + 1, k))
            end do
         end do
      end do
   end subroutine

   subroutine grid_interp_z(nx, ny, nz1, nz2, source, target, koffset) bind(c)
      integer(c_int), intent(in), value :: nx, ny, nz1, nz2, koffset
      real(c_double), intent(in)        :: source(nx, ny, nz1)
      real(c_double), intent(inout)     :: target(nx, ny, nz2)

      integer :: i, j, k

      do k = 1, nz1 - 1
         do j = 1, ny
            do i = 1, nx
               target(i, j, k + koffset) = 0.5_c_double * (source(i, j, k) + source(i, j, k + 1))
            end do
         end do
      end do
   end subroutine

   subroutine grid_interp_xy(nx1, ny1, nx2, ny2, nz, source, target, ioffset, joffset) bind(c)
      integer(c_int), intent(in), value :: nx1, ny1, nx2, ny2, nz, ioffset, joffset
      real(c_double), intent(in)        :: source(nx1, ny1, nz)
      real(c_double), intent(inout)     :: target(nx2, ny2, nz)

      integer :: i, j, k

      do k = 1, nz
         do j = 1, ny1 - 1
            do i = 1, nx1 - 1
               target(i + ioffset, j + joffset, k) = 0.25_c_double * (source(i, j, k) + source(i + 1, j, k) &
                  + source(i, j + 1, k) + source(i + 1, j + 1, k))
            end do
         end do
      end do
   end subroutine

end module
