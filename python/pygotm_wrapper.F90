module pygotm
   use iso_c_binding, only: c_ptr, c_int, c_double, c_char, c_loc, c_f_pointer, C_NULL_CHAR, C_NULL_PTR

   use turbulence, only: init_turbulence, post_init_turbulence, do_turbulence, tke, tkeo, eps, L, num, nuh, clean_turbulence
   use mtridiagonal, only: init_tridiagonal, clean_tridiagonal
   use yaml_settings

   implicit none

   logical, save :: initialized = .false.

contains

   subroutine initialize(nlev, nml_file, ptke, ptkeo, peps, pL, pnum, pnuh) bind(c)
      integer(c_int), value,          intent(in)  :: nlev
      character(kind=c_char), target, intent(in)  :: nml_file(*)
      type(c_ptr),                    intent(out) :: ptke, ptkeo, peps, pL, pnum, pnuh

      type (type_settings) :: branch
      integer, parameter :: iunit = 60
      integer :: namelen
      character(len=256), pointer :: pnml_file

      call c_f_pointer(c_loc(nml_file), pnml_file)

      if (initialized) then
         call clean_turbulence()
         call clean_tridiagonal()
      end if
      call init_turbulence(branch)
      namelen = index(pnml_file, C_NULL_CHAR) - 1
      if (namelen > 0) call init_turbulence(iunit, pnml_file(:namelen))
      call init_tridiagonal(nlev)
      call post_init_turbulence(nlev)

      ptke = c_loc(tke)
      ptkeo = c_loc(tkeo)
      peps = c_loc(eps)
      pL   = c_loc(L)
      pnum = c_loc(num)
      pnuh = c_loc(nuh)

      initialized = .true.
   end subroutine

   subroutine calculate(nlev, dt, h, D, u_taus, u_taub, z0s, z0b, NN, SS) bind(c)
      integer(c_int), value, intent(in) :: nlev
      real(c_double), value, intent(in) :: dt, D, u_taus, u_taub, z0s, z0b
      real(c_double),        intent(in) :: h(0:nlev), SS(0:nlev), NN(0:nlev)
      call do_turbulence(nlev, dt, D, u_taus, u_taub, z0s, z0b, h, NN, SS)
   end subroutine

   subroutine calculate_3d(nx, ny, nz, istart, istop, jstart, jstop, dt, mask, h3d, D, u_taus, u_taub, z0s, z0b, NN, SS, &
         tke3d, tkeo3d, eps3d, L3d, num3d, nuh3d) bind(c)
      integer(c_int), intent(in), value                        :: nx, ny, nz, istart, istop, jstart, jstop
      real(c_double), intent(in), value                        :: dt
      integer(c_int), intent(in),    dimension(nx, ny)         :: mask
      real(c_double), intent(in),    dimension(nx, ny)         :: D, u_taus, u_taub, z0s, z0b
      real(c_double), intent(in),    dimension(nx, ny, nz)     :: h3d
      real(c_double), intent(in),    dimension(nx, ny, nz + 1) :: SS, NN
      real(c_double), intent(inout), dimension(nx, ny, nz + 1) :: tke3d, tkeo3d, eps3d, L3d, num3d, nuh3d

      real(c_double) :: h(0:nz), NN_loc(0:nz), SS_loc(0:nz)

      integer i, j
      do j = jstart, jstop
         do i = istart, istop
            if (mask(i, j) == 1) then
               tke(:) = tke3d(i, j, :)
               tkeo(:) = tkeo3d(i, j, :)
               eps(:) = eps3d(i, j, :)
               L(:) = L3d(i, j, :)
               num(:) = num3d(i, j, :)
               nuh(:) = nuh3d(i, j, :)
               h(1:nz) = h3d(i, j, :)
               NN_loc(:) = NN(i, j, :)
               SS_loc(:) = SS(i, j, :)
               call do_turbulence(nz, dt, D(i, j), u_taus(i, j), u_taub(i, j), z0s(i, j), z0b(i, j), h, NN_loc, SS_loc)
               tke3d(i, j, :) = tke(:)
               tkeo3d(i, j, :) = tkeo(:)
               eps3d(i, j, :) = eps(:)
               L3d(i, j, :) = L(:)
               num3d(i, j, :) = num(:)
               nuh3d(i, j, :) = nuh(:)
            end if
         end do
      end do
   end subroutine

   subroutine diff(nlev, dt, cnpar, posconc, h, Bcup, Bcdw, Yup, Ydw, nuY, Lsour, Qsour, Taur, Yobs, Y) bind(c)
      integer(c_int), value,             intent(in) :: nlev, posconc, Bcup, Bcdw
      real(c_double), value,             intent(in) :: dt, cnpar, Yup, Ydw
      real(c_double), dimension(0:nlev), intent(in) :: h, nuY, Lsour, Qsour, Taur, Yobs, Y
      call diff_center(nlev, dt, cnpar, posconc, h, Bcup, Bcdw, Yup, Ydw, nuY, Lsour, Qsour, Taur, Yobs, Y)
   end subroutine

end module
