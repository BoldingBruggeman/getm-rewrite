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
      character(len=256), pointer :: pnml_file

      call c_f_pointer(c_loc(nml_file), pnml_file)

      if (initialized) then
         call clean_turbulence()
         call clean_tridiagonal()
      end if

      call init_turbulence(branch)
      if (pnml_file /= '') call init_turbulence(iunit, pnml_file)
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

   subroutine diff(nlev, dt, cnpar, posconc, h, Bcup, Bcdw, Yup, Ydw, nuY, Lsour, Qsour, Taur, Yobs, Y) bind(c)
      integer(c_int), value,             intent(in) :: nlev, posconc, Bcup, Bcdw
      real(c_double), value,             intent(in) :: dt, cnpar, Yup, Ydw
      real(c_double), dimension(0:nlev), intent(in) :: h, nuY, Lsour, Qsour, Taur, Yobs, Y
      call diff_center(nlev, dt, cnpar, posconc, h, Bcup, Bcdw, Yup, Ydw, nuY, Lsour, Qsour, Taur, Yobs, Y)
   end subroutine

end module
