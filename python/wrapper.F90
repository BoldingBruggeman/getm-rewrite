module pygetm

   use iso_c_binding, only: c_ptr, c_int, c_double, c_char, c_loc, c_f_pointer
   use iso_fortran_env, only: real64

   use getm_domain, only: type_getm_domain, type_getm_grid
   use getm_operators, only: type_advection
   use memory_manager

   implicit none

contains

   function domain_create(imin, imax, jmin, jmax, kmin, kmax, halox, haloy, haloz) result(pdomain) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: domain_create
      integer(c_int), intent(in), value :: imin, imax, jmin, jmax, kmin, kmax
      integer(c_int), intent(out)       :: halox, haloy, haloz
      type(c_ptr)                       :: pdomain

      type (type_getm_domain), pointer :: domain
      integer                          :: stat 

      allocate(domain)
      call domain%configure(imin=imin,imax=imax,jmin=jmin,jmax=jmax,kmin=kmin,kmax=kmax)
      pdomain = c_loc(domain)
      domain%domain_type = 1
      halox = imin - lbound(domain%T%c1, 1)
      haloy = jmin - lbound(domain%T%c2, 1)
      haloz = 0
   end function

   function domain_get_grid(pdomain, grid_type) result(pgrid) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: domain_get_grid
      type(c_ptr),            intent(in), value :: pdomain
      character(kind=c_char), intent(in), value :: grid_type
      type(c_ptr)                               :: pgrid

      type (type_getm_domain), pointer :: domain
      type (type_getm_grid),   pointer :: grid

      call c_f_pointer(pdomain, domain)
      select case (grid_type)
      case ('T'); grid => domain%T
      case ('U'); grid => domain%U
      case ('V'); grid => domain%V
      case ('X'); grid => domain%X
      end select
      pgrid = c_loc(grid)
   end function

   subroutine grid_get_arrays(pgrid, pc1, pc2, px, py, pdx, pdy, plon, plat, pdlon, pdlat, parea, pinv_area, pH, pD, pmask) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: grid_get_arrays
      type(c_ptr), intent(in), value :: pgrid
      type(c_ptr), intent(out)       :: pc1, pc2, px, py, pdx, pdy, plon, plat, pdlon, pdlat, parea, pinv_area, pH, pD, pmask

      type (type_getm_grid), pointer :: grid

      call c_f_pointer(pgrid, grid)

      pc1 = c_loc(grid%c1)
      pc2 = c_loc(grid%c2)
      px = c_loc(grid%x)
      py = c_loc(grid%y)
      pdx = c_loc(grid%dx)
      pdy = c_loc(grid%dy)
      plon = c_loc(grid%lon)
      plat = c_loc(grid%lat)
      pdlon = c_loc(grid%dlon)
      pdlat = c_loc(grid%dlat)
      parea = c_loc(grid%area)
      pinv_area = c_loc(grid%inv_area)
      pH = c_loc(grid%H)
      pD = c_loc(grid%D)
      pmask = c_loc(grid%mask)
   end subroutine

   subroutine domain_initialize(pdomain) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: domain_initialize
      type(c_ptr), intent(in), value :: pdomain

      type (type_getm_domain), pointer :: domain

      call c_f_pointer(pdomain, domain)
      call domain%initialize()
   end subroutine

   function advection_create() result(padvection) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: advection_create
      type(c_ptr) :: padvection

      type (type_advection), pointer :: advection

      allocate(advection)
      padvection = c_loc(advection)
   end function

   subroutine advection_calculate(padvection, scheme, pdomain,  pu, pv, timestep, pvar) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: advection_calculate
      integer(c_int), intent(in), value :: scheme
      real(c_double), intent(in), value :: timestep
      type(c_ptr),    intent(in), value :: padvection, pdomain, pu, pv, pvar

      type (type_advection),    pointer                 :: advection
      type (type_getm_domain),  pointer                 :: domain
      real(real64), contiguous, pointer, dimension(:,:) :: u, v, var

      call c_f_pointer(padvection, advection)
      call c_f_pointer(pdomain, domain)
      call c_f_pointer(pu, u, domain%T%u(1:2) - domain%T%l(1:2) + 1)
      call c_f_pointer(pv, v, domain%T%u(1:2) - domain%T%l(1:2) + 1)
      call c_f_pointer(pvar, var, domain%T%u(1:2) - domain%T%l(1:2) + 1)
     call advection%advection_calculate_2d(scheme, domain%U, u, domain%V, v, timestep, domain%T,var)
   end subroutine

end module