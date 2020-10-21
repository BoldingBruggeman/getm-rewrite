module pylo

   use iso_c_binding, only: c_ptr, c_int, c_double, c_char, c_loc, c_f_pointer
   use iso_fortran_env, only: real64

   use getm_domain, only: type_getm_domain, type_getm_grid
   use getm_operators, only: type_advection
   use memory_manager

   implicit none

contains

   function domain_create(imin, imax, jmin, jmax, kmin, kmax, halo) result(pdomain) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: domain_create
      integer(c_int), intent(in), value :: imin, imax, jmin, jmax, kmin, kmax, halo
      type(c_ptr)                       :: pdomain

      type (type_getm_domain), pointer :: domain
      integer                          :: stat 

      allocate(domain)
      call domain%T%configure(imin=imin,imax=imax,jmin=jmin,jmax=jmax,kmin=kmin,kmax=kmax,halo=(/halo,halo,0/))
      call mm_s('c1', domain%T%c1, domain%T%l(1), domain%T%u(1), def=0._real64,stat=stat)
      call mm_s('c2', domain%T%c2, domain%T%l(2), domain%T%u(2), def=0._real64,stat=stat)
      call mm_s('H', domain%T%H, domain%T%l(1:2), domain%T%u(1:2), def=-99._real64, stat=stat)
      call mm_s('mask', domain%T%mask, domain%T%l(1:2), domain%T%u(1:2), def=0, stat=stat)
      pdomain = c_loc(domain)
      domain%domain_type = 1
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

   subroutine grid_get_arrays(pgrid, pc1, pc2, pH, pmask) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: grid_get_arrays
      type(c_ptr), intent(in), value :: pgrid
      type(c_ptr), intent(out)       :: pc1, pc2, pH, pmask

      type (type_getm_grid), pointer :: grid

      call c_f_pointer(pgrid, grid)

      pc1 = c_loc(grid%c1)
      pc2 = c_loc(grid%c2)
      pH = c_loc(grid%H)
      pmask = c_loc(grid%mask)
   end subroutine

   subroutine domain_initialize(pdomain) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: domain_initialize
      type(c_ptr), intent(in), value :: pdomain

      type (type_getm_domain), pointer :: domain

      call c_f_pointer(pdomain, domain)
      call domain%configure()
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