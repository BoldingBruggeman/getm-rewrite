module pygetm

   use iso_c_binding, only: c_ptr, c_int, c_double, c_char, c_loc, c_f_pointer, C_NULL_CHAR, C_NULL_PTR
   use iso_fortran_env, only: real64

   use getm_domain, only: type_getm_domain, type_getm_grid
   use getm_operators, only: type_advection
   use getm_sealevel, only: type_getm_sealevel
   use getm_pressure, only: type_getm_pressure
   use getm_momentum, only: type_getm_momentum
   use memory_manager

   implicit none

contains

   function domain_create(imin, imax, jmin, jmax, kmin, kmax, halox, haloy, haloz) result(pdomain) bind(c)
      integer(c_int), intent(in), value :: imin, imax, jmin, jmax, kmin, kmax
      integer(c_int), intent(out)       :: halox, haloy, haloz
      type(c_ptr)                       :: pdomain

      type (type_getm_domain), pointer :: domain
      integer                          :: stat 

      allocate(domain)
      domain%have_metrics = .true.
      call domain%configure(imin=imin,imax=imax,jmin=jmin,jmax=jmax,kmin=kmin,kmax=kmax)
      pdomain = c_loc(domain)
      domain%domain_type = 1
      halox = imin - lbound(domain%T%c1, 1)
      haloy = jmin - lbound(domain%T%c2, 1)
      haloz = 0
   end function

   subroutine domain_finalize(pdomain) bind(c)
      type(c_ptr), intent(in), value :: pdomain

      type (type_getm_domain), pointer :: domain

      call c_f_pointer(pdomain, domain)
      deallocate(domain)
   end subroutine

   subroutine domain_initialize_open_boundaries(pdomain, nwb, nnb, neb, nsb, nbdyp) bind(c)
      type(c_ptr),    intent(in), value :: pdomain
      integer(c_int), intent(in), value :: nwb, nnb, neb, nsb, nbdyp

      type (type_getm_domain), pointer :: domain

      call c_f_pointer(pdomain, domain)
      call domain%initialize_open_boundaries(nwb=nwb, nnb=nnb, neb=neb, nsb=nsb, nbdyp=nbdyp)
   end subroutine

   function domain_get_grid(pdomain, grid_type) result(pgrid) bind(c)
      type(c_ptr),         intent(in), value :: pdomain
      integer(kind=c_int), intent(in), value :: grid_type
      type(c_ptr)                            :: pgrid

      type (type_getm_domain),   pointer :: domain
      type (type_getm_grid),     pointer :: grid
      type (type_getm_pressure), pointer :: pressure

      call c_f_pointer(pdomain, domain)
      select case (grid_type)
      case (1); grid => domain%T
      case (2); grid => domain%U
      case (3); grid => domain%V
      case (4); grid => domain%X
      case default; grid => null()
      end select
      pgrid = c_loc(grid)
   end function

   subroutine get_array(source_type, obj, name, grid_type, data_type, p) bind(c)
      integer(c_int),  value,         intent(in)  :: source_type
      type(c_ptr), value,             intent(in)  :: obj
      character(kind=c_char), target, intent(in)  :: name(*)
      integer(c_int),                 intent(out) :: grid_type, data_type
      type(c_ptr),                    intent(out) :: p

      type (type_getm_grid),     pointer :: grid
      type (type_getm_momentum), pointer :: momentum
      type (type_getm_pressure), pointer :: pressure
      character(len=10),         pointer :: pname

      call c_f_pointer(c_loc(name), pname)

      p = C_NULL_PTR
      grid_type = 1   ! TGRID (1: TGRID, 2: UGRID, 3: VGRID, 4: XGRID)
      data_type = 0   ! double (use 1 for integer)
      select case (source_type)
      case (0)
         call c_f_pointer(obj, grid)
         grid_type = grid%grid_type
         select case (pname(:index(pname, C_NULL_CHAR) - 1))
         case ('c1'); p = c_loc(grid%c1)
         case ('c2'); p = c_loc(grid%c2)
         case ('x'); p = c_loc(grid%x)
         case ('y'); p = c_loc(grid%y)
         case ('dx'); p = c_loc(grid%dx)
         case ('dy'); p = c_loc(grid%dy)
         case ('lon'); p = c_loc(grid%lon)
         case ('lat'); p = c_loc(grid%lat)
         case ('dlon'); p = c_loc(grid%dlon)
         case ('dlat'); p = c_loc(grid%dlat)
         case ('area'); p = c_loc(grid%area)
         case ('iarea'); p = c_loc(grid%iarea)
         case ('H'); p = c_loc(grid%H)
         case ('D'); p = c_loc(grid%D)
         case ('mask'); p = c_loc(grid%mask); data_type = 1
         case ('z'); p = c_loc(grid%z)
         case ('zo'); p = c_loc(grid%zo)
         case ('cor'); p = c_loc(grid%cor)
         end select
      case (1)
         call c_f_pointer(obj, momentum)
         select case (pname(:index(pname, C_NULL_CHAR) - 1))
         case ('U');   p = c_loc(momentum%U); grid_type = 2
         case ('V');   p = c_loc(momentum%V); grid_type = 3
         case ('fU');  p = c_loc(momentum%fU); grid_type = 2
         case ('fV');  p = c_loc(momentum%fV); grid_type = 3
         case ('advU');  p = c_loc(momentum%advU); grid_type = 2
         case ('advV');  p = c_loc(momentum%advV); grid_type = 3
         end select
      case (2)
         call c_f_pointer(obj, pressure)
         select case (pname(:index(pname, C_NULL_CHAR) - 1))
         case ('dpdx'); p = c_loc(pressure%dpdx); grid_type = 2
         case ('dpdy'); p = c_loc(pressure%dpdy); grid_type = 3
         case default; p = C_NULL_PTR
         end select
      end select
   end subroutine

   subroutine domain_initialize(pdomain, runtype, maxdt) bind(c)
      type(c_ptr),    intent(in), value :: pdomain
      integer(c_int), intent(in), value :: runtype
      real(c_double), intent(out)       :: maxdt

      type (type_getm_domain), pointer :: domain

      call c_f_pointer(pdomain, domain)
      call domain%initialize(runtype)
      maxdt = domain%maxdt
   end subroutine

   subroutine domain_update_depths(pdomain) bind(c)
      type(c_ptr), intent(in), value :: pdomain

      type (type_getm_domain), pointer :: domain

      call c_f_pointer(pdomain, domain)
      call domain%update_depths()
   end subroutine

   function advection_create(scheme, ptgrid) result(padvection) bind(c)
      integer(c_int), intent(in), value :: scheme
      type(c_ptr),    intent(in), value :: ptgrid
      type(c_ptr) :: padvection

      type (type_getm_grid), pointer :: tgrid
      type (type_advection), pointer :: advection

      call c_f_pointer(ptgrid, tgrid)
      allocate(advection)
      call advection%initialize(scheme, tgrid)
      padvection = c_loc(advection)
   end function

   subroutine advection_2d_start(padvection, ptgrid) bind(c)
      type(c_ptr),    intent(in), value :: padvection, ptgrid

      type (type_advection), pointer :: advection
      type (type_getm_grid), pointer :: tgrid

      call c_f_pointer(padvection, advection)
      call c_f_pointer(ptgrid, tgrid)
      advection%D(:,:) = tgrid%D
   end subroutine

   subroutine advection_2d_calculate(direction, padvection, ptgrid, pugrid,  pu, timestep, pvar) bind(c)
      integer(c_int), intent(in), value :: direction
      real(c_double), intent(in), value :: timestep
      type(c_ptr),    intent(in), value :: padvection, ptgrid, pugrid, pu, pvar

      type (type_advection),    pointer                 :: advection
      type (type_getm_grid),  pointer                   :: tgrid, ugrid
      real(real64), contiguous, pointer, dimension(:,:) :: u, var

      call c_f_pointer(padvection, advection)
      call c_f_pointer(ptgrid, tgrid)
      call c_f_pointer(pugrid, ugrid)
      call c_f_pointer(pu, u, tgrid%u(1:2) - tgrid%l(1:2) + 1)
      call c_f_pointer(pvar, var, tgrid%u(1:2) - tgrid%l(1:2) + 1)
      select case (direction)
         case (1)
            call advection%advection_calculate_u2d(ugrid, u, timestep, tgrid, var)
         case (2)
            call advection%advection_calculate_v2d(ugrid, u, timestep, tgrid, var)
      end select
   end subroutine

   function momentum_create(runtype, pdomain, padvection, advection_scheme, apply_bottom_friction) result(pmomentum) bind(c)
      integer(c_int), intent(in), value :: runtype
      type(c_ptr),    intent(in), value :: pdomain
      type(c_ptr),    intent(in), value :: padvection
      integer(c_int), intent(in), value :: advection_scheme
      integer(c_int), intent(in), value :: apply_bottom_friction
      type(c_ptr) :: pmomentum

      type (type_getm_domain),   pointer :: domain
      type (type_advection),     pointer :: advection
      type (type_getm_momentum), pointer :: momentum

      call c_f_pointer(pdomain, domain)
      call c_f_pointer(padvection, advection)
      allocate(momentum)
      call momentum%configure()
      momentum%advection_scheme = advection_scheme
      momentum%apply_bottom_friction = (apply_bottom_friction == 1)
      call momentum%initialize(runtype, domain, advection)
      pmomentum = c_loc(momentum)
   end function

   subroutine momentum_uv_momentum_2d(pmomentum, runtype, timestep, ptausx, ptausy, pdpdx, pdpdy) bind(c)
      type(c_ptr),    intent(in), value :: pmomentum
      integer(c_int), intent(in), value :: runtype
      real(c_double), intent(in), value :: timestep
      type(c_ptr),    intent(in), value :: ptausx, ptausy, pdpdx, pdpdy

      type (type_getm_momentum), pointer :: momentum
      real(real64), contiguous, pointer, dimension(:,:) :: tausx, tausy, dpdx, dpdy

      call c_f_pointer(pmomentum, momentum)
      call c_f_pointer(ptausx, tausx, momentum%domain%U%u(1:2) - momentum%domain%U%l(1:2) + 1)
      call c_f_pointer(ptausy, tausy, momentum%domain%V%u(1:2) - momentum%domain%V%l(1:2) + 1)
      call c_f_pointer(pdpdx, dpdx, momentum%domain%T%u(1:2) - momentum%domain%T%l(1:2) + 1)
      call c_f_pointer(pdpdy, dpdy, momentum%domain%T%u(1:2) - momentum%domain%T%l(1:2) + 1)
      call momentum%uv_momentum_2d(runtype, timestep, tausx, tausy, dpdx, dpdy)
   end subroutine

   function pressure_create(runtype, pdomain) result(ppressure) bind(c)
      integer(c_int), intent(in), value :: runtype
      type(c_ptr),    intent(in), value :: pdomain
      type(c_ptr) :: ppressure

      type (type_getm_domain),   pointer :: domain
      type (type_getm_pressure), pointer :: pressure

      call c_f_pointer(pdomain, domain)
      allocate(pressure)
      call pressure%configure()
      call pressure%initialize(runtype, domain)
      ppressure = c_loc(pressure)
   end function

   subroutine pressure_surface(ppressure, pz, psp) bind(c)
      type(c_ptr), intent(in), value :: ppressure
      type(c_ptr), intent(in), value :: pz, psp

      type (type_getm_pressure), pointer :: pressure
      real(real64), contiguous, pointer, dimension(:,:) :: z, sp

      call c_f_pointer(ppressure, pressure)
      call c_f_pointer(pz, z, pressure%domain%T%u(1:2) - pressure%domain%T%l(1:2) + 1)
      call c_f_pointer(psp, sp, pressure%domain%T%u(1:2) - pressure%domain%T%l(1:2) + 1)
      call pressure%surface(z, sp)
   end subroutine

   function sealevel_create(pdomain) result(psealevel) bind(c)
      type(c_ptr), intent(in), value :: pdomain
      type(c_ptr) :: psealevel

      type (type_getm_domain),   pointer :: domain
      type (type_getm_sealevel), pointer :: sealevel

      call c_f_pointer(pdomain, domain)
      allocate(sealevel)
      call sealevel%configure()
      call sealevel%initialize(domain)
      psealevel = c_loc(sealevel)
   end function

   subroutine sealevel_update(psealevel, timestep, pU, pV) bind(c)
      type(c_ptr), intent(in),    value :: psealevel
      real(c_double), intent(in), value :: timestep
      type(c_ptr), intent(in),    value :: pU, pV

      type (type_getm_sealevel), pointer :: sealevel
      real(real64), contiguous, pointer, dimension(:,:) :: U, V

      call c_f_pointer(psealevel, sealevel)
      call c_f_pointer(pU, U, sealevel%domain%T%u(1:2) - sealevel%domain%T%l(1:2) + 1)
      call c_f_pointer(pV, V, sealevel%domain%T%u(1:2) - sealevel%domain%T%l(1:2) + 1)
      call sealevel%t(timestep, U, V)
      call sealevel%uvx()
   end subroutine

   subroutine sealevel_update_uvx(psealevel) bind(c)
      type(c_ptr), intent(in),    value :: psealevel

      type (type_getm_sealevel), pointer :: sealevel

      call c_f_pointer(psealevel, sealevel)
      call sealevel%uvx()
   end subroutine

end module
