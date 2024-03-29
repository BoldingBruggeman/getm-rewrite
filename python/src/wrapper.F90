module pygetm

   use iso_c_binding, only: c_ptr, c_int, c_double, c_char, c_loc, c_f_pointer, c_associated, C_NULL_CHAR, C_NULL_PTR
   use iso_fortran_env, only: real64

   use getm_domain, only: type_getm_domain, type_getm_grid, g
   use getm_operators, only: type_advection, type_vertical_diffusion
   use getm_sealevel, only: type_getm_sealevel
   use getm_pressure, only: type_getm_pressure
   use getm_momentum, only: type_getm_momentum, kappa, rho0
   use memory_manager

   implicit none

   type array
      real(c_double), allocatable :: rdata(:)
      integer(c_int), allocatable :: idata(:)
   end type

contains

   subroutine c_allocate_array(n, data_type, ptype, pdata) bind(c)
      integer(c_int), intent(in), value :: n, data_type
      type(c_ptr),    intent(out)       :: ptype, pdata

      type(array), pointer :: p

      allocate(p)
      ptype = c_loc(p)
      if (data_type == 0) then
         allocate(p%rdata(n))
         pdata = c_loc(p%rdata)
      else
         allocate(p%idata(n))
         pdata = c_loc(p%idata)
      end if
   end subroutine

   subroutine c_deallocate_array(ptype) bind(c)
      type(c_ptr), intent(in), value :: ptype

      type(array), pointer :: p

      call c_f_pointer(ptype, p)
      deallocate(p)
   end subroutine

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
      call domain%cleanup()
      deallocate(domain)
   end subroutine

   function domain_get_grid(pdomain, grid_type, imin, imax, jmin, jmax, kmin, kmax, halox, haloy, haloz) result(pgrid) bind(c)
      type(c_ptr),         intent(in), value :: pdomain
      integer(kind=c_int), intent(in), value :: grid_type
      integer(kind=c_int), intent(in), value :: imin, imax, jmin, jmax, kmin, kmax, halox, haloy, haloz
      type(c_ptr)                            :: pgrid

      type (type_getm_domain),   pointer :: domain
      type (type_getm_grid),     pointer :: grid
      type (type_getm_pressure), pointer :: pressure
      integer                            :: halo(3)

      call c_f_pointer(pdomain, domain)
      select case (grid_type)
      case (1); grid => domain%T
      case (2); grid => domain%U
      case (3); grid => domain%V
      case (4); grid => domain%X
      case default
         allocate(grid)
         halo = (/halox, haloy, haloz/)
         call grid%configure(grid_type=grid_type, imin=imin, imax=imax, &
            jmin=jmin, jmax=jmax, kmin=kmin, kmax=kmax, halo=halo)
      end select
      pgrid = c_loc(grid)
   end function

   subroutine grid_finalize(pgrid) bind(c)
      type(c_ptr), intent(in), value :: pgrid

      type (type_getm_grid), pointer :: grid

      call c_f_pointer(pgrid, grid)
      deallocate(grid)
   end subroutine

   subroutine domain_do_vertical(pdomain, timestep) bind(c)
      type(c_ptr),    intent(in), value :: pdomain
      real(c_double), intent(in), value :: timestep

      type (type_getm_domain), pointer :: domain

      call c_f_pointer(pdomain, domain)
      call domain%do_vertical(timestep)
   end subroutine

   subroutine get_array(source_type, obj, name, grid_type, sub_type, data_type, p) bind(c)
      integer(c_int),  value,         intent(in)  :: source_type
      type(c_ptr), value,             intent(in)  :: obj
      character(kind=c_char), target, intent(in)  :: name(*)
      integer(c_int),                 intent(out) :: grid_type, sub_type, data_type
      type(c_ptr),                    intent(out) :: p

      type (type_getm_grid),     pointer :: grid
      type (type_getm_momentum), pointer :: momentum
      character(len=10),         pointer :: pname

      integer, parameter :: subtype_depth_explicit = 1
      integer, parameter :: subtype_depth_explicit_interfaces = 2

      call c_f_pointer(c_loc(name), pname)

      p = C_NULL_PTR
      grid_type = 1   ! TGRID (1: TGRID, 2: UGRID, 3: VGRID, 4: XGRID)
      sub_type = 0    ! on-grid: 0, on boundary points: 1
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
         case ('idx'); p = c_loc(grid%idx)
         case ('idy'); p = c_loc(grid%idy)
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
         case ('zio'); p = c_loc(grid%zio)
         case ('zin'); p = c_loc(grid%zin)
         case ('cor'); p = c_loc(grid%cor)
         case ('z0b'); p = c_loc(grid%z0b)
         case ('z0b_min'); p = c_loc(grid%z0b_min)
         case ('hn'); p = c_loc(grid%hn); sub_type = subtype_depth_explicit
         case ('ho'); p = c_loc(grid%ho); sub_type = subtype_depth_explicit
         case ('zc'); p = c_loc(grid%zc); sub_type = subtype_depth_explicit
         case ('zf'); p = c_loc(grid%zf); sub_type = subtype_depth_explicit_interfaces
         case ('alpha'); p = c_loc(grid%alpha)
         end select
      case (1)
         call c_f_pointer(obj, momentum)
         select case (pname(:index(pname, C_NULL_CHAR) - 1))
         case ('U');   p = c_loc(momentum%U); grid_type = 2
         case ('V');   p = c_loc(momentum%V); grid_type = 3
         case ('fU');  p = c_loc(momentum%fU); grid_type = 3
         case ('fV');  p = c_loc(momentum%fV); grid_type = 2
         case ('fpk');  if (allocated(momentum%fpk)) p = c_loc(momentum%fpk); grid_type = 3; sub_type = subtype_depth_explicit
         case ('fqk');  if (allocated(momentum%fqk)) p = c_loc(momentum%fqk); grid_type = 2; sub_type = subtype_depth_explicit
         case ('advU');  p = c_loc(momentum%advU); grid_type = 2
         case ('advV');  p = c_loc(momentum%advV); grid_type = 3
         case ('diffU');  p = c_loc(momentum%diffU); grid_type = 2
         case ('diffV');  p = c_loc(momentum%diffV); grid_type = 3
         case ('dampU');  p = c_loc(momentum%dampU); grid_type = 2
         case ('dampV');  p = c_loc(momentum%dampV); grid_type = 3
         case ('u1');   p = c_loc(momentum%u1); grid_type = 2
         case ('v1');   p = c_loc(momentum%v1); grid_type = 3
         case ('uk');   if (allocated(momentum%uk)) p = c_loc(momentum%uk); grid_type = 2; sub_type = subtype_depth_explicit
         case ('vk');   if (allocated(momentum%vk)) p = c_loc(momentum%vk); grid_type = 3; sub_type = subtype_depth_explicit
         case ('ru');   p = c_loc(momentum%ru); grid_type = 2
         case ('rv');   p = c_loc(momentum%rv); grid_type = 3
         case ('rru');   p = c_loc(momentum%rru); grid_type = 2
         case ('rrv');   p = c_loc(momentum%rrv); grid_type = 3
         case ('ustar2_s');   if (allocated(momentum%ustar2_s)) p = c_loc(momentum%ustar2_s)
         case ('ustar2_b');   if (allocated(momentum%ustar2_b))  p = c_loc(momentum%ustar2_b)
         case ('pk');   if (allocated(momentum%pk)) p = c_loc(momentum%pk); grid_type = 2; sub_type = subtype_depth_explicit
         case ('qk');   if (allocated(momentum%qk)) p = c_loc(momentum%qk); grid_type = 3; sub_type = subtype_depth_explicit
         case ('ww');   if (allocated(momentum%ww)) p = c_loc(momentum%ww); sub_type = subtype_depth_explicit_interfaces
         case ('advpk'); if (allocated(momentum%advpk)) p = c_loc(momentum%advpk); grid_type = 2; sub_type = subtype_depth_explicit
         case ('advqk'); if (allocated(momentum%advqk)) p = c_loc(momentum%advqk); grid_type = 3; sub_type = subtype_depth_explicit
         case ('diffpk'); if (allocated(momentum%diffpk)) p = c_loc(momentum%diffpk); grid_type = 2 &
            ;sub_type = subtype_depth_explicit
         case ('diffqk'); if (allocated(momentum%diffqk)) p = c_loc(momentum%diffqk); grid_type = 3 &
            ;sub_type = subtype_depth_explicit
         case ('Ui');   p = c_loc(momentum%Ui); grid_type = 2
         case ('Vi');   p = c_loc(momentum%Vi); grid_type = 3
         case ('SS');   if (allocated(momentum%SS)) p = c_loc(momentum%SS); sub_type = subtype_depth_explicit_interfaces
         case ('SxB');   p = c_loc(momentum%SxB); grid_type = 2
         case ('SyB');   p = c_loc(momentum%SyB); grid_type = 3
         case ('SxA');   p = c_loc(momentum%SxA); grid_type = 2
         case ('SyA');   p = c_loc(momentum%SyA); grid_type = 3
         case ('SxD');   p = c_loc(momentum%SxD); grid_type = 2
         case ('SyD');   p = c_loc(momentum%SyD); grid_type = 3
         case ('SxF');   p = c_loc(momentum%SxF); grid_type = 2
         case ('SyF');   p = c_loc(momentum%SyF); grid_type = 3
         end select
      end select
   end subroutine

   subroutine domain_initialize(pdomain, runtype, Dmin, method_vertical_coordinates, ddl, ddu, Dgamma, gamma_surf, maxdt) bind(c)
      type(c_ptr),    intent(in), value :: pdomain
      integer(c_int), intent(in), value :: runtype
      real(c_double), intent(in), value :: Dmin
      integer(c_int), intent(in), value :: method_vertical_coordinates
      real(c_double), intent(in), value :: ddl, ddu, Dgamma
      integer(c_int), intent(in), value :: gamma_surf
      real(c_double), intent(out)       :: maxdt

      type (type_getm_domain), pointer :: domain

      call c_f_pointer(pdomain, domain)
      domain%Dmin = Dmin
      domain%ddl = ddl
      domain%ddu = ddu
      domain%Dgamma = Dgamma
      domain%gamma_surf = gamma_surf /= 0
      domain%method_vertical_coordinates = method_vertical_coordinates
      call domain%initialize(runtype)
      maxdt = domain%maxdt
   end subroutine

   function vertical_diffusion_create(ptgrid) result(pdiffusion) bind(c)
      type(c_ptr),    intent(in), value :: ptgrid
      type(c_ptr) :: pdiffusion

      type (type_getm_grid), pointer :: tgrid
      type (type_vertical_diffusion), pointer :: diffusion

      call c_f_pointer(ptgrid, tgrid)
      allocate(diffusion)
      call diffusion%initialize(tgrid)
      pdiffusion = c_loc(diffusion)
   end function

   subroutine vertical_diffusion_finalize(pdiffusion) bind(c)
      type(c_ptr), intent(in), value :: pdiffusion

      type (type_vertical_diffusion), pointer :: diffusion

      call c_f_pointer(pdiffusion, diffusion)
      deallocate(diffusion)
   end subroutine

   subroutine c_vertical_diffusion_prepare(pdiffusion, nx, ny, nz, molecular, nuh, timestep, cnpar, mask, ho, hn) bind(c)
      type(c_ptr),    intent(in), value :: pdiffusion
      integer(c_int), intent(in), value :: nx, ny, nz
      real(c_double), intent(in), value :: molecular, timestep, cnpar
      integer(c_int), intent(in) :: mask(nx, ny)
      real(c_double), intent(in) :: nuh(nx, ny, 0:nz), ho(nx, ny, nz), hn(nx, ny, nz)

      type (type_vertical_diffusion), pointer :: diffusion
      real(real64), contiguous, pointer, dimension(:,:,:) :: ea2, ea4

      call c_f_pointer(pdiffusion, diffusion)
      call diffusion%prepare(timestep, cnpar, mask, ho, hn, molecular, nuh(:, :, 1:nz-1))
   end subroutine

   subroutine c_vertical_diffusion_apply(pdiffusion, nx, ny, nz, mask, ho, hn, var, pea2, pea4) bind(c)
      type(c_ptr),    intent(in), value :: pdiffusion
      integer(c_int), intent(in), value :: nx, ny, nz
      integer(c_int), intent(in) :: mask(nx, ny)
      real(c_double), intent(in) :: ho(nx, ny, nz), hn(nx, ny, nz)
      real(c_double), intent(inout) :: var(nx, ny, nz)
      type(c_ptr),    intent(in), value :: pea2, pea4

      type (type_vertical_diffusion), pointer :: diffusion
      real(real64), contiguous, pointer, dimension(:,:,:) :: ea2, ea4

      call c_f_pointer(pdiffusion, diffusion)
      if (c_associated(pea2)) call c_f_pointer(pea2, ea2, (/nx, ny, nz/))
      if (c_associated(pea4)) call c_f_pointer(pea4, ea4, (/nx, ny, nz/))
      if (c_associated(pea2) .and. c_associated(pea4)) then
         call diffusion%apply(mask, ho, hn, var, ea2=ea2, ea4=ea4)
      elseif (c_associated(pea2)) then
         call diffusion%apply(mask, ho, hn, var, ea2=ea2)
      elseif (c_associated(pea4)) then
         call diffusion%apply(mask, ho, hn, var, ea4=ea4)
      else
         call diffusion%apply(mask, ho, hn, var)
      end if
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
      !pD = c_loc(advection%D)
      !phn = c_loc(advection%hn)
   end function

   subroutine advection_finalize(padvection) bind(c)
      type(c_ptr), intent(in), value :: padvection

      type (type_advection), pointer :: advection

      call c_f_pointer(padvection, advection)
      deallocate(advection)
   end subroutine

   subroutine advection_uv_calculate(direction, nk, padvection, ptgrid, pugrid, pu, pAh, timestep, ph, phu, pvar) bind(c)
      integer(c_int), intent(in), value :: direction, nk
      type(c_ptr),    intent(in), value :: pAh
      real(c_double), intent(in), value :: timestep
      type(c_ptr),    intent(in), value :: padvection, ptgrid, pugrid, pu, ph, phu, pvar

      type (type_advection),    pointer                   :: advection
      type (type_getm_grid),  pointer                     :: tgrid, ugrid
      real(real64), contiguous, pointer, dimension(:,:)   :: Ah
      real(real64), contiguous, pointer, dimension(:,:,:) :: u, h, hu, var
      integer                                             :: k
      logical                                             :: apply_diffusion

      call c_f_pointer(padvection, advection)
      if (.not. allocated(advection%op)) return
      call c_f_pointer(ptgrid, tgrid)
      call c_f_pointer(pugrid, ugrid)
      apply_diffusion = c_associated(pAh)
      if (apply_diffusion) then
         call c_f_pointer(pAh, Ah, (/ugrid%u(1) - ugrid%l(1) + 1, ugrid%u(2) - ugrid%l(2) + 1/))
      else
         call c_f_pointer(pu, Ah, (/ugrid%u(1) - ugrid%l(1) + 1, ugrid%u(2) - ugrid%l(2) + 1/))
      end if
      call c_f_pointer(pu, u, (/ugrid%u(1) - ugrid%l(1) + 1, ugrid%u(2) - ugrid%l(2) + 1, nk/))
      call c_f_pointer(ph, h, (/tgrid%u(1) - tgrid%l(1) + 1, tgrid%u(2) - tgrid%l(2) + 1, nk/))
      call c_f_pointer(phu, hu, (/ugrid%u(1) - ugrid%l(1) + 1, ugrid%u(2) - ugrid%l(2) + 1, nk/))
      call c_f_pointer(pvar, var, (/tgrid%u(1) - tgrid%l(1) + 1, tgrid%u(2) - tgrid%l(2) + 1, nk/))
      select case (direction)
         case (1)
            do k = 1, nk
               call advection%op%u2d(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax,tgrid%halo, &
                       ugrid%mask,ugrid%idx,ugrid%dy,hu(:,:,k),u(:,:,k), &
                       tgrid%mask,tgrid%iarea,apply_diffusion,Ah,timestep,h(:,:,k),var(:,:,k))
            end do
         case (2)
            do k = 1, nk
               call advection%op%v2d(tgrid%imin,tgrid%imax,tgrid%jmin,tgrid%jmax,tgrid%halo, &
                       ugrid%mask,ugrid%dx,ugrid%idy,hu(:,:,k),u(:,:,k), &
                       tgrid%mask,tgrid%iarea,apply_diffusion,Ah,timestep,h(:,:,k),var(:,:,k))
            end do
      end select
   end subroutine

   subroutine advection_w_calculate(padvection, ptgrid, pw, pw_var, timestep, ph, pvar) bind(c)
      real(c_double), intent(in), value :: timestep
      type(c_ptr),    intent(in), value :: padvection, ptgrid, pw, pw_var, ph, pvar

      type (type_advection),    pointer                   :: advection
      type (type_getm_grid),    pointer                   :: tgrid
      real(real64), contiguous, pointer, dimension(:,:,:) :: w, w_var, h, var

      call c_f_pointer(padvection, advection)
      if (.not. allocated(advection%op)) return
      call c_f_pointer(ptgrid, tgrid)
      call c_f_pointer(pw, w, (/tgrid%u(1) - tgrid%l(1) + 1, tgrid%u(2) - tgrid%l(2) + 1, tgrid%kmax + 1/))
      call c_f_pointer(pw_var, w_var, (/tgrid%u(1) - tgrid%l(1) + 1, tgrid%u(2) - tgrid%l(2) + 1, tgrid%kmax + 1/))
      call c_f_pointer(ph, h, tgrid%u - tgrid%l + 1)
      call c_f_pointer(pvar, var, tgrid%u - tgrid%l + 1)
      call advection%op%w3d(tgrid%imin, tgrid%imax, tgrid%jmin, tgrid%jmax, tgrid%kmax, tgrid%halo, &
                    w, w_var, tgrid%mask, timestep, h, var)
   end subroutine

   function momentum_create(runtype, pdomain, Am0, cnpar, coriolis_scheme) result(pmomentum) bind(c)
      integer(c_int), intent(in), value :: runtype, coriolis_scheme
      type(c_ptr),    intent(in), value :: pdomain
      real(c_double), intent(in), value :: Am0, cnpar
      type(c_ptr) :: pmomentum

      type (type_getm_domain),   pointer :: domain
      type (type_getm_momentum), pointer :: momentum

      call c_f_pointer(pdomain, domain)
      allocate(momentum)
      call momentum%configure()
      momentum%advection_scheme = 0
      momentum%Am0 = Am0
      momentum%cnpar = cnpar
      momentum%coriolis_scheme = coriolis_scheme
      call momentum%initialize(runtype, domain)
      pmomentum = c_loc(momentum)
   end function

   subroutine momentum_finalize(pmomentum) bind(c)
      type(c_ptr), intent(in), value :: pmomentum

      type (type_getm_momentum), pointer :: momentum

      call c_f_pointer(pmomentum, momentum)
      deallocate(momentum)
   end subroutine

   subroutine momentum_diffusion_driver(pmomentum, nk, ph, phx, pu, pv, pdiffu, pdiffv) bind(c)
      integer(c_int), intent(in), value :: nk
      type(c_ptr),    intent(in), value :: pmomentum, ph, phx, pu, pv, pdiffu, pdiffv

      type (type_getm_momentum), pointer :: momentum
      real(real64), contiguous, pointer, dimension(:,:,:) :: h, hx, u, v, diffu, diffv
      integer :: k

      call c_f_pointer(pmomentum, momentum)
      associate(TG => momentum%domain%T, UG => momentum%domain%U, VG => momentum%domain%V, XG => momentum%domain%X)
      call c_f_pointer(ph, h,   (/TG%u(1) - TG%l(1) + 1, TG%u(2) - TG%l(2) + 1, nk/))
      call c_f_pointer(phx, hx, (/XG%u(1) - XG%l(1) + 1, XG%u(2) - XG%l(2) + 1, nk/))
      call c_f_pointer(pu, u,   (/UG%u(1) - UG%l(1) + 1, UG%u(2) - UG%l(2) + 1, nk/))
      call c_f_pointer(pv, v,   (/VG%u(1) - VG%l(1) + 1, VG%u(2) - VG%l(2) + 1, nk/))
      call c_f_pointer(pdiffu, diffu, (/UG%u(1) - UG%l(1) + 1, UG%u(2) - UG%l(2) + 1, nk/))
      call c_f_pointer(pdiffv, diffv, (/VG%u(1) - VG%l(1) + 1, VG%u(2) - VG%l(2) + 1, nk/))
      end associate
      do k = 1, nk
         call momentum%diffusion_driver(h(:,:,k), hx(:,:,k), u(:,:,k), v(:,:,k), diffu(:,:,k), diffv(:,:,k))
      end do
   end subroutine

   subroutine c_thickness2center_depth(nx, ny, nz, istart, istop, jstart, jstop, mask, h, out) bind(c)
      integer(c_int), intent(in), value :: nx, ny, nz, istart, istop, jstart, jstop
      integer(c_int), intent(in)        :: mask(nx, ny)
      real(c_double), intent(in)        :: h(nx, ny, nz)
      real(c_double), intent(inout)     :: out(nx, ny, nz)

      integer :: k

      where (mask /= 0) out(:,:,nz) = 0.5_c_double * h(:,:,nz)
      do k=nz-1,1,-1
         where (mask /= 0) out(:,:,k) = out(:,:,k+1) + 0.5_c_double * (h(:,:,k) + h(:,:,k+1))
      end do
   end subroutine

   subroutine c_thickness2vertical_coordinates(nx, ny, nz, mask, bottom_depth, h, zc, zf) bind(c)
      integer(c_int), intent(in), value :: nx, ny, nz
      integer(c_int), intent(in)        :: mask(nx, ny)
      real(c_double), intent(in)        :: bottom_depth(nx, ny), h(nx, ny, nz)
      real(c_double), intent(inout)     :: zc(nx, ny, nz), zf(nx, ny, 0:nz)

      integer :: k

      where (mask /= 0) zc(:,:,1) = -bottom_depth(:,:) + 0.5_c_double * h(:,:,1)
      do k=2,nz
         where (mask /= 0) zc(:,:,k) = zc(:,:,k-1) + 0.5_c_double * (h(:,:,k-1) + h(:,:,k))
      end do

      where (mask /= 0) zf(:,:,0) = -bottom_depth(:,:)
      do k=1,nz
         where (mask /= 0) zf(:,:,k) = zf(:,:,k-1) + h(:,:,k)
      end do
   end subroutine

   subroutine c_alpha(n, D, Dmin, Dcrit, mask, alpha) bind(c)
      integer(c_int), intent(in), value :: n
      integer(c_int), intent(in)        :: mask(n)
      real(c_double), intent(in)        :: D(n)
      real(c_double), intent(in), value :: Dmin, Dcrit
      real(c_double), intent(inout)     :: alpha(n)
      where (mask /= 0) alpha = max(0._c_double, min(1._c_double, (D - Dmin) / (Dcrit - Dmin)))
   end subroutine

   subroutine c_clip_z(n, z, H, Dmin, mask) bind(c)
      integer(c_int), intent(in), value :: n
      real(c_double), intent(inout)     :: z(n)
      real(c_double), intent(in)        :: H(n)
      real(c_double), intent(in), value :: Dmin
      integer(c_int), intent(in)        :: mask(n)
      where (mask /= 0) z = max(-H + Dmin, z)
   end subroutine

   subroutine c_vertical_advection_to_sources(nx, ny, nz, halo, mask, c, w, h, s) bind(c)
      ! first-order upstream-biased advection, e.g., for integrating FABM sinking/floating into source term
      integer(c_int), intent(in), value :: nx, ny, nz, halo
      integer(c_int), intent(in)        :: mask(nx, ny, nz)
      real(c_double), intent(in)        :: c(nx, ny, nz), w(nx, ny, nz), h(nx, ny, nz)
      real(c_double), intent(inout)     :: s(nx, ny, nz)

      logical :: active
      integer :: i, j, k
      real(c_double) :: local_w, flux
      real(c_double) :: upward = -1.0_c_double  ! -1 for surface-to-bottom ordering!

      active = .false.
      outer: do k=1,nz
         do j=1+halo,ny-halo
            active = any(w(1+halo:nx-halo,j,k) /= 0.0_c_double)  ! Note FABM guarantees w is 0 in masked points
            if (active) exit outer
         end do
      end do outer
      if (.not. active) return

      do k=1,nz-1
         do j=1+halo,ny-halo
            do i=1+halo,nx-halo
               if (mask(i,j,k) == 1 .and. mask(i,j,k+1) == 1) then
                  local_w = upward * 0.5_c_double * (w(i,j,k) + w(i,j,k+1))
                  if (local_w > 0.0_c_double) then
                     ! Towards greater k
                     flux = local_w * c(i,j,k)
                  else
                     ! Towards smaller k
                     flux = local_w * c(i,j,k+1)
                  end if
                  s(i,j,k)   = s(i,j,k) - flux / h(i,j,k)
                  s(i,j,k+1) = s(i,j,k+1) + flux / h(i,j,k+1)
               end if
            end do
         end do
      end do
   end subroutine

   SUBROUTINE c_multiply_add(n, tgt, add, scale_factor) bind(c)
      integer(c_int), value, intent(in) :: n
      real(c_double), intent(inout) :: tgt(n)
      real(c_double), intent(in) :: add(n)
      real(c_double), value, intent(in) :: scale_factor
      tgt = tgt + scale_factor * add
   END SUBROUTINE

   SUBROUTINE c_advance_surface_elevation(nx, ny, halox, haloy, mask, dyu, dxv, iarea, z, U, V, fwf, dt) bind(c)
      integer(c_int), value, intent(in) :: nx, ny
      integer(c_int), value, intent(in) :: halox
      integer(c_int), value, intent(in) :: haloy
      integer(c_int), intent(in) :: mask(nx, ny)
      real(c_double), intent(inout) :: z(nx, ny)
      real(c_double), intent(in) :: dyu(nx, ny)
      real(c_double), intent(in) :: dxv(nx, ny)
      real(c_double), intent(in) :: iarea(nx, ny)
      real(c_double), intent(in) :: U(nx, ny)
      real(c_double), intent(in) :: V(nx, ny)
      real(c_double), intent(in) :: fwf(nx, ny)
      real(c_double), value, intent(in) :: dt

      integer :: i, j

      do j=1+haloy,ny-haloy
         do i=1+halox,nx-halox
            if (mask(i,j) == 1) then
               z(i,j) = z(i,j) & ! [GETM Scientific Report: eq. 4.28]
                           + dt * ((  U(i-1,j  ) * dyu(i-1,j) - U(i,j) * dyu(i,j)  &
                                    + V(i  ,j-1) * dxv(i,j-1) - V(i,j) * dxv(i,j)) &
                                   * iarea(i,j) &
                                   + fwf(i,j))
            end if
         end do
      end do
   END SUBROUTINE

   SUBROUTINE c_surface_pressure_gradient(nx, ny, imin, imax, jmin, jmax, umask, vmask, idxu, idyv, &
         z, sp, H, D, Dmin, dpdx, dpdy) bind(c)
      integer(c_int), value, intent(in) :: nx, ny
      integer(c_int), value, intent(in) :: imin, imax, jmin, jmax
      integer(c_int), intent(in) :: umask(nx, ny)
      integer(c_int), intent(in) :: vmask(nx, ny)
      real(c_double), intent(in) :: idxu(nx, ny)
      real(c_double), intent(in) :: idyv(nx, ny)
      real(c_double), intent(in) :: z(nx, ny)
      real(c_double), intent(in) :: sp(nx, ny)
      real(c_double), intent(in) :: H(nx, ny)
      real(c_double), intent(in) :: D(nx, ny)
      real(c_double), value, intent(in) :: Dmin
      real(c_double), intent(inout) :: dpdx(nx, ny)
      real(c_double), intent(inout) :: dpdy(nx, ny)

      integer :: i, j
      real(real64) :: zp, zm
   
      real(real64), parameter :: gammai = 1._real64/(g * rho0)

      do j = jmin, jmax
         do i = imin, imax
            if (umask(i,j) == 1) then
               zp = max(z(i+1,j), -H(i  ,j)+min(Dmin,D(i+1,j)))
               zm = max(z(i  ,j), -H(i+1,j)+min(Dmin,D(i  ,j)))
               dpdx(i,j) = (zp - zm + (sp(i+1,j)-sp(i,j))*gammai) * idxu(i,j)
            end if
         end do
      end do

      do j = jmin, jmax
         do i = imin, imax
            if (vmask(i,j) == 1) then
               zp = max(z(i,j+1), -H(i  ,j)+min(Dmin,D(i,j+1)))
               zm = max(z(i,j  ), -H(i,j+1)+min(Dmin,D(i,j  )))
               dpdy(i,j) = (zp - zm + (sp(i,j+1)-sp(i,j))*gammai) * idyv(i,j)
            end if
         end do
      end do
   END SUBROUTINE c_surface_pressure_gradient

end module
