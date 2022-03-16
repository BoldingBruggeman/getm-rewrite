module pygetm

   use iso_c_binding, only: c_ptr, c_int, c_double, c_char, c_loc, c_f_pointer, c_associated, C_NULL_CHAR, C_NULL_PTR
   use iso_fortran_env, only: real64

   use getm_domain, only: type_getm_domain, type_getm_grid
   use getm_operators, only: type_advection, type_vertical_diffusion
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

   subroutine domain_initialize_open_boundaries(pdomain, nbdyp, nwb, nnb, neb, nsb, bdy_info, bdy_i, bdy_j) bind(c)
      type(c_ptr),    intent(in), value :: pdomain
      integer(c_int), intent(in), value :: nwb, nnb, neb, nsb, nbdyp
      integer(c_int), intent(in)        :: bdy_info(nwb + nnb + neb + nsb, 6), bdy_i(nbdyp), bdy_j(nbdyp)

      type (type_getm_domain), pointer :: domain
      integer :: nbdy

      call c_f_pointer(pdomain, domain)
      call domain%initialize_open_boundaries(nbdyp=nbdyp, nwb=nwb, nnb=nnb, neb=neb, nsb=nsb)

      nbdy = 0
      if (nwb > 0) then
         domain%wi(:)  = bdy_info(nbdy + 1:nbdy + nwb, 1) + 1
         domain%wfj(:) = bdy_info(nbdy + 1:nbdy + nwb, 2) + 1
         domain%wlj(:) = bdy_info(nbdy + 1:nbdy + nwb, 3)
         nbdy = nbdy + nwb
      end if
      if (nnb > 0) then
         domain%nj(:) = bdy_info(nbdy + 1:nbdy + nnb, 1) + 1
         domain%nfi(:) = bdy_info(nbdy + 1:nbdy + nnb, 2) + 1
         domain%nli(:) = bdy_info(nbdy + 1:nbdy + nnb, 3)
         nbdy = nbdy + nnb
      end if
      if (neb > 0) then
         domain%ei(:) = bdy_info(nbdy + 1:nbdy + neb, 1) + 1
         domain%efj(:) = bdy_info(nbdy + 1:nbdy + neb, 2) + 1
         domain%elj(:) = bdy_info(nbdy + 1:nbdy + neb, 3)
         nbdy = nbdy + neb
      end if
      if (nsb > 0) then
         domain%sj(:) = bdy_info(nbdy + 1:nbdy + nsb, 1) + 1
         domain%sfi(:) = bdy_info(nbdy + 1:nbdy + nsb, 2) + 1
         domain%sli(:) = bdy_info(nbdy + 1:nbdy + nsb, 3)
         nbdy = nbdy + nsb
      end if
      if (nbdy > 0) then
         domain%bdy_2d_type(:) = bdy_info(:, 4)
         domain%bdy_3d_type(:) = bdy_info(:, 5)
         domain%bdy_index(:) = bdy_info(:, 6) + 1
         domain%bdy_map(:, 1) = bdy_i + 1
         domain%bdy_map(:, 2) = bdy_j + 1
      end if
   end subroutine

   function domain_get_grid(pdomain, grid_type, imin, imax, jmin, jmax, kmin, kmax, halox, haloy, haloz) result(pgrid) bind(c)
      type(c_ptr),         intent(in), value :: pdomain
      integer(kind=c_int), intent(in), value :: grid_type
      integer(kind=c_int), intent(in), value :: imin, imax, jmin, jmax, kmin, kmax, halox, haloy, haloz
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
      case default
         allocate(grid)
         call grid%configure(grid_type=grid_type, imin=imin, imax=imax, &
            jmin=jmin, jmax=jmax, kmin=kmin, kmax=kmax, halo=(/halox, haloy, haloz/))
      end select
      pgrid = c_loc(grid)
   end function

   subroutine domain_do_vertical(pdomain) bind(c)
      type(c_ptr), intent(in), value :: pdomain

      type (type_getm_domain), pointer :: domain

      call c_f_pointer(pdomain, domain)
      call domain%do_vertical(1._c_double)
   end subroutine

   subroutine get_array(source_type, obj, name, grid_type, sub_type, data_type, p) bind(c)
      integer(c_int),  value,         intent(in)  :: source_type
      type(c_ptr), value,             intent(in)  :: obj
      character(kind=c_char), target, intent(in)  :: name(*)
      integer(c_int),                 intent(out) :: grid_type, sub_type, data_type
      type(c_ptr),                    intent(out) :: p

      type (type_getm_grid),     pointer :: grid
      type (type_getm_momentum), pointer :: momentum
      type (type_getm_pressure), pointer :: pressure
      type (type_getm_sealevel), pointer :: sealevel
      character(len=10),         pointer :: pname

      integer, parameter :: subtype_boundary = 1
      integer, parameter :: subtype_depth_explicit = 2
      integer, parameter :: subtype_depth_explicit_interfaces = 3

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
         case ('z0b'); p = c_loc(grid%z0b)
         case ('z0b_min'); p = c_loc(grid%z0b_min)
         case ('hn'); p = c_loc(grid%hn); sub_type = subtype_depth_explicit
         case ('ho'); p = c_loc(grid%ho); sub_type = subtype_depth_explicit
         case ('zc'); p = c_loc(grid%zc); sub_type = subtype_depth_explicit
         case ('zf'); p = c_loc(grid%zf); sub_type = subtype_depth_explicit_interfaces
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
         case ('u1');   p = c_loc(momentum%u1); grid_type = 2
         case ('v1');   p = c_loc(momentum%v1); grid_type = 3
         case ('bdyu');   if (allocated(momentum%bdyu)) p = c_loc(momentum%bdyu); grid_type = 2; sub_type = subtype_boundary
         case ('bdyv');   if (allocated(momentum%bdyv)) p = c_loc(momentum%bdyv); grid_type = 3; sub_type = subtype_boundary
         case ('uk');   if (allocated(momentum%uk)) p = c_loc(momentum%uk); grid_type = 2; sub_type = subtype_depth_explicit
         case ('vk');   if (allocated(momentum%vk)) p = c_loc(momentum%vk); grid_type = 3; sub_type = subtype_depth_explicit
         case ('ru');   p = c_loc(momentum%ru); grid_type = 2
         case ('rv');   p = c_loc(momentum%rv); grid_type = 3
         case ('rru');   p = c_loc(momentum%rru); grid_type = 2
         case ('rrv');   p = c_loc(momentum%rrv); grid_type = 3
         case ('pk');   if (allocated(momentum%pk)) p = c_loc(momentum%pk); grid_type = 2; sub_type = subtype_depth_explicit
         case ('qk');   if (allocated(momentum%qk)) p = c_loc(momentum%qk); grid_type = 3; sub_type = subtype_depth_explicit
         case ('ww');   if (allocated(momentum%ww)) p = c_loc(momentum%ww); sub_type = subtype_depth_explicit_interfaces
         case ('advpk'); if (allocated(momentum%advpk)) p = c_loc(momentum%advpk); grid_type = 2; sub_type = subtype_depth_explicit
         case ('advqk'); if (allocated(momentum%advqk)) p = c_loc(momentum%advqk); grid_type = 3; sub_type = subtype_depth_explicit
         end select
      case (2)
         call c_f_pointer(obj, pressure)
         select case (pname(:index(pname, C_NULL_CHAR) - 1))
         case ('dpdx'); p = c_loc(pressure%dpdx); grid_type = 2
         case ('dpdy'); p = c_loc(pressure%dpdy); grid_type = 3
         case default; p = C_NULL_PTR
         end select
      case (3)
         call c_f_pointer(obj, sealevel)
         select case (pname(:index(pname, C_NULL_CHAR) - 1))
         case ('zbdy'); if (allocated(sealevel%zbdy)) p = c_loc(sealevel%zbdy); sub_type = subtype_boundary
         case default; p = C_NULL_PTR
         end select
      end select
   end subroutine

   subroutine grid_interp_x(pgrid, psource, ptarget, ioffset, n) bind(c)
      type(c_ptr),    intent(in), value :: pgrid, psource, ptarget
      integer(c_int), intent(in), value :: ioffset, n

      type (type_getm_grid),    pointer :: grid
      real(real64), contiguous, pointer :: source(:,:,:), target(:,:,:)
      integer :: i, j, k

      call c_f_pointer(pgrid, grid)
      call c_f_pointer(psource, source, (/grid%u(1) - grid%l(1) + 1, grid%u(2) - grid%l(2) + 1, n/))
      call c_f_pointer(ptarget, target, (/grid%u(1) - grid%l(1) + 1, grid%u(2) - grid%l(2) + 1, n/))
      do k = 1, n
         do j = 1, size(source, 2)
            do i = 1, size(source, 1) - 1
               target(i + ioffset, j, k) = 0.5_real64 * (source(i, j, k) + source(i + 1, j, k))
            end do
         end do
      end do
   end subroutine

   subroutine grid_interp_y(pgrid, psource, ptarget, joffset, n) bind(c)
      type(c_ptr),    intent(in), value :: pgrid, psource, ptarget
      integer(c_int), intent(in), value :: joffset, n

      type (type_getm_grid),    pointer :: grid
      real(real64), contiguous, pointer :: source(:,:,:), target(:,:,:)
      integer :: i, j, k

      call c_f_pointer(pgrid, grid)
      call c_f_pointer(psource, source, (/grid%u(1) - grid%l(1) + 1, grid%u(2) - grid%l(2) + 1, n/))
      call c_f_pointer(ptarget, target, (/grid%u(1) - grid%l(1) + 1, grid%u(2) - grid%l(2) + 1, n/))
      do k = 1, n
         do j = 1, size(source, 2) - 1
            do i = 1, size(source, 1)
               target(i, j + joffset, k) = 0.5_real64 * (source(i, j, k) + source(i, j + 1, k))
            end do
         end do
      end do
   end subroutine

   subroutine grid_interp_xy(psource_grid, psource, ptarget_grid, ptarget, ioffset, joffset, n) bind(c)
      type(c_ptr),    intent(in), value :: psource_grid, psource, ptarget_grid, ptarget
      integer(c_int), intent(in), value :: ioffset, joffset, n

      type (type_getm_grid),    pointer :: source_grid, target_grid
      real(real64), contiguous, pointer :: source(:,:,:), target(:,:,:)
      integer :: i, j, k

      call c_f_pointer(psource_grid, source_grid)
      call c_f_pointer(ptarget_grid, target_grid)
      call c_f_pointer(psource, source, (/source_grid%u(1) - source_grid%l(1) + 1, source_grid%u(2) - source_grid%l(2) + 1, n/))
      call c_f_pointer(ptarget, target, (/target_grid%u(1) - target_grid%l(1) + 1, target_grid%u(2) - target_grid%l(2) + 1, n/))
      do k = 1, n
         do j = 1, size(source, 2) - 1
            do i = 1, size(source, 1) - 1
               target(i + ioffset, j + joffset, k) = 0.25_real64 * (source(i, j, k) + source(i + 1, j, k) &
                  + source(i, j + 1, k) + source(i + 1, j + 1, k))
            end do
         end do
      end do
   end subroutine

   subroutine domain_initialize(pdomain, runtype, Dmin, maxdt) bind(c)
      type(c_ptr),    intent(in), value :: pdomain
      integer(c_int), intent(in), value :: runtype
      real(c_double), intent(in), value :: Dmin
      real(c_double), intent(out)       :: maxdt

      type (type_getm_domain), pointer :: domain

      call c_f_pointer(pdomain, domain)
      domain%Dmin = Dmin
      call domain%initialize(runtype)
      maxdt = domain%maxdt
   end subroutine

   subroutine domain_update_depths(pdomain) bind(c)
      type(c_ptr), intent(in), value :: pdomain

      type (type_getm_domain), pointer :: domain

      call c_f_pointer(pdomain, domain)
      call domain%update_depths()
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

   subroutine vertical_diffusion_calculate(pdiffusion, ptgrid, molecular, pnuh, timestep, cnpar, pvar, pea2, pea4) bind(c)
      real(c_double), intent(in), value :: timestep, cnpar, molecular
      type(c_ptr),    intent(in), value :: pdiffusion, ptgrid, pnuh, pvar, pea2, pea4

      type (type_getm_grid),          pointer :: tgrid
      type (type_vertical_diffusion), pointer :: diffusion
      real(real64), contiguous, pointer, dimension(:,:,:) :: nuh, var, ea2, ea4

      call c_f_pointer(pdiffusion, diffusion)
      call c_f_pointer(ptgrid, tgrid)
      call c_f_pointer(pnuh, nuh, (/tgrid%u(1) - tgrid%l(1) + 1, tgrid%u(2) - tgrid%l(2) + 1, tgrid%u(3) - tgrid%l(3) + 2/))
      call c_f_pointer(pvar, var, tgrid%u - tgrid%l + 1)
      if (c_associated(pea2)) call c_f_pointer(pea2, ea2, tgrid%u - tgrid%l + 1)
      if (c_associated(pea4)) call c_f_pointer(pea4, ea4, tgrid%u - tgrid%l + 1)
      if (c_associated(pea2) .and. c_associated(pea4)) then
         call diffusion%calculate(timestep, cnpar, tgrid%mask, tgrid%ho, tgrid%hn, molecular, nuh(:, :, 2:size(nuh,3)-1), var, &
            ea2=ea2, ea4=ea4)
      elseif (c_associated(pea2)) then
         call diffusion%calculate(timestep, cnpar, tgrid%mask, tgrid%ho, tgrid%hn, molecular, nuh(:, :, 2:size(nuh,3)-1), var, &
            ea2=ea2)
      elseif (c_associated(pea4)) then
         call diffusion%calculate(timestep, cnpar, tgrid%mask, tgrid%ho, tgrid%hn, molecular, nuh(:, :, 2:size(nuh,3)-1), var, &
            ea4=ea4)
      else
         call diffusion%calculate(timestep, cnpar, tgrid%mask, tgrid%ho, tgrid%hn, molecular, nuh(:, :, 2:size(nuh,3)-1), var)
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

   subroutine advection_2d_calculate(direction, padvection, ptgrid, pugrid, pu, timestep, pD, pvar) bind(c)
      integer(c_int), intent(in), value :: direction
      real(c_double), intent(in), value :: timestep
      type(c_ptr),    intent(in), value :: padvection, ptgrid, pugrid, pu, pD, pvar

      type (type_advection),    pointer                 :: advection
      type (type_getm_grid),  pointer                   :: tgrid, ugrid
      real(real64), contiguous, pointer, dimension(:,:) :: u, D, var

      call c_f_pointer(padvection, advection)
      if (.not. allocated(advection%op)) return
      call c_f_pointer(ptgrid, tgrid)
      call c_f_pointer(pugrid, ugrid)
      call c_f_pointer(pu, u, ugrid%u(1:2) - ugrid%l(1:2) + 1)
      call c_f_pointer(pD, D, tgrid%u(1:2) - tgrid%l(1:2) + 1)
      call c_f_pointer(pvar, var, tgrid%u(1:2) - tgrid%l(1:2) + 1)
      select case (direction)
         case (1)
            call advection%advection_calculate_u2d(ugrid, u, timestep, tgrid, D, var)
         case (2)
            call advection%advection_calculate_v2d(ugrid, u, timestep, tgrid, D, var)
      end select
   end subroutine

   subroutine advection_w_calculate(padvection, ptgrid, pw, timestep, ph, pvar) bind(c)
      real(c_double), intent(in), value :: timestep
      type(c_ptr),    intent(in), value :: padvection, ptgrid, pw, ph, pvar

      type (type_advection),    pointer                   :: advection
      type (type_getm_grid),    pointer                   :: tgrid
      real(real64), contiguous, pointer, dimension(:,:,:) :: w, h, var

      call c_f_pointer(padvection, advection)
      if (.not. allocated(advection%op)) return
      call c_f_pointer(ptgrid, tgrid)
      call c_f_pointer(pw, w, (/tgrid%u(1) - tgrid%l(1) + 1, tgrid%u(2) - tgrid%l(2) + 1, tgrid%kmax + 1/))
      call c_f_pointer(ph, h, tgrid%u - tgrid%l + 1)
      call c_f_pointer(pvar, var, tgrid%u - tgrid%l + 1)
      call advection%advection_calculate_w3d(w, timestep, tgrid, h, var)
   end subroutine

   function momentum_create(runtype, pdomain, apply_bottom_friction) result(pmomentum) bind(c)
      integer(c_int), intent(in), value :: runtype
      type(c_ptr),    intent(in), value :: pdomain
      integer(c_int), intent(in), value :: apply_bottom_friction
      type(c_ptr) :: pmomentum

      type (type_getm_domain),   pointer :: domain
      type (type_getm_momentum), pointer :: momentum

      call c_f_pointer(pdomain, domain)
      allocate(momentum)
      call momentum%configure()
      momentum%advection_scheme = 0
      momentum%apply_bottom_friction = (apply_bottom_friction == 1)
      momentum%apply_diffusion = .false.
      call momentum%initialize(runtype, domain)
      allocate(momentum%vertical_diffusion)
      call momentum%vertical_diffusion%initialize(domain%T)
      pmomentum = c_loc(momentum)
   end function

   subroutine momentum_u_2d(direction, pmomentum, timestep, ptausx, pdpdx) bind(c)
      integer(c_int), intent(in), value :: direction
      type(c_ptr),    intent(in), value :: pmomentum
      real(c_double), intent(in), value :: timestep
      type(c_ptr),    intent(in), value :: ptausx, pdpdx

      type (type_getm_momentum), pointer :: momentum
      real(real64), contiguous, pointer, dimension(:,:) :: tausx, dpdx

      call c_f_pointer(pmomentum, momentum)
      call c_f_pointer(ptausx, tausx, momentum%domain%U%u(1:2) - momentum%domain%U%l(1:2) + 1)
      call c_f_pointer(pdpdx, dpdx, momentum%domain%T%u(1:2) - momentum%domain%T%l(1:2) + 1)
      select case (direction)
         case (1)
            call momentum%u_2d(timestep,tausx,dpdx)
         case (2)
            call momentum%v_2d(timestep,tausx,dpdx)
      end select
   end subroutine

   subroutine momentum_u_3d(direction, pmomentum, timestep, ptausx, pdpdx, pidpdx, pviscosity) bind(c)
      integer(c_int), intent(in), value :: direction
      type(c_ptr),    intent(in), value :: pmomentum
      real(c_double), intent(in), value :: timestep
      type(c_ptr),    intent(in), value :: ptausx, pdpdx, pidpdx, pviscosity

      type (type_getm_momentum), pointer :: momentum
      real(real64), contiguous, pointer, dimension(:,:) :: tausx, dpdx
      real(real64), contiguous, pointer, dimension(:,:,:) :: idpdx, viscosity

      call c_f_pointer(pmomentum, momentum)
      call c_f_pointer(ptausx, tausx, momentum%domain%U%u(1:2) - momentum%domain%U%l(1:2) + 1)
      call c_f_pointer(pdpdx, dpdx, momentum%domain%U%u(1:2) - momentum%domain%U%l(1:2) + 1)
      call c_f_pointer(pidpdx, idpdx, momentum%domain%U%u - momentum%domain%U%l + 1)
      call c_f_pointer(pviscosity, viscosity, (/momentum%domain%T%u(1) - momentum%domain%T%l(1) + 1, &
         momentum%domain%T%u(2) - momentum%domain%T%l(2) + 1, momentum%domain%T%u(3) - momentum%domain%T%l(3) + 2/))
      select case (direction)
         case (1)
            call momentum%pk_3d(timestep,tausx,dpdx,idpdx,viscosity)
         case (2)
            call momentum%qk_3d(timestep,tausx,dpdx,idpdx,viscosity)
      end select
   end subroutine

   subroutine momentum_w_3d(pmomentum, timestep) bind(c)
      type(c_ptr),    intent(in), value :: pmomentum
      real(c_double), intent(in), value :: timestep

      type (type_getm_momentum), pointer :: momentum

      call c_f_pointer(pmomentum, momentum)
      call momentum%w_momentum_3d(timestep)
   end subroutine

   subroutine momentum_uv_coriolis(direction, pmomentum) bind(c)
      integer(c_int), intent(in), value :: direction
      type(c_ptr),    intent(in), value :: pmomentum

      type (type_getm_momentum), pointer :: momentum

      call c_f_pointer(pmomentum, momentum)
      select case (direction)
         case (1)
            call momentum%coriolis_fu()
         case (2)
            call momentum%coriolis_fv()
      end select
   end subroutine

   subroutine momentum_uv_coriolis_3d(direction, pmomentum) bind(c)
      integer(c_int), intent(in), value :: direction
      type(c_ptr),    intent(in), value :: pmomentum

      type (type_getm_momentum), pointer :: momentum

      call c_f_pointer(pmomentum, momentum)
      select case (direction)
         case (1)
            call momentum%coriolis_fpk()
         case (2)
            call momentum%coriolis_fqk()
      end select
   end subroutine

   subroutine momentum_bottom_friction_2d(pmomentum, runtype) bind(c)
      type(c_ptr),    intent(in), value :: pmomentum
      integer(c_int), intent(in), value :: runtype

      type (type_getm_momentum), pointer :: momentum

      call c_f_pointer(pmomentum, momentum)
      call momentum%bottom_friction_2d(runtype)
   end subroutine

   subroutine momentum_bottom_friction_3d(pmomentum) bind(c)
      type(c_ptr),    intent(in), value :: pmomentum

      type (type_getm_momentum), pointer :: momentum

      call c_f_pointer(pmomentum, momentum)
      call momentum%bottom_friction_3d()
   end subroutine

   subroutine momentum_shear_frequency(pmomentum, pviscosity) bind(c)
      type(c_ptr),    intent(in), value :: pmomentum
      type(c_ptr),    intent(in), value :: pviscosity

      type (type_getm_momentum), pointer :: momentum
      real(real64), contiguous, pointer, dimension(:,:,:) :: viscosity

      call c_f_pointer(pmomentum, momentum)
      call c_f_pointer(pviscosity, viscosity, momentum%domain%T%u - momentum%domain%T%l + 1)
      call momentum%shear_frequency(viscosity)
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

   subroutine sealevel_boundaries(psealevel, pmomentum, timestep) bind(c)
      type(c_ptr), intent(in), value :: psealevel, pmomentum
      real(c_double), intent(in), value :: timestep

      type (type_getm_sealevel), pointer :: sealevel
      type (type_getm_momentum), pointer :: momentum

      call c_f_pointer(psealevel, sealevel)
      call c_f_pointer(pmomentum, momentum)
      call sealevel%boundaries(timestep, momentum%U, momentum%V, momentum%bdyu, momentum%bdyv)
   end subroutine

end module
