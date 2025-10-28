module pygetm

   use iso_c_binding, only: c_ptr, c_int, c_double, c_char, c_loc, c_f_pointer, c_associated, C_NULL_CHAR, C_NULL_PTR
   use iso_fortran_env, only: real64

   use getm_domain, only: type_getm_grid
   use getm_operators, only: type_advection, type_vertical_diffusion

   implicit none

   real(c_double), parameter, public :: rho0 = 1025._c_double
   real(c_double), parameter, public :: kappa = 0.4_c_double
   real(c_double), parameter, public :: g = 9.81_c_double

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

   function create_grid(imin, imax, jmin, jmax, kmin, kmax, halox, haloy) result(pgrid) bind(c)
      integer(kind=c_int), intent(in), value :: imin, imax, jmin, jmax, kmin, kmax, halox, haloy
      type(c_ptr)                            :: pgrid

      type (type_getm_grid), pointer :: grid
      integer                        :: halo(3)

      allocate(grid)
      halo = (/halox, haloy, 0/)
      call grid%configure(imin=imin, imax=imax, jmin=jmin, jmax=jmax, kmin=kmin, kmax=kmax, halo=halo)
      pgrid = c_loc(grid)
   end function

   subroutine grid_finalize(pgrid) bind(c)
      type(c_ptr), intent(in), value :: pgrid

      type (type_getm_grid), pointer :: grid

      call c_f_pointer(pgrid, grid)
      deallocate(grid)
   end subroutine

   subroutine get_array(source_type, obj, name, grid_type, sub_type, data_type, p) bind(c)
      integer(c_int),  value,         intent(in)  :: source_type
      type(c_ptr), value,             intent(in)  :: obj
      character(kind=c_char), target, intent(in)  :: name(*)
      integer(c_int),                 intent(out) :: grid_type, sub_type, data_type
      type(c_ptr),                    intent(out) :: p

      type (type_getm_grid), pointer :: grid
      character(len=10),     pointer :: pname

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
      end select
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

   subroutine c_momentum_diffusion(pugrid, pvgrid, puugrid, puvgrid, pvugrid, pvvgrid, nk, phuu, phuv, phvu, phvv, pu, pv, Am0, &
         pdiffu, pdiffv) bind(c)
      integer(c_int), intent(in), value :: nk
      type(c_ptr),    intent(in), value :: pugrid, pvgrid, puugrid, puvgrid, pvugrid, pvvgrid, phuu, phuv, phvu, phvv, pu, pv
      real(c_double), intent(in), value :: Am0
      type(c_ptr),    intent(in), value :: pdiffu, pdiffv

      type (type_getm_grid),      pointer                   :: UG, VG, UUG, UVG, VUG, VVG
      real(c_double), contiguous, pointer, dimension(:,:,:) :: huu, huv, hvu, hvv, u, v, diffu, diffv
      integer :: k

      call c_f_pointer(pugrid, UG)
      call c_f_pointer(pvgrid, VG)
      call c_f_pointer(puugrid, UUG)
      call c_f_pointer(puvgrid, UVG)
      call c_f_pointer(pvugrid, VUG)
      call c_f_pointer(pvvgrid, VVG)
      call c_f_pointer(phuu, huu, (/UUG%u(1) - UUG%l(1) + 1, UUG%u(2) - UUG%l(2) + 1, nk/))
      call c_f_pointer(phuv, huv, (/UVG%u(1) - UVG%l(1) + 1, UVG%u(2) - UVG%l(2) + 1, nk/))
      call c_f_pointer(phvu, hvu, (/VUG%u(1) - VUG%l(1) + 1, VUG%u(2) - VUG%l(2) + 1, nk/))
      call c_f_pointer(phvv, hvv, (/VVG%u(1) - VVG%l(1) + 1, VVG%u(2) - VVG%l(2) + 1, nk/))
      call c_f_pointer(pu, u,   (/UG%u(1) - UG%l(1) + 1, UG%u(2) - UG%l(2) + 1, nk/))
      call c_f_pointer(pv, v,   (/VG%u(1) - VG%l(1) + 1, VG%u(2) - VG%l(2) + 1, nk/))
      call c_f_pointer(pdiffu, diffu, (/UG%u(1) - UG%l(1) + 1, UG%u(2) - UG%l(2) + 1, nk/))
      call c_f_pointer(pdiffv, diffv, (/VG%u(1) - VG%l(1) + 1, VG%u(2) - VG%l(2) + 1, nk/))
      do k = 1, nk
         call horizontal_momentum_diffusion(UG, VG, UUG, UVG, VUG, VVG, huu(:,:,k), huv(:,:,k), hvu(:,:,k), hvv(:,:,k), &
            u(:,:,k), v(:,:,k), Am0, diffu(:,:,k), diffv(:,:,k))
      end do
   end subroutine

   subroutine horizontal_momentum_diffusion(UG, VG, UUG, UVG, VUG, VVG, huu, huv, hvu, hvv, u, v, Am0, diffu, diffv)
      type (type_getm_grid), intent(in) :: UG, VG, UUG, UVG, VUG, VVG
      real(c_double), dimension(:,:), intent(in) :: huu(UUG%l(1):,UUG%l(2):)
      real(c_double), dimension(:,:), intent(in) :: huv(UVG%l(1):,UVG%l(2):)
      real(c_double), dimension(:,:), intent(in) :: hvu(VUG%l(1):,VUG%l(2):)
      real(c_double), dimension(:,:), intent(in) :: hvv(VVG%l(1):,VVG%l(2):)
      real(c_double), dimension(:,:), intent(in) :: u(UG%l(1):,UG%l(2):)
      real(c_double), dimension(:,:), intent(in) :: v(VG%l(1):,VG%l(2):)
      real(c_double),                 intent(in) :: Am0
      real(c_double), dimension(:,:), intent(inout) :: diffu(UG%l(1):,UG%l(2):)
      real(c_double), dimension(:,:), intent(inout) :: diffv(VG%l(1):,VG%l(2):)

      integer :: i,j
      real(c_double), allocatable :: flux(:,:)

      allocate(flux(UG%l(1):UG%u(1), UG%l(2):UG%u(2)))

      ! Central for dx(2*Am*dx(U/DU))
      do j=UUG%jmin,UUG%jmax
         do i=UUG%imin-1,UUG%imax ! shear defined on T-points
            flux(i,j)=0._real64
            if (UUG%mask(i,j) /= 0) then
               flux(i,j) = 2._real64 * Am0 * UUG%dy(i,j) * huu(i,j) * (u(i+1,j) - u(i,j)) * UUG%idx(i,j)
            end if
         end do
      end do
      do j=UG%jmin,UG%jmax
         do i=UG%imin,UG%imax ! diffu defined on U-points
            diffu(i,j)=0._real64
            if (UG%mask(i,j) == 1) then
               diffu(i,j) = (flux(i,j) - flux(i-1,j)) * UG%iarea(i,j)
            end if
         end do
      end do

      ! Central for dy(Am*(dy(U/DU)+dx(V/DV)))
      do j=UVG%jmin-1,UVG%jmax ! work2d defined on X-points
         do i=UVG%imin,UVG%imax
            flux(i,j)=0._real64
            if (UVG%mask(i,j) /= 0) then
               flux(i,j) = Am0 * UVG%dx(i,j) * huv(i,j) * ((u(i,j+1) - u(i,j)) * UVG%idy(i,j) + (v(i+1,j) - v(i,j)) * UVG%idx(i,j))
            end if
         end do
      end do
      do j=UG%jmin,UG%jmax ! diffu defined on U-points
         do i=UG%imin,UG%imax
            if (UG%mask(i,j) == 1 .or. UG%mask(i,j) == 2) then
               diffu(i,j) = diffu(i,j) + (flux(i,j) - flux(i,j-1)) * UG%iarea(i,j)
            end if
         end do
      end do

      ! Central for dx(Am*(dy(U/DU)+dx(V/DV)))
      do j=VUG%jmin,VUG%jmax
         do i=VUG%imin-1,VUG%imax ! work2d defined on X-points
            flux(i,j)=0._real64
            if (VUG%mask(i,j) /= 0) then
               flux(i,j) = Am0 * VUG%dy(i,j) * hvu(i,j) * ((u(i,j+1) - u(i,j)) * VUG%idy(i,j) + (v(i+1,j) - v(i,j)) * VUG%idx(i,j))
            end if
         end do
      end do
      do j=VG%jmin,VG%jmax ! diffv defined on V-points
         do i=VG%imin,VG%imax ! diffv defined on V-points
            diffv(i,j)=0._real64
            if (VG%mask(i,j) == 1) then
               diffv(i,j) = (flux(i,j) - flux(i-1,j)) * VG%iarea(i,j)
            end if
         end do
      end do

      ! Central for dy(2*Am*dy(V/DV))
      do j=VVG%jmin-1,VVG%jmax ! work2d defined on T-points
         do i=VVG%imin,VVG%imax
            flux(i,j)=0._real64
            if (VVG%mask(i,j) /= 0) then
               flux(i,j) = 2._real64 * Am0 * VVG%dx(i,j) * hvv(i,j) * (v(i,j+1) - v(i,j)) * VVG%idy(i,j)
            end if
         end do
      end do
      do j=VG%jmin,VG%jmax ! diffv defined on V-points
         do i=VG%imin,VG%imax
            if (VG%mask(i,j) == 1) then
               diffv(i,j) = diffv(i,j) + (flux(i,j) - flux(i,j-1)) * VG%iarea(i,j)
            end if
         end do
      end do

   end subroutine horizontal_momentum_diffusion

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
      where (mask == 1) alpha = max(0._c_double, min(1._c_double, (D - Dmin) / (Dcrit - Dmin)))
   end subroutine

   subroutine c_elevation2depth(n, z, H, Dmin, mask, D) bind(c)
      integer(c_int), intent(in), value :: n
      real(c_double), intent(in)        :: z(n)
      real(c_double), intent(in)        :: H(n)
      real(c_double), intent(in), value :: Dmin
      integer(c_int), intent(in)        :: mask(n)
      real(c_double), intent(inout)     :: D(n)
      where (mask /= 0) D = max(H + z, Dmin)
   end subroutine

   subroutine c_vertical_advection_to_sources(nx, ny, nz, halox, haloy, mask, c, w, h, s) bind(c)
      ! first-order upstream-biased advection, e.g., for integrating FABM sinking/floating into source term
      integer(c_int), intent(in), value :: nx, ny, nz, halox, haloy
      integer(c_int), intent(in)        :: mask(nx, ny, nz)
      real(c_double), intent(in)        :: c(nx, ny, nz), w(nx, ny, nz), h(nx, ny, nz)
      real(c_double), intent(inout)     :: s(nx, ny, nz)

      logical :: active
      integer :: i, j, k
      real(c_double) :: local_w, flux
      real(c_double) :: upward = -1.0_c_double  ! -1 for surface-to-bottom ordering!

      active = .false.
      outer: do k=1,nz
         do j=1+haloy,ny-haloy
            active = any(w(1+halox:nx-halox,j,k) /= 0.0_c_double)  ! Note FABM guarantees w is 0 in masked points
            if (active) exit outer
         end do
      end do outer
      if (.not. active) return

      do k=1,nz-1
         do j=1+haloy,ny-haloy
            do i=1+halox,nx-halox
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

   subroutine c_update_gvc(nx, ny, nz, dsigma, dbeta, Dgamma, kk, D, mask, h) bind(c)
      integer(c_int), value, intent(in) :: nx, ny, nz
      real(c_double), value, intent(in) :: dsigma
      real(c_double), intent(in) :: dbeta(nz)
      real(c_double), value, intent(in) :: Dgamma
      integer(c_int), value, intent(in) :: kk
      real(c_double), intent(in) :: D(nx, ny)
      integer(c_int), intent(in) :: mask(nx, ny)
      real(c_double), intent(inout) :: h(nx, ny, nz)

      real(c_double), allocatable :: alpha(:, :)
      integer :: k

      allocate(alpha(nx, ny))

      ! The aim is to give the reference layer (surface or bottom) a constant
      ! thickness of Dgamma / nz = Dgamma * dsigma
      ! The final calculation of the reference layer thickness blends
      ! fractional thicknesses dsigma and dbeta, giving a thickness in m of
      !   (alpha * dsigma + (1 - alpha) * dbeta) * D
      ! Setting this equal to Dgamma * dsigma and rearranging, we obtain
      !   alpha = (Dgamma / D * dsigma - dbeta) / (dsigma - dbeta)
      ! If we additionally reduce the target thickness to D * dsigma when
      ! the column height drops below Dgamma, we effectively substitute
      ! min(Dgamma, D) for Dgamma. That leads to:
      !   alpha = (min(Dgamma / D, 1.0) * dsigma - dbeta) / (dsigma - dbeta)
      alpha = (min(Dgamma / D, 1.0_c_double) * dsigma - dbeta(kk)) / (dsigma - dbeta(kk))

      do k = 1, nz
         ! Blend equal thicknesses (dsigma) with zoomed thicknesses (dbeta)
         ! alpha * self.dsigma + (1.0 - alpha) * self.dbeta
         ! Then multiply with water depth to obtain layer thicknesses in m
         where (mask /= 0) h(:, :, k) = D * (dbeta(k) + alpha * (dsigma - dbeta(k)))
      end do
   end subroutine

end module
