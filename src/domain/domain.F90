! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!>  This module provides all variables related to the bathymetry and
!>  model grid. The public subroutine $init\_domain()$ is called once
!>  and upon successful completion the bathymetry has been read and
!>  optionally modified, the calculation masks have been setup and
!>  all grid related variables have been initialised.\newline
!>  The $domain$-module depends on another module doing the actual
!>  reading of variables from files. This is provided through the
!>  generic subroutine $read\_topo\_file$. This subroutine takes two
!>  parameters - 1) a fileformat and 2) a filename. Adding a new
!>  input file format is thus straight forward and can be done
!>  without any changes to $domain$.
!>  Public variables defined in this module is used through out the
!>  code.
!>  @note
!>  Hej kurt
!>  @endnote
!>  @bug
!>  Hej kurt
!>  @endbug
!>  @warning
!>  Hej kurt
!>  @endwarning
!>  @todo
!>  Hej kurt
!>  @endtodo

MODULE getm_domain

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use grid_module
   use memory_manager
   use logging
   use field_manager

   IMPLICIT NONE

   PRIVATE  ! Private scope by default

!  Module constants
   integer, parameter :: halo(3)=(/2,2,0/)
   integer, parameter :: cartesian=1, spherical=2, curvilinear=3
      !! Fortran precision
   real(real64), parameter, public :: g = 9.81_real64
      !! Gravity
   real(real64), parameter :: rearth_default = 6378815._real64
   real(real64), parameter :: rearth = 6378815._real64
      !! Earth radius
   real(real64), parameter :: pi = atan(1._real64)*4._real64
      !! \( \pi \)
   real(real64), parameter :: omega = 2._real64*pi/86164._real64
      !! Earth frequency
   real(real64), parameter :: deg2rad = pi/180._real64
      !! degrees to radians conversion
   real(real64), parameter :: Hland = -10._real64

!  Module types and variables
   integer :: vel_depth_method = 0

   ENUM, BIND(C)
      ENUMERATOR :: TGRID=1
      ENUMERATOR :: UGRID=2
      ENUMERATOR :: VGRID=3
      ENUMERATOR :: XGRID=4
   END ENUM

   type, public :: type_grid_config
      !! author: Karsten Bolding
      !! version: v0.1
      !!
      !! Grid configuration

      character(len=256) :: bathymetry = "bathymetry.nc"
   end type type_grid_config

   type, public, extends(type_3d_grid) :: type_getm_grid
      !! author: Karsten Bolding
      !! version: v0.1
      !!
      !! The GETM calculation grid

      class(type_logging), pointer :: logs => null()

      integer :: grid_type=-1
        !!  T, U, V or X grid
      integer :: iextr=-1,jextr=-1
        !!  global grid extend
      integer :: ioff=0,joff=0
        !!  offset in subdomain decomposition
      integer :: ilg=-1,ihg=-1,jlg=-1,jhg=-1
        !!  global index range
      integer :: ill=-1,ihl=-1,jll=-1,jhl=-1
        !!  local index range

      integer :: dim_2d_ids(2)
      integer :: dim_3d_ids(3)
      real(real64), dimension(:), allocatable :: c1,c2,c3
        !! cordinate variables from NetCDF file
      real(real64), dimension(:,:), allocatable :: lon,lat
        !! longitude and latitude of grid
      real(real64), dimension(:,:), allocatable :: dlon,dlat
        !! grid spacing in degrees
      real(real64), dimension(:,:), allocatable :: x,y
        !! cartesian cordinates off grid - e.g. UTM
      real(real64), dimension(:,:), allocatable :: dx,dy
        !! grid spacing in meters
      real(real64), dimension(:,:), allocatable :: cor
        !! Coriolis term
      real(real64), dimension(:,:), allocatable :: area,inv_area
        !! grid area in m2 and inverse area
      real(real64), dimension(:,:), allocatable :: H
        !! undisturbed water depth
      real(real64), dimension(:,:), allocatable :: z
        !! sea surface elevation
      real(real64), dimension(:,:), allocatable :: zo
        !! previous time step
      integer, dimension(:,:), allocatable :: mask
        !! mask=0 -> land
      real(real64), dimension(:,:), allocatable :: D
        !! total water depth - time varying
      real(real64), dimension(:,:), allocatable :: zin
        !! elevation at T-points baroclinic time step
      real(real64), dimension(:,:), allocatable :: zio
        !! previous timstep
      real(real64), dimension(:,:,:), allocatable :: hn
        !! layer heights - new time step
      real(real64), dimension(:,:,:), allocatable :: ho
        !! layer heights - old time step
      real(real64), dimension(:,:,:), allocatable :: zf, zc
        !! depth to grid faces and grid centers
      real(real64), dimension(:,:), allocatable :: alpha
        !! drying factor
      logical :: is_initialized = .false.

      contains

      procedure :: configure => grid_configure
      procedure :: print_info => grid_print_info
      procedure :: print_mask => grid_print_mask
      procedure :: report => grid_report
   end type type_getm_grid

   TYPE, public :: type_getm_domain
      class(type_logging), pointer :: logs => null()
      class(type_field_manager), pointer :: fm => null()

      integer :: domain_type
         !! Cartesian, spherical or curvi-linear
         !! Infer from NetCDF bathymetry file
      logical :: domain_ready = .false.

      real(real64) :: maxdt=huge(1._real64)

      character(len=32) :: layout='Arakawa C - NE index convention'
      TYPE(type_getm_grid) :: T, U, V, X
         !! C-grid
         !! T, U, V, X grid - Tracer, U-momentum, V-momentum and corners
      integer ::id_dim_x, id_dim_y, id_dim_z, id_dim_time
         !! dimension ids for the central points
      integer ::id_dim_xi, id_dim_yi, id_dim_zi
         !! dimension ids for the interface points
      integer ::id_dim_xx, id_dim_yx
         !! dimension ids for the X points

      real(real64) :: Dmin=1._real64
      real(real64) :: Dcrit=2._real64
      real(real64) :: Dgamma=1._real64
      real(real64) :: Dmax=1._real64
      logical :: gamma_surf
      real(real64) :: ddl=-1._real64, ddu=-1._real64
      real(real64) :: lat0=-999._real64
      integer :: method_vertical_coordinates=1
      integer :: nwb=0,nnb=0,neb=0,nsb=0
         !! # of western, northern, eastern and southern boundaries
      integer :: nbdy=0
         !! # of open boundaries
      integer :: nbdyp=0
         !! # of open boundary points
      integer, allocatable :: wi(:),wfj(:),wlj(:)
         !! indices for western boundaries
      integer, allocatable :: nj(:),nfi(:),nli(:)
         !! indices for northen boundaries
      integer, allocatable :: ei(:),efj(:),elj(:)
         !! indices for eastern boundaries
      integer, allocatable :: sj(:),sfi(:),sli(:)
         !! indices for southern boundaries
      integer, allocatable :: bdy_2d_type(:)
      integer, allocatable :: bdy_3d_type(:)
      integer, allocatable :: bdy_index(:)
      integer, allocatable :: bdy_map(:,:)

      contains

      procedure :: configure => domain_configure
      procedure :: initialize => domain_initialize
      procedure :: report => domain_report
      procedure :: metrics => metrics
      procedure :: register => register
      procedure :: uvx_depths => uvx_depths
      procedure :: mirror_bdy_2d => mirror_bdy_2d
      procedure :: mirror_bdy_3d => mirror_bdy_3d
      generic   :: mirror_bdys => mirror_bdy_2d, mirror_bdy_3d
      procedure :: cfl_check => cfl_check
      procedure :: depth_update => depth_update
      procedure :: start_3d => start_3d
      procedure :: init_vertical => init_vertical
      procedure :: do_vertical => do_vertical
      procedure :: cleanup => domain_cleanup
   end type type_getm_domain

   INTERFACE

      module subroutine metrics(self)
         class(type_getm_domain), intent(inout) :: self
      end subroutine metrics

      module subroutine register(self,runtype)
         class(type_getm_domain), intent(inout) :: self
         integer, intent(in) :: runtype
      end subroutine register

      module subroutine uvx_depths(self)
         class(type_getm_domain), intent(inout) :: self
      end subroutine uvx_depths

      module subroutine cfl_check(self)
         class(type_getm_domain), intent(inout) :: self
      end subroutine cfl_check

      module subroutine depth_update(self)
         class(type_getm_domain), intent(inout) :: self
      end subroutine depth_update

      module subroutine mirror_bdy_2d(self,grid,f)
         class(type_getm_domain), intent(inout) :: self
         class(type_getm_grid), intent(in) :: grid
         real(real64), dimension(:,:), intent(inout) :: f(grid%l(1):,grid%l(2):)
      end subroutine mirror_bdy_2d

      module subroutine mirror_bdy_3d(self,grid,f)
         class(type_getm_domain), intent(inout) :: self
         class(type_getm_grid), intent(in) :: grid
         real(real64), dimension(:,:,:), intent(inout) :: f(grid%l(1):,grid%l(2):,grid%l(3):)
      end subroutine mirror_bdy_3d

!      module subroutine init_vertical(self,z,zo)
      module subroutine init_vertical(self)
         class(type_getm_domain), intent(inout) :: self
!         real(real64), dimension(:,:), intent(in) :: z
!         real(real64), dimension(:,:), intent(in) :: z,zo
      end subroutine init_vertical

!      module subroutine init_vertical(self,z,zo)
      module subroutine do_vertical(self,dt)
         class(type_getm_domain), intent(inout) :: self
         real(real64), intent(in)  :: dt
!         real(real64), dimension(:,:), intent(in) :: z,zo
      end subroutine do_vertical

      module subroutine grid_configure(self,logs,imin,imax,jmin,jmax,kmin,kmax,halo)
         class(type_getm_grid), intent(inout) :: self
         class(type_logging), intent(in), target, optional :: logs
         integer, intent(in), optional :: imin,imax
         integer, intent(in), optional :: jmin,jmax
         integer, intent(in), optional :: kmin,kmax
         integer, intent(in), dimension(3), optional :: halo
      end subroutine grid_configure

      module subroutine grid_report(self,logs,unit,header)
         class(type_getm_grid), intent(inout) :: self
         class(type_logging), intent(in) :: logs
         integer, intent(in) :: unit
         character(len=*), intent(in) :: header
      end subroutine grid_report

      module subroutine grid_print_info(self,logs)
         class(type_getm_grid), intent(in) :: self
         class(type_logging), intent(in) :: logs
      end subroutine grid_print_info

      module subroutine grid_print_mask(self,unit)
         class(type_getm_grid), intent(in) :: self
         integer, intent(in) :: unit
      end subroutine grid_print_mask

      module subroutine deallocate_grid_variables(self)
         class(type_getm_grid), intent(inout) :: self
      end subroutine deallocate_grid_variables
   END INTERFACE

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

SUBROUTINE domain_configure(self,imin,imax,jmin,jmax,kmin,kmax,nwb,nnb,neb,nsb,logs,fm)

   !! Configure the type_grid

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_domain), intent(inout) :: self
   integer, intent(in) :: imin,imax,jmin,jmax,kmin,kmax
   integer, intent(in), optional :: nwb,nnb,neb,nsb
   class(type_logging), intent(in), target, optional :: logs
   TYPE(type_field_manager), intent(inout), target, optional :: fm

!  Local constants

!  Local variables
   integer :: stat
!-----------------------------------------------------------------------
   if (present(logs)) then
      self%logs => logs
      call self%logs%info('domain_configure()',level=1)
      call self%logs%info(trim(self%layout),level=2)
   end if
   if (present(fm)) then
      self%fm => fm
   end if

   ! grid types must be set - not part of API for now
   call self%T%configure(logs,imin=imin,imax=imax,jmin=jmin,jmax=jmax,kmin=kmin,kmax=kmax,halo=halo)
   self%T%grid_type=TGRID
   call self%U%configure(logs,imin=imin,imax=imax,jmin=jmin,jmax=jmax,kmin=kmin,kmax=kmax,halo=self%T%halo)
   self%U%grid_type=UGRID
   call self%V%configure(logs,imin=imin,imax=imax,jmin=jmin,jmax=jmax,kmin=kmin,kmax=kmax,halo=self%T%halo)
   self%V%grid_type=VGRID
   call self%X%configure(logs,imin=imin-1,imax=imax,jmin=jmin-1,jmax=jmax,kmin=kmin,kmax=kmax,halo=self%T%halo)
   self%X%grid_type=XGRID

   if (present(nwb)) then
      self%nwb=nwb
      if (self%nwb > 0) then
         self%nbdy=self%nbdy+self%nwb
         allocate(self%wi(self%nwb),stat=stat)
         allocate(self%wfj,self%wlj,mold=self%wi,stat=stat)
      end if
   end if
   if (present(nnb)) then
      self%nnb=nnb
      if (self%nnb > 0) then
         self%nbdy=self%nbdy+self%nnb
         allocate(self%nj(self%nnb),stat=stat)
         allocate(self%nfi,self%nli,mold=self%nj,stat=stat)
      end if
   end if
   if (present(neb)) then
      self%neb=neb
      if (self%neb > 0) then
         self%nbdy=self%nbdy+self%neb
         allocate(self%ei(self%neb),stat=stat)
         allocate(self%efj,self%elj,mold=self%ei,stat=stat)
      end if
   end if
   if (present(nsb)) then
      self%nsb=nsb
      if (self%nsb > 0) then
         self%nbdy=self%nbdy+self%nsb
         allocate(self%sj(self%nsb),stat=stat)
         allocate(self%sfi,self%sli,mold=self%sj,stat=stat)
      end if
   end if
   if (self%nbdy > 0) then
      allocate(self%bdy_2d_type(self%nbdy),stat=stat)
      allocate(self%bdy_3d_type(self%nbdy),stat=stat)
   end if

!  set local and global grid extend
!KB   self%T%ill=self%T%l(1); self%T%ihl=self%T%u(1)
!KB   self%T%jll=self%T%l(2); self%T%jhl=self%T%u(2)

   self%domain_ready = self%T%grid_ready .and. self%U%grid_ready .and. &
                       self%V%grid_ready .and. self%X%grid_ready
END SUBROUTINE domain_configure

!-----------------------------------------------------------------------------

SUBROUTINE domain_initialize(self,runtype)

   !! Configure the type_grid

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_domain), intent(inout) :: self
   integer, intent(in) :: runtype

!  Local constants

!  Local variables
!-----------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('domain_initialize()',level=1)

   call boundary_bookkeeping(self)
   call self%mirror_bdys(self%T,self%T%H)
   call self%metrics()
   call self%uvx_depths()
   call self%register(runtype)
   where (self%T%mask > 0)
      self%T%z = 0._real64
      self%T%zo = 0._real64
      self%T%zin = 0._real64
      self%T%zio = 0._real64
   end where
   where (self%U%mask > 0)
      self%U%z = 0._real64
      self%U%zo = 0._real64
      self%U%zin = 0._real64
      self%U%zio = 0._real64
   end where
   where (self%V%mask > 0)
      self%V%z = 0._real64
      self%V%zo = 0._real64
      self%V%zin = 0._real64
      self%V%zio = 0._real64
   end where
   call self%cfl_check()
   if (self%T%kmax > 1) then
      call self%init_vertical()
      call self%do_vertical(1._real64) ! KB
   end if

   self%domain_ready = .true.

!>  @todo
!>  various modifications - boundaries, mask and bathymetry
!>  part_domain
!>  @endtodo

END SUBROUTINE domain_initialize

!---------------------------------------------------------------------------

SUBROUTINE boundary_bookkeeping(self)
   !! open boundary  book keeping

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_domain), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: stat
   integer :: i,j,il,jl,k,l,n
!-----------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('boundary_bookkeeping()',level=2)
   self%nbdyp=0
   do n=1,self%nwb
      self%nbdyp=self%nbdyp+self%wlj(n)-self%wfj(n)+1
   end do
   do n=1,self%nnb
      self%nbdyp=self%nbdyp+self%nli(n)-self%nfi(n)+1
   end do
   do n=1,self%neb
      self%nbdyp=self%nbdyp+self%elj(n)-self%efj(n)+1
   end do
   do n=1,self%nsb
      self%nbdyp=self%nbdyp+self%sli(n)-self%sfi(n)+1
   end do
   if (self%nbdy > 0) then
      allocate(self%bdy_index(self%nbdyp),stat=stat)
      allocate(self%bdy_map(self%nbdyp,2),stat=stat)
   end if
   k=1
   l=1
   do n=1,self%nwb
      self%bdy_index(l) = k
      l=l+1
      i=self%wi(n)
      do j=self%wfj(n),self%wlj(n)
         il=i-self%T%ioff; jl=j-self%T%joff
         self%bdy_map(k,1)=i
         self%bdy_map(k,2)=j
         self%T%mask(il,jl)=2
         k = k+1
      end do
   end do
   do n=1,self%nnb
      self%bdy_index(l)=k
      l=l+1
      j=self%nj(n)
      do i=self%nfi(n),self%nli(n)
         il=i-self%T%ioff; jl=j-self%T%joff
         self%bdy_map(k,1)=i
         self%bdy_map(k,2)=j
         self%T%mask(il,jl)=2
         k = k+1
      end do
   end do
   do n=1,self%neb
      self%bdy_index(l)=k
      l=l+1
      i=self%ei(n)
      do j=self%efj(n),self%elj(n)
         il=i-self%T%ioff; jl=j-self%T%joff
         self%bdy_map(k,1)=i
         self%bdy_map(k,2)=j
         self%T%mask(il,jl)=2
         k=k+1
      end do
   end do
   do n=1,self%nsb
      self%bdy_index(l)=k
      l=l+1
      j=self%sj(n)
      do i=self%sfi(n),self%sli(n)
         il=i-self%T%ioff; jl=j-self%T%joff
         self%bdy_map(k,1)=i
         self%bdy_map(k,2)=j
         self%T%mask(il,jl)=2
         k=k+1
      end do
   end do
END SUBROUTINE boundary_bookkeeping

!-----------------------------------------------------------------------------

SUBROUTINE deallocate_variables(self)
   !! Allocate all domain related variables

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_domain), intent(inout) :: self

!  Local constants

!  Local variables
!-----------------------------------------------------------------------------
#ifndef _STATIC_
   if (associated(self%logs)) call self%logs%info('deallocate_variables()',level=3)
   call deallocate_grid_variables(self%T)
   call deallocate_grid_variables(self%U)
   call deallocate_grid_variables(self%V)
   call deallocate_grid_variables(self%X)
#endif
END SUBROUTINE deallocate_variables

!-----------------------------------------------------------------------------

SUBROUTINE domain_report(self)

   !! Report the configured domain

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_domain), intent(inout) :: self

!  Local variables
   integer :: gridunit
!-----------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('domain_report()',level=1)

   if (.not. self%domain_ready) then
      write(*,*) 'domain is not - fully - configured'
   else
      open(newunit=gridunit,file='tgrid.dat')
      call self%T%report(self%logs,gridunit,'T-grid info: ')
      close(gridunit)
      open(newunit=gridunit,file='ugrid.dat')
      call self%U%report(self%logs,gridunit,'U-grid info: ')
      close(gridunit)
      open(newunit=gridunit,file='vgrid.dat')
      call self%V%report(self%logs,gridunit,'V-grid info: ')
      close(gridunit)
      open(newunit=gridunit,file='xgrid.dat')
      call self%X%report(self%logs,gridunit,'X-grid info: ')
      close(gridunit)
   end if
END SUBROUTINE domain_report

!---------------------------------------------------------------------------

SUBROUTINE start_3d(self)
   !! Prepare for 3D calculations

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_domain), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: i,j
!-----------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('start_3d()',level=2)
   TGrid: associate( TG => self%T )
   TG%zio=TG%zin
   TG%zin=TG%z

   UGrid: associate( UG => self%U )
   UG%zio=UG%zin
   do j=UG%l(2),UG%u(2)
      do i=UG%l(1),UG%u(1)-1
         if (UG%mask(i,j) > 0) then
            UG%zin(i,j)=0.25_real64*(TG%zio(i,j)+TG%zio(i+1,j)+TG%zin(i,j)+TG%zin(i+1,j))
            UG%zin(i,j)=max(UG%zin(i,j),-UG%H(i,j)+self%Dmin)
         end if
      end do
   end do
   UG%ho=UG%hn
   end associate UGrid

   VGrid: associate( VG => self%V )
   VG%zio=VG%zin
   do j=VG%l(2),VG%u(2)-1
      do i=VG%l(1),VG%u(1)
         if (VG%mask(i,j) > 0) then
            VG%zin(i,j)=0.25_real64*(TG%zio(i,j)+TG%zio(i,j+1)+TG%zin(i,j)+TG%zin(i,j+1))
            VG%zin(i,j)=max(VG%zin(i,j),-VG%H(i,j)+self%Dmin)
         end if
      end do
   end do
   self%V%ho=self%V%hn
   end associate VGrid
   end associate TGrid
END SUBROUTINE start_3d

!---------------------------------------------------------------------------

SUBROUTINE domain_cleanup(self)
   !! Clean-up the domain module

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_domain), intent(inout) :: self

!  Local constants

!  Local variables
!-----------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('domain_cleanup()',level=2)
   call deallocate_variables(self)
END SUBROUTINE domain_cleanup

!---------------------------------------------------------------------------

END MODULE getm_domain
