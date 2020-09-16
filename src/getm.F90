! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!> In this module, the temperature equation is processed by
!> reading in the namelist {\tt temp} and initialising the temperature field
!> (this is done in {\tt init\_temperature}),
!> and calculating the advection-diffusion-equation, which includes
!> penetrating short-wave radiation as source term (see {\tt do\_temperature}).
!> |url|

!> @note
!>    1. Use __mold__ to allocate arrays
!> @endnote

MODULE getm_model

   !! Description:
   !!   < Say what this module contains >
   !!
   !! Current Code Owner: < Name of person responsible for this code >
   !!
   !! Code Description:
   !!   Language: Fortran 90.
   !!   This code is written to JULES coding standards v1.
   !! [egon](|url|/page/science/getm_physics.html)

!   use iso_c_binding
!   use iso_fortran_env
   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use datetime_module, only: datetime, timedelta, clock, strptime
   use logging
   use field_manager
   use getm_domain
   use getm_operators
   use getm_airsea
   use getm_physics
   use getm_dynamics
   use getm_input
   use getm_bathymetry
   use getm_output

   IMPLICIT NONE

   PRIVATE  ! Private scope by default

!  Module constants

!  Module types and variables
   type, public :: type_getm_model
      !! author: Karsten Bolding
      !! version: v0.1
      !!
      !! GETM model type

      character(len=256) :: name = "GETM"
      TYPE(datetime) :: sim_start, sim_stop
      TYPE(type_logging) :: logs
      TYPE(type_field_manager) :: fm
      TYPE(type_getm_domain) :: domain
      TYPE(type_vertical_diffusion) :: vertical_diffusion
      TYPE(type_advection) :: advection
      TYPE(type_getm_airsea) :: airsea
      TYPE(type_bathymetry) :: bathymetry
      TYPE(type_getm_physics) :: physics
      TYPE(type_getm_dynamics) :: dynamics
      TYPE(type_getm_output) :: output

#ifdef _STATIC_
      integer, parameter :: imin=_IMIN_, imax=_IMAX_
      integer, parameter :: jmin=_JMIN_, jmax=_JMAX_
      integer, parameter :: kmin=_KMIN_, kmax=_KMAX_
#else
#ifdef _PARALLEL_
      integer, parameter :: imin=1,imax=3,jmin=1,jmax=6,kmin=0,kmax=25,halo=2
#endif
#endif
      integer :: runtype
      integer :: mode_split
      integer :: info_frequency
      real(real64) :: timestep

      contains

      procedure :: settings => getm_settings
      procedure :: configure => getm_configure
      procedure :: initialize => getm_initialize
      procedure :: integrate => getm_integrate
      procedure :: finalize => getm_finalize

   end type type_getm_model

CONTAINS

!---------------------------------------------------------------------------

SUBROUTINE getm_settings(self)

   !! Configure all component of the model

   IMPLICIT NONE

! Subroutine arguments
   class(type_getm_model) :: self

!  Local constants
   character(len=256), parameter :: time_format='%Y-%m-%d %H:%M:%S'

!  Local variables
   TYPE(timedelta) :: timediff
!---------------------------------------------------------------------------

!write(*,*) 'This file was compiled by ', compiler_version(), ' using the options ', compiler_options()
   self%logs%prepend = ''
   call self%logs%info('getm_settings()')

#if 0
   self%bathymetry%depth%f = '/data/kb/getm-setups/box_spherical/topo.nc'
   self%bathymetry%depth%f = '/data/kb/getm-setups/box_cartesian/topo.nc'
   self%bathymetry%depth%f = '/data/kb/getm-setups/sylt/topo.nc'
   self%bathymetry%depth%f = '/data/kb/getm-setups/seamount/topo.nc'
   self%bathymetry%depth%f = '/data/kb/getm-setups/sylt/topo.nc'
   self%bathymetry%depth%f = '/data/kb/getm-setups/NorthSea/topo.nc'
#else
   if (command_argument_count() .ne. 1 ) then
      write(*,*) 'ERROR, must be called like:'
      STOP ' getm_exe <name of bathymetry file>'
   end if
   call get_command_argument(1,self%bathymetry%depth%f)
   self%bathymetry%depth%v = 'bathymetry'
#endif
   self%runtype = 1
   self%info_frequency = 3
   self%mode_split = 20
   self%mode_split = 100000
   self%timestep = 3600._real64
   self%timestep = 10._real64
   self%sim_start = strptime("2020-01-01 00:00:00", time_format)
   self%sim_stop = strptime("2020-01-01 12:00:00", time_format)
   self%sim_stop = strptime("2020-01-01 00:02:00", time_format)

   self%domain%ddl = 1.0_real64
   self%domain%ddu = 2.0_real64
   self%domain%Dmin = 0.5_real64

   self%airsea%taux0 = 0.001_real64
   self%airsea%tauy0 = 0.001_real64

   return
END SUBROUTINE getm_settings

!---------------------------------------------------------------------------

SUBROUTINE getm_configure(self)

   !! Configure all component of the model

   IMPLICIT NONE

! Subroutine arguments
   class(type_getm_model) :: self

! Local constants

! Local variables
   integer :: length
   integer :: kmin=1,kmax=20
!-----------------------------------------------------------------------------
   self%logs%prepend = ''
   call self%logs%info('getm_configure()')

   call self%domain%T%configure(self%logs,kmin=kmin,kmax=kmax,halo=(/0, 0, 0/))

#ifdef _STATIC_
   call self%logs%debug('_STATIC_ compilation')
   call self%domain%configure()
#else
   call self%logs%debug('_DYNAMIC_ compilation')
#ifdef _PARALLEL_
   call self%domain%configure(self%logs,imin=imin,imax=imax,jmin=jmin,jmax=jmax)
#else
   call self%bathymetry%initialize(self%logs,self%domain%T,self%domain%domain_type)
   call self%domain%configure(self%logs,self%fm)
#endif
#endif

   ! now we can configure the field_manager
   !
   length = size(self%bathymetry%pgrid%c1)
   call self%fm%register_dimension(trim(self%bathymetry%c1%v),length, &
                                   global_length=length,offset=0,newid=self%domain%id_dim_x)
   call self%fm%register_dimension(trim(self%bathymetry%c1%v)//'i',length, &
                                   global_length=length,offset=0,newid=self%domain%id_dim_xi)
   length = size(self%bathymetry%pgrid%c2)
   call self%fm%register_dimension(trim(self%bathymetry%c2%v),length, &
                                   global_length=length,offset=0,newid=self%domain%id_dim_y)
   call self%fm%register_dimension(trim(self%bathymetry%c2%v)//'i',length, &
                                   global_length=length,offset=0,newid=self%domain%id_dim_yi)
   call self%fm%register_dimension('z',kmax,global_length=kmax,offset=0,newid=self%domain%id_dim_z)
   call self%fm%register_dimension('zi',kmax+1,global_length=kmax+1,offset=-1,newid=self%domain%id_dim_zi)
   call self%fm%register_dimension('time',id=id_dim_time)

   self%domain%T%dim_2d_ids = (/ self%domain%id_dim_x, self%domain%id_dim_y /)
   self%domain%T%dim_3d_ids = (/ self%domain%id_dim_x, self%domain%id_dim_y, self%domain%id_dim_z /)
   self%domain%U%dim_2d_ids = (/ self%domain%id_dim_xi, self%domain%id_dim_y /)
   self%domain%U%dim_3d_ids = (/ self%domain%id_dim_xi, self%domain%id_dim_y, self%domain%id_dim_z /)
   self%domain%V%dim_2d_ids = (/ self%domain%id_dim_x, self%domain%id_dim_yi /)
   self%domain%V%dim_3d_ids = (/ self%domain%id_dim_x, self%domain%id_dim_yi, self%domain%id_dim_z /)
   call self%fm%initialize(append_by_default=(/id_dim_time/))

   call self%airsea%configure(self%logs,self%fm)
   call self%physics%configure(self%logs,self%fm)
   call self%dynamics%configure(self%logs,self%fm)
   call self%output%configure(self%logs)

   call self%logs%info('done')
   return
END SUBROUTINE getm_configure

!---------------------------------------------------------------------------

SUBROUTINE getm_initialize(self)

   !! Initialize all components of the model

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_model) :: self

! Local constants

! Local variables

!-----------------------------------------------------------------------------
   call self%logs%info('getm_initialize()')
   call self%domain%initialize()
   call self%domain%report()
   call self%airsea%initialize(self%domain)
   call self%physics%initialize(self%domain,self%advection,self%vertical_diffusion)
   call self%dynamics%initialize(self%domain)
   call self%domain%depth_update()
   call self%output%initialize(self%fm)
   call self%logs%info('done')
!KB   call self%fm%list()
   return
END SUBROUTINE getm_initialize

!---------------------------------------------------------------------------

SUBROUTINE getm_integrate(self)

   !! Configure all component of the model

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_model) :: self

! Local constants

! Local variables
  integer :: n
  TYPE(datetime) :: sim_time
  TYPE(timedelta) :: dt, remain
  integer :: seconds, milliseconds

!KB
  real(real64), allocatable :: nuh(:,:,:)

!-----------------------------------------------------------------------------
!KB   Momentum: associate( momentum => self%dynamics%momentum )
   if (self%sim_start .gt. self%sim_stop) return
   call self%logs%info('getm_integrate()',level=0)

   n = 0
   sim_time = self%sim_start
   seconds = int(self%timestep)
   milliseconds = int(1000*(self%timestep-seconds))
   dt = timedelta(seconds = nint(self%timestep), milliseconds = 0)
   remain = self%sim_stop - sim_time

   do while (remain%total_seconds() .gt. 0._real64)
!   do while (remain%total_seconds() .gt. 0._real64 .and. n .lt. 50)
      n = n+1
      sim_time = sim_time + dt
      self%logs%global_info_silence = mod(n,self%info_frequency) .ne. 0
      call self%logs%info(sim_time%isoformat(),level=1)
      call self%output%prepare_output(sim_time,n)

      ! call do_input(n)
      call self%airsea%update(n)
      ! total surface stress at T-points
      ! taus(i,j)=rho0i*sqrt( tausx(i,j)**2 + tausy(i,j)**2 )

      ! This is the GETM 2D call order - not carved in stone
      ! call uv_advect(U,V,DU,DV)
      ! uv_diffusion(An_method,U,V,D,DU,DV) ! Has to be called after uv_advect.
      ! call momentum(loop,tausx,tausy,airp)
      ! Ui and Vi are updated in momentum_2d.F90
      ! if (have_boundaries) call update_2d_bdy(loop,bdy2d_ramp)
      ! call sealevel(loop)
      ! call depth_update()

      ! diffusion/advectio is still missing

      call self%dynamics%momentum%advection_2d(self%timestep)
!KB      call momentum%advection_2d(self%timestep)
      call self%dynamics%pressure%surface(self%domain%T%z,self%airsea%sp)
      call self%dynamics%momentum%do_2d(self%timestep,self%dynamics%pressure%dpdx,self%airsea%taux, &
                                                      self%dynamics%pressure%dpdy,self%airsea%tauy)
      call self%dynamics%sealevel%update(self%timestep,self%dynamics%momentum%U,self%dynamics%momentum%V)
      call self%domain%depth_update()

      if (self%runtype > 1 .and. mod(n,self%mode_split) == 0) then ! 3D calculations
         ! This is the GETM 3D call order - not carved in stone
         ! call start_macro()
         ! huo=hun; hvo=hvn
         ! call structure_friction_3d
         ! if (ufirst) then
         !    call uu_momentum_3d(n,bdy3d)
         !    call vv_momentum_3d(n,bdy3d)
         !    ufirst=.false.
         ! else
         !    call vv_momentum_3d(n,bdy3d)
         !    call uu_momentum_3d(n,bdy3d)
         !    ufirst=.true.
         ! end if
         ! call coordinates(.false.)
         ! call ww_momentum_3d()
         !!!! call uv_advect_3d()
         !!!! call uv_diffusion_3d()  ! Must be called after uv_advect_3d
         ! call stresses_3d()
         !!!! call gotm()
         !!!! if (calc_temp) call do_temperature(n)
         !!!! if (calc_salt) call do_salinity(n)
         ! call do_eqstate()
         ! call slow_bottom_friction()
         !!!! call uv_advect(Uint,Vint,Dun,Dvn)
         !!!! call uv_diffusion(0,Uint,Vint,Dn,Dun,Dvn) ! Has to be called after uv_advect.
         ! call slow_terms()
         ! call stop_macro()

         call self%domain%start_3d()
! moved to momentum_3d
!KB         self%dynamics%momentum%Uint=self%dynamics%momentum%Uint/self%mode_split
!KB         self%dynamics%momentum%Vint=self%dynamics%momentum%Vint/self%mode_split
         call self%dynamics%pressure%surface(self%domain%T%sseo,self%airsea%sp)
         call self%dynamics%momentum%do_3d(self%mode_split)
         call self%dynamics%momentum%vel_3d()
         call self%domain%do_vertical()
!KB         call self%dynamics%momentum%do_w(self%mode_split*self%timestep)
         call self%dynamics%momentum%stresses()

         if (self%runtype > 3) then
            xSalinity: associate( salinity => self%physics%salinity )
            xTemperature: associate( temperature => self%physics%temperature )
            xDensity: associate( density => self%physics%density )
            xMomentum: associate( momentum => self%dynamics%momentum )

            call salinity%calculate(self%timestep,momentum%pk,momentum%qk,nuh)
            call temperature%calculate(self%timestep,momentum%pk,momentum%qk,nuh)
            call density%density(salinity%S,temperature%T)
            call density%buoyancy()

            end associate xMomentum
            end associate xDensity
            end associate xTemperature
            end associate xSalinity
         end if

         call self%dynamics%momentum%slow_bottom_friction()
         !
         !
         call self%dynamics%momentum%slow_terms(self%dynamics%pressure%idpdx,self%dynamics%pressure%idpdy)
         ! reset variables
         self%dynamics%momentum%Uio=self%dynamics%momentum%Ui; self%dynamics%momentum%Ui=0._real64
         self%dynamics%momentum%Vio=self%dynamics%momentum%Vi; self%dynamics%momentum%Vi=0._real64
      end if

      call self%output%do_output(sim_time)
      remain = self%sim_stop - sim_time
   end do
   call self%logs%info('done',level=0)
!KB   end associate Momentum
   return
END SUBROUTINE getm_integrate

!---------------------------------------------------------------------------

SUBROUTINE getm_finalize(self)

   !! Finalize all components of the model

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_model) :: self

! Local constants

! Local variables
!-----------------------------------------------------------------------------
   call self%logs%info('getm_finalize()')
   call self%logs%info('done')
   return
END SUBROUTINE getm_finalize

!---------------------------------------------------------------------------

END MODULE getm_model
