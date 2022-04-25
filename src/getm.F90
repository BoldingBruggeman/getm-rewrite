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
!KB   use getm_parallel
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
!KB      TYPE(type_getm_parallel) :: parallel
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
      integer :: imin,imax,jmin,jmax,kmin,kmax
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
   integer :: adv_scheme=1
!---------------------------------------------------------------------------

!write(*,*) 'This file was compiled by ', compiler_version(), ' using the options ', compiler_options()
   self%logs%prepend = ''
   call self%logs%info('getm_settings()')

   if (command_argument_count() .ne. 1 ) then
      write(*,*) 'ERROR, must be called like:'
      STOP ' getm_exe <name of bathymetry file>'
   end if
   call get_command_argument(1,self%bathymetry%depth%f)
   self%bathymetry%depth%v = 'bathymetry'

! defines to support different box_* simulations
!#define OPEN_BDY
#define WITH_SURFACE_FORCING

   self%runtype = 4
   self%info_frequency = 10
   self%mode_split = 100000
   self%mode_split = 20
   self%timestep = 10._real64
   self%sim_start = strptime("2020-01-01 00:00:00", time_format)
   self%sim_stop  = strptime("2020-01-02 00:00:00", time_format)
   self%sim_stop  = strptime("2020-01-02 01:00:00", time_format)
   self%sim_stop  = strptime("2020-01-01 01:00:00", time_format)
   self%sim_stop  = strptime("2020-01-01 00:10:00", time_format)

   self%imin = 1
   self%imax = 100
   self%jmin = 1
   self%jmax = 30
   self%kmin = 1
   self%kmax = 20

   self%domain%domain_type = 1
   self%domain%method_vertical_coordinates = 1
   self%domain%Dmin = 0.5_real64
   self%domain%Dcrit = 2.0_real64
   self%domain%ddu = 1.0_real64
   self%domain%ddl = 1.0_real64
   self%domain%ddu = 0.0_real64
   self%domain%ddl = 0.0_real64
   self%domain%lat0 = 45._real64
   self%domain%lat0 = 0._real64

#ifdef WITH_SURFACE_FORCING
   self%airsea%taux0 = 0.00_real64
   self%airsea%tauy0 = 0.01_real64
   self%airsea%swr0 =  25._real64
   self%airsea%shf0 = -25._real64
#else
   self%airsea%taux0 = 0.0_real64
   self%airsea%tauy0 = 0.0_real64
#endif

adv_scheme=0
   self%physics%mixing%method_mixing = 2 !KB constant
   self%physics%mixing%num0 = 1.e-4_real64
   self%physics%mixing%nuh0 = 1.e-4_real64
   self%physics%temperature%advection_scheme = adv_scheme
   self%physics%salinity%advection_scheme = adv_scheme

   self%dynamics%pressure%method_internal_pressure = 1
   self%dynamics%momentum%Am0 = 1.e-4_real64
self%dynamics%momentum%apply_bottom_friction = .false.
   self%dynamics%momentum%advection_scheme = adv_scheme
self%dynamics%momentum%apply_diffusion = .false.
   self%dynamics%momentum%store_advection = .true.
   self%dynamics%momentum%store_diffusion = .true.
   self%dynamics%momentum%store_damping = .true.
   self%dynamics%momentum%store_slowterms = .true.

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
!-----------------------------------------------------------------------------
  self%logs%prepend = ''
  call self%logs%info('getm_configure()',level=0)

#ifdef _STATIC_
   call self%logs%debug('_STATIC_ compilation')
#else
   call self%logs%debug('_DYNAMIC_ compilation')
#endif

#ifdef OPEN_BDY
!KB   call self%domain%configure(self%imin,self%imax,self%jmin,self%jmax,self%kmin,self%kmax,nwb=1,nnb=0,neb=1,nsb=0,logs=self%logs,fm=self%fm)
   call self%domain%configure(self%imin,self%imax,self%jmin,self%jmax,self%kmin,self%kmax,nwb=1,nnb=0,neb=1,nsb=0,logs=self%logs,fm=self%fm)
#else
   call self%domain%configure(self%imin,self%imax,self%jmin,self%jmax,self%kmin,self%kmax,nwb=0,nnb=0,neb=0,nsb=0,logs=self%logs,fm=self%fm)
#endif

!KB   call self%parallel%configure(self%domain%T)

   call self%bathymetry%initialize(self%logs,self%domain%T,self%domain%domain_type)
#ifdef OPEN_BDY
   boundaries: block
   integer :: n=0
   integer :: bdytype=3
   if (self%domain%nwb > 0) then
      self%domain%wi(1)=3;  self%domain%wfj(1)=3; self%domain%wlj(1)=28
      n=n+1
      self%domain%bdy_2d_type(n)=bdytype
   end if
   if (self%domain%nnb > 0) then
      self%domain%nj(1)=28;  self%domain%nfi(1)=3; self%domain%nli(1)=98
      n=n+1
      self%domain%bdy_2d_type(n)=bdytype
   end if
   if (self%domain%neb > 0) then
      self%domain%ei(1)=98; self%domain%efj(1)=3; self%domain%elj(1)=28
      n=n+1
      self%domain%bdy_2d_type(n)=bdytype
   end if
   if (self%domain%nsb > 0) then
      self%domain%sj(1)=3;  self%domain%sfi(1)=3; self%domain%sli(1)=98
      n=n+1
      self%domain%bdy_2d_type(n)=bdytype
   end if
   end block boundaries
#endif

   ! now we can configure the field_manager
   !
   length = self%imax-self%imin+1
   call self%fm%register_dimension('x',length, &
                                   global_length=length,offset=0,newid=self%domain%id_dim_x)
   call self%fm%register_dimension('xi'//'i',length, &
                                   global_length=length,offset=0,newid=self%domain%id_dim_xi)

   length = self%jmax-self%jmin+1
   call self%fm%register_dimension('y',length, &
                                   global_length=length,offset=0,newid=self%domain%id_dim_y)
   call self%fm%register_dimension('yi'//'i',length, &
                                   global_length=length,offset=0,newid=self%domain%id_dim_yi)

   call self%fm%register_dimension('z',self%kmax,global_length=self%kmax,offset=0,newid=self%domain%id_dim_z)
   call self%fm%register_dimension('zi',self%kmax+1,global_length=self%kmax+1,offset=-1,newid=self%domain%id_dim_zi)
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
   call self%output%configure(self%logs,self%fm)
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
   call self%domain%initialize(self%runtype)
   call self%domain%report()
   call self%vertical_diffusion%initialize(self%domain%T)
   call self%airsea%initialize(self%domain)
   call self%physics%initialize(self%runtype,self%domain,self%advection,self%vertical_diffusion)
   call self%dynamics%initialize(self%runtype,self%domain,self%advection,self%vertical_diffusion)
   call self%domain%depth_update()
   call self%output%initialize()
!KB   call self%fm%list()
END SUBROUTINE getm_initialize

!---------------------------------------------------------------------------

SUBROUTINE getm_integrate(self)

   !! Configure all component of the model

   use gsw_mod_toolbox, only: gsw_pt_from_CT
   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_model) :: self

! Local constants

! Local variables
  integer :: n
  TYPE(datetime) :: sim_time
  TYPE(timedelta) :: dt, remain
  real(real64) :: dte,dti
  integer :: seconds, milliseconds
!-----------------------------------------------------------------------------
   xOutput: associate( output => self%output )
   xDomain: associate( domain => self%domain )
   xAirsea: associate( airsea => self%airsea )
   xPHYSICS: associate(mixing=>self%physics%mixing, salinity=>self%physics%salinity, &
                       temperature=>self%physics%temperature, radiation=>self%physics%radiation, &
                       density=>self%physics%density)
   xDYNAMICS: associate(momentum=>self%dynamics%momentum, pressure=>self%dynamics%pressure, &
                        sealevel=>self%dynamics%sealevel)

   if (self%sim_start .gt. self%sim_stop) return
   call self%logs%info('getm_integrate()',level=0)

   n = 0
   sim_time = self%sim_start
   seconds = int(self%timestep)
   milliseconds = int(1000*(self%timestep-seconds))
   dt = timedelta(seconds = nint(self%timestep), milliseconds = 0)
   remain = self%sim_stop - sim_time

   dte=self%timestep
   dti = self%mode_split*self%timestep
   call self%output%do_output(self%sim_start)
   do while (remain%total_seconds() .gt. 0._real64)
      n = n+1
      sim_time = sim_time + dt
      self%logs%global_info_silence = mod(n,self%info_frequency) .ne. 0
      call self%logs%info(sim_time%isoformat(),level=1)

      call self%output%prepare_output(sim_time,n)

      ! call do_input(n)
      call airsea%update(n)

      ! 2D barotropic
      call pressure%surface(domain%T%z,airsea%sp)
      call momentum%uv_momentum_2d(self%runtype,dte,airsea%taux,airsea%tauy,pressure%dpdx,pressure%dpdy)
      call momentum%velocities_2d()
#ifdef OPEN_BDY
      sealevel%zbdy=0.001
#endif
      call sealevel%boundaries(dte,momentum%U,momentum%V,momentum%bdyu,momentum%bdyv)
      call sealevel%update(dte,momentum%U,momentum%V)
      call domain%depth_update()

      ! 3D barotropic
      if (self%runtype > 1 .and. mod(n,self%mode_split) == 0) then
         momentum%Ui=momentum%Ui/self%mode_split
         momentum%Vi=momentum%Vi/self%mode_split

         call self%domain%start_3d()
         call domain%do_vertical(dti)
         call pressure%surface(domain%T%zio,airsea%sp)
         if (self%runtype > 3) call pressure%internal(density%buoy,momentum%SxB,momentum%SyB)
         call momentum%uvw_momentum_3d(dti,airsea%taux,airsea%tauy,pressure%dpdx,pressure%dpdy, &
                                       pressure%idpdx,pressure%idpdy,mixing%num)
         call momentum%stresses(airsea%taux,airsea%tauy)

         call mixing%calculate(dti,momentum%ustar2_s,momentum%ustar2_b,momentum%SS,density%NN)

         ! 3D baroclinic
         if (self%runtype > 3) then
            call radiation%calculate(airsea%swr,airsea%albedo)
            call temperature%calculate(dti,momentum%pk,momentum%qk,mixing%nuh,radiation%rad,airsea%shf)
            call salinity%calculate(dti,momentum%pk,momentum%qk,mixing%nuh)
            temperature%sst=gsw_pt_from_CT(salinity%S(:,:,domain%T%kmax),temperature%T(:,:,domain%T%kmax))
            call density%density(salinity%S,temperature%T)
            call density%buoyancy()
            call density%brunt_vaisala(salinity%S,temperature%T)
         end if

         ! MUST be the last routines to call in the 3D loop
         ! update the slow terms with contribution from 3D advection, diffusion, buoyancy and friction
         call momentum%slow_momentum_terms(dte)
      end if

      call self%output%do_output(sim_time)
      remain = self%sim_stop - sim_time
   end do
   end associate xDYNAMICS
   end associate xPHYSICS
   end associate xAirsea
   end associate xDomain
   end associate xOutput
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
   call self%domain%cleanup()
END SUBROUTINE getm_finalize

!---------------------------------------------------------------------------

END MODULE getm_model
