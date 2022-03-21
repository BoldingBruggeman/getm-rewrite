! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!> In this module, the temperature equation is processed by
!> reading in the namelist {\tt temp} and initialising the temperature field
!> (this is done in {\tt init\_temperature}),
!> and calculating the advection-diffusion-equation, which includes
!> penetrating short-wave radiation as source term (see {\tt do\_temperature}).

!> @note
!> loop boundaries in uuadvgrid, uvadvgrid, vvadvgrid and vvadvgrid
!> @nednote

MODULE getm_momentum

   !! Description:
   !!   < Say what this module contains >
   !!
   !! Current Code Owner: < Name of person responsible for this code >
   !!
   !! Code Description:
   !!   Language: Fortran 90.
   !!   This code is written to JULES coding standards v1.

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use memory_manager
   use logging
   use field_manager
   use getm_domain
   use getm_operators

   IMPLICIT NONE

   PRIVATE  ! Private scope by default

!  Module constants
   real(real64), parameter :: rho0 = 1025._real64
      !! Reference density
   real(real64), parameter :: kappa = 0.4_real64
      !! constant
   real(real64), parameter :: avmmol = 0.001_real64 !KB
      !! constant

!  Module types and variables
   type, public :: type_getm_momentum

      class(type_logging), pointer :: logs => null()
      class(type_field_manager), pointer :: fm => null()
      class(type_getm_domain), pointer :: domain
      class(type_advection), pointer :: advection => null()
      class(type_vertical_diffusion), pointer :: vertical_diffusion => null()

#ifdef _STATIC_
      real(real64), dimension(E2DFIELD) :: U, V
      real(real64), dimension(I3DFIELD) :: pk, qk, w
#else
      real(real64), dimension(:,:), allocatable :: U,V
      real(real64), dimension(:,:), allocatable :: Ui,Vi
      real(real64), dimension(:,:), allocatable :: Uio,Vio
      real(real64), dimension(:,:), allocatable :: Ua,Va
      real(real64), dimension(:,:), allocatable :: fU,fV
      real(real64), dimension(:,:), allocatable :: advU,advV
      real(real64), dimension(:,:), allocatable :: diffu1,diffv1
      real(real64), dimension(:,:), allocatable :: dampU,dampV
      real(real64), dimension(:,:), allocatable :: u1,v1
      real(real64), dimension(:,:), allocatable :: SxA,SyA ! Slow advection
      real(real64), dimension(:,:), allocatable :: SxB,SyB ! Slow internal pressure
      real(real64), dimension(:,:), allocatable :: SxD,SyD ! Slow diffusion
      real(real64), dimension(:,:), allocatable :: SxF,SyF ! Slow friction
      real(real64), dimension(:,:), allocatable :: rru,rrv
      real(real64), dimension(:,:), allocatable :: ru,rv
      real(real64), dimension(:,:), allocatable :: Am,An
      real(real64), dimension(:,:,:), allocatable :: pk,qk,ww
      real(real64), dimension(:,:,:), allocatable :: pka,qka
      real(real64), dimension(:,:,:), allocatable :: fpk,fqk
      real(real64), dimension(:,:,:), allocatable :: advpk,advqk
      real(real64), dimension(:,:,:), allocatable :: diffuk,diffvk
      real(real64), dimension(:,:,:), allocatable :: uk,vk
      real(real64), dimension(:,:), allocatable :: taus,taub
      real(real64), dimension(:,:), allocatable :: taubx,tauby
      real(real64), dimension(:,:,:), allocatable :: SS
      real(real64), dimension(:), allocatable :: bdyu,bdyv
      ! help variables
      real(real64), dimension(:,:,:), allocatable :: num,ea2,ea4
      real(real64), dimension(:,:), allocatable :: work2d
      type(type_getm_grid) :: uuadvgrid,uvadvgrid,vuadvgrid,vvadvgrid
      logical :: apply_bottom_friction=.true.
      integer :: advection_scheme=1
      logical :: apply_diffusion=.true.
      integer :: coriolis_scheme=1
      real(real64) :: molecular=0._real64
      real(real64) :: cnpar=1._real64
      integer :: An_method=0
      real(real64) :: Am0=0._real64
      real(real64) :: An0=0._real64
      logical :: store_advection=.false.
      logical :: store_diffusion=.true.
      logical :: store_damping=.true.
      logical :: store_slowterms=.false.
      logical :: ufirst = .false.
#endif

      contains

      procedure :: configure => momentum_configure
      procedure :: initialize => momentum_initialize
      procedure :: register => momentum_register
      procedure :: initialize_2d => uv_initialize_2d
      procedure :: uv_momentum_2d => uv_momentum_2d
      procedure :: bottom_friction_2d => bottom_friction_2d
      procedure :: u_2d => u_2d
      procedure :: v_2d => v_2d
      procedure :: coriolis_fu => coriolis_fu
      procedure :: coriolis_fv => coriolis_fv
      procedure :: uv_advection_2d => uv_advection_2d
      procedure :: uv_diffusion_2d => uv_diffusion_2d
      procedure :: velocities_2d => velocities_2d
      procedure :: initialize_3d => uv_initialize_3d
      procedure :: uvw_momentum_3d => uvw_momentum_3d
      procedure :: bottom_friction_3d => bottom_friction_3d
      procedure :: w_momentum_3d => w_momentum_3d
      procedure :: pk_3d => pk_3d
      procedure :: qk_3d => qk_3d
      procedure :: coriolis_fpk => coriolis_fpk
      procedure :: coriolis_fqk => coriolis_fqk
      procedure :: uv_advection_3d => uv_advection_3d
      procedure :: uv_diffusion_3d => uv_diffusion_3d
      procedure :: velocities_3d => velocities_3d
      procedure :: slow_momentum_terms => slow_momentum_terms
      procedure :: slow_advection => slow_advection
      procedure :: slow_diffusion => slow_diffusion
      procedure :: slow_bottom_friction => slow_bottom_friction
      procedure :: shear_frequency => shear_frequency
      procedure :: stresses => stresses
      procedure :: list => momentum_list
      procedure :: set => momentum_set

   end type type_getm_momentum

   INTERFACE
      MODULE SUBROUTINE momentum_register(self,runtype)
         class(type_getm_momentum), intent(inout) :: self
         integer, intent(in) :: runtype
      END SUBROUTINE momentum_register

      ! 2D momentum
      MODULE SUBROUTINE uv_initialize_2d(self)
         class(type_getm_momentum), intent(inout) :: self
      END SUBROUTINE uv_initialize_2d

      MODULE SUBROUTINE uv_momentum_2d(self,runtype,dt,tausx,tausy,dpdx,dpdy)
         class(type_getm_momentum), intent(inout) :: self
         integer, intent(in) :: runtype
         real(real64), intent(in) :: dt
            !! timestep [s]
#define _U2_ self%domain%U%l(1):,self%domain%U%l(2):
         real(real64), intent(in) :: tausx(_U2_)
            !! surface stress
#undef _U2_
#define _V2_ self%domain%V%l(1):,self%domain%V%l(2):
         real(real64), intent(in) :: tausy(_V2_)
            !! surface stress
#undef _V2_
#define _U2_ self%domain%U%l(1):,self%domain%U%l(2):
         real(real64), intent(in) :: dpdx(_U2_)
            !! surface pressure gradient - including air pressure
#undef _U2_
#define _V2_ self%domain%V%l(1):,self%domain%V%l(2):
         real(real64), intent(in) :: dpdy(_V2_)
            !! surface pressure gradient - including air pressure
#undef _V2_
      END SUBROUTINE uv_momentum_2d

      MODULE SUBROUTINE u_2d(self,dt,tausx,dpdx)
         class(type_getm_momentum), intent(inout) :: self
         real(real64), intent(in) :: dt
            !! timestep [s]
#define _U2_ self%domain%T%l(1):,self%domain%T%l(2):
         real(real64), intent(in) :: tausx(_U2_)
            !! surface stresses
         real(real64), intent(in) :: dpdx(_U2_)
            !! surface pressure gradient - including air pressure
#undef _U2_
      END SUBROUTINE u_2d

      MODULE SUBROUTINE v_2d(self,dt,tausy,dpdy)
         class(type_getm_momentum), intent(inout) :: self
         real(real64), intent(in) :: dt
            !! timestep [s]
#define _V2_ self%domain%T%l(1):,self%domain%T%l(2):
         real(real64), intent(in) :: tausy(_V2_)
            !! surface stresses
         real(real64), intent(in) :: dpdy(_V2_)
            !! surface pressure gradient - including air pressure
#undef _V2_
   END SUBROUTINE v_2d

      ! 3D routines
      MODULE SUBROUTINE uv_initialize_3d(self)
         class(type_getm_momentum), intent(inout) :: self
      END SUBROUTINE uv_initialize_3d
      MODULE SUBROUTINE uvw_momentum_3d(self,dt,tausx,tausy,dpdx,dpdy,idpdx,idpdy,viscosity)
         !! Solve the 3D momemtum equations
         class(type_getm_momentum), intent(inout) :: self
         real(real64), intent(in) :: dt
            !! timestep [s]
#define _U2_ self%domain%U%l(1):,self%domain%U%l(2):
   real(real64), intent(in) :: tausx(_U2_)
     !! surface stress - x
#undef _U2_
#define _V2_ self%domain%V%l(1):,self%domain%V%l(2):
   real(real64), intent(in) :: tausy(_V2_)
     !! surface stress - y
#undef _V2_
#define _U2_ self%domain%U%l(1):,self%domain%U%l(2):
   real(real64), intent(in) :: dpdx(_U2_)
     !! surface pressure (including air pressure) - x-gradient
#undef _U2_
#define _V2_ self%domain%V%l(1):,self%domain%V%l(2):
   real(real64), intent(in) :: dpdy(_V2_)
     !! surface pressure (including air pressure) - y-gradient
#undef _V2_
#define _U3_ self%domain%U%l(1):,self%domain%U%l(2):,self%domain%U%l(3):
   real(real64), intent(in) :: idpdx(_U3_)
     !! internal pressure - x-gradient
#undef _U3_
#define _V3_ self%domain%V%l(1):,self%domain%V%l(2):,self%domain%V%l(3):
   real(real64), intent(in) :: idpdy(_V3_)
     !! internal pressure - y-gradient
#undef _V3_
#define _T3_ self%domain%T%l(1):,self%domain%T%l(2):,self%domain%T%l(3):
   real(real64), intent(in) :: viscosity(_T3_)
     !! viscosity
#undef _T3_
      END SUBROUTINE uvw_momentum_3d
      MODULE SUBROUTINE w_momentum_3d(self,dt)
         class(type_getm_momentum), intent(inout) :: self
         real(real64), intent(in) :: dt
            !! timestep [s]
   END SUBROUTINE w_momentum_3d


   MODULE SUBROUTINE pk_3d(self,dt,tausx,dpdx,idpdx,viscosity)
!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
      !! GETM momentum type
   real(real64), intent(in) :: dt
      !! timestep [s]
#define _U2_ self%domain%U%l(1):,self%domain%U%l(2):
   real(real64), intent(in) :: tausx(_U2_)
      !! surface stress in local x-direction
   real(real64), intent(in) :: dpdx(_U2_)
      !! surface pressure gradient - including air pressure
#undef _U2_
#define _U3_ self%domain%U%l(1):,self%domain%U%l(2):,self%domain%U%l(3):
   real(real64), intent(in) :: idpdx(_U3_)
      !! internal pressure gradient
   real(real64), intent(in) :: viscosity(_U3_)
      !! viscosity
#undef _U3_
      END SUBROUTINE pk_3d


   MODULE SUBROUTINE qk_3d(self,dt,tausy,dpdy,idpdy,viscosity)
!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
      !! GETM momentum type
   real(real64), intent(in) :: dt
      !! timestep [s]
#define _V2_ self%domain%V%l(1):,self%domain%V%l(2):
   real(real64), intent(in) :: tausy(_V2_)
      !! surface stress in local y-direction
   real(real64), intent(in) :: dpdy(_V2_)
      !! surface pressure gradient - including air pressure
#undef _V2_
#define _V3_ self%domain%V%l(1):,self%domain%V%l(2):,self%domain%V%l(3):
   real(real64), intent(in) :: idpdy(_V3_)
      !! internal pressure gradient
   real(real64), intent(in) :: viscosity(_V3_)
      !! viscosity
#undef _V3_
      END SUBROUTINE qk_3d

      ! friction
      MODULE SUBROUTINE bottom_friction_2d(self,runtype)
         class(type_getm_momentum), intent(inout) :: self
         integer, intent(in) :: runtype
      END SUBROUTINE bottom_friction_2d
      MODULE SUBROUTINE bottom_friction_3d(self)
         class(type_getm_momentum), intent(inout) :: self
      END SUBROUTINE bottom_friction_3d
      MODULE SUBROUTINE slow_bottom_friction(self)
         class(type_getm_momentum), intent(inout) :: self
      END SUBROUTINE slow_bottom_friction

      ! Coriolis
      MODULE SUBROUTINE coriolis_fu(self)
         class(type_getm_momentum), intent(inout) :: self
      END SUBROUTINE coriolis_fu
      MODULE SUBROUTINE coriolis_fv(self)
         class(type_getm_momentum), intent(inout) :: self
      END SUBROUTINE coriolis_fv
      MODULE SUBROUTINE coriolis_fpk(self)
         class(type_getm_momentum), intent(inout) :: self
      END SUBROUTINE coriolis_fpk
      MODULE SUBROUTINE coriolis_fqk(self)
         class(type_getm_momentum), intent(inout) :: self
      END SUBROUTINE coriolis_fqk

      ! advection
      MODULE SUBROUTINE uv_advection_2d(self,dt)
         class(type_getm_momentum), intent(inout) :: self
         real(real64), intent(in) :: dt
            !! timestep [s]
      END SUBROUTINE uv_advection_2d
      MODULE SUBROUTINE slow_advection(self,dt)
         class(type_getm_momentum), intent(inout) :: self
         real(real64), intent(in) :: dt
            !! timestep [s]
      END SUBROUTINE slow_advection
      MODULE SUBROUTINE uv_advection_3d(self,dt)
         class(type_getm_momentum), intent(inout) :: self
         real(real64), intent(in) :: dt
            !! timestep [s]
      END SUBROUTINE uv_advection_3d

      ! diffusion
      MODULE SUBROUTINE uv_diffusion_2d(self,dt)
         class(type_getm_momentum), intent(inout) :: self
         real(real64), intent(in) :: dt
            !! timestep [s]
      END SUBROUTINE uv_diffusion_2d
      MODULE SUBROUTINE slow_diffusion(self)
         class(type_getm_momentum), intent(inout) :: self
      END SUBROUTINE slow_diffusion
      MODULE SUBROUTINE uv_diffusion_3d(self,dt)
         class(type_getm_momentum), intent(inout) :: self
         real(real64), intent(in) :: dt
            !! timestep [s]
      END SUBROUTINE uv_diffusion_3d

      MODULE SUBROUTINE slow_momentum_terms(self,dt)
         class(type_getm_momentum), intent(inout) :: self
         real(real64), intent(in) :: dt
      END SUBROUTINE slow_momentum_terms

      MODULE SUBROUTINE velocities_2d(self)
         class(type_getm_momentum), intent(inout) :: self
      END SUBROUTINE velocities_2d

      MODULE SUBROUTINE velocities_3d(self)
         class(type_getm_momentum), intent(inout) :: self
      END SUBROUTINE velocities_3d

      MODULE SUBROUTINE shear_frequency(self,num)
         class(type_getm_momentum), intent(inout) :: self
#define _T3_ self%domain%T%l(1):,self%domain%T%l(2):,self%domain%T%l(3):
         real(real64), intent(in) :: num(_T3_)
#undef _T3_
      END SUBROUTINE shear_frequency

      MODULE SUBROUTINE stresses(self,tausx,tausy)
         class(type_getm_momentum), intent(inout) :: self
#define _U2_ self%domain%U%l(1):,self%domain%U%l(2):
         real(real64), intent(in) :: tausx(_U2_)
            !! surface stresses - x-direction
#undef _U2_
#define _V2_ self%domain%V%l(1):,self%domain%V%l(2):
         real(real64), intent(in) :: tausy(_V2_)
            !! surface stresses - y-direction
#undef _V2_
      END SUBROUTINE stresses

   END INTERFACE

CONTAINS

!---------------------------------------------------------------------------

SUBROUTINE momentum_configure(self,logs,fm)

   !! Configure the components belonging to the dynamics

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
   class(type_logging), intent(in), target, optional :: logs
   TYPE(type_field_manager), intent(inout), target, optional :: fm

!  Local constants

!  Local variables
!---------------------------------------------------------------------------
   if (present(logs)) then
      self%logs => logs
      call self%logs%info('momentum_configuration()',level=2)
   end if
   if (present(fm)) then
      self%fm => fm
   end if
END SUBROUTINE momentum_configure

!---------------------------------------------------------------------------

SUBROUTINE momentum_initialize(self,runtype,domain,advection,vertical_diffusion)

   !! Initialize all dynamical components

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
   integer, intent(in) :: runtype
   TYPE(type_getm_domain), intent(inout), target :: domain
   class(type_advection), intent(in), target, optional :: advection
   class(type_vertical_diffusion), intent(in), target, optional :: vertical_diffusion

!  Local constants

!  Local variables
   integer :: i,j,k
   integer :: stat
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('momentum_initialize()',level=2)
   self%domain => domain
   if (present(advection)) then
      self%advection => advection
   end if
   if (present(vertical_diffusion)) then
      self%vertical_diffusion => vertical_diffusion
   end if
#ifndef _STATIC_
   XGrid: associate( XG => self%domain%X )
   TGrid: associate( TG => self%domain%T )
   UGrid: associate( UG => self%domain%U )
   VGrid: associate( VG => self%domain%V )
   call self%initialize_2d()
   if (runtype > 1) then
      call self%initialize_3d()
   end if
   call mm_s('work2d',self%work2d,self%U,def=0._real64,stat=stat)
#endif
   if (associated(self%fm)) then
      call self%register(runtype)
   end if
   if (self%domain%nbdyp > 0) then
      allocate(self%bdyu(self%domain%nbdyp),stat=stat)
      allocate(self%bdyv(self%domain%nbdyp),stat=stat)
      self%bdyu=0._real64
      self%bdyv=0._real64
   end if
   end associate VGrid
   end associate UGrid
   end associate TGrid
   end associate XGrid
END SUBROUTINE momentum_initialize

!---------------------------------------------------------------------------

SUBROUTINE momentum_list(self)
   !! List the configurable parameters in type_getm_momentum

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self

!  Local constants

!  Local variables
!---------------------------------------------------------------------------
   write(*,*) 'Momentum settings:'
   write(*,*) 'advection_scheme=    ',self%advection_scheme
   write(*,*) 'cnpar=               ',self%cnpar
   write(*,*) 'molecular viscosity= ',self%molecular
END SUBROUTINE momentum_list

!---------------------------------------------------------------------------

SUBROUTINE momentum_set(self,adv,cnpar,mol)
   !! List the configurable parameters in type_getm_momentum

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
   integer, intent(in), optional :: adv
   real(real64), intent(in), optional :: cnpar
   real(real64), intent(in), optional :: mol

!  Local constants

!  Local variables
!---------------------------------------------------------------------------
   if (present(adv)) then
      self%advection_scheme=adv
      write(*,*) 'setting advection_scheme to:',self%advection_scheme
   end if
   if (present(cnpar)) then
      self%cnpar=cnpar
      write(*,*) 'setting cnpar to:',self%cnpar
   end if
   if (present(cnpar)) then
      self%molecular=mol
      write(*,*) 'setting molecular to:',self%molecular
   end if
END SUBROUTINE momentum_set

!---------------------------------------------------------------------------

END MODULE getm_momentum
