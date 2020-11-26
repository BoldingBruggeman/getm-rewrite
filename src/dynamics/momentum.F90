! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!> In this module, the temperature equation is processed by
!> reading in the namelist {\tt temp} and initialising the temperature field
!> (this is done in {\tt init\_temperature}),
!> and calculating the advection-diffusion-equation, which includes
!> penetrating short-wave radiation as source term (see {\tt do\_temperature}).

!> @note
!> loop boundaries in uadvgrid and vadvgrid
!> @nednote

#ifdef _STATIC_
#include "dimensions.h"
#endif

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
   real(real64), parameter :: g = 9.81_real64
      !! Gravity
   real(real64), parameter :: kappa = 0.4_real64
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
      real(real64), dimension(:,:), allocatable :: u1,v1
      real(real64), dimension(:,:), allocatable :: SxA,SyA ! Slow advection
      real(real64), dimension(:,:), allocatable :: SxB,SyB ! Slow internal pressure
      real(real64), dimension(:,:), allocatable :: SxD,SyD ! Slow diffusion
      real(real64), dimension(:,:), allocatable :: SxF,SyF ! Slow friction
      real(real64), dimension(:,:), allocatable :: ru,rv
      real(real64), dimension(:,:,:), allocatable :: pk,qk,ww
      real(real64), dimension(:,:,:), allocatable :: pka,qka
      real(real64), dimension(:,:,:), allocatable :: fpk,fqk
      real(real64), dimension(:,:,:), allocatable :: advpk,advqk
      real(real64), dimension(:,:,:), allocatable :: diffuk,diffvk
      real(real64), dimension(:,:,:), allocatable :: uk,vk
      real(real64), dimension(:,:), allocatable :: taub,taubx,tauby
      real(real64), dimension(:,:), allocatable :: rru,rrv
      real(real64), dimension(:,:), allocatable :: zub,zvb
      ! help variables
      real(real64), dimension(:,:,:), allocatable :: num,ea2,ea4
      type(type_getm_grid) :: uadvgrid,vadvgrid
      integer :: advection_scheme=1
      real(real64) :: molecular=0._real64
      real(real64) :: cnpar=1._real64
      real(real64) :: Am=0.0001_real64 ! KB
      integer :: An_method=0
#endif

      contains

      procedure :: configure => momentum_configure
      procedure :: list => momentum_list
      procedure :: set => momentum_set
      procedure :: initialize => momentum_initialize
      procedure :: register => momentum_register
      procedure :: initialize_2d => uv_initialize_2d
      procedure :: uv_momentum_2d => uv_momentum_2d
      procedure :: uv_advection_2d => uv_advection_2d
      procedure :: uivi_advection_2d => uivi_advection_2d
      procedure :: uv_diffusion_2d => uv_diffusion_2d
      procedure :: uivi_diffusion_2d => uivi_diffusion_2d
      procedure :: initialize_3d => uv_initialize_3d
      procedure :: uv_momentum_3d => uv_momentum_3d
      procedure :: w_momentum_3d => w_momentum_3d
      procedure :: advection_3d => uv_advection_3d
      procedure :: uv_diffusion_3d => uv_diffusion_3d
      procedure :: vel_2d => velocities_2d
      procedure :: vel_3d => velocities_3d
      procedure :: stresses => stresses
      procedure :: shear => velocity_shear_frequency
      procedure :: slow_terms => slow_terms
      procedure :: slow_bottom_friction => slow_bottom_friction

   end type type_getm_momentum

   INTERFACE
      ! 2D routines
      MODULE SUBROUTINE uv_initialize_2d(self)
         class(type_getm_momentum), intent(inout) :: self
      END SUBROUTINE uv_initialize_2d
      MODULE SUBROUTINE uv_momentum_2d(self,dt,tausx,tausy,dpdx,dpdy)
         class(type_getm_momentum), intent(inout) :: self
         real(real64), intent(in) :: dt
            !! timestep [s]
#define _T2_ self%domain%T%l(1):,self%domain%T%l(2):
         real(real64), intent(in) :: tausx(_T2_),tausy(_T2_)
            !! surface stresses
         real(real64), intent(in) :: dpdx(_T2_), dpdy(_T2_)
            !! surface pressure gradient - including air pressure
#undef _T2_
      END SUBROUTINE uv_momentum_2d

      MODULE SUBROUTINE uv_advection_2d(self,dt)
         class(type_getm_momentum), intent(inout) :: self
         real(real64), intent(in) :: dt
            !! timestep [s]
      END SUBROUTINE uv_advection_2d

      MODULE SUBROUTINE uivi_advection_2d(self,dt)
         class(type_getm_momentum), intent(inout) :: self
         real(real64), intent(in) :: dt
            !! timestep [s]
      END SUBROUTINE uivi_advection_2d

      MODULE SUBROUTINE uv_diffusion_2d(self)
         class(type_getm_momentum), intent(inout) :: self
      END SUBROUTINE uv_diffusion_2d

      MODULE SUBROUTINE uivi_diffusion_2d(self)
         class(type_getm_momentum), intent(inout) :: self
      END SUBROUTINE uivi_diffusion_2d

      ! 3D routines
      MODULE SUBROUTINE uv_initialize_3d(self)
         class(type_getm_momentum), intent(inout) :: self
      END SUBROUTINE uv_initialize_3d
      MODULE SUBROUTINE uv_momentum_3d(self,mode_split,dt,tausx,tausy,dpdx,dpdy,idpdx,idpdy,viscosity)
         !! Solve the 3D momemtum equations
         class(type_getm_momentum), intent(inout) :: self
         real(real64), intent(in) :: dt
            !! timestep [s]
         integer, intent(in) :: mode_split
#define _T2_ self%domain%T%l(1):,self%domain%T%l(2):
   real(real64), intent(in) :: tausx(_T2_)
     !! surface stress - x
   real(real64), intent(in) :: tausy(_T2_)
     !! surface stress - y
   real(real64), intent(in) :: dpdx(_T2_)
     !! surface pressure (including air pressure) - x-gradient
   real(real64), intent(in) :: dpdy(_T2_)
     !! surface pressure (including air pressure) - y-gradient
#undef _T2_
#define _T3_ self%domain%T%l(1):,self%domain%T%l(2):,self%domain%T%l(3):
   real(real64), intent(in) :: idpdx(_T3_)
     !! internal pressure - x-gradient
   real(real64), intent(in) :: idpdy(_T3_)
     !! internal pressure - y-gradient
   real(real64), intent(in) :: viscosity(_T3_)
     !! viscosity
#undef _T3_
      END SUBROUTINE uv_momentum_3d
      MODULE SUBROUTINE w_momentum_3d(self,dt)
         class(type_getm_momentum), intent(inout) :: self
         real(real64), intent(in) :: dt
            !! timestep [s]
      END SUBROUTINE w_momentum_3d
      MODULE SUBROUTINE uv_advection_3d(self,dt)
         class(type_getm_momentum), intent(inout) :: self
         real(real64), intent(in) :: dt
            !! timestep [s]
      END SUBROUTINE uv_advection_3d
      MODULE SUBROUTINE uv_diffusion_3d(self)
         class(type_getm_momentum), intent(inout) :: self
      END SUBROUTINE uv_diffusion_3d

      MODULE SUBROUTINE slow_terms(self,idpdx,idpdy)
#define _T3_ self%domain%T%l(1):,self%domain%T%l(2):,self%domain%T%l(3):
         class(type_getm_momentum), intent(inout) :: self
         real(real64), intent(in) :: idpdx(_T3_)
         real(real64), intent(in) :: idpdy(_T3_)
#undef _T3_
      END SUBROUTINE slow_terms
      MODULE SUBROUTINE slow_bottom_friction(self)
         class(type_getm_momentum), intent(inout) :: self
      END SUBROUTINE slow_bottom_friction

      MODULE SUBROUTINE velocities_2d(self)
         class(type_getm_momentum), intent(inout) :: self
      END SUBROUTINE velocities_2d

      MODULE SUBROUTINE velocities_3d(self)
         class(type_getm_momentum), intent(inout) :: self
      END SUBROUTINE velocities_3d

      MODULE SUBROUTINE velocity_shear_frequency(self,num,SS)
         class(type_getm_momentum), intent(inout) :: self
#define _T3_ self%domain%T%l(1):,self%domain%T%l(2):,self%domain%T%l(3):
         real(real64), intent(in) :: num(_T3_)
         real(real64), intent(inout) :: SS(_T3_)
#undef _T3_
      END SUBROUTINE velocity_shear_frequency

      MODULE SUBROUTINE stresses(self)
         class(type_getm_momentum), intent(inout) :: self
      END SUBROUTINE stresses

      MODULE SUBROUTINE momentum_register(self)
         class(type_getm_momentum), intent(inout) :: self
      END SUBROUTINE momentum_register

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

SUBROUTINE momentum_initialize(self,domain,advection,vertical_diffusion)

   !! Initialize all dynamical components

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
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
   call self%initialize_3d()

!KB should move to velocities
   call mm_s('u1',self%u1,self%U,def=0._real64,stat=stat)
   call mm_s('v1',self%v1,self%V,def=0._real64,stat=stat)
   call mm_s('uk',self%uk,self%pk,def=0._real64,stat=stat)
   call mm_s('vk',self%vk,self%qk,def=0._real64,stat=stat)

!KB some of these should move to momentum_3d
   call mm_s('taub',self%taub,self%domain%T%l(1:2),self%domain%T%u(1:2),def=0._real64,stat=stat)
   call mm_s('taubx',self%taubx,self%U,def=0._real64,stat=stat)
   call mm_s('tauby',self%tauby,self%V,def=0._real64,stat=stat)
   call mm_s('rru',self%rru,self%U,def=0._real64,stat=stat)
   call mm_s('rrv',self%rrv,self%V,def=0._real64,stat=stat)
   call mm_s('zub',self%zub,self%U,def=0._real64,stat=stat)
   call mm_s('zvb',self%zvb,self%V,def=0._real64,stat=stat)
#endif
   if (associated(self%fm)) then
      call self%register()
   end if

   end associate VGrid
   end associate UGrid
   end associate TGrid
   end associate XGrid
END SUBROUTINE momentum_initialize

!---------------------------------------------------------------------------

END MODULE getm_momentum

!---------------------------------------------------------------------------

