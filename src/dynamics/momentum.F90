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

      class(type_logging), pointer :: logs
      class(type_field_manager), pointer :: fm
      class(type_getm_domain), pointer :: domain
      class(type_advection), pointer :: advection => null()
      class(type_vertical_diffusion), pointer :: vertical_diffusion => null()

#ifdef _STATIC_
      real(real64), dimension(E2DFIELD) :: U, V
      real(real64), dimension(I3DFIELD) :: pk, qk, w
#else
      real(real64), dimension(:,:), allocatable :: U, V
      real(real64), dimension(:,:), allocatable :: Ui, Vi
      real(real64), dimension(:,:), allocatable :: Uio, Vio
      real(real64), dimension(:,:), allocatable :: Uadv,Vadv
      real(real64), dimension(:,:), allocatable :: u1, v1
      real(real64), dimension(:,:), allocatable :: fU, fV
      real(real64), dimension(:,:), allocatable :: Slru, Slrv ! Should go away
      real(real64), dimension(:,:), allocatable :: UEx, VEx ! SxA, SyA, SxB, SyB, SxD, SyD, SxF, SyF
      real(real64), dimension(:,:), allocatable :: SxA, SyA ! Slow advection
      real(real64), dimension(:,:), allocatable :: SxB, SyB ! Slow internal pressure
      real(real64), dimension(:,:), allocatable :: SxD, SyD ! Slow diffusion
      real(real64), dimension(:,:), allocatable :: SxF, SyF ! Slow friction
      real(real64), dimension(:,:), allocatable :: SlUx, SlVx ! Should go away
      real(real64), dimension(:,:), allocatable :: ru, rv
      real(real64), dimension(:,:,:), allocatable :: pk, qk, ww
      real(real64), dimension(:,:,:), allocatable :: pkadv, qkadv
      real(real64), dimension(:,:,:), allocatable :: uk, vk
      real(real64), dimension(:,:,:), allocatable :: uuEx,vvEx ! 0
      real(real64), dimension(:,:), allocatable :: taub,taubx, tauby
      real(real64), dimension(:,:), allocatable :: rru,rrv
      real(real64), dimension(:,:), allocatable :: zub,zvb
      type(type_getm_grid) :: uadvgrid,vadvgrid
      integer :: advection_scheme=1
      real(real64) :: molecular=0._real64
      real(real64) :: cnpar
#endif

      contains

      procedure :: configuration => momentum_configuration
      procedure :: initialize => momentum_initialize
      procedure :: register => momentum_register
      procedure :: uv_momentum_2d => uv_momentum_2d
      procedure :: advection_2d => uv_advection_2d
      procedure :: uv_momentum_3d => uv_momentum_3d
      procedure :: w_momentum_3d => w_momentum_3d
      procedure :: advection_3d => uv_advection_3d
      procedure :: vel_2d => velocities_2d
      procedure :: vel_3d => velocities_3d
      procedure :: stresses => stresses
      procedure :: shear => velocity_shear_frequency
      procedure :: slow_terms => slow_terms
      procedure :: slow_bottom_friction => slow_bottom_friction

   end type type_getm_momentum

   INTERFACE
      module subroutine uv_momentum_2d(self,dt,tausx,tausy,dpdx,dpdy)
         class(type_getm_momentum), intent(inout) :: self
         real(real64), intent(in) :: dt
            !! timestep [s]
#define _T_ self%domain%T%l(1):,self%domain%T%l(2):
         real(real64), dimension(:,:), intent(in) :: tausx(_T_),tausy(_T_)
            !! surface stresses
         real(real64), dimension(:,:), intent(in) :: dpdx(_T_), dpdy(_T_)
            !! surface pressure gradient - including air pressure
#undef _T_
      end subroutine uv_momentum_2d

      module subroutine uv_advection_2d(self,dt)
         class(type_getm_momentum), intent(inout) :: self
         real(real64), intent(in) :: dt
            !! timestep [s]
      end subroutine uv_advection_2d

      module subroutine uv_momentum_3d(self,mode_split,dt,tausx,tausy,dpdx,dpdy,idpdx,idpdy,viscosity)
         !! Solve the 3D momemtum equations
         class(type_getm_momentum), intent(inout) :: self
         real(real64), intent(in) :: dt
            !! timestep [s]
         integer, intent(in) :: mode_split
         real(real64), dimension(:,:), intent(in) :: tausx,tausy
           !! surface stresses
         real(real64), dimension(:,:), intent(in) :: dpdx,dpdy
           !! surface pressure gradients - including air pressure
         real(real64), dimension(:,:,:), intent(in) :: idpdx,idpdy
           !! internal pressure gradients
         real(real64), dimension(:,:,:), intent(in) :: viscosity
           !! viscosity
      end subroutine uv_momentum_3d

      module subroutine w_momentum_3d(self,dt)
         class(type_getm_momentum), intent(inout) :: self
         real(real64), intent(in) :: dt
            !! timestep [s]
      end subroutine w_momentum_3d

      module subroutine uv_advection_3d(self,dt)
         class(type_getm_momentum), intent(inout) :: self
         real(real64), intent(in) :: dt
            !! timestep [s]
      end subroutine uv_advection_3d

      module subroutine velocities_2d(self)
         class(type_getm_momentum), intent(inout) :: self
      end subroutine velocities_2d

      module subroutine velocities_3d(self)
         class(type_getm_momentum), intent(inout) :: self
      end subroutine velocities_3d

      module subroutine velocity_shear_frequency(self,num,SS)
         class(type_getm_momentum), intent(inout) :: self
         real(real64), dimension(:,:,:), intent(in) :: num
         real(real64), dimension(:,:,:), intent(inout) :: SS
      end subroutine velocity_shear_frequency

      module subroutine stresses(self)
         class(type_getm_momentum), intent(inout) :: self
      end subroutine stresses

      module subroutine momentum_register(self)
         class(type_getm_momentum), intent(inout) :: self
      end subroutine momentum_register

      module subroutine slow_terms(self,idpdx,idpdy)
         class(type_getm_momentum), intent(inout) :: self
         real(real64), dimension(:,:,:), intent(in) :: idpdx
         real(real64), dimension(:,:,:), intent(in) :: idpdy
      end subroutine slow_terms

      module subroutine slow_bottom_friction(self)
         class(type_getm_momentum), intent(inout) :: self
      end subroutine slow_bottom_friction
   END INTERFACE

CONTAINS

!---------------------------------------------------------------------------

SUBROUTINE momentum_configuration(self,logs,fm)

   !! Configure the components belonging to the dynamics

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
   class(type_logging), intent(in), target :: logs
   TYPE(type_field_manager), intent(inout), target :: fm

!  Local constants

!  Local variables
!---------------------------------------------------------------------------
   self%logs => logs
   call self%logs%info('momentum_configuration()',level=2)
   self%fm => fm
   return
END SUBROUTINE momentum_configuration

!---------------------------------------------------------------------------

SUBROUTINE momentum_initialize(self,domain,advection,vertical_diffusion)

   !! Initialize all dynamical components

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
   TYPE(type_getm_domain), intent(inout), target :: domain
   class(type_advection), intent(in), optional, target :: advection
   class(type_vertical_diffusion), intent(in), optional, target :: vertical_diffusion

!  Local constants

!  Local variables
   integer :: i,j,k
   integer :: stat
!---------------------------------------------------------------------------
   call self%logs%info('momentum_initialize()',level=2)
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

   call mm_s('U',self%U,UG%l(1:2),UG%u(1:2),def=0._real64,stat=stat)
   call mm_s('V',self%V,VG%l(1:2),VG%u(1:2),def=0._real64,stat=stat)
   call mm_s('Ui',self%Ui,self%U,def=0._real64,stat=stat)
   call mm_s('Vi',self%Vi,self%V,def=0._real64,stat=stat)
   call mm_s('Uio',self%Uio,self%U,def=0._real64,stat=stat)
   call mm_s('Vio',self%Vio,self%V,def=0._real64,stat=stat)
   call mm_s('u1',self%u1,self%U,def=0._real64,stat=stat)
   call mm_s('v1',self%v1,self%V,def=0._real64,stat=stat)
   call mm_s('fU',self%fU,self%domain%T%l(1:2),self%domain%T%u(1:2),def=0._real64,stat=stat)
   call mm_s('fV',self%fV,self%domain%T%l(1:2),self%domain%T%u(1:2),def=0._real64,stat=stat)
   call mm_s('SxA',self%SxA,self%U,def=0._real64,stat=stat)
   call mm_s('SyA',self%SyA,self%V,def=0._real64,stat=stat)
   call mm_s('SxB',self%SxB,self%U,def=0._real64,stat=stat)
   call mm_s('SyB',self%SyB,self%V,def=0._real64,stat=stat)
   call mm_s('SxD',self%SxD,self%U,def=0._real64,stat=stat)
   call mm_s('SyD',self%SyD,self%V,def=0._real64,stat=stat)
   call mm_s('SxF',self%SxF,self%U,def=0._real64,stat=stat)
   call mm_s('SyF',self%SyF,self%V,def=0._real64,stat=stat)
   call mm_s('UEx',self%UEx,self%domain%T%l(1:2),self%domain%T%u(1:2),def=0._real64,stat=stat)
   call mm_s('VEx',self%VEx,self%domain%T%l(1:2),self%domain%T%u(1:2),def=0._real64,stat=stat)
   call mm_s('SlUx',self%SlUx,self%domain%T%l(1:2),self%domain%T%u(1:2),def=0._real64,stat=stat)
   call mm_s('SlVx',self%SlVx,self%domain%T%l(1:2),self%domain%T%u(1:2),def=0._real64,stat=stat)
   call mm_s('Slru',self%Slru,self%domain%T%l(1:2),self%domain%T%u(1:2),def=0._real64,stat=stat)
   call mm_s('Slrv',self%Slrv,self%domain%T%l(1:2),self%domain%T%u(1:2),def=0._real64,stat=stat)
   call mm_s('ru',self%ru,self%domain%T%l(1:2),self%domain%T%u(1:2),def=0._real64,stat=stat)
   call mm_s('rv',self%rv,self%domain%T%l(1:2),self%domain%T%u(1:2),def=0._real64,stat=stat)
   call mm_s('taub',self%taub,self%domain%T%l(1:2),self%domain%T%u(1:2),def=0._real64,stat=stat)
   call mm_s('taubx',self%taubx,self%U,def=0._real64,stat=stat)
   call mm_s('tauby',self%tauby,self%V,def=0._real64,stat=stat)
   call mm_s('rru',self%rru,self%U,def=0._real64,stat=stat)
   call mm_s('rrv',self%rrv,self%V,def=0._real64,stat=stat)
   call mm_s('zub',self%zub,self%U,def=0._real64,stat=stat)
   call mm_s('zvb',self%zvb,self%V,def=0._real64,stat=stat)
   call mm_s('pk',self%pk,UG%l(1:3),UG%u(1:3),def=0._real64,stat=stat)
   call mm_s('qk',self%qk,VG%l(1:3),VG%u(1:3),def=0._real64,stat=stat)
   call mm_s('uk',self%uk,self%pk,def=0._real64,stat=stat)
   call mm_s('vk',self%vk,self%qk,def=0._real64,stat=stat)
   call mm_s('ww',self%ww,TG%l(1:3),TG%u(1:3),def=0._real64,stat=stat)
   call mm_s('uuEx',self%uuEx,UG%l(1:3),UG%u(1:3),def=0._real64,stat=stat)
   call mm_s('vvEx',self%vvEx,VG%l(1:3),VG%u(1:3),def=0._real64,stat=stat)
!   call mm_s('uuEx',self%uuEx,self%pk,def=0._real64,stat=stat)
!   call mm_s('vvEx',self%vvEx,self%qk,def=0._real64,stat=stat)
#endif
   if (associated(self%fm)) then
      call self%register()
   end if

   ! Grids for U and V advection - updates of time varying fields in advection calling routine
   call mm_s('uadvmask',self%uadvgrid%mask,TG%mask,def=0,stat=stat)
   call mm_s('uadvdx',self%uadvgrid%dx,TG%dx,def=0._real64,stat=stat)
   call mm_s('uadvdy',self%uadvgrid%dy,TG%dy,def=0._real64,stat=stat)
   call mm_s('uadvD',self%uadvgrid%D,TG%D,def=0._real64,stat=stat)
   call mm_s('uadvhn',self%uadvgrid%hn,TG%hn,def=0._real64,stat=stat)

   call mm_s('vadvmask',self%vadvgrid%mask,TG%mask,def=0,stat=stat)
   call mm_s('vadvdx',self%vadvgrid%dx,TG%dx,def=0._real64,stat=stat)
   call mm_s('vadvdy',self%vadvgrid%dy,TG%dy,def=0._real64,stat=stat)
   call mm_s('vadvD',self%vadvgrid%D,TG%D,def=0._real64,stat=stat)
   call mm_s('vadvhn',self%vadvgrid%hn,TG%hn,def=0._real64,stat=stat)

   call mm_s('Uadv',self%Uadv,self%U,def=0._real64,stat=stat)
   call mm_s('Vadv',self%Vadv,self%V,def=0._real64,stat=stat)
   call mm_s('pkadv',self%pkadv,self%pk,def=0._real64,stat=stat)
   call mm_s('qkadv',self%qkadv,self%qk,def=0._real64,stat=stat)
   do j=UG%jmin,UG%jmax
      do i=UG%imin,UG%imax-1 !KB - note
         self%uadvgrid%mask(i,j) = TG%mask(i+1,j) ! check this
         self%uadvgrid%dx(i,j) = TG%dx(i+1,j)
         self%uadvgrid%dy(i,j) = TG%dy(i+1,j)
         self%vadvgrid%dx(i,j) = XG%dx(i,j)
         self%vadvgrid%dy(i,j) = XG%dy(i,j)
         self%Uadv(i,j) = 0.5_real64*(self%U(i,j) + self%U(i+1,j))
         self%Vadv(i,j) = 0.5_real64*(self%V(i,j) + self%V(i+1,j))
      end do
   end do
   do j=UG%jmin,UG%jmax-1 !KB - note
      do i=UG%imin,UG%imax
         self%uadvgrid%dx(i,j) = XG%dx(i,j)
         self%uadvgrid%dy(i,j) = XG%dy(i,j)
         self%vadvgrid%mask(i,j) = TG%mask(i,j+1) ! check this
         self%vadvgrid%dx(i,j) = TG%dx(i,j+1)
         self%vadvgrid%dy(i,j) = TG%dy(i,j+1)
         self%Uadv(i,j) = 0.5_real64*(self%U(i,j) + self%U(i,j+1))
         self%Vadv(i,j) = 0.5_real64*(self%V(i,j) + self%V(i,j+1))
      end do
   end do
   end associate VGrid
   end associate UGrid
   end associate TGrid
   end associate XGrid
   return
END SUBROUTINE momentum_initialize

END MODULE getm_momentum

!---------------------------------------------------------------------------

