! Copyright (C) 2020 Bolding & Bruggeman

!> In this module, the temperature equation is processed by
!> reading in the namelist {\tt temp} and initialising the temperature field
!> (this is done in {\tt init\_temperature}),
!> and calculating the advection-diffusion-equation, which includes
!> penetrating short-wave radiation as source term (see {\tt do\_temperature}).

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

   IMPLICIT NONE

   PRIVATE  ! Private scope by default

!  Module constants
   real(real64), parameter :: rho_0 = 1025._real64
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

#ifdef _STATIC_
      real(real64), dimension(E2DFIELD) :: U, V
      real(real64), dimension(I3DFIELD) :: uu, vv, w
#else
      real(real64), dimension(:,:), allocatable :: U, V
      real(real64), dimension(:,:), allocatable :: Uint, Vint
      real(real64), dimension(:,:), allocatable :: Uinto, Vinto
      real(real64), dimension(:,:), allocatable :: vel2dx, vel2dy
      real(real64), dimension(:,:), allocatable :: fU, fV
      real(real64), dimension(:,:), allocatable :: Slru, Slrv
      real(real64), dimension(:,:), allocatable :: UEx, VEx
      real(real64), dimension(:,:), allocatable :: SlUx, SlVx
      real(real64), dimension(:,:), allocatable :: ru, rv
      real(real64), dimension(:,:), allocatable :: dry_u, dry_v
      real(real64), dimension(:,:,:), allocatable :: uu, vv, ww
      real(real64), dimension(:,:,:), allocatable :: vel3dx, vel3dy
      real(real64), dimension(:,:,:), allocatable :: uuEx,vvEx
      real(real64), dimension(:,:), allocatable :: taub,taubx, tauby
      real(real64), dimension(:,:), allocatable :: rru,rrv
      real(real64), dimension(:,:), allocatable :: zub,zvb
#endif

      contains

      procedure :: configuration => momentum_configuration
      procedure :: initialize => momentum_initialize
      procedure :: x_2d => momentum_x_2d
      procedure :: y_2d => momentum_y_2d
!KB      procedure :: do_2d => momentum_2d
      procedure :: do_3d => momentum_3d
      procedure :: do_w => momentum_z_3d
      procedure :: vel_2d => velocities_2d
      procedure :: vel_3d => velocities_3d
      procedure :: stresses => stresses
      procedure :: shear => velocity_shear_frequency
      procedure :: register => momentum_register
      procedure :: slow_terms => slow_terms
      procedure :: slow_bottom_friction => slow_bottom_friction

   end type type_getm_momentum

   INTERFACE
      module subroutine momentum_x_2d(self,dt,dpdx,taus)
         class(type_getm_momentum), intent(inout) :: self
         real(real64), intent(in) :: dt
            !! timestep [s]
         real(real64), dimension(:,:), intent(in) :: dpdx
            !! surface pressure gradient - including air pressure
         real(real64), dimension(:,:), intent(in) :: taus
            !! surface stress in X-direction
      end subroutine momentum_x_2d
      module subroutine momentum_y_2d(self,dt,dpdy,taus)
         class(type_getm_momentum), intent(inout) :: self
         real(real64), intent(in) :: dt
            !! timestep [s]
         real(real64), dimension(:,:), intent(in) :: dpdy
            !! surface pressure gradient - including air pressure
         real(real64), dimension(:,:), intent(in) :: taus
            !! surface stress in Y-direction
      end subroutine momentum_y_2d
      module subroutine momentum_x_3d(self,dt,dpdx,taus)
         class(type_getm_momentum), intent(inout) :: self
         real(real64), intent(in) :: dt
            !! timestep [s]
         real(real64), dimension(:,:), intent(in) :: dpdx
            !! surface pressure gradient - including air pressure
         real(real64), dimension(:,:), intent(in) :: taus
            !! surface stress in X-direction
      end subroutine momentum_x_3d
      module subroutine momentum_y_3d(self,dt,dpdy,taus)
         class(type_getm_momentum), intent(inout) :: self
         real(real64), intent(in) :: dt
            !! timestep [s]
         real(real64), dimension(:,:), intent(in) :: dpdy
            !! surface pressure gradient - including air pressure
         real(real64), dimension(:,:), intent(in) :: taus
            !! surface stress in Y-direction
      end subroutine momentum_y_3d
      module subroutine momentum_z_3d(self,dt)
         class(type_getm_momentum), intent(inout) :: self
         real(real64), intent(in) :: dt
            !! timestep [s]
      end subroutine momentum_z_3d
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
      module subroutine slow_bottom_friction(self)
         class(type_getm_momentum), intent(inout) :: self
      end subroutine slow_bottom_friction
      module subroutine slow_terms(self,idpdx,idpdy)
         class(type_getm_momentum), intent(inout) :: self
         real(real64), dimension(:,:,:), intent(in) :: idpdx
         real(real64), dimension(:,:,:), intent(in) :: idpdy
      end subroutine slow_terms
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

SUBROUTINE momentum_initialize(self,domain)

   !! Initialize all dynamical components

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
   TYPE(type_getm_domain), intent(inout), target :: domain

!  Local constants

!  Local variables
   integer :: imin,imax,jmin,jmax,kmax
   integer :: stat
!---------------------------------------------------------------------------
   call self%logs%info('momentum_initialize()',level=2)
   self%domain => domain
#if 0
   imin = grid_dims%imin; imax = grid_dims%imax
   jmin = grid_dims%jmin; jmax = grid_dims%jmax
   kmax = grid_dims%kmax
#endif
#ifndef _STATIC_
   call mm_s('U',self%U,self%domain%U%l(1:2),self%domain%U%u(1:2),def=0._real64,stat=stat)
   call mm_s('V',self%V,self%domain%V%l(1:2),self%domain%V%u(1:2),def=0._real64,stat=stat)
   call mm_s('Uint',self%Uint,self%U,def=0._real64,stat=stat)
   call mm_s('Vint',self%Vint,self%V,def=0._real64,stat=stat)
   call mm_s('Uinto',self%Uinto,self%U,def=0._real64,stat=stat)
   call mm_s('Vinto',self%Vinto,self%V,def=0._real64,stat=stat)
   call mm_s('vel2dx',self%vel2dx,self%U,def=0._real64,stat=stat)
   call mm_s('vel2dy',self%vel2dy,self%V,def=0._real64,stat=stat)
   call mm_s('fU',self%fU,self%domain%T%l(1:2),self%domain%T%u(1:2),def=0._real64,stat=stat)
   call mm_s('fV',self%fV,self%domain%T%l(1:2),self%domain%T%u(1:2),def=0._real64,stat=stat)
   call mm_s('UEx',self%UEx,self%domain%T%l(1:2),self%domain%T%u(1:2),def=0._real64,stat=stat)
   call mm_s('VEx',self%VEx,self%domain%T%l(1:2),self%domain%T%u(1:2),def=0._real64,stat=stat)
   call mm_s('SlUx',self%SlUx,self%domain%T%l(1:2),self%domain%T%u(1:2),def=0._real64,stat=stat)
   call mm_s('SlVx',self%SlVx,self%domain%T%l(1:2),self%domain%T%u(1:2),def=0._real64,stat=stat)
   call mm_s('Slru',self%Slru,self%domain%T%l(1:2),self%domain%T%u(1:2),def=0._real64,stat=stat)
   call mm_s('Slrv',self%Slrv,self%domain%T%l(1:2),self%domain%T%u(1:2),def=0._real64,stat=stat)
   call mm_s('ru',self%ru,self%domain%T%l(1:2),self%domain%T%u(1:2),def=0._real64,stat=stat)
   call mm_s('rv',self%rv,self%domain%T%l(1:2),self%domain%T%u(1:2),def=0._real64,stat=stat)
   call mm_s('dry_u',self%dry_u,self%domain%T%l(1:2),self%domain%T%u(1:2),def=0._real64,stat=stat)
   call mm_s('dry_v',self%dry_v,self%domain%T%l(1:2),self%domain%T%u(1:2),def=0._real64,stat=stat)
   call mm_s('taub',self%taub,self%domain%T%l(1:2),self%domain%T%u(1:2),def=0._real64,stat=stat)
   call mm_s('taubx',self%taubx,self%U,def=0._real64,stat=stat)
   call mm_s('tauby',self%tauby,self%V,def=0._real64,stat=stat)
   call mm_s('rru',self%rru,self%U,def=0._real64,stat=stat)
   call mm_s('rrv',self%rrv,self%V,def=0._real64,stat=stat)
   call mm_s('uu',self%uu,self%domain%U%l(1:3),self%domain%U%u(1:3),def=0._real64,stat=stat)
   call mm_s('vv',self%vv,self%domain%V%l(1:3),self%domain%V%u(1:3),def=0._real64,stat=stat)
   call mm_s('ww',self%ww,self%domain%T%l(1:3),self%domain%T%u(1:3),def=0._real64,stat=stat)
   call mm_s('uuEx',self%uuEx,self%domain%U%l(1:3),self%domain%U%u(1:3),def=0._real64,stat=stat)
   call mm_s('vvEx',self%vvEx,self%domain%V%l(1:3),self%domain%V%u(1:3),def=0._real64,stat=stat)
!   call mm_s('uuEx',self%uuEx,self%uu,def=0._real64,stat=stat)
!   call mm_s('vvEx',self%vvEx,self%vv,def=0._real64,stat=stat)
#endif
   if (associated(self%fm)) then
      call self%register()
   end if
   return
END SUBROUTINE momentum_initialize

!---------------------------------------------------------------------------

SUBROUTINE momentum_2d(self,dt,tausx,dpdx,tausy,dpdy)
   !! Solve the 2D momemtum equations

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
   real(real64), intent(in) :: dt
   real(real64), intent(in), dimension(:,:) :: tausx,dpdx,tausy,dpdy

!  Local constants

!  Local variables
   logical, save :: ufirst=.false. ! should likely not have save
   !---------------------------------------------------------------------------
   if(ufirst) then
      call self%x_2d(dt,tausx,dpdx)
      call self%y_2d(dt,tausy,dpdy)
      ufirst = .false.
   else
      call self%y_2d(dt,tausy,dpdy)
      call self%x_2d(dt,tausx,dpdx)
      ufirst = .true.
   end if
   return
END SUBROUTINE momentum_2d

!---------------------------------------------------------------------------

SUBROUTINE momentum_3d(self,mode_split)
   !! Solve the 3D momemtum equations

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
   integer, intent(in) :: mode_split

!  Local constants

!  Local variables
   logical, save :: ufirst=.false. ! should likely not have save
!---------------------------------------------------------------------------
   self%Uint=self%Uint/mode_split
   self%Vint=self%Vint/mode_split

   if(ufirst) then
!      call self%x_3d(dt,tausx,dpdx)
!      call self%y_3d(dt,tausy,dpdy)
      ufirst = .false.
   else
 !     call self%y_3d(dt,tausy,dpdy)
 !     call self%x_3d(dt,tausx,dpdx)
      ufirst = .true.
   end if
   return
END SUBROUTINE momentum_3d

!---------------------------------------------------------------------------

END MODULE getm_momentum
