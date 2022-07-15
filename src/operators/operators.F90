! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

MODULE getm_operators

   !! Description:
   !!   Provides a public type - type_vertical_diffusion -
   !!   providing a generic routine for calculating the vertical
   !!   diffusion for any variable provided as a subroutine
   !!   argument.
   !!   https://www.tek-tips.com/viewthread.cfm?qid=1726162
   !!   https://pdfs.semanticscholar.org/0798/fa452cda22b0b501cf1388a021931efe1686.pdf
   !!
   !! Current Code Owner: Karsten Bolding
   !!
   !! Code Description:
   !!   Language: Fortran 90.
   !!   Solves the - vertical - diffusion equation.

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   USE getm_domain
   use advection_base

   IMPLICIT NONE

   PRIVATE  ! Private scope by default

!  Module constants

!  Module types and variables
   type, public :: type_advection
      !! author: Karsten Bolding
      !! version: v0.1
      !!

      class (type_advection_base), allocatable :: op
      real(real64), allocatable, private :: flux(:,:), QU(:,:)
      !real(real64), allocatable :: D(:,:)
      !real(real64), allocatable :: hn(:,:,:)

      real(real64) :: flux_time
      real(real64) :: adv_time

      contains

      procedure :: initialize => advection_initialize

      procedure :: advection_calculate_2d
      procedure :: advection_calculate_3d
      generic   :: calculate => advection_calculate_2d, advection_calculate_3d

   end type type_advection

   INTERFACE
      module subroutine advection_initialize(self, scheme, tgrid)
         class(type_advection), intent(inout) :: self
         integer, intent(in) :: scheme
         type(type_getm_grid), intent(in) :: tgrid
      end subroutine advection_initialize

      module subroutine advection_calculate_2d(self, ugrid, u, vgrid, v, apply_diffusion, Ah, dt, tgrid, f)
         class(type_advection), intent(inout) :: self
         type(type_getm_grid), intent(in) :: ugrid, vgrid
         real(real64), intent(in) :: u(:,:), v(:,:)
         logical, intent(in) :: apply_diffusion
         real(real64), intent(in) :: Ah(:,:)
         real(real64), intent(in) :: dt
         type(type_getm_grid), intent(inout) :: tgrid
         real(real64), intent(inout) :: f(:,:)
      end subroutine advection_calculate_2d

      module subroutine advection_calculate_3d(self, ugrid, u, vgrid, v, apply_diffusion, Ah, dt, tgrid, f)
         class(type_advection), intent(inout) :: self
         type(type_getm_grid), intent(in) :: ugrid, vgrid
         real(real64), intent(in) :: u(:,:,:), v(:,:,:)
         logical, intent(in) :: apply_diffusion
         real(real64), intent(in) :: Ah(:,:)
         real(real64), intent(in) :: dt
         type(type_getm_grid), intent(inout) :: tgrid
         real(real64), intent(inout) :: f(:,:,:)
      end subroutine advection_calculate_3d

   END INTERFACE

   type, public :: type_vertical_diffusion
      !! author: Karsten Bolding
      !! version: v0.1
      !!
      real(real64), private, dimension(:,:,:), allocatable :: auxo, auxn
      real(real64), private, dimension(:,:,:), allocatable :: a1,a2,a3,a4

      integer, private :: imin,imax,jmin,jmax,kmin,kmax,halo(3)=(/2,2,0/)

      real(real64) :: matrix_time
      real(real64) :: tridiag_time

      contains

      procedure :: initialize_field => vertical_diffusion_initialize_field
      procedure :: initialize_grid => vertical_diffusion_initialize_grid
      generic   :: initialize => initialize_field, initialize_grid
      procedure :: calculate => vertical_diffusion_calculate

   end type type_vertical_diffusion

   INTERFACE
      module subroutine vertical_diffusion_initialize_grid(self,grid)
         class(type_vertical_diffusion), intent(inout) :: self
         type(type_getm_grid), intent(in) :: grid
            !! to get loop boundaries
      end subroutine vertical_diffusion_initialize_grid

      module subroutine vertical_diffusion_initialize_field(self,f,halo)
         class(type_vertical_diffusion), intent(inout) :: self
         real(real64), dimension(:,:,:), intent(in) :: f
            !! variable to be diffused
         integer, dimension(3), intent(in), optional :: halo
      end subroutine vertical_diffusion_initialize_field

      module subroutine vertical_diffusion_calculate(self,dt,cnpar,mask,dzo,dzn,molecular,nuh,var,ea2,ea4)
         class(type_vertical_diffusion), intent(inout) :: self
         real(real64), intent(in) :: dt
         real(real64), intent(in) :: cnpar
#define _T2_ self%imin-self%halo(1):,self%jmin-self%halo(2):
         integer, intent(in) :: mask(_T2_)
#undef _T2_
#define _T3_ self%imin-self%halo(1):,self%jmin-self%halo(2):,self%kmin:
         real(real64), intent(in) :: dzo(_T3_)
         real(real64), intent(in) :: dzn(_T3_)
         real(real64), intent(in) :: molecular
         real(real64), intent(in) :: nuh(_T3_)
         real(real64), intent(inout) :: var(_T3_)
         real(real64), intent(in), optional :: ea2(_T3_)
         real(real64), intent(in), optional :: ea4(_T3_)
#undef _T3_
      end subroutine vertical_diffusion_calculate
   END INTERFACE

!---------------------------------------------------------------------------

END MODULE getm_operators

