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

   IMPLICIT NONE

   PRIVATE  ! Private scope by default

!  Module constants

!  Module types and variables
   type, public :: type_advection
      !! author: Karsten Bolding
      !! version: v0.1
      !!

      real(real64), allocatable, private :: flux(:,:), QU(:,:)

      real(real64) :: flux_time
      real(real64) :: adv_time

      contains

      procedure :: initialize => advection_initialize

      procedure :: advection_calculate_2d
      procedure :: advection_calculate_3d
      generic   :: calculate => advection_calculate_2d, advection_calculate_3d

   end type type_advection

   INTERFACE
      module subroutine advection_initialize(self,scheme)
         class(type_advection), intent(inout) :: self
         integer, intent(in) :: scheme
      end subroutine advection_initialize

      module subroutine advection_calculate_2d(self, scheme, ugrid, u, vgrid, v, dt, tgrid, f)
         class(type_advection), intent(inout) :: self
         integer, intent(in) :: scheme
         type(type_getm_grid), intent(in) :: ugrid, vgrid
         real(real64), intent(in) :: u(:,:), v(:,:)
         real(real64), intent(in) :: dt
         type(type_getm_grid), intent(inout) :: tgrid
         real(real64), intent(inout) :: f(:,:)
      end subroutine advection_calculate_2d

      module subroutine advection_calculate_3d(self, scheme, ugrid, u, vgrid, v, dt, tgrid, f)
         class(type_advection), intent(inout) :: self
         integer, intent(in) :: scheme
         type(type_getm_grid), intent(in) :: ugrid, vgrid
         real(real64), intent(in) :: u(:,:,:), v(:,:,:)
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

      integer, private :: imin,imax,jmin,jmax,kmin,kmax

      real(real64) :: matrix_time
      real(real64) :: tridiag_time

      contains

      procedure :: initialize => vertical_diffusion_initialize
      procedure :: calculate => vertical_diffusion_calculate

   end type type_vertical_diffusion

   INTERFACE
      module subroutine vertical_diffusion_initialize(self,var)
         class(type_vertical_diffusion), intent(inout) :: self
         real(real64), dimension(:,:,:), intent(in) :: var
            !! variable to be diffused
      end subroutine vertical_diffusion_initialize

      module subroutine vertical_diffusion_calculate(self,mask,dzo,dzn,dt,cnpar,molecular,nuh,var,ea2,ea4)
         class(type_vertical_diffusion), intent(inout) :: self
         integer, dimension(:,:), intent(in) :: mask
         real(real64), dimension(:,:,:), intent(in) :: dzo,dzn
         real(real64), intent(in) :: dt
         real(real64), intent(in) :: cnpar
         real(real64), intent(in) :: molecular
         real(real64), dimension(:,:,:), intent(in) :: nuh
         real(real64), dimension(:,:,:), intent(inout) :: var
         real(real64), dimension(:,:,:), intent(in), optional :: ea2
         real(real64), dimension(:,:,:), intent(in), optional :: ea4
      end subroutine vertical_diffusion_calculate
   END INTERFACE

!---------------------------------------------------------------------------

END MODULE getm_operators

