! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

SUBMODULE (getm_operators) advection_smod

!KB   use getm_domain, only: type_getm_domain

!-----------------------------------------------------------------------------

INTERFACE
   module subroutine upstream_initialize()
!      class(type_advection), intent(inout) :: self
   end subroutine upstream_initialize
   module subroutine upstream_calculate(self,domain,dt,var)
      class(type_advection), intent(inout) :: self
      class(type_getm_domain), intent(in) :: domain
      real(real64), intent(in) :: dt
      real(real64), dimension(:,:,:), intent(inout) :: var
   end subroutine upstream_calculate
#if 0
   module subroutine init_gvc(self)
      class(type_getm_domain), intent(inout) :: self
   end subroutine init_gvc
   module subroutine do_gvc(self,dt)
      class(type_getm_domain), intent(inout) :: self
      real(real64), intent(in) :: dt
   end subroutine do_gvc
#endif
END INTERFACE

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

module SUBROUTINE advection_initialize(self,scheme)

   !! Initialize the salinity field

   IMPLICIT NONE

   ! Subroutine arguments
   class(type_advection), intent(inout) :: self
   integer, intent(in) :: scheme
!   real(real64), dimension(:,:,:), intent(in) :: var

!  Local constants

!  Local variables
!---------------------------------------------------------------------------
   select case (scheme)
      case (1)
         call upstream_initialize()
   end select
   return
END SUBROUTINE advection_initialize

!---------------------------------------------------------------------------

module SUBROUTINE advection_calculate(self,scheme,domain,dt,var)

   IMPLICIT NONE

   ! Subroutine arguments
   class(type_advection), intent(inout) :: self
   integer, intent(in) :: scheme
   class(type_getm_domain), intent(in) :: domain
   real(real64), intent(in) :: dt
   real(real64), dimension(:,:,:), intent(inout) :: var

!  Local constants

!  Local variables
   integer :: i,j,k
!---------------------------------------------------------------------------
   select case (scheme)
      case (1)
!         call do_upstream(grid,dt,var)
   end select
   return
END SUBROUTINE advection_calculate

!---------------------------------------------------------------------------

END SUBMODULE advection_smod
