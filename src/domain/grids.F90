! Copyright (C) 2020 Bolding & Bruggeman

!>  This module provides all variables related to the bathymetry and

SUBMODULE (getm_domain) getm_grids_smod

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

module SUBROUTINE grid_configure(self,logs,imin,imax,jmin,jmax,kmin,kmax,halo)

   !! Configure the type_grid

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_grid), intent(inout) :: self
   class(type_logging), intent(in) :: logs
   integer, intent(in), optional :: imin,imax
   integer, intent(in), optional :: jmin,jmax
   integer, intent(in), optional :: kmin,kmax
   integer, intent(in), dimension(3), optional :: halo
      !! grid dim in case of dynamic memory allocation

!  Local constants

!  Local variables
!-----------------------------------------------------------------------
   if (self%is_initialized) return
   call logs%info('grid_configure()',level=2)

   call self%type_3d_grid%create(imin=imin,imax=imax,jmin=jmin,jmax=jmax,kmin=kmin,kmax=kmax,halo=halo)

!KB   self%is_initialized = .true.
   return
END SUBROUTINE grid_configure

!-----------------------------------------------------------------------------

module SUBROUTINE grid_report(self,logs,unit,header)

   !! Report the configured domain

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_grid), intent(inout) :: self
   class(type_logging), intent(in) :: logs
   integer, intent(in) :: unit
   character(len=*), intent(in) :: header

!  Local constants

!  Local variables
!-----------------------------------------------------------------------
   call logs%info('grid_report()',level=2)

   call logs%costum(unit,trim(header))
   call self%print(unit)
   call logs%costum(unit,'mask')
   call self%print_mask(unit)
   call logs%costum(unit,'')
   return
END SUBROUTINE grid_report

!-----------------------------------------------------------------------------

module SUBROUTINE grid_print_info(self,logs)
  IMPLICIT NONE

! Subroutine arguments
   class(type_getm_grid), intent(in) :: self
   class(type_logging), intent(in) :: logs

!  Local constants

!  Local variables
   integer :: i,j
!-----------------------------------------------------------------------------
   call self%type_3d_grid%print(1)
   return
END SUBROUTINE grid_print_info

!-----------------------------------------------------------------------------

module SUBROUTINE grid_print_mask(self,unit)

  IMPLICIT NONE

! Subroutine arguments
   class(type_getm_grid), intent(in) :: self
   integer, intent(in) :: unit

!  Local constants

!  Local variables
   integer :: i,j
!-----------------------------------------------------------------------------
#if 0
   do j=jmax+HALO,jmin-HALO,-1
   do j=self%u(2),self%l(2),-1
      write(unit,'(*(i1))') (mask(i,j), i=self%l(1),self%u(1))
   end do
#else
   do j=self%jmax,self%jmin,-1
      write(unit,'(*(i1))') (self%mask(i,j), i=self%imin,self%imax)
   end do
#endif
   return
END SUBROUTINE grid_print_mask

!-----------------------------------------------------------------------------

module SUBROUTINE allocate_grid_variables(self)
   !! Allocate all domain related variables

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_grid), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: stat
!-----------------------------------------------------------------------------
#ifndef _STATIC_
   call mm_s('H',self%H,self%l(1:2),self%u(1:2),def=-10._real64,stat=stat)
   call mm_s('x',self%x,self%H,def=-9999._real64,stat=stat)
   call mm_s('y',self%y,self%H,def=-9999._real64,stat=stat)
   call mm_s('dx',self%dx,self%H,def=-9999._real64,stat=stat)
   call mm_s('dy',self%dy,self%H,def=-9999._real64,stat=stat)
   call mm_s('lon',self%lon,self%H,def=-9999._real64,stat=stat)
   call mm_s('lat',self%lat,self%H,def=-9999._real64,stat=stat)
   call mm_s('dlon',self%dlon,self%H,def=-9999._real64,stat=stat)
   call mm_s('dlat',self%dlat,self%H,def=-9999._real64,stat=stat)
   call mm_s('area',self%area,self%H,def=-9999._real64,stat=stat)
   call mm_s('inv_area',self%inv_area,self%H,def=-9999._real64,stat=stat)
   call mm_s('cor',self%cor,self%H,def=-9999._real64,stat=stat)
   call mm_s('mask',self%mask,self%l(1:2),self%u(1:2),def=0,stat=stat)
   call mm_s('z',self%z,self%H,def=-9999._real64,stat=stat)
   call mm_s('zo',self%zo,self%H,def=-9999._real64,stat=stat)
   call mm_s('D',self%D,self%H,def=-9999._real64,stat=stat)
   call mm_s('ssen',self%ssen,self%H,def=-9999._real64,stat=stat)
   call mm_s('sseo',self%sseo,self%H,def=-9999._real64,stat=stat)
   call mm_s('hn',self%hn,self%l(1:3),self%u(1:3),def=-9999._real64,stat=stat)
   call mm_s('ho',self%ho,self%hn,def=-9999._real64,stat=stat)
#if 0
   call mm_s('zf',self%zf,self%hn,def=-9999._real64,stat=stat)
   call mm_s('zc',self%zc,self%S%hn,def=-9999._real64,stat=stat)
#endif
#endif
   return
END SUBROUTINE allocate_grid_variables

!---------------------------------------------------------------------------

END SUBMODULE getm_grids_smod
