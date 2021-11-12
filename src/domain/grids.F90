! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!>  This module provides all variables related to the bathymetry and

SUBMODULE (getm_domain) getm_grids_smod

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

module SUBROUTINE grid_configure(self,logs,grid_type,imin,imax,jmin,jmax,kmin,kmax,halo)

   !! Configure the type_grid

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_grid), intent(inout) :: self
   class(type_logging), intent(in), target, optional :: logs
   integer, intent(in), optional :: grid_type
   integer, intent(in), optional :: imin,imax
   integer, intent(in), optional :: jmin,jmax
   integer, intent(in), optional :: kmin,kmax
   integer, intent(in), dimension(3), optional :: halo
      !! grid dim in case of dynamic memory allocation

!  Local constants

!  Local variables
!-----------------------------------------------------------------------
   if (self%is_initialized) return
   if (present(logs)) then
      self%logs => logs
      call logs%info('grid_configure()',level=2)
   end if

   call self%type_3d_grid%create(imin=imin,imax=imax,jmin=jmin,jmax=jmax,kmin=kmin,kmax=kmax,halo=halo)
   call allocate_grid_variables(self)

   if (present(grid_type)) self%grid_type=grid_type

   self%is_initialized = .true.
END SUBROUTINE grid_configure

!-----------------------------------------------------------------------------

module SUBROUTINE set_mask(self,logs,il,jl,ih,jh,val)

   !! Manually masking an area of the grid

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_grid), intent(inout) :: self
   class(type_logging), intent(in) :: logs
   integer, intent(in) :: il,jl,ih,jh,val

!  Local constants

!  Local variables
   integer :: i,j
!-----------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('set_mask()',level=2)
   self%mask(il-self%ioff:ih-self%ioff,jl-self%joff:jh-self%joff) = val
END SUBROUTINE set_mask

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
   if (associated(self%logs)) then
      call logs%info('grid_report()',level=2)
      call logs%costum(unit,trim(header))
      call self%print(unit)
      call logs%costum(unit,'mask')
      call self%print_mask(unit)
      call logs%costum(unit,'')
   end if
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
!KB   if (associated(self%logs)) call self%type_3d_grid%print(1)
   call self%type_3d_grid%print(6)
END SUBROUTINE grid_print_info

!-----------------------------------------------------------------------------

module SUBROUTINE grid_print_mask(self,unit,print_halos)

  IMPLICIT NONE

! Subroutine arguments
   class(type_getm_grid), intent(in) :: self
   integer, intent(in) :: unit
   logical, intent(in), optional :: print_halos

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
      if (mod(j,10) == 0) then
         write(unit,'(a,*(i1))') '* ',(self%mask(i,j), i=self%imin,self%imax)
      else
         write(unit,'(a,*(i1))') '  ',(self%mask(i,j), i=self%imin,self%imax)
      endif
   end do
#endif
END SUBROUTINE grid_print_mask

!-----------------------------------------------------------------------------

!module SUBROUTINE allocate_grid_variables(self)
SUBROUTINE allocate_grid_variables(self)
   !! Allocate all domain related variables

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_grid), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: stat, l(3)
!-----------------------------------------------------------------------------
#ifndef _STATIC_
   call mm_s('H',self%H,self%l(1:2),self%u(1:2),def=-10._real64,stat=stat)
   call mm_s('z0b',self%z0b,self%l(1:2),self%u(1:2),def=0._real64,stat=stat)
   call mm_s('z0b_min',self%z0b_min,self%l(1:2),self%u(1:2),def=0._real64,stat=stat)
   call mm_s('c1',self%c1,self%l(1),self%u(1),def=0._real64,stat=stat)
   call mm_s('c2',self%c2,self%l(2),self%u(2),def=0._real64,stat=stat)
   call mm_s('lon',self%lon,self%H,def=-9999._real64,stat=stat)
   call mm_s('lat',self%lat,self%H,def=-9999._real64,stat=stat)
   call mm_s('x',self%x,self%H,def=-9999._real64,stat=stat)
   call mm_s('y',self%y,self%H,def=-9999._real64,stat=stat)
   call mm_s('mask',self%mask,self%l(1:2),self%u(1:2),def=0,stat=stat)
   call mm_s('dlon',self%dlon,self%H,def=-9999._real64,stat=stat)
   call mm_s('dlat',self%dlat,self%H,def=-9999._real64,stat=stat)
   call mm_s('dx',self%dx,self%H,def=-9999._real64,stat=stat)
   call mm_s('dy',self%dy,self%H,def=-9999._real64,stat=stat)
   call mm_s('idx',self%idx,self%H,def=-9999._real64,stat=stat)
   call mm_s('idy',self%idy,self%H,def=-9999._real64,stat=stat)
   call mm_s('area',self%area,self%H,def=-9999._real64,stat=stat)
   call mm_s('iarea',self%iarea,self%H,def=-9999._real64,stat=stat)
   call mm_s('cor',self%cor,self%H,def=0._real64,stat=stat)
   call mm_s('z',self%z,self%H,def=-9999._real64,stat=stat)
   call mm_s('zo',self%zo,self%H,def=-9999._real64,stat=stat)
   call mm_s('D',self%D,self%H,def=-9999._real64,stat=stat)
   call mm_s('zin',self%zin,self%H,def=-9999._real64,stat=stat)
   call mm_s('zio',self%zio,self%H,def=-9999._real64,stat=stat)
   call mm_s('hn',self%hn,self%l(1:3),self%u(1:3),def=-9999._real64,stat=stat)
   call mm_s('ho',self%ho,self%hn,def=-9999._real64,stat=stat)
   l = self%l+(/0,0,-1/)
   call mm_s('zf',self%zf,l,self%u,def=-9999._real64,stat=stat)
   call mm_s('zc',self%zc,self%hn,def=-9999._real64,stat=stat)
   call mm_s('alpha',self%alpha,self%H,def=1._real64,stat=stat)
#endif
END SUBROUTINE allocate_grid_variables

!-----------------------------------------------------------------------------

module SUBROUTINE deallocate_grid_variables(self)
   !! Allocate all domain related variables

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_grid), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: stat
!-----------------------------------------------------------------------------
#ifndef _STATIC_
   deallocate(self%c1)
   deallocate(self%c2)
   deallocate(self%lon)
   deallocate(self%lat)
   deallocate(self%x)
   deallocate(self%y)
   deallocate(self%H)
   deallocate(self%z0b)
   deallocate(self%z0b_min)
   deallocate(self%mask)
   deallocate(self%dlon)
   deallocate(self%dlat)
   deallocate(self%dx)
   deallocate(self%dy)
   deallocate(self%idx)
   deallocate(self%idy)
   deallocate(self%area)
   deallocate(self%iarea)
   deallocate(self%cor)
   deallocate(self%z)
   deallocate(self%zo)
   deallocate(self%D)
   deallocate(self%zin)
   deallocate(self%zio)
   deallocate(self%hn)
   deallocate(self%ho)
   deallocate(self%zf)
   deallocate(self%zc)
   deallocate(self%alpha)
#endif
END SUBROUTINE deallocate_grid_variables

!---------------------------------------------------------------------------

END SUBMODULE getm_grids_smod
