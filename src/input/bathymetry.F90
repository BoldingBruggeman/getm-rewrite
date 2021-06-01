! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!#define _STATIC_

!>  This module provides all input from external files

MODULE getm_bathymetry

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use memory_manager
   use logging
   use input_module
   use getm_domain

   IMPLICIT NONE

!-----------------------------------------------------------------------------

   PRIVATE  ! Private scope by default

! Module constants

! Module types and variables

!   type, extends(type_input_data), public :: type_input_bathymetry

!-----------------------------------------------------------------------------
   type, public :: type_bathymetry

      TYPE(type_netcdf_input) :: H
      TYPE(type_netcdf_input) :: area
      TYPE(type_netcdf_input) :: dx
      TYPE(type_netcdf_input) :: dy
      TYPE(type_netcdf_input) :: mask
      TYPE(type_netcdf_input) :: cor
      TYPE(type_netcdf_input) :: Hu
      TYPE(type_netcdf_input) :: areau
      TYPE(type_netcdf_input) :: dxu
      TYPE(type_netcdf_input) :: dyu
      TYPE(type_netcdf_input) :: masku
      TYPE(type_netcdf_input) :: coru
      TYPE(type_netcdf_input) :: Hv
      TYPE(type_netcdf_input) :: areav
      TYPE(type_netcdf_input) :: dxv
      TYPE(type_netcdf_input) :: dyv
      TYPE(type_netcdf_input) :: maskv
      TYPE(type_netcdf_input) :: corv
      TYPE(type_netcdf_input) :: Hx
      TYPE(type_netcdf_input) :: areax
      TYPE(type_netcdf_input) :: dxx
      TYPE(type_netcdf_input) :: dyx
      TYPE(type_netcdf_input) :: maskx
      TYPE(type_netcdf_input) :: corx

   CONTAINS

      procedure :: initialize => initialize_bathymetry

   end type type_bathymetry

CONTAINS

!-----------------------------------------------------------------------------
SUBROUTINE initialize_bathymetry(self,logs,domain,domain_type)
   !! Read the bathymetry from an external file and sets the mask

   IMPLICIT NONE

! Subroutine arguments
   class(type_bathymetry), intent(inout) :: self
   class(type_logging), intent(in) :: logs
      !! grid size in case of dynamic memory allocation
   class(type_getm_domain), intent(inout), target :: domain
      !! The domain type to hold bathymetry and grid sizes
   integer, intent(inout) :: domain_type

! Local constants

! Local variables
  integer :: stat
  integer :: ncid=-1
  logical :: verbose=.true.
!-----------------------------------------------------------------------------
   call logs%info('initialize_bathymetry()',level=1)

!  open the file with bathymetri information
   call self%H%open(ncid)

!  reading T-grid
   TGrid: associate( TG => domain%T )
   self%H%v = 'H'
   call self%H%initialize(ncid)
   self%H%p2dreal64 => TG%H(TG%imin:TG%imax,TG%jmin:TG%jmax)
   call self%H%get()
   if (verbose) call self%H%print_info()

   self%mask%v = 'mask'
   call self%mask%initialize(ncid)
   self%mask%p2dint32 => TG%mask(TG%imin:TG%imax,TG%jmin:TG%jmax)
   call self%mask%get()
   if (verbose) call self%mask%print_info()

   self%dx%v = 'dx'
   call self%dx%initialize(ncid)
   self%dx%p2dreal64 => TG%dx(TG%imin:TG%imax,TG%jmin:TG%jmax)
   call self%dx%get()
   if (verbose) call self%dx%print_info()

   self%dy%v = 'dy'
   call self%dy%initialize(ncid)
   self%dy%p2dreal64 => TG%dy(TG%imin:TG%imax,TG%jmin:TG%jmax)
   call self%dy%get()
   if (verbose) call self%dy%print_info()

   self%area%v = 'area'
   call self%area%initialize(ncid)
   self%area%p2dreal64 => TG%area(TG%imin:TG%imax,TG%jmin:TG%jmax)
   call self%area%get()
   if (verbose) call self%area%print_info()

   self%cor%v = 'cor'
   call self%cor%initialize(ncid)
   self%cor%p2dreal64 => TG%cor(TG%imin:TG%imax,TG%jmin:TG%jmax)
   call self%cor%get()
   if (verbose) call self%cor%print_info()
   end associate TGrid

!  reading U-grid
   UGrid: associate( UG => domain%U )
   self%Hu%v = 'Hu'
   call self%Hu%initialize(ncid)
   self%Hu%p2dreal64 => UG%H(UG%imin:UG%imax,UG%jmin:UG%jmax)
   call self%Hu%get()
   if (verbose) call self%Hu%print_info()

   self%masku%v = 'masku'
   call self%masku%initialize(ncid)
   self%masku%p2dint32 => UG%mask(UG%imin:UG%imax,UG%jmin:UG%jmax)
   call self%masku%get()
   if (verbose) call self%masku%print_info()

   self%dxu%v = 'dxu'
   call self%dxu%initialize(ncid)
   self%dxu%p2dreal64 => UG%dx(UG%imin:UG%imax,UG%jmin:UG%jmax)
   call self%dxu%get()
   if (verbose) call self%dxu%print_info()

   self%dyu%v = 'dyu'
   call self%dyu%initialize(ncid)
   self%dyu%p2dreal64 => UG%dy(UG%imin:UG%imax,UG%jmin:UG%jmax)
   call self%dyu%get()
   if (verbose) call self%dyu%print_info()

   self%areau%v = 'areau'
   call self%areau%initialize(ncid)
   self%areau%p2dreal64 => UG%area(UG%imin:UG%imax,UG%jmin:UG%jmax)
   call self%areau%get()
   if (verbose) call self%areau%print_info()

   self%coru%v = 'cor'
   call self%coru%initialize(ncid)
   self%coru%p2dreal64 => UG%cor(UG%imin:UG%imax,UG%jmin:UG%jmax)
   call self%coru%get()
   if (verbose) call self%coru%print_info()
   end associate UGrid

!  reading V-grid
   VGrid: associate( VG => domain%V )
   self%Hv%v = 'Hv'
   call self%Hv%initialize(ncid)
   self%Hv%p2dreal64 => VG%H(VG%imin:VG%imax,VG%jmin:VG%jmax)
   call self%Hv%get()
   if (verbose) call self%Hv%print_info()

   self%maskv%v = 'maskv'
   call self%maskv%initialize(ncid)
   self%maskv%p2dint32 => VG%mask(VG%imin:VG%imax,VG%jmin:VG%jmax)
   call self%maskv%get()
   if (verbose) call self%maskv%print_info()

   self%dxv%v = 'dxv'
   call self%dxv%initialize(ncid)
   self%dxv%p2dreal64 => VG%dx(VG%imin:VG%imax,VG%jmin:VG%jmax)
   call self%dxv%get()
   if (verbose) call self%dxv%print_info()

   self%dyv%v = 'dyv'
   call self%dyv%initialize(ncid)
   self%dyv%p2dreal64 => VG%dy(VG%imin:VG%imax,VG%jmin:VG%jmax)
   call self%dyv%get()
   if (verbose) call self%dyv%print_info()

   self%areav%v = 'areav'
   call self%areav%initialize(ncid)
   self%areav%p2dreal64 => VG%area(VG%imin:VG%imax,VG%jmin:VG%jmax)
   call self%areav%get()
   if (verbose) call self%areav%print_info()

   self%corv%v = 'cor'
   call self%corv%initialize(ncid)
   self%corv%p2dreal64 => VG%cor(VG%imin:VG%imax,VG%jmin:VG%jmax)
   call self%corv%get()
   if (verbose) call self%corv%print_info()
   end associate VGrid

!  reading X-grid
   XGrid: associate( XG => domain%X )
   self%Hx%v = 'Hx'
   call self%Hx%initialize(ncid)
   self%Hx%p2dreal64 => XG%H(XG%imin:XG%imax,XG%jmin:XG%jmax)
   call self%Hx%get()
   if (verbose) call self%Hx%print_info()

   self%maskx%v = 'maskx'
   call self%maskx%initialize(ncid)
   self%maskx%p2dint32 => XG%mask(XG%imin:XG%imax,XG%jmin:XG%jmax)
   call self%maskx%get()
   if (verbose) call self%maskx%print_info()

   self%dxx%v = 'dxx'
   call self%dxx%initialize(ncid)
   self%dxx%p2dreal64 => XG%dx(XG%imin:XG%imax,XG%jmin:XG%jmax)
   call self%dxx%get()
   if (verbose) call self%dxx%print_info()

   self%dyx%v = 'dyx'
   call self%dyx%initialize(ncid)
   self%dyx%p2dreal64 => XG%dy(XG%imin:XG%imax,XG%jmin:XG%jmax)
   call self%dyx%get()
   if (verbose) call self%dyx%print_info()

   self%areax%v = 'areax'
   call self%areax%initialize(ncid)
   self%areax%p2dreal64 => XG%area(XG%imin:XG%imax,XG%jmin:XG%jmax)
   call self%areax%get()
   if (verbose) call self%areax%print_info()

   self%corx%v = 'cor'
   call self%corx%initialize(ncid)
   self%corx%p2dreal64 => XG%cor(XG%imin:XG%imax,XG%jmin:XG%jmax)
   call self%corx%get()
   if (verbose) call self%corx%print_info()
   end associate XGrid

   call self%H%close()

   domain%have_metrics=.true.

END SUBROUTINE initialize_bathymetry

!-----------------------------------------------------------------------------

END MODULE getm_bathymetry
