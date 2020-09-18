! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

MODULE getm_output

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use datetime_module
   use julian_days, only: init_epoch, jd=>julian_day, cal=>calendar_date
   use logging
   use field_manager
   use output_manager_core, only:output_manager_host=>host, type_output_manager_host=>type_host
   use output_manager

   IMPLICIT NONE

!-----------------------------------------------------------------------------

   PRIVATE  ! Private scope by default

!  Module constants
!   integer, parameter :: rjd=2400000
     !! Modified Julian Day - 0h Nov 16, 1858

!  Module types and variables
!   type(datetime) :: epoch
      !! used as reference time for Julian Days calculations

   type, public, extends(type_output_manager_host) :: type_getm_output
      class(type_logging), pointer :: logs => null()
      class(type_field_manager), pointer :: fm
   contains
      procedure :: configure => configure_output
      procedure :: initialize => initialize_output
      procedure :: prepare_output => prepare_output
      procedure :: do_output => do_output
      procedure :: julian_day => getm_julian_day
      procedure :: calendar_date => getm_calendar_date
   end type type_getm_output

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

SUBROUTINE configure_output(self,logs)
   !! Configure the getm_manager

   IMPLICIT NONE

! Subroutine arguments
   class(type_getm_output), intent(inout) :: self
   class(type_logging), intent(in), optional, target :: logs

! Local constants

! Local variables
   integer :: n
!-----------------------------------------------------------------------------
   if (present(logs)) self%logs => logs
   call self%logs%info('output_configure()',level=1)
   call init_epoch()
   if (associated(self%logs)) call self%logs%info('done',level=1)
   return
END SUBROUTINE configure_output

!-----------------------------------------------------------------------------

SUBROUTINE initialize_output(self,fm)
   !! Initialize the getm_manager

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_output), intent(inout) :: self
   class(type_field_manager), intent(in), target :: fm

!  Local constants

!  Local variables
   character(len=128) :: title='getm'
!-----------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('output_initialize()',level=1)
   self%fm => fm
   allocate(type_getm_output::output_manager_host)
!KB   call self%fm%list()
   call output_manager_init(self%fm,title)
   if (associated(self%logs)) call self%logs%info('done',level=1)
   return
END SUBROUTINE initialize_output

!-----------------------------------------------------------------------------

SUBROUTINE prepare_output(self,t,n)
   !! Save data to NetCDF files

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_output), intent(inout) :: self
   type(datetime), intent(in) :: t
   integer, intent(in) :: n

!  Local constants

!  Local variables
   integer :: julday
   integer :: secs, microsecs
!-----------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('do_output()',level=2)
   call self%julian_day(t%getYear(),t%getMonth(),t%getDay(),julday)
   secs = 3600*t%getHour()+60*t%getMinute()+t%getSecond()
   microsecs = 1000*t%getmilliSecond()
   call output_manager_prepare_save(julday,secs,microsecs,n)
   return
END SUBROUTINE prepare_output

!-----------------------------------------------------------------------------

SUBROUTINE do_output(self,t)
   !! Save data to NetCDF files

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_output), intent(inout) :: self
   type(datetime), intent(in) :: t

!  Local constants

!  Local variables
   integer :: julday
   integer :: secs, microsecs
!-----------------------------------------------------------------------------
   if (associated(self%logs)) then
      call self%logs%info('do_output()',level=2)
   end if
   call self%julian_day(t%getYear(),t%getMonth(),t%getDay(),julday)
   secs = 3600*t%getHour()+60*t%getMinute()+t%getSecond()
   microsecs = 1000*t%getmilliSecond()
   call output_manager_save(julday,secs,microsecs)
   return
END SUBROUTINE do_output

!-----------------------------------------------------------------------------

SUBROUTINE getm_julian_day(self,yyyy,mm,dd,julian)
   !! Return the julian day based on year, month and day

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_output), intent(in) :: self
   integer, intent(in) :: yyyy,mm,dd
   integer, intent(out)  :: julian
!-----------------------------------------------------------------------------
   call jd(yyyy,mm,dd,julian)
   return
END SUBROUTINE getm_julian_day

!-----------------------------------------------------------------------------

SUBROUTINE getm_calendar_date(self,julian,yyyy,mm,dd)
   !! Return year, month and day for the julian day

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_output), intent(in) :: self
   integer, intent(in)  :: julian
   integer, intent(out) :: yyyy,mm,dd
!-----------------------------------------------------------------------------
   call cal(julian,yyyy,mm,dd)
   return
END SUBROUTINE getm_calendar_date

#if 0
!-----------------------------------------------------------------------------

SUBROUTINE close_output(self)

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_output), intent(inout) :: self

!  Local constants

!  Local variables
!-----------------------------------------------------------------------------
   return
END SUBROUTINE close_output
#endif

!-----------------------------------------------------------------------------

END MODULE getm_output
