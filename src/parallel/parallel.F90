! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

MODULE getm_parallel

   !! Description:
   !!   < Say what this module contains >
   !!
   !! Current Code Owner: < Name of person responsible for this code >
   !!
   !! Code Description:
   !!   Language: Fortran 90.
   !!   This code is written to JULES coding standards v1.

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use getm_domain, only: type_getm_grid
#ifdef _PARALLEL_
!#define _F08_
#ifdef _F08_
   use mpi_f08
#else
   use mpi
#endif
#endif

   IMPLICIT NONE

   PRIVATE  ! Private scope by default

!  Module constants

!  Module types and variables
   type, public :: type_getm_parallel

#ifdef _PARALLEL_
      integer :: x_line,x_lines,y_line,y_lines
      integer :: halo_line,halo_square
      integer :: xz_slice,xz_slices
      integer :: yz_slice,yz_slices
      integer :: z_column
      integer :: halo_columns
      integer :: x_size,y_size,z_size
      integer :: xy_size,xz_size,yz_size,xyz_size
#if 0
      TYPE(MPI_Request) :: requests(16)
      TYPE(MPI_Status), allocatable :: status_array(:)
!KB      TYPE(MPI_Status) :: status_array(16)
#else
      integer :: requests(16)
      integer :: status_array(MPI_STATUS_SIZE,16)
      integer :: comm=MPI_COMM_NULL
      integer :: nprocs
      integer :: myid
      CHARACTER(LEN=MPI_MAX_PROCESSOR_NAME) :: pname
#endif
#endif

      contains

      procedure :: configure => parallel_configure
      procedure :: initialize => parallel_initialize
      procedure :: barrier => parallel_barrier
      procedure :: finalize => parallel_finalize
      procedure :: update_halos_xy
      procedure :: update_halos_xyz
      generic :: update_halos => update_halos_xy, update_halos_xyz
      procedure :: wait_halo_updates

   end type type_getm_parallel

!---------------------------------------------------------------------------

CONTAINS

!---------------------------------------------------------------------------

SUBROUTINE parallel_configure(self,grid)

   !! Configure the components belonging to the dynamics

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_parallel), intent(inout) :: self
   type(type_getm_grid), intent(in) :: grid

!  Local constants

!  Local variables
   integer :: imin,imax,jmin,jmax,kmin,kmax
   integer :: n
   integer :: stat
!-----------------------------------------------------------------------
#ifdef _PARALLEL_
#ifdef _F08_
   call MPI_init(stat=stat)
#else
   call MPI_init(stat)
   ! check on error
   self%comm = MPI_COMM_WORLD
   call MPI_COMM_SIZE(self%comm,self%nprocs,stat)
write(*,*) stat
   call MPI_COMM_RANK(self%comm,self%myid,stat)
write(*,*) stat
   call MPI_GET_PROCESSOR_NAME(self%pname,n,stat)
write(*,*) stat
   write(*,*) 'nprocs= ',self%nprocs
   write(*,*) 'myid=   ',self%myid,trim(self%pname)
#endif

   call grid%type_3d_grid%print(6)
   call data_types(self)
#endif
END SUBROUTINE parallel_configure

!---------------------------------------------------------------------------

SUBROUTINE parallel_initialize(self)

   !! Initialize all dynamical components

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_parallel), intent(inout) :: self

!  Local constants

!  Local variables
!---------------------------------------------------------------------------
END SUBROUTINE parallel_initialize

!---------------------------------------------------------------------------

SUBROUTINE parallel_barrier(self)

   !! Initialize all dynamical components

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_parallel), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: stat
!---------------------------------------------------------------------------
#ifdef _PARALLEL_
   call MPI_BARRIER(self%comm,stat)
#endif
END SUBROUTINE parallel_barrier

!---------------------------------------------------------------------------

SUBROUTINE parallel_finalize(self)

   !! Initialize all dynamical components

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_parallel), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: stat
!---------------------------------------------------------------------------
#ifdef _PARALLEL_
   call MPI_Finalize(stat)
#endif
END SUBROUTINE parallel_finalize

!---------------------------------------------------------------------------

SUBROUTINE data_types(self)

   !! Initialize all dynamical components

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_parallel), intent(inout) :: self

!  Local constants

!  Local variables
!---------------------------------------------------------------------------
#ifdef _PARALLEL_
write(*,*) 'data_types()'
#if 0
USE mpi_f08
MPI_Type_create_subarray(ndims, array_of_sizes, array_of_subsizes,
        array_of_starts, order, oldtype, newtype, ierror)
    INTEGER, INTENT(IN) :: ndims, array_of_sizes(ndims),
    array_of_subsizes(ndims), array_of_starts(ndims), order
    TYPE(MPI_Datatype), INTENT(IN) :: oldtype
    TYPE(MPI_Datatype), INTENT(OUT) :: newtype
    INTEGER, OPTIONAL, INTENT(OUT) :: ierror
#endif

#endif
END SUBROUTINE data_types

!---------------------------------------------------------------------------

SUBROUTINE update_halos_xy(self,f)

   !! Initialize all dynamical components

   IMPLICIT NONE
   real(real64) :: f(:,:)

!  Subroutine arguments
   class(type_getm_parallel), intent(inout) :: self

!  Local constants

!  Local variables
!---------------------------------------------------------------------------
#ifdef _PARALLEL_
write(*,*) 'update_halos_xy()'
#endif
END SUBROUTINE update_halos_xy

!---------------------------------------------------------------------------

SUBROUTINE update_halos_xyz(self,f)

   !! Initialize all dynamical components

   IMPLICIT NONE
   real(real64) :: f(:,:,:)

!  Subroutine arguments
   class(type_getm_parallel), intent(inout) :: self

!  Local constants

!  Local variables
!---------------------------------------------------------------------------
#ifdef _PARALLEL_
write(*,*) 'update_halos_xyz()'
#endif
END SUBROUTINE update_halos_xyz

!---------------------------------------------------------------------------

SUBROUTINE wait_halo_updates(self)

   !! Initialize all dynamical components

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_parallel), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: stat
   integer :: N=16
!---------------------------------------------------------------------------
#ifdef _PARALLEL_
!KB   call MPI_Waitall(N,self%requests,self%status_array,stat)
   call MPI_Waitall(N,self%requests,self%status_array,stat)
#else
#endif
END SUBROUTINE wait_halo_updates

!---------------------------------------------------------------------------

END MODULE getm_parallel
