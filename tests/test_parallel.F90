! Copyright (C) 2020 Bolding & Bruggeman

PROGRAM test_parallel

   USE, INTRINSIC :: ISO_FORTRAN_ENV
#ifdef _PARALLEL_
!   use mpi_f08
!   use mpi_f08
#endif
   use logging
   use memory_manager
   use getm_domain, only: type_getm_grid
   use getm_parallel
   IMPLICIT NONE

!  Local constants

!  Local variables
   TYPE(type_logging) :: logs
   TYPE(type_getm_grid) :: grid
   TYPE(type_getm_parallel) :: parallel
   integer, parameter :: imin=1,imax=70,jmin=1,jmax=20,kmin=0,kmax=23
   real(real64), allocatable :: v2(:,:), v3(:,:,:)
   integer :: i, j,stat
!-----------------------------------------------------------------------------
   logs%prepend = ''
   call logs%info('test_parallel:')
   call grid%configure(imin=imin,imax=imax,jmin=jmin,jmax=jmax,kmin=kmin,kmax=kmax,halo=(/2,2,0/),logs=logs)

   call parallel%configure(grid)
   call mm_s('v2',v2,grid%l(1:2),grid%u(1:2),stat=stat)
   call mm_s('v3',v3,grid%l(1:3),grid%u(1:3),stat=stat)
   call parallel%update_halos(v2)
   call parallel%wait_halo_updates()
   call parallel%update_halos(v3)
   call parallel%wait_halo_updates()
   call parallel%barrier()
   call parallel%finalize()
END PROGRAM test_parallel
