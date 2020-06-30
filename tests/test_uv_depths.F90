! Copyright (C) 2020 Bolding & Bruggeman

PROGRAM test_uv_depths
   !! Testing calculation of depth at U and V points

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use logging
   use getm_domain
   IMPLICIT NONE

!  Local constants

!  Local variables
   TYPE(type_logging) :: logs
   TYPE(type_getm_domain) :: domain
   integer, parameter :: imin=1, imax=70, jmin=1, jmax=20
   integer :: i, j
!-----------------------------------------------------------------------------
   logs%prepend = ''
   call logs%info('testing test_uv_depths():')
!   write(*,*) ('testing test_uv_depths():')
   domain%T%imin = imin
   domain%T%imax = imax
   domain%T%jmin = jmin
   domain%T%jmax = jmax

!   call domain%print(logs)
   call domain%T%configure(logs,imin=imin,imax=imax,jmin=jmin,jmax=jmax,kmin=-1,kmax=-1)
   call domain%configure(logs)

   call domain%initialize()

   call random_number(domain%T%H)
   where (domain%T%H > 0.10_real64)
      domain%T%mask = 1
   else where
      domain%T%mask = 0
   end where

   where (domain%T%mask == 1)
      domain%T%H = 5._real64 + 20._real64*domain%T%H
   else where
      domain%T%H = -10._real64
   end where

   call logs%info('H:')
   call domain%T%print_info(logs)
   call domain%T%print_mask(6)

   call domain%uv_depths()

   call logs%info('HU:')
   call domain%U%print_mask(6)
!   write(*,*) domain%U%H

   call logs%info('HV:')
   call domain%V%print_mask(6)
!   write(*,*) domain%V%H

END PROGRAM test_uv_depths