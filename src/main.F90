! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

PROGRAM getm

   !! Lake Ocean Model in 3d - LOM3d
   !! A re-write of GETM using modern Fortran methods.

   use getm_model
   IMPLICIT NONE

!  Local constants

!  Local variables
   TYPE(type_getm_model) :: themodel
!-----------------------------------------------------------------------------

   call themodel%settings()
   call themodel%configure()
   call themodel%initialize()
   call themodel%integrate()
   call themodel%finalize()

END PROGRAM getm

! https://gist.github.com/n-s-k/522f2669979ed6d0582b8e80cf6c95fd
