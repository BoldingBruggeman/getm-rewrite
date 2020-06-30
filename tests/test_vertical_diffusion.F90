! Copyright (C) 2020 Bolding & Bruggeman

PROGRAM test_vertical_diffusion
   !! Testing calculation of depth at U and V points

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   use getm_operators, only: type_vertical_diffusion
   IMPLICIT NONE

!  Local constants
   integer, parameter :: imin=1, imax=100, jmin=1, jmax=100, kmin=0, kmax=100
   real(real64), parameter :: tmax=10._real64
   integer, parameter :: Nmax=100

!  Local variables
   type(type_vertical_diffusion) :: vertical_diffusion
   integer :: mask(imin:imax,jmin:jmax)
   real(real64) :: z(imin:imax,jmin:jmax,kmin:kmax)
   real(real64) :: dz(imin:imax,jmin:jmax,kmin:kmax)
   real(real64) :: dt, cnpar, avmol
   real(real64) :: nuh(imin:imax,jmin:jmax,kmin:kmax)
   real(real64) :: var(imin:imax,jmin:jmax,kmin:kmax)

   integer :: i,j,k,n
!-----------------------------------------------------------------------------

   mask = 1
   dt=tmax/Nmax
   cnpar=0.5_real64
   avmol=1._real64
   nuh = 0._real64

   ! initial condition
   i=1; j=1
   z(:,:,0)= -kmax*1._real64
   dz(:,:,:) = 1._real64
   var(:,:,0)= 1._real64
   do k=1,kmax
     z(:,:,k)= z(:,:,k-1)+dz(:,:,k)
     var(:,:,k)=tanh(-(z(i,j,k)+50._real64))
   end do
   var(:,:,kmax)= -1._real64

   call vertical_diffusion%initialize(var)

!   write(*,*) lbound(var,1), ubound(var,1)                                                                                                                                          
!   write(*,*) lbound(var,2), ubound(var,2)                                                                                                                                          
!   write(*,*) lbound(var,3), ubound(var,3)                                                                                                                                          

   do k=1,kmax
      write(105,*) k,var(i,j,k)
   end do

   do n=1,Nmax
!   do n=1,1
      call vertical_diffusion%calculate(mask,dz,dt,cnpar,avmol,nuh,var)
         write(100,*) n,var(i,j,55)
   end do
   do k=1,kmax
      write(110,*) k,var(i,j,k)
   end do

   write(*,*) vertical_diffusion%matrix_time/Nmax,vertical_diffusion%tridiag_time/Nmax

END PROGRAM test_vertical_diffusion
