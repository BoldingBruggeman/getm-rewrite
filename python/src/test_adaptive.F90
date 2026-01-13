program test_adaptive

   use, intrinsic :: iso_fortran_env
   use c_adaptive

   implicit NONE

   integer nx, ny, nz
   integer halox, haloy
   integer, allocatable, dimension (:,:) :: mask
   real(real64), allocatable, dimension (:,:) :: H, D, zo
   real(real64), allocatable, dimension (:,:,:) :: ho
   real(real64), allocatable, dimension (:,:,:) :: NN, SS
   real(real64), allocatable, dimension (:,:,:) :: nu
   real(real64) :: decay
   integer :: hpow
   real(real64) :: chsurf, hsurf
   real(real64) :: chmidd, hmidd
   real(real64) :: chbott, hbott
   real(real64) :: cneigh, rneigh
   real(real64) :: cNN, drho
   real(real64) :: cSS, dvel
   real(real64) :: chmin, hmin
   real(real64) :: dt

   real(real64) :: csigma = 0.01_real64


   integer :: i, j, k, n, nmax=10

   nx=100
   ny=1
   nz=20
   halox=2
   haloy=2
   decay = 0.1_real64
   hpow = 3
   chsurf = -0.5_real64; hsurf = 0.5_real64
   chmidd = -0.2_real64; hmidd = 4._real64
   chbott = -0.3_real64; hbott = -0.25_real64
   cneigh = 0.1_real64; rneigh = 0.25_real64
   cNN = -1._real64; drho = 0.3_real64
   cSS = -1._real64; dvel = 0.1_real64
   chmin = 1._real64; hmin = 0.5_real64
   dt = 600._real64
   nmax=10

#define _A2_ -halox+1:nx+halox,-haloy+1:ny+haloy
   allocate(mask(_A2_))
   allocate(H(_A2_))
   allocate(D(_A2_))
   allocate(zo(_A2_))
#undef _A2_
#define _A3_ -halox+1:nx+halox,-haloy+1:ny+haloy,nz
   allocate(ho(_A3_))
   allocate(NN(_A3_))
   allocate(SS(_A3_))
   allocate(nu(_A3_))
#undef _A3_

    NN = 1000.1_real64

   mask=0
   mask(1:nx,1:ny) = 1
   where (mask > 0) H=10.
   where (mask > 0) D=10.1
   where (mask > 0) zo=-0.1
   do k=1,nz
      do j=1,ny
         do i=1,nx
            if (mask(i,j) < 1) cycle

            ho(i,j,k)=1._real64
         end do
      end do
   end do

   do n=1, nmax
      write(*,*) n, ' of ', nmax
      do k=1,nz
         where (mask > 0) nu(:,:,k)=csigma
      end do
      call c_update_adaptive(nx, ny, nz, halox, haloy, &
                             mask, H, D, zo, ho, &
                             NN, SS, &
                             nu, &
                             decay, hpow,  &
                             chsurf, hsurf,  &
                             chmidd, hmidd, &
                             chbott, hbott, &
                             cneigh, rneigh, &
                             cNN, drho, &
                             cSS, dvel, &
                             chmin, hmin, &
                             dt)
   end do

   !write(*,*) mask(1:nx,1)
   !write(*,*) nu(1:nx,1,:)

end program test_adaptive
