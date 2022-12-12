module c_internal_pressure

   use iso_c_binding, only: c_ptr, c_int, c_double, c_char, c_loc, c_f_pointer, c_associated, C_NULL_CHAR, C_NULL_PTR

   implicit none

   contains

   subroutine c_blumberg_mellor(nx, ny, nz, imin, imax, jmin, jmax, umask, vmask, idxu, idyv, hu, hv, zf, &
         buoy, idpdx, idpdy) bind(c)
      integer(c_int), value, intent(in)    :: nx, ny, nz, imin, imax, jmin, jmax
      integer(c_int),        intent(in)    :: umask(nx,ny), vmask(nx,ny)
      real(c_double),        intent(in)    :: idxu(nx,ny), idyv(nx,ny)
      real(c_double),        intent(in)    :: hu(nx,ny,nz), hv(nx,ny,nz), zf(nx,ny,0:nz)
      real(c_double),        intent(in)    :: buoy(nx,ny,nz)
      real(c_double),        intent(inout) :: idpdx(nx,ny,nz), idpdy(nx,ny,nz)

      integer :: i,j,k
      real(c_double) :: grdl,grdu,buoyl,buoyu,prgr,dxz,dyz

      do j=jmin,jmax
         do i=imin,imax
            if (umask(i,j) > 0) then
               grdl=(buoy(i+1,j,nz)-buoy(i,j,nz))*idxu(i,j)
               buoyl=0.5*(buoy(i+1,j,nz)+buoy(i,j,nz))
               prgr=grdl*0.5_c_double*hu(i,j,nz)
               idpdx(i,j,nz)=hu(i,j,nz)*prgr
               do k=nz-1,1,-1
                  grdu=grdl
                  grdl=(buoy(i+1,j,k)-buoy(i,j,k))*idxu(i,j)
                  buoyu=buoyl
                  buoyl=0.5_c_double*(buoy(i+1,j,k)+buoy(i,j,k))
                  dxz=(zf(i+1,j,k)-zf(i,j,k))*idxu(i,j)
                  prgr=prgr+0.5_c_double*(grdu+grdl)*0.5_c_double*(hu(i,j,k)+hu(i,j,k+1))-dxz*(buoyu-buoyl)
                  idpdx(i,j,k)=hu(i,j,k)*prgr
               end do
            end if
         end do
      end do

      do j=jmin,jmax
         do i=imin,imax
            if (vmask(i,j) > 0) then
               grdl=(buoy(i,j+1,nz)-buoy(i,j,nz))*idyv(i,j)
               buoyl=0.5_c_double*(buoy(i,j+1,nz)+buoy(i,j,nz))
               prgr=grdl*0.5_c_double*hv(i,j,nz)
               idpdy(i,j,nz)=hv(i,j,nz)*prgr
               do k=nz-1,1,-1
                  grdu=grdl
                  grdl=(buoy(i,j+1,k)-buoy(i,j,k))*idyv(i,j)
                  buoyu=buoyl
                  buoyl=0.5_c_double*(buoy(i,j+1,k)+buoy(i,j,k))
                  dyz=(zf(i,j+1,k)-zf(i,j,k))*idyv(i,j)
                  prgr=prgr+0.5_c_double*(grdu+grdl)*0.5_c_double*(hv(i,j,k)+hv(i,j,k+1))-dyz*(buoyu-buoyl)
                  idpdy(i,j,k)=hv(i,j,k)*prgr
               end do
            end if
         end do
      end do
   end subroutine c_blumberg_mellor


   subroutine c_shchepetkin_mcwilliams(nx, ny, nz, imin, imax, jmin, jmax, mask, umask, vmask, idxu, idyv, h, z, zc, &
         buoy, idpdx, idpdy) bind(c)
      integer(c_int), value, intent(in)    :: nx, ny, nz, imin, imax, jmin, jmax
      integer(c_int),        intent(in)    :: mask(nx,ny), umask(nx,ny), vmask(nx,ny)
      real(c_double),        intent(in)    :: idxu(nx,ny), idyv(nx,ny)
      real(c_double),        intent(in)    :: h(nx,ny,nz), z(nx,ny), zc(nx,ny,nz)
      real(c_double),        intent(in)    :: buoy(nx,ny,nz)
      real(c_double),        intent(inout) :: idpdx(nx,ny,nz), idpdy(nx,ny,nz)

      real(c_double), parameter :: eps=1.e-10_c_double
      integer :: i,j,k
      real(c_double) :: cff
      real(c_double), parameter :: x=1._c_double/12._c_double
      real(c_double) :: FC
      real(c_double), allocatable :: dR(:,:,:), dZ(:,:,:), P(:,:,:), dZx(:,:), dRx(:,:)
      allocate(dR(nx,ny,0:nz), dZ(nx,ny,0:nz), P(nx,ny,nz), dZx(nx,ny), dRx(nx,ny))

      do j=jmin,jmax+1
         do i=imin,imax+1
            if (mask(i,j) > 0) then
               do k=nz-1,1,-1
                  dR(i,j,k)=buoy(i,j,k+1)-buoy(i,j,k)
                  dZ(i,j,k)=zc(i,j,k+1)-zc(i,j,k)
               end do
               dR(i,j,nz)=dR(i,j,nz-1)
               dZ(i,j,nz)=dZ(i,j,nz-1)
               dR(i,j,0)=dR(i,j,1)
               dZ(i,j,0)=dZ(i,j,1)

               do k=nz,1,-1
                  cff=2._c_double*dR(i,j,k)*dR(i,j,k-1)
                  if (cff > eps) then
                     dR(i,j,k)=cff/(dR(i,j,k)+dR(i,j,k-1))
                  else
                     dR(i,j,k)=0._c_double
                  end if
                  dZ(i,j,k)=2._c_double*dZ(i,j,k)*dZ(i,j,k-1)/(dZ(i,j,k)+dZ(i,j,k-1))
               end do

               if (nz > 1) then
                  cff=0.5_c_double*(buoy(i,j,nz)-buoy(i,j,nz-1))*0.5_c_double*h(i,j,nz) &
                      /(zc(i,j,nz)-zc(i,j,nz-1))
               else
                  cff=0.0_c_double
               end if
               P(i,j,nz)=(buoy(i,j,nz)+cff)*0.5_c_double*h(i,j,nz)
               do k=nz-1,1,-1
                  P(i,j,k)=P(i,j,k+1)+0.5_c_double*((buoy(i,j,k+1)+buoy(i,j,k)) &
                       *(zc(i,j,k+1)-zc(i,j,k))-0.2_c_double*((dR(i,j,k+1)-dR(i,j,k)) &
                       *(zc(i,j,k+1)-zc(i,j,k)-x*(dZ(i,j,k+1)+dZ(i,j,k))) &
                       -(dZ(i,j,k+1)-dZ(i,j,k))*(buoy(i,j,k+1)-buoy(i,j,k) &
                       -x*(dR(i,j,k+1)+dR(i,j,k)))))
               end do
            end if
         end do
      end do

      do k=nz,1,-1
         do j=jmin,jmax
            do i=imin,imax+2
               if (umask(i-1,j) > 0) then
                  dZx(i,j)=zc(i,j,k)-zc(i-1,j,k)
                  dRx(i,j)=buoy(i,j,k)-buoy(i-1,j,k)
               else
                  dZx(i,j)=0._c_double
                  dRx(i,j)=0._c_double
               end if
            end do
         end do

         do j=jmin,jmax
            do i=imin,imax+1
               cff=2._c_double*dZx(i,j)*dZx(i+1,j)
               if (cff > eps) then
                  dZx(i,j)=cff/(dZx(i,j)+dZx(i+1,j))
               else
                  dZx(i,j)=0._c_double
               end if
               cff=2._c_double*dRx(i,j)*dRx(i+1,j)
               if (cff > eps) then
                  dRx(i,j)=cff/(dRx(i,j)+dRx(i+1,j))
               else
                  dRx(i,j)=0._c_double
               end if
            end do
         end do

         do j=jmin,jmax
            do i=imin,imax
               if (umask(i,j) == 1) then
                  FC=0.5_c_double*((buoy(i+1,j,k)+buoy(i,j,k))*(zc(i+1,j,k)-zc(i,j,k)) &
#ifndef _STD_JACOBIAN_
                    -0.2_c_double*((dRx(i+1,j)-dRx(i,j)) &
                                *(zc(i+1,j,k)-zc(i,j,k)-x*(dZx(i+1,j)+dZx(i,j))) &
                                -(dZx(i+1,j)-dZx(i,j)) &
                                *( buoy(i+1,j,k)- buoy(i,j,k)-x*(dRx(i+1,j)+dRx(i,j)))) &
#endif
                                )
                  idpdx(i,j,k)=0.5_c_double*(h(i,j,k)+h(i+1,j,k))*idxu(i,j) &
                                    *(P(i+1,j,k)-P(i,j,k)+FC &
                                    -(z(i+1,j)-z(i,j))*0.5_c_double*(buoy(i+1,j,nz)+buoy(i,j,nz)))
               end if
            end do
         end do
      end do

      do k=nz,1,-1
         do j=jmin,jmax+2
            do i=imin,imax
               if (vmask(i,j-1) > 0) then
                  dZx(i,j)=zc(i,j,k)-zc(i,j-1,k)
                  dRx(i,j)=buoy(i,j,k)-buoy(i,j-1,k)
               else
                  dZx(i,j)=0._c_double
                  dRx(i,j)=0._c_double
               end if
            end do
         end do

         do j=jmin,jmax+1
            do i=imin,imax
               cff=2._c_double*dZx(i,j)*dZx(i,j+1)
               if (cff > eps) then
                  dZx(i,j)=cff/(dZx(i,j)+dZx(i,j+1))
               else
                  dZx(i,j)=0._c_double
               end if
               cff=2._c_double*dRx(i,j)*dRx(i,j+1)
               if (cff > eps) then
                  dRx(i,j)=cff/(dRx(i,j)+dRx(i,j+1))
               else
                  dRx(i,j)=0._c_double
               end if
            end do
         end do

         do j=jmin,jmax
            do i=imin,imax
               if (vmask(i,j) == 1) then
                  FC=0.5_c_double*((buoy(i,j+1,k)+buoy(i,j,k))*(zc(i,j+1,k)-zc(i,j,k)) &
#ifndef _STD_JACOBIAN_
                    -0.2_c_double*((dRx(i,j+1)-dRx(i,j)) &
                                *(zc(i,j+1,k)-zc(i,j,k)-x*(dZx(i,j+1)+dZx(i,j))) &
                                -(dZx(i,j+1)-dZx(i,j)) &
                                *(buoy(i,j+1,k)-buoy(i,j,k)-x*(dRx(i,j+1)+dRx(i,j)))) &
#endif
                                )
                  idpdy(i,j,k)=0.5_c_double*(h(i,j,k)+h(i,j+1,k))*idyv(i,j) &
                                    *(P(i,j+1,k)-P(i,j,k)+FC &
                                    -(z(i,j+1)-z(i,j))*0.5_c_double*(buoy(i,j+1,nz)+buoy(i,j,nz)))
               end if
            end do
         end do
      end do
   end subroutine c_shchepetkin_mcwilliams

end module
