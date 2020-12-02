! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!!{!./code/dynamics/ip_blumberg_mellor.md!}

SUBMODULE (getm_pressure : pressure_internal_smod) blumberg_mellor_smod

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

MODULE SUBROUTINE blumberg_mellor(self,buoy)

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_pressure), intent(inout) :: self
#define _T3_ self%domain%T%l(1):,self%domain%T%l(2):,self%domain%T%l(3):
   real(real64), intent(in) :: buoy(:,:,:)
#undef _T3_

!  Local constants

!  Local variables
   integer :: i,j,k
   real(real64) :: dxm1,dym1
   real(real64) :: grdl,grdu,buoyl,buoyu,prgr,dxz,dyz
!-----------------------------------------------------------------------
   TGrid: associate( TG => self%domain%T )
   UGrid: associate( UG => self%domain%U )
   do j=UG%jmin,UG%jmax
      do i=UG%imin,UG%imax
         self%idpdx(i,j,0) = 0._real64
         if (UG%mask(i,j) > 0) then
            dxm1=1._real64/UG%dx(i,j)
            grdl=(buoy(i+1,j,UG%kmax)-buoy(i,j,UG%kmax))*dxm1
            buoyl=0.5*(buoy(i+1,j,UG%kmax)+buoy(i,j,UG%kmax))
            prgr=grdl*0.5_real64*UG%hn(i,j,UG%kmax)
            self%idpdx(i,j,UG%kmax)=UG%hn(i,j,UG%kmax)*prgr
            do k=UG%kmax-1,1,-1
               grdu=grdl
               grdl=(buoy(i+1,j,k)-buoy(i,j,k))*dxm1
               buoyu=buoyl
               buoyl=0.5_real64*(buoy(i+1,j,k)+buoy(i,j,k))
write(*,*) i,j,k,dxm1
               dxz=(TG%zf(i+1,j,k)-TG%zf(i,j,k))*dxm1
               prgr=prgr+0.5_real64*(grdu+grdl)*0.5_real64*(UG%hn(i,j,k)+UG%hn(i,j,k+1))-dxz*(buoyu-buoyl)
               self%idpdx(i,j,k)=UG%hn(i,j,k)*prgr
            end do
         end if
      end do
   end do
   end associate UGrid

   VGrid: associate( VG => self%domain%V )
   do j=VG%jmin,VG%jmax
      do i=VG%imin,VG%imax
         self%idpdy(i,j,0) = 0._real64
         if (VG%mask(i,j) > 0) then
            dym1 = 1._real64/VG%dy(i,j)
            grdl=(buoy(i,j+1,VG%kmax)-buoy(i,j,VG%kmax))*dym1
            buoyl=0.5_real64*(buoy(i,j+1,VG%kmax)+buoy(i,j,VG%kmax))
            prgr=grdl*0.5_real64*VG%hn(i,j,VG%kmax)
            self%idpdy(i,j,VG%kmax)=VG%hn(i,j,VG%kmax)*prgr
            do k=VG%kmax-1,1,-1
               grdu=grdl
               grdl=(buoy(i,j+1,k)-buoy(i,j,k))*dym1
               buoyu=buoyl
               buoyl=0.5_real64*(buoy(i,j+1,k)+buoy(i,j,k))
               dyz=(TG%zf(i,j+1,k)-TG%zf(i,j,k))*dym1
               prgr=prgr+0.5_real64*(grdu+grdl)*0.5_real64*(VG%hn(i,j,k)+VG%hn(i,j,k+1))-dyz*(buoyu-buoyl)
               self%idpdy(i,j,k)=VG%hn(i,j,k)*prgr
            end do
         end if
      end do
   end do
   end associate VGrid
   end associate TGrid
   end subroutine blumberg_mellor

!---------------------------------------------------------------------------

END SUBMODULE blumberg_mellor_smod
