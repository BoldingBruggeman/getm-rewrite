! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!KB#define _STD_JACOBIAN_

!!{!./code/dynamics/ip_shchepetkin_mcwilliams.md!}

SUBMODULE (getm_pressure : pressure_internal_smod) shchepetkin_mcwilliams_smod

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   MODULE SUBROUTINE init_shchepetkin_mcwilliams(self)
!
   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_pressure), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: stat, l(3)
!-----------------------------------------------------------------------
   l = self%domain%T%l+(/0,0,-1/)
   call mm_s('dR',self%dR,l,self%domain%T%u,def=0._real64,stat=stat)
   call mm_s('dZ',self%dZ,self%dR,def=0._real64,stat=stat)
   call mm_s('P',self%P,self%domain%T%l,self%domain%T%u,def=0._real64,stat=stat)
   call mm_s('dZx',self%dZx,self%domain%T%l(1:2),self%domain%T%u(1:2),def=0._real64,stat=stat)
   call mm_s('dRx',self%dRx,self%dZx,def=0._real64,stat=stat)
   end subroutine init_shchepetkin_mcwilliams

!-----------------------------------------------------------------------------

   MODULE SUBROUTINE shchepetkin_mcwilliams(self,buoy)
!
   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_pressure), intent(inout) :: self
#define _T3_ self%domain%T%l(1):,self%domain%T%l(2):,self%domain%T%l(3):
   real(real64), intent(in) :: buoy(_T3_)
#undef _T3_

!  Local constants

!  Local variables
   real(real64), parameter :: eps=1.e-10_real64
   integer :: i,j,k
#ifdef _TIMING_
   real(real64) :: P_time,U_time,V_time
#endif
!-----------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('shchepetkin_mcwilliams()',level=3)
   TGrid: associate( TG => self%domain%T )

   Pblock: block
   real(real64) :: cff
   real(real64), parameter :: x=1._real64/12._real64
#ifdef _TIMING_
   real(real64) :: P_start, P_stop
   call cpu_time(P_start)
#endif
   do j=TG%jmin,TG%jmax+1
      do i=TG%imin,TG%imax+1
         if (TG%mask(i,j) > 0) then
            do k=TG%kmax-1,1,-1
               self%dR(i,j,k)=buoy(i,j,k+1)-buoy(i,j,k)
               self%dZ(i,j,k)=TG%zc(i,j,k+1)-TG%zc(i,j,k)
            end do
            self%dR(i,j,TG%kmax)=self%dR(i,j,TG%kmax-1)
            self%dZ(i,j,TG%kmax)=self%dZ(i,j,TG%kmax-1)
            self%dR(i,j,0)=self%dR(i,j,1)
            self%dZ(i,j,0)=self%dZ(i,j,1)

            do k=TG%kmax,1,-1
               cff=2._real64*self%dR(i,j,k)*self%dR(i,j,k-1)
               if (cff > eps) then
                  self%dR(i,j,k)=cff/(self%dR(i,j,k)+self%dR(i,j,k-1))
               else
                  self%dR(i,j,k)=0._real64
               end if
               self%dZ(i,j,k)=2._real64*self%dZ(i,j,k)*self%dZ(i,j,k-1)/(self%dZ(i,j,k)+self%dZ(i,j,k-1))
            end do

            if (TG%kmax > 1) then
               cff=0.5_real64*(buoy(i,j,TG%kmax)-buoy(i,j,TG%kmax-1))*0.5_real64*TG%hn(i,j,TG%kmax) &
                   /(TG%zc(i,j,TG%kmax)-TG%zc(i,j,TG%kmax-1))
            else
               cff=0.0_real64
            end if
            self%P(i,j,TG%kmax)=(buoy(i,j,TG%kmax)+cff)*0.5_real64*TG%hn(i,j,TG%kmax)
            do k=TG%kmax-1,1,-1
               self%P(i,j,k)=self%P(i,j,k+1)+0.5_real64*((buoy(i,j,k+1)+buoy(i,j,k)) &
                    *(TG%zc(i,j,k+1)-TG%zc(i,j,k))-0.2_real64*((self%dR(i,j,k+1)-self%dR(i,j,k)) &
                    *(TG%zc(i,j,k+1)-TG%zc(i,j,k)-x*(self%dZ(i,j,k+1)+self%dZ(i,j,k))) &
                    -(self%dZ(i,j,k+1)-self%dZ(i,j,k))*(buoy(i,j,k+1)-buoy(i,j,k) &
                    -x*(self%dR(i,j,k+1)+self%dR(i,j,k)))))
            end do
         end if
      end do
   end do
#ifdef _TIMING_
   call cpu_time(P_stop)
   P_time=P_time+P_stop-P_start
#endif
   end block Pblock

   Ublock: block
   real(real64) :: cff
   real(real64), parameter :: x=1._real64/12._real64
   real(real64) :: FC
#ifdef _TIMING_
   real(real64) :: U_start, U_stop
   call cpu_time(U_start)
#endif
   do k=TG%kmax,1,-1
      UGrid: associate( UG => self%domain%U )
      do j=TG%jmin,TG%jmax
         do i=TG%imin,TG%imax+2
            if (UG%mask(i-1,j) > 0) then
               self%dZx(i,j)=TG%zc(i,j,k)-TG%zc(i-1,j,k)
               self%dRx(i,j)=buoy(i,j,k)-buoy(i-1,j,k)
            else
               self%dZx(i,j)=0._real64
               self%dRx(i,j)=0._real64
            end if
         end do
      end do

      do j=TG%jmin,TG%jmax
         do i=TG%imin,TG%imax+1
            cff=2._real64*self%dZx(i,j)*self%dZx(i+1,j)
            if (cff > eps) then
               self%dZx(i,j)=cff/(self%dZx(i,j)+self%dZx(i+1,j))
            else
               self%dZx(i,j)=0._real64
            end if
            cff=2._real64*self%dRx(i,j)*self%dRx(i+1,j)
            if (cff > eps) then
               self%dRx(i,j)=cff/(self%dRx(i,j)+self%dRx(i+1,j))
            else
               self%dRx(i,j)=0._real64
            end if
         end do
      end do

      do j=UG%jmin,UG%jmax
         do i=UG%imin,UG%imax
            if (UG%mask(i,j) == 1) then
               FC=0.5_real64*((buoy(i+1,j,k)+buoy(i,j,k))*(TG%zc(i+1,j,k)-TG%zc(i,j,k)) &
#ifndef _STD_JACOBIAN_
                 -0.2_real64*((self%dRx(i+1,j)-self%dRx(i,j)) &
                             *(TG%zc(i+1,j,k)-TG%zc(i,j,k)-x*(self%dZx(i+1,j)+self%dZx(i,j))) &
                             -(self%dZx(i+1,j)-self%dZx(i,j)) &
                             *( buoy(i+1,j,k)- buoy(i,j,k)-x*(self%dRx(i+1,j)+self%dRx(i,j)))) &
#endif
                             )
               self%idpdx(i,j,k)=0.5_real64*(TG%hn(i,j,k)+TG%hn(i+1,j,k))*UG%idx(i,j) &
                                 *(self%P(i+1,j,k)-self%P(i,j,k)+FC &
                                 -(TG%zin(i+1,j)-TG%zin(i,j))*0.5_real64*(buoy(i+1,j,TG%kmax)+buoy(i,j,TG%kmax)))
            end if
         end do
      end do
      end associate UGrid
   end do
#ifdef _TIMING_
   call cpu_time(U_stop)
   U_time=U_time+U_stop-U_start
#endif
   end block Ublock

   Vblock: block
   real(real64) :: cff,cff1,cff2
   real(real64), parameter :: x=1._real64/12._real64
   real(real64) :: FC
#ifdef _TIMING_
   real(real64) :: V_start, V_stop
   call cpu_time(V_start)
#endif
   do k=TG%kmax,1,-1
      VGrid: associate( VG => self%domain%V )
      do j=TG%jmin,TG%jmax+2
         do i=TG%imin,TG%imax
            if (VG%mask(i,j-1) > 0) then
               self%dZx(i,j)=TG%zc(i,j,k)-TG%zc(i,j-1,k)
               self%dRx(i,j)=buoy(i,j,k)-buoy(i,j-1,k)
            else
               self%dZx(i,j)=0._real64
               self%dRx(i,j)=0._real64
            end if
         end do
      end do

      do j=TG%jmin,TG%jmax+1
         do i=TG%imin,TG%imax
            cff=2._real64*self%dZx(i,j)*self%dZx(i,j+1)
            if (cff > eps) then
               self%dZx(i,j)=cff/(self%dZx(i,j)+self%dZx(i,j+1))
            else
               self%dZx(i,j)=0._real64
            end if
            cff=2._real64*self%dRx(i,j)*self%dRx(i,j+1)
            if (cff > eps) then
               self%dRx(i,j)=cff/(self%dRx(i,j)+self%dRx(i,j+1))
            else
               self%dRx(i,j)=0._real64
            end if
         end do
      end do

      do j=VG%jmin,VG%jmax
         do i=VG%imin,VG%imax
            if (VG%mask(i,j) == 1) then
               FC=0.5_real64*((buoy(i,j+1,k)+buoy(i,j,k))*(TG%zc(i,j+1,k)-TG%zc(i,j,k)) &
#ifndef _STD_JACOBIAN_
                 -0.2_real64*((self%dRx(i,j+1)-self%dRx(i,j)) &
                             *(TG%zc(i,j+1,k)-TG%zc(i,j,k)-x*(self%dZx(i,j+1)+self%dZx(i,j))) &
                             -(self%dZx(i,j+1)-self%dZx(i,j)) &
                             *(buoy(i,j+1,k)-buoy(i,j,k)-x*(self%dRx(i,j+1)+self%dRx(i,j)))) &
#endif
                             )
               self%idpdy(i,j,k)=0.5_real64*(TG%hn(i,j,k)+TG%hn(i,j+1,k))*VG%idy(i,j) &
                                 *(self%P(i,j+1,k)-self%P(i,j,k)+FC &
                                 -(TG%zin(i,j+1)-TG%zin(i,j))*0.5_real64*(buoy(i,j+1,TG%kmax)+buoy(i,j,TG%kmax)))
            end if
         end do
      end do
      end associate VGrid
   end do
#ifdef _TIMING_
   call cpu_time(V_stop)
   V_time=V_time+V_stop-V_start
#endif
   end block Vblock
   end associate TGrid
#ifdef _TIMING_
   write(34,*) P_time,U_time,V_time
#endif
   end subroutine shchepetkin_mcwilliams

!---------------------------------------------------------------------------

END SUBMODULE shchepetkin_mcwilliams_smod
