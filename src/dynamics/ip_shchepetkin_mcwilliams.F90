! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

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
   integer :: stat
!-----------------------------------------------------------------------
   call mm_s('dR',self%dR,self%domain%T%l+(/0,0,-1/),self%domain%T%u,def=0._real64,stat=stat)
   call mm_s('dZ',self%dZ,self%dR,def=0._real64,stat=stat)
   call mm_s('P',self%P,self%domain%T%l,self%domain%T%u,def=0._real64,stat=stat)
   call mm_s('FC',self%FC,self%domain%T%l(1:2),self%domain%T%u(1:2),def=0._real64,stat=stat)
   call mm_s('dZx',self%dZx,self%FC,def=0._real64,stat=stat)
   call mm_s('dRx',self%dRx,self%FC,def=0._real64,stat=stat)
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
   real(real64), parameter :: eps=1.e-10
   real(real64), parameter :: OneFifth = 0.2
   integer :: i,j,k
   real(real64) :: P_time,U_time,V_time
!-----------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('shchepetkin_mcwilliams()',level=3)
   TGrid: associate( TG => self%domain%T )
#if 0
write(*,*) lbound(buoy)
write(*,*) ubound(buoy)
write(*,*) lbound(self%dR)
write(*,*) ubound(self%dR)
write(*,*) lbound(self%FC)
write(*,*) ubound(self%FC)
stop 'aaa'
#endif

   Pblock: block
   real(real64) :: P_start, P_stop
   real(real64) :: cff
   real(real64), parameter :: x=1._real64/12._real64
   call cpu_time(P_start)
   do j=TG%jmin,TG%jmax+1
      do i=TG%imin,TG%imax+1
         if (TG%mask(i,j) > 0) then
            do k=1,TG%kmax-1
               self%dR(i,j,k)=buoy(i,j,k+1)-buoy(i,j,k)
               self%dZ(i,j,k)=TG%zc(i,j,k+1)-TG%zc(i,j,k)
            end do
            self%dR(i,j,TG%kmax)=self%dR(i,j,TG%kmax-1)
            self%dZ(i,j,TG%kmax)=self%dZ(i,j,TG%kmax-1)
            self%dR(i,j,0)=self%dR(i,j,1)
            self%dZ(i,j,0)=self%dZ(i,j,1)

            do k=TG%kmax,1,-1
               cff=2.0*self%dR(i,j,k)*self%dR(i,j,k-1)
               if (cff > eps) then
                  self%dR(i,j,k)=cff/(self%dR(i,j,k)+self%dR(i,j,k-1))
               else
                  self%dR(i,j,k)=0._real64
               end if
               self%dZ(i,j,k)=2._real64*self%dZ(i,j,k)*self%dZ(i,j,k-1)/(self%dZ(i,j,k)+self%dZ(i,j,k-1))
            end do

            cff=0.5_real64*(buoy(i,j,TG%kmax)-buoy(i,j,TG%kmax-1)) &
               *0.5_real64*TG%hn(i,j,TG%kmax)/(TG%zc(i,j,TG%kmax)-TG%zc(i,j,TG%kmax-1))

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
   call cpu_time(P_stop)
   P_time=P_time+P_stop-P_start
   end block Pblock

   Ublock: block
   real(real64) :: U_start, U_stop
   real(real64) :: cff,cff1,cff2
   real(real64), parameter :: x=1._real64/12._real64
   call cpu_time(U_start)
   do k=TG%kmax,1,-1
      do j=TG%jmin,TG%jmax
         do i=TG%imin,TG%imax+2
            self%dZx(i,j)=TG%zc(i,j,k)-TG%zc(i-1,j,k)
            self%dRx(i,j)=buoy(i,j,k)-buoy(i-1,j,k)
#if 1
if ((i < 5 .or. i > 96) .and. j == 15) then
if (self%domain%U%mask(i,j) > 0) then
write(*,*) i,j,self%domain%U%mask(i,j),TG%zc(i,j,k)
write(*,*) TG%zc(i,j,k)-TG%zc(i-1,j,k),buoy(i,j,k)-buoy(i-1,j,k)
end if
end if
#endif
         end do
      end do
stop 'shchepetkin_mcwilliams'

      do j=TG%jmin,TG%jmax
         do i=TG%imin,TG%imax+1
            cff=2.0*self%dZx(i,j)*self%dZx(i+1,j)
            if (cff > eps) then
               cff1=1.0/(self%dZx(i,j)+self%dZx(i+1,j))
               self%dZx(i,j)=cff*cff1
            else
               self%dZx(i,j)=0._real64
            end if
            cff1=2.0*self%dRx(i,j)*self%dRx(i+1,j)
            if (cff1 > eps) then
               cff2=1.0/(self%dRx(i,j)+self%dRx(i+1,j))
               self%dRx(i,j)=cff1*cff2
            else
               self%dRx(i,j)=0._real64
            end if
         end do
      end do

      UGrid: associate( UG => self%domain%U )
      do j=UG%jmin,UG%jmax
         do i=UG%imin,UG%imax
            if (UG%mask(i,j) > 0) then
               self%FC(i,j) = 0.5*((buoy(i+1,j,k)+buoy(i,j,k))*(TG%zc(i+1,j,k)-TG%zc(i,j,k)) &
#if 1
                             -0.2_real64*((self%dRx(i+1,j)-self%dRx(i,j)) &
                                     *(TG%zc(i+1,j,k)-TG%zc(i,j,k)-x*(self%dZx(i+1,j)+self%dZx(i,j))) &
                                     -(self%dZx(i+1,j)-self%dZx(i,j)) &
                                     *( buoy(i+1,j,k)- buoy(i,j,k)-x*(self%dRx(i+1,j)+self%dRx(i,j)))) &
#endif
                                  )
               self%idpdx(i,j,k)=UG%hn(i,j,k)/UG%dx(i,j) &
                                 *(self%P(i+1,j,k)-self%P(i,j,k)+self%FC(i,j) &
                                 -(UG%zio(i+1,j)-UG%zio(i,j))*0.5_real64*(buoy(i+1,j,TG%kmax)+buoy(i,j,TG%kmax)))
            end if
         end do
      end do
      end associate UGrid
   end do
   call cpu_time(U_stop)
   U_time=U_time+U_stop-U_start
   end block Ublock

   Vblock: block
   real(real64) :: V_start, V_stop
   real(real64) :: cff,cff1,cff2
   real(real64), parameter :: x=1._real64/12._real64
   call cpu_time(V_start)
   do k=TG%kmax,1,-1
      do j=TG%jmin,TG%jmax+2
         do i=TG%imin,TG%imax
            self%dZx(i,j)=TG%zc(i,j,k)-TG%zc(i,j-1,k)
            self%dRx(i,j)=buoy(i,j,k)-buoy(i,j-1,k)
         end do
      end do

      do j=TG%jmin,TG%jmax+1
         do i=TG%imin,TG%imax
            cff=2._real64*self%dZx(i,j)*self%dZx(i,j+1)
            if (cff > eps) then
               cff1=1._real64/(self%dZx(i,j)+self%dZx(i,j+1))
               self%dZx(i,j)=cff*cff1
            else
               self%dZx(i,j)=0._real64
            end if
            cff1=2._real64*self%dRx(i,j)*self%dRx(i,j+1)
            if (cff1 > eps) then
               cff2=1._real64/(self%dRx(i,j)+self%dRx(i,j+1))
               self%dRx(i,j)=cff1*cff2
            else
               self%dRx(i,j)=0._real64
            end if
         end do
      end do

      VGrid: associate( VG => self%domain%V )
      do j=VG%jmin,VG%jmax
         do i=VG%imin,VG%imax
            if (VG%mask(i,j) > 0) then
               self%FC(i,j) = 0.5_real64*((buoy(i,j+1,k)+buoy(i,j,k))*(TG%zc(i,j+1,k)-TG%zc(i,j,k)) &
#if 1
                             -0.2_real64*((self%dRx(i,j+1)-self%dRx(i,j)) &
                                     *(TG%zc(i,j+1,k)-TG%zc(i,j,k)-x*(self%dZx(i,j+1)+self%dZx(i,j))) &
                                     -(self%dZx(i,j+1)-self%dZx(i,j)) &
                                     *(buoy(i,j+1,k)-buoy(i,j,k)-x*(self%dRx(i,j+1)+self%dRx(i,j)))) &
#endif
                                  )
               self%idpdy(i,j,k)=VG%hn(i,j,k)/VG%dy(i,j) &
                                 *(self%P(i,j+1,k)-self%P(i,j,k)+self%FC(i,j) &
                                 -(VG%zio(i,j+1)-VG%zio(i,j))*0.5_real64*(buoy(i,j+1,TG%kmax)+buoy(i,j,TG%kmax)))
            end if
         end do
      end do
      end associate VGrid
   end do
   call cpu_time(V_stop)
   V_time=V_time+V_stop-V_start
   end block Vblock
   end associate TGrid
   write(34,*) P_time,U_time,V_time
   end subroutine shchepetkin_mcwilliams

!---------------------------------------------------------------------------

END SUBMODULE shchepetkin_mcwilliams_smod
