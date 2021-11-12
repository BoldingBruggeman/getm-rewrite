! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard
!! @note
!! check with GETM
!! @endnote

SUBMODULE (getm_momentum) bottom_friction_smod

CONTAINS

!---------------------------------------------------------------------------

MODULE SUBROUTINE bottom_friction_2d(self,runtype)
   !!

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
   integer, intent(in) :: runtype

!  Local constants

!  Local variables
   integer :: i,j,k
   real(real64) :: hh,ustar
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('bottom_friction_2d()',level=3)
   ! x-direction
   UGrid: associate( UG => self%domain%U )
   do j=UG%jmin,UG%jmax
      do i=UG%imin,UG%imax
         if (UG%mask(i,j) > 0) then
            hh=max(self%domain%Dmin,UG%D(i,j))
            self%rru(i,j)=(kappa/log((UG%z0b(i,j)+0.5_real64*hh)/UG%z0b(i,j)))**2
            self%work2d(i,j)=0.25_real64*(self%v1(i,j)+self%v1(i+1,j)+self%v1(i,j-1)+self%v1(i+1,j-1))
         end if
      end do
   end do
   if (runtype .eq. 1) then
      do j=UG%jmin,UG%jmax
         do i=UG%imin,UG%imax
            if (UG%mask(i,j) > 0) then
               ustar=sqrt(self%rru(i,j)*(self%u1(i,j)**2+self%work2d(i,j)**2))
               hh=max(self%domain%Dmin,UG%D(i,j))
               UG%z0b(i,j)=min(hh,UG%z0b_min(i,j)+0.1_real64*avmmol/max(avmmol,ustar))
               self%rru(i,j)=(kappa/log((UG%z0b(i,j)+0.5_real64*hh)/UG%z0b(i,j)))**2
            end if
         end do
      end do
   end if
   do j=UG%jmin,UG%jmax
      do i=UG%imin,UG%imax
         if (UG%mask(i,j) > 0) then
            self%ru(i,j)=self%rru(i,j)*sqrt(self%u1(i,j)**2+self%work2d(i,j)**2)
         end if
      end do
   end do
   end associate UGrid

   ! y-direction
   VGrid: associate( VG => self%domain%V )
   do j=VG%jmin,VG%jmax
      do i=VG%imin,VG%imax
         if (VG%mask(i,j) > 0) then
            hh=max(self%domain%Dmin,VG%D(i,j))
            self%rrv(i,j)=(kappa/log((VG%z0b(i,j)+0.5_real64*hh)/VG%z0b(i,j)))**2
            self%work2d(i,j)=0.25_real64*(self%u1(i,j)+self%u1(i+1,j)+self%u1(i,j-1)+self%u1(i+1,j-1))
         end if
      end do
   end do
   if (runtype .eq. 1) then
      do j=VG%jmin,VG%jmax
         do i=VG%imin,VG%imax
            if (VG%mask(i,j) > 0) then
               ustar=sqrt(self%rrv(i,j)*(self%work2d(i,j)**2+self%v1(i,j)**2))
               hh=max(self%domain%Dmin,VG%D(i,j))
               VG%z0b(i,j)=min(hh,VG%z0b_min(i,j)+0.1_real64*avmmol/max(avmmol,ustar))
               self%rrv(i,j)=(kappa/log((VG%z0b(i,j)+0.5_real64*hh)/VG%z0b(i,j)))**2
            end if
         end do
      end do
   end if
   do j=VG%jmin,VG%jmax
      do i=VG%imin,VG%imax
         if (VG%mask(i,j) > 0) then
            self%rv(i,j)=self%rrv(i,j)*sqrt(self%work2d(i,j)**2+self%v1(i,j)**2)
         end if
      end do
   end do
   end associate VGrid
END SUBROUTINE bottom_friction_2d

!---------------------------------------------------------------------------

MODULE SUBROUTINE bottom_friction_3d(self)
   !!
   !! @note
   !! is it assumed the velocity profile is resolved - as there is no effort to calculate
   !! 'would be' velocity in half grid distance from bottom
   !! @endnote

   IMPLICIT NONE
!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: i,j,k
   real(real64) :: hh,ustar,r
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('bottom_friction_3d()',level=3)
   UGrid: associate( UG => self%domain%U )
   k=UG%kmin
   do j=UG%jmin,UG%jmax
      do i=UG%imin,UG%imax
         self%work2d(i,j)=0._real64
         if (UG%mask(i,j) > 0) then
            hh=max(self%domain%Dmin/UG%kmax,UG%hn(i,j,k))
            r=(kappa/log((UG%z0b(i,j)+0.5*hh)/UG%z0b(i,j)))**2 ! GETM online report - (127)
            self%work2d(i,j)=0.25_real64*(self%vk(i,j,k)+self%vk(i+1,j,k)+self%vk(i,j-1,k)+self%vk(i+1,j-1,k))
#if 0
            ustar=sqrt(r*(uuloc(i,j)**2+uvloc(i,j)**2))
            UG%z0b(i,j)=min(hh,UG%z0b_min(i,j)+0.1_real64*avmmol/max(avmmol,ustar))
            r=(kappa/log((UG%z0b(i,j)+0.5_real64*hh)/UG%z0b(i,j)))**2
#endif
            self%rru(i,j)=r*sqrt(self%uk(i,j,k)**2+self%work2d(i,j)**2) ! GETM online report - (126)
         end if
      end do
   end do
   end associate UGrid

   VGrid: associate( VG => self%domain%V )
   k=VG%kmin
   do j=VG%jmin,VG%jmax
      do i=VG%imin,VG%imax
         self%work2d(i,j)=0._real64
         if (VG%mask(i,j) > 0) then
            hh=max(self%domain%Dmin/VG%kmax,VG%hn(i,j,k))
            r=(kappa/log((VG%z0b(i,j)+0.5_real64*hh)/VG%z0b(i,j)))**2
            self%work2d(i,j)=0.25_real64*(self%uk(i,j,k)+self%uk(i-1,j,k)+self%uk(i,j+1,k)+self%uk(i-1,j+1,k))
#if 0
            ustar=sqrt(r*(vuloc(i,j)**2+vvloc(i,j)**2))
            VG%z0b(i,j)=min(hh,VG%z0b_min(i,j)+0.1_real64*avmmol/max(avmmol,ustar))
            r=(kappa/log((VG%z0b(i,j)+0.5*hh)/VG%z0b(i,j)))**2
#endif
            self%rrv(i,j)=r*sqrt(self%work2d(i,j)**2+self%vk(i,j,k)**2)
         end if
      end do
   end do
   end associate VGrid
END SUBROUTINE bottom_friction_3d

!---------------------------------------------------------------------------

MODULE SUBROUTINE slow_bottom_friction(self)
   !! Slow terms

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: stat
   integer :: i,j,k
   real(real64) :: hh,uloc,vloc
   real(real64), allocatable :: Ui(:,:),Vi(:,:)
!---------------------------------------------------------------------------
   if(associated(self%logs)) call self%logs%info('slow_bottom_friction()',level=3)

!KB   allocate(ruu, mold=self%U, stat=stat)
!KB   allocate(rvv, mold=self%V, stat=stat)
!KB   allocate(ru, mold=self%U, stat=stat)
!KB   allocate(rv, mold=self%V, stat=stat)

   allocate(Ui, mold=self%U, stat=stat)
   allocate(Vi, mold=self%V, stat=stat)

   UGrid: associate( UG => self%domain%U )
   VGrid: associate( VG => self%domain%V )

   do j=UG%l(2),UG%u(2)
      do i=UG%l(1),UG%u(1)
         if (UG%mask(i,j) .ge. 1) then
            Ui(i,j)=self%Uio(i,j)/(UG%zio(i,j)+UG%H(i,j))
         else
            Ui(i,j)=0._real64
         end if
      end do
   end do

   do j=VG%l(2),VG%u(2)
      do i=VG%l(1),VG%u(1)
         if (VG%mask(i,j) .ge. 1) then
            Vi(i,j)=self%Vio(i,j)/(VG%zio(i,j)+VG%H(i,j))
         else
            Vi(i,j)=0._real64
         end if
      end do
   end do

   do j=UG%l(2),UG%u(2)
      do i=UG%l(1),UG%u(1)
         self%ru(i,j)=0._real64
         self%work2d(i,j)=0._real64
         if (UG%mask(i,j) .ge. 1) then
            hh=max(self%domain%Dmin,UG%zin(i,j)+UG%H(i,j))
            self%rru(i,j)=(kappa/log((UG%z0b(i,j)+0.5_real64*HH)/UG%z0b(i,j)))**2
            vloc = 0.25_real64*(Vi(i,j)+Vi(i+1,j)+Vi(i,j-1)+Vi(i+1,j-1))
            self%ru(i,j) = self%rru(i,j)*sqrt(Ui(i,j)**2+vloc**2)
         end if
      end do
   end do

   do j=VG%l(2),VG%u(2)
      do i=VG%l(1),VG%u(1)
         self%rv(i,j)=0._real64
         self%work2d(i,j)=0._real64
         if (VG%mask(i,j) .ge. 1) then
            ! why not VG%D !KB
            hh=max(self%domain%Dmin,VG%zin(i,j)+VG%H(i,j))
            self%rrv(i,j)=(kappa/log((VG%z0b(i,j)+0.5_real64*hh)/VG%z0b(i,j)))**2
            uloc = 0.25_real64*(Ui(i,j)+Ui(i+1,j)+Ui(i,j-1)+Ui(i+1,j-1))
            self%rv(i,j) = self%rrv(i,j)*sqrt(uloc**2+Vi(i,j)**2)
         end if
      end do
   end do

   k=UG%kmin
   do j=UG%jmin,UG%jmax
      do i=UG%imin,UG%imax
         self%SxF(i,j)= 0._real64
         if (UG%mask(i,j) .ge. 1) then
!KB            if (UG%kmax .gt. 1) then
               self%SxF(i,j)=-self%ru(i,j)*self%Ui(i,j)/(0.5_real64*(UG%zio(i,j)+UG%zin(i,j))+UG%H(i,j)) &
                             +self%rru(i,j)*self%pk(i,j,k)/(0.5_real64*(UG%ho(i,j,k)+UG%hn(i,j,k)))
!KB            end if
         end if
      end do
   end do

   k=VG%kmin
   do j=VG%jmin,VG%jmax
      do i=VG%imin,VG%imax
         self%SyF(i,j)= 0._real64
         if (VG%mask(i,j) .ge. 1) then
!KB            if (VG%kmax .gt. 1) then
               self%SyF(i,j)=-self%rv(i,j)*self%Vi(i,j)/(0.5_real64*(VG%zio(i,j)+VG%zin(i,j))+VG%H(i,j)) &
                             +self%rrv(i,j)*self%qk(i,j,k)/(0.5_real64*(VG%ho(i,j,k)+VG%hn(i,j,k)))
!KB            end if
         end if
      end do
   end do
   end associate VGrid
   end associate UGrid
END SUBROUTINE slow_bottom_friction

!---------------------------------------------------------------------------

END SUBMODULE bottom_friction_smod

#if 0
#ifdef STRUCTURE_FRICTION
               do k=1,kmax
                  Slru(i,j)=Slru(i,j)+pk(i,j,k)*_HALF_*(sf(i,j,k)+sf(i+1,j,k))
               end do
#endif
#ifdef STRUCTURE_FRICTION
               do k=1,kmax
                  Slrv(i,j)=Slrv(i,j)+pk(i,j,k)*_HALF_*(sf(i,j,k)+sf(i+1,j,k))
               end do
#endif
#endif

#if 0
MODULE PROCEDURE slow_bottom_friction

   IMPLICIT NONE

!  Local constants

!  Local variables
   integer :: i,j,k
   real(real64) :: uloc,vloc,HH
   real(real64), dimension(:,:), allocatable :: ruu,rvv,Ui,Vi
   integer :: stat
!---------------------------------------------------------------------------
   if(associated(self%logs)) call self%logs%info('slow_bottom_friction()',level=2)

   allocate(ruu, mold=self%U, stat=stat)
   allocate(rvv, mold=self%V, stat=stat)
   allocate(Ui, mold=self%U, stat=stat)
   allocate(Vi, mold=self%V, stat=stat)
   self%zub = 0.01_real64
   self%zvb = 0.01_real64

   UGrid: associate( UG => self%domain%U )
   VGrid: associate( VG => self%domain%V )

   do j=UG%l(2),UG%u(2)
      do i=UG%l(1),UG%u(1)
         if (UG%mask(i,j) .ge. 1) then
            Ui(i,j)=self%Uio(i,j)/(UG%zio(i,j)+UG%H(i,j))
         else
            Ui(i,j)=0._real64
         end if
      end do
   end do

   do j=VG%l(2),VG%u(2)
      do i=VG%l(1),VG%u(1)
         if (VG%mask(i,j) .ge. 1) then
            Vi(i,j)=self%Vio(i,j)/(VG%zio(i,j)+VG%H(i,j))
         else
            Vi(i,j)=0._real64
         end if
      end do
   end do

   do j=UG%l(2),UG%u(2)
      do i=UG%l(1),UG%u(1)
         if (UG%mask(i,j) .ge. 1) then
            HH=max(UG%zin(i,j)+UG%H(i,j),self%domain%Dmin)
            ruu(i,j)=(self%zub(i,j)+0.5_real64*HH)/self%zub(i,j)
#if 0
            if (ruu(i,j) .le. 1._real64) then
               LEVEL1 i,j,UG%ssuo(i,j)
               LEVEL1 'Bottom xfriction coefficient infinite.'
               FATAL 'slow_bottom_friction()'
!KB               call getm_error('slow_bottom_friction()','Bottom xfriction coefficient infinite.')
               STOP ! Just if getm_error doesnt actually halt
            end if
#endif
            ruu(i,j)=(kappa/log(ruu(i,j)))**2
            vloc = 0.25_real64*(Vi(i,j)+Vi(i+1,j)+Vi(i,j-1)+Vi(i+1,j-1))
            self%ru(i,j) = ruu(i,j)*sqrt(Ui(i,j)**2+vloc**2)
         else
            self%ru(i,j)=0._real64
         end if
      end do
   end do

   do j=VG%l(2),VG%u(2)
      do i=VG%l(1),VG%u(1)
         if (VG%mask(i,j) .ge. 1) then
            ! why not VG%D !KB
            HH=max(VG%zin(i,j)+VG%H(i,j),self%domain%Dmin)
            ruu(i,j)=(self%zvb(i,j)+0.5_real64*HH)/self%zvb(i,j)
#if 0
            if (ruu(i,j) .le. 1._real64) then
               LEVEL1 i,j,UG%ssuo(i,j)
               LEVEL1 'Bottom xfriction coefficient infinite.'
               FATAL 'slow_bottom_friction()'
!KB               call getm_error('slow_bottom_friction()','Bottom xfriction coefficient infinite.')
               STOP ! Just if getm_error doesnt actually halt
            end if
#endif
            ruu(i,j)=(kappa/log(ruu(i,j)))**2
            uloc = 0.25_real64*(Ui(i,j)+Ui(i+1,j)+Ui(i,j-1)+Ui(i+1,j-1))
            self%rv(i,j) = ruu(i,j)*sqrt(uloc**2+Vi(i,j)**2)
         else
            self%rv(i,j)=0._real64
         end if
      end do
   end do

   end associate VGrid
   end associate UGrid
END PROCEDURE slow_bottom_friction

#endif

#if 0
   do j=UG%jmin,UG%jmax
      do i=UG%imin,UG%imax
         if (UG%mask(i,j) .ge. 1) then
            !k=kumin(i,j)
            k=1
            if (UG%kmax .gt. 1) then
!               self%Slru(i,j)=-self%ru(i,j)*self%Ui(i,j) &
               self%SxF(i,j)=-self%ru(i,j)*self%Ui(i,j) &
                             /(0.5_real64*(UG%zio(i,j)+UG%zin(i,j))+UG%H(i,j)) &
                             +self%rru(i,j)*self%pk(i,j,k) &
                             /(0.5_real64*(UG%ho(i,j,k)+UG%hn(i,j,k)))
#ifdef STRUCTURE_FRICTION
               do k=1,kmax
                  Slru(i,j)=Slru(i,j)+pk(i,j,k)*_HALF_*(sf(i,j,k)+sf(i+1,j,k))
               end do
#endif
            else
               self%SxF(i,j)= 0._real64
            end if
         end if
      end do
   end do
   end associate UGrid

      do i=VG%imin,VG%imax
         if (VG%mask(i,j) .ge. 1) then
            !k=kumin(i,j)
            k=1
            if (VG%kmax .gt. 1) then
!               self%Slrv(i,j)=-self%rv(i,j)*self%Vi(i,j) &
               self%SyF(i,j)=-self%rv(i,j)*self%Vi(i,j) &
                             /(0.5_real64*(VG%zio(i,j)+VG%zin(i,j))+VG%H(i,j)) &
                             +self%rrv(i,j)*self%qk(i,j,k) &
                             /(0.5_real64*(VG%ho(i,j,k)+VG%hn(i,j,k)))
#ifdef STRUCTURE_FRICTION
               do k=1,kmax
                  Slrv(i,j)=Slrv(i,j)+pk(i,j,k)*_HALF_*(sf(i,j,k)+sf(i+1,j,k))
               end do
#endif
            else
               self%SyF(i,j)= 0._real64
            end if
         end if
      end do
   end do
#endif
