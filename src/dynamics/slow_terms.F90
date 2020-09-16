! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard
!! @note
!! check with GETM
!! @endnote

SUBMODULE (getm_momentum) slow_terms_smod

CONTAINS

!---------------------------------------------------------------------------

module SUBROUTINE slow_terms(self,idpdx,idpdy)
   !! Slow terms

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
   real(real64), dimension(:,:,:), intent(in) :: idpdx
   real(real64), dimension(:,:,:), intent(in) :: idpdy

!  Local constants

!  Local variables
   integer :: i,j,k
!KB   real(real64) :: vertsum,ip_fac
!---------------------------------------------------------------------------
   call self%logs%info('slow_terms()',level=2)

   UGrid: associate( UG => self%domain%U )
   do j=UG%jmin,UG%jmax
      do i=UG%imin,UG%imax
         if (UG%mask(i,j) .ge. 1) then
            self%SlUx(i,j)=-self%UEx(i,j)+SUM(self%uuEx(i,j,1:))-SUM(idpdx(i,j,1:))
            self%SxA(i,j)=0._real64 ! KB
            self%SxB(i,j)=-SUM(idpdx(i,j,1:))
            self%SxD(i,j)=SUM(self%uuEx(i,j,1:))
         end if
      end do
   end do

   do j=UG%jmin,UG%jmax
      do i=UG%imin,UG%imax
         if (UG%mask(i,j) .ge. 1) then
            !k=kumin(i,j)
            k=1
            if (UG%kmax .gt. 1) then
!               self%Slru(i,j)=-self%ru(i,j)*self%Ui(i,j) &
               self%SxF(i,j)=-self%ru(i,j)*self%Ui(i,j) &
                             /(0.5_real64*(UG%sseo(i,j)+UG%ssen(i,j))+UG%H(i,j)) &
                             +self%rru(i,j)*self%pk(i,j,k) &
                             /(0.5_real64*(UG%ho(i,j,k)+UG%hn(i,j,k)))
#ifdef STRUCTURE_FRICTION
               do k=1,kmax
                  Slru(i,j)=Slru(i,j)+pk(i,j,k)*_HALF_*(sf(i,j,k)+sf(i+1,j,k))
               end do
#endif
            else
               self%Slru(i,j)= 0._real64
            end if
         end if
      end do
   end do
   end associate UGrid

   VGrid: associate( VG => self%domain%V )
   do j=VG%jmin,VG%jmax
      do i=VG%imin,VG%imax
         if (VG%mask(i,j) .ge. 1) then
            self%SlVx(i,j)=-self%VEx(i,j)+SUM(self%vvEx(i,j,1:))-SUM(idpdy(i,j,1:))
            self%SyA(i,j)=0._real64 ! KB
            self%SyB(i,j)=-SUM(idpdy(i,j,1:))
            self%SyD(i,j)=SUM(self%vvEx(i,j,1:))
         end if
      end do
   end do

   do j=VG%jmin,VG%jmax
      do i=VG%imin,VG%imax
         if (VG%mask(i,j) .ge. 1) then
            !k=kumin(i,j)
            k=1
            if (VG%kmax .gt. 1) then
!               self%Slrv(i,j)=-self%rv(i,j)*self%Vi(i,j) &
               self%SyF(i,j)=-self%rv(i,j)*self%Vi(i,j) &
                             /(0.5_real64*(VG%sseo(i,j)+VG%ssen(i,j))+VG%H(i,j)) &
                             +self%rrv(i,j)*self%qk(i,j,k) &
                             /(0.5_real64*(VG%ho(i,j,k)+VG%hn(i,j,k)))
#ifdef STRUCTURE_FRICTION
               do k=1,kmax
                  Slrv(i,j)=Slrv(i,j)+pk(i,j,k)*_HALF_*(sf(i,j,k)+sf(i+1,j,k))
               end do
#endif
            else
               self%Slrv(i,j)= 0._real64
            end if
         end if
      end do
   end do
   end associate VGrid

END SUBROUTINE slow_terms

!---------------------------------------------------------------------------

module SUBROUTINE slow_bottom_friction(self)

   !!

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: i,j,k
   real(real64) :: uloc,vloc,HH
   real(real64), dimension(:,:), allocatable :: ruu,rvv,Ui,Vi
   integer :: stat
!---------------------------------------------------------------------------
   call self%logs%info('slow_bottom_friction()',level=2)

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
            Ui(i,j)=self%Uio(i,j)/(UG%sseo(i,j)+UG%H(i,j))
         else
            Ui(i,j)=0._real64
         end if
      end do
   end do

   do j=VG%l(2),VG%u(2)
      do i=VG%l(1),VG%u(1)
         if (VG%mask(i,j) .ge. 1) then
            Vi(i,j)=self%Vio(i,j)/(VG%sseo(i,j)+VG%H(i,j))
         else
            Vi(i,j)=0._real64
         end if
      end do
   end do

   do j=UG%l(2),UG%u(2)
      do i=UG%l(1),UG%u(1)
         if (UG%mask(i,j) .ge. 1) then
            HH=max(UG%ssen(i,j)+UG%H(i,j),self%domain%Dmin)
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
            HH=max(VG%ssen(i,j)+VG%H(i,j),self%domain%Dmin)
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

END SUBROUTINE slow_bottom_friction

!---------------------------------------------------------------------------

END SUBMODULE slow_terms_smod
