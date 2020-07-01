! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!! @note
!! ddl and ddu == 0
!! @endnote

SUBMODULE (getm_domain : vertical_coordinates_smod) vertical_gvc_smod

!  Module types and variables
   real(real64), dimension(:,:,:), allocatable :: gga

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

module SUBROUTINE init_gvc(self)
   !! A wrapper for vertical coordinate calculations

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_domain), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: stat
   real(real64), dimension(:), allocatable :: ga
   real(real64), dimension(:), allocatable :: be
   real(real64), dimension(:), allocatable :: sig
   real(real64) :: HH, alpha
   integer :: i,j,k,kk
!-----------------------------------------------------------------------------
   call self%logs%info('init_vertical_gvc()',level=3)
   allocate(be(0:self%T%kmax),stat=stat)     ! dimensionless beta-coordinate
   if (stat /= 0) STOP 'coordinates: Error allocating (be)'
   allocate(ga(0:self%T%kmax),stat=stat)
   if (stat /= 0) stop 'coordinates: Error allocating (ga)'
   allocate(sig(0:self%T%kmax),stat=stat)    ! dimensionless sigma-coordinate
   if (stat /= 0) STOP 'coordinates: Error allocating (sig)'
   do k=0,self%T%kmax
      ga(k) = k
   end do
!   call mm_s('gga',gga,self%T%l(1:3),self%T%u(1:3),stat=stat)
   allocate(gga(self%T%imax,self%T%jmax,0:self%T%kmax),stat=stat)    ! dimensionless sigma-coordinate
   if (stat /= 0) stop 'coordinates: Error allocating memory (gga)'
   be(0)=  -1._real64
   sig(0)= -1._real64
   if (self%ddu .lt. 0._real64) self%ddu=0._real64
   if (self%ddl .lt. 0._real64) self%ddl=0._real64
   !KB
!   if (self%ddu .le. 0._real64) self%ddu=0.1_real64
!   if (self%ddl .le. 0._real64) self%ddl=0.2_real64

   do k=1,self%T%kmax
      be(k)=tanh((self%ddl+self%ddu)*k/float(self%T%kmax)-self%ddl)+tanh(self%ddl)
!      be(k)=be(k)/(tanh(self%ddl)+tanh(self%ddu))-1._real64
      sig(k)=k/float(self%T%kmax)-1._real64
   end do
   if (self%gamma_surf) then
      kk=self%T%kmax
   else
      kk=1
   end if
   do j=self%T%l(2),self%T%u(2)
      do i=self%T%l(1),self%T%u(1)
         HH=max(self%T%sseo(i,j)+self%T%H(i,j),self%Dmin)
         alpha=min(&
                  ((be(kk)-be(kk-1))-self%Dgamma/HH &
                   *(sig(kk)-sig(kk-1))) &
                   /((be(kk)-be(kk-1))-(sig(kk)-sig(kk-1))),1._real64)
         gga(i,j,0)=-1._real64
         do k=1,self%T%kmax
            gga(i,j,k)=alpha*sig(k)+(1.-alpha)*be(k)
!            if (gga(i,j,k) .lt. gga(i,j,k-1)) then
!               STDERR kk,(be(kk)-be(kk-1)),(sig(kk)-sig(kk-1))
!               STDERR Dgamma,HH
!               STDERR alpha
!               STDERR k-1,gga(i,j,k-1),be(k-1),sig(k-1)
!               STDERR k,gga(i,j,k),be(k),sig(k)
!               stop 'coordinates'
!            end if
         end do
      end do
   end do

   write(*,*) self%ddl,self%ddu
   write(*,*) tanh(0._real64)
   write(*,*) be
   write(*,*) sig
!   write(*,*) gga
!   stop
!  Here, the initial layer distribution is calculated.
#define TG self%T
   do k=1,TG%kmax
      do j=TG%l(2),TG%u(2)
         do i=TG%l(1),TG%u(1)
            if (TG%mask(i,j) > 0) then
            HH=max(TG%sseo(i,j)+TG%H(i,j),self%Dmin)
            write(*,*) HH
            TG%hn(i,j,k)=HH*(gga(i,j,k)-gga(i,j,k-1))
            end if
         end do
      end do
   end do
#undef TG

#define UG self%U
   do k=1,UG%kmax
      do j=UG%l(2),UG%u(2)
         do i=UG%l(1),UG%u(1)-1
            if (UG%mask(i,j) > 0) then
            HH=max(UG%sseo(i,j)+UG%H(i,j),self%Dmin)
            UG%ho(i,j,k)=HH*0.5*            &
             (gga(i,j,k)-gga(i,j,k-1)+gga(i+1,j,k)-gga(i+1,j,k-1))
            UG%hn(i,j,k)=UG%ho(i,j,k)
            end if
         end do
      end do
   end do
#undef UG

#define VG self%V
   do k=1,VG%kmax
      do j=VG%l(2),VG%u(2)-1
         do i=VG%l(1),VG%u(1)
            if (VG%mask(i,j) > 0) then
            HH=max(VG%sseo(i,j)+VG%H(i,j),self%Dmin)
            VG%ho(i,j,k)=HH*0.5*            &
             (gga(i,j,k)-gga(i,j,k-1)+gga(i,j+1,k)-gga(i,j+1,k-1))
            VG%hn(i,j,k)=VG%ho(i,j,k)
            end if
         end do
      end do
   end do
#undef VG

   deallocate(be)
   deallocate(ga)
   deallocate(sig)
   return
END SUBROUTINE init_gvc

!---------------------------------------------------------------------------

module SUBROUTINE do_gvc(self,dt)
   !! A wrapper for vertical coordinate calculations

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_domain), intent(inout) :: self
   real(real64), intent(in) :: dt

!  Local constants

!  Local variables
   integer :: i,j,k
   integer :: stat
   real(real64) :: cord_relax !!!!!!!!KB
   real(real64) :: HH
   real(real64) :: zz
   real(real64) :: r
!-----------------------------------------------------------------------------
   call self%logs%info('do_gvc()',level=3)
   write(*,*) 'aaaaaa'
   return

   !! why not ho=hn as sseo=ssen
   !! why not use total water depth %T%D
!KB   self%T%ho=self%T%hn
   do j=self%T%l(2),self%T%u(2)
      do i=self%T%l(1),self%T%u(1)
         if (self%T%mask(i,j) > 0) then
            r=cord_relax/dt*self%T%H(i,j)/self%Dmax
            HH=self%T%ssen(i,j)+self%T%H(i,j)
            if (HH .lt. self%Dgamma) then
               do k=1,self%T%kmax
                  self%T%ho(i,j,k)=self%T%hn(i,j,k)
                  self%T%hn(i,j,k)=HH/self%T%kmax
               end do
            else
               zz=-self%T%H(i,j)
               do k=1,self%T%kmax-1
                  self%T%ho(i,j,k)=self%T%hn(i,j,k)
                  self%T%hn(i,j,k)=(self%T%ho(i,j,k)*r+HH*(gga(i,j,k)-gga(i,j,k-1)))/(r+1.)
                  zz=zz+self%T%hn(i,j,k)
               end do
               self%T%ho(i,j,self%T%kmax)=self%T%hn(i,j,self%T%kmax)
               self%T%hn(i,j,self%T%kmax)=self%T%ssen(i,j)-zz
            end if
         end if
      end do
   end do


   !! if gga, ssen and H are updated in halo zones - extend to all domain
   !! what about mask
!KB   self%U%ho=self%U%hn
   do j=self%U%l(2),self%U%u(2)
      do i=self%U%l(1),self%U%u(1)-1
         if (self%U%mask(i,j) > 0) then
!KBK         if (au(i,j) .gt. 0) then
         r=cord_relax/dt*self%U%H(i,j)/self%Dmax
         zz=-self%U%H(i,j)
         HH=self%U%ssen(i,j)+self%U%H(i,j)
         do k=1,self%U%kmax-1
            self%U%ho(i,j,k)=self%U%hn(i,j,k)
            self%U%hn(i,j,k)=(self%U%ho(i,j,k)*r+HH*0.5*(gga(i,j,k)-gga(i,j,k-1) &
                      +gga(i+1,j,k)-gga(i+1,j,k-1)))/(r+1.)
            zz=zz+self%U%hn(i,j,k)
         end do
         self%U%ho(i,j,self%U%kmax)=self%U%hn(i,j,self%U%kmax)
         self%U%hn(i,j,self%U%kmax)=self%U%ssen(i,j)-zz
!KBK         end if
         end if
      end do
   end do

   !! if gga, ssen and H are updated in halo zones - extend to all domain
!KB   self%V%ho=self%V%hn
   do j=self%V%l(2),self%V%u(2)-1
      do i=self%V%l(1),self%V%u(1)
         if (self%V%mask(i,j) > 0) then
!KBK         if (av(i,j).gt.0) then
         r=cord_relax/dt*self%V%H(i,j)/self%Dmax
         zz=-self%V%H(i,j)
         HH=self%V%ssen(i,j)+self%V%H(i,j)
         do k=1,self%V%kmax-1
            self%V%ho(i,j,k)=self%V%hn(i,j,k)
            self%V%hn(i,j,k)=(self%V%ho(i,j,k)*r+HH*0.5*(gga(i,j,k)-gga(i,j,k-1) &
                      +gga(i,j+1,k)-gga(i,j+1,k-1)))/(r+1.)
            zz=zz+self%V%hn(i,j,k)
         end do
         self%V%ho(i,j,self%V%kmax)=self%V%hn(i,j,self%V%kmax)
         self%V%hn(i,j,self%V%kmax)=self%V%ssen(i,j)-zz
!KBK         end if
         end if
      end do
   end do
   return
END SUBROUTINE do_gvc

!---------------------------------------------------------------------------

END SUBMODULE vertical_gvc_smod
