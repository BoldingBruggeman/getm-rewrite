! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!! @note
!! ddl and ddu == 0
!! zc and zf could be calculated in here - careful about definitions and signs of zc and zf
!! @endnote

SUBMODULE (getm_domain : vertical_coordinates_smod) vertical_gvc_smod

!  Module types and variables
   real(real64), dimension(:,:,:), allocatable :: gga

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

MODULE SUBROUTINE init_gvc(self)
   !! Initializing the general vertical coordinates

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
   if (associated(self%logs)) call self%logs%info('init_vertical_gvc()',level=3)

   allocate(be(0:self%T%kmax),stat=stat)     ! dimensionless beta-coordinate
   if (stat /= 0) STOP 'coordinates: Error allocating (be)'
   allocate(ga(0:self%T%kmax),stat=stat)
   if (stat /= 0) stop 'coordinates: Error allocating (ga)'
   allocate(sig(0:self%T%kmax),stat=stat)    ! dimensionless sigma-coordinate
   if (stat /= 0) STOP 'coordinates: Error allocating (sig)'
   do k=0,self%T%kmax
      ga(k) = k
   end do
   call mm_s('gga',gga,self%T%l+(/0,0,-1/),self%T%u,stat=stat)
   if (stat /= 0) stop 'coordinates: Error allocating memory (gga)'
   be(0)=  -1._real64
   sig(0)= -1._real64
   if (self%ddu < 0._real64) self%ddu=0._real64
   if (self%ddl < 0._real64) self%ddl=0._real64

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
   TGrid: associate( TG => self%T )
   do j=TG%l(2),TG%u(2)
      do i=TG%l(1),TG%u(1)
         HH=max(TG%zio(i,j)+TG%H(i,j),self%Dmin)
         alpha=min(((be(kk)-be(kk-1))-self%Dgamma/HH*(sig(kk)-sig(kk-1))) &
                  /((be(kk)-be(kk-1))-(sig(kk)-sig(kk-1))),1._real64)
         gga(i,j,0)=-1._real64
         do k=1,TG%kmax
            gga(i,j,k)=alpha*sig(k)+(1._real64-alpha)*be(k)
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

!  Here, the initial layer distribution is calculated.
   do k=1,TG%kmax
      do j=TG%l(2),TG%u(2)
         do i=TG%l(1),TG%u(1)
            if (TG%mask(i,j) > 0) then
            HH=max(TG%zio(i,j)+TG%H(i,j),self%Dmin)
            TG%hn(i,j,k)=HH*(gga(i,j,k)-gga(i,j,k-1))
            end if
         end do
      end do
   end do
   end associate TGrid

   UGrid: associate( UG => self%U )
   do k=1,UG%kmax
      do j=UG%l(2),UG%u(2)
         do i=UG%l(1),UG%u(1)-1
            if (UG%mask(i,j) > 0) then
            HH=max(UG%zio(i,j)+UG%H(i,j),self%Dmin)
            UG%ho(i,j,k)=HH*0.5*(gga(i,j,k)-gga(i,j,k-1)+gga(i+1,j,k)-gga(i+1,j,k-1))
            UG%hn(i,j,k)=UG%ho(i,j,k)
            end if
         end do
      end do
   end do
   end associate UGrid

   VGrid: associate( VG => self%V )
   do k=1,VG%kmax
      do j=VG%l(2),VG%u(2)-1
         do i=VG%l(1),VG%u(1)
            if (VG%mask(i,j) > 0) then
            HH=max(VG%zio(i,j)+VG%H(i,j),self%Dmin)
            VG%ho(i,j,k)=HH*0.5*(gga(i,j,k)-gga(i,j,k-1)+gga(i,j+1,k)-gga(i,j+1,k-1))
            VG%hn(i,j,k)=VG%ho(i,j,k)
            end if
         end do
      end do
   end do
   end associate VGrid
END SUBROUTINE init_gvc

!---------------------------------------------------------------------------

MODULE SUBROUTINE do_gvc(self,dt)
   !! A wrapper for vertical coordinate calculations

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_domain), intent(inout) :: self
   real(real64), intent(in) :: dt

!  Local constants

!  Local variables
   real(real64) :: cord_relax=0.0_real64
   real(real64) :: HH
   real(real64) :: zz
   real(real64) :: r
   integer :: i,j,k
   real(real64) :: Lstart,Lstop
!-----------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('do_gvc()',level=3)
write(*,*) 'do_gvc - has not been tested yet'

   !! why not ho=hn as zio=zin
   !! why not use total water depth %T%D
!KB   self%T%ho=self%T%hn
   TGrid: associate( TG => self%T )
   call cpu_time(Lstart)
write(*,*) cord_relax,dt,self%Dmax
stop
   do j=TG%l(2),TG%u(2)
      do i=TG%l(1),TG%u(1)
         if (TG%mask(i,j) > 0) then
            r=cord_relax/dt*TG%H(i,j)/self%Dmax
            HH=TG%zin(i,j)+TG%H(i,j)
            if (HH .lt. self%Dgamma) then
               do k=1,TG%kmax
                  TG%ho(i,j,k)=TG%hn(i,j,k)
                  TG%hn(i,j,k)=HH/TG%kmax
               end do
            else
               zz=-TG%H(i,j)
               do k=1,TG%kmax-1
                  TG%ho(i,j,k)=TG%hn(i,j,k)
                  TG%hn(i,j,k)=(TG%ho(i,j,k)*r+HH*(gga(i,j,k)-gga(i,j,k-1)))/(r+1.)
                  zz=zz+TG%hn(i,j,k)
               end do
               TG%ho(i,j,TG%kmax)=TG%hn(i,j,TG%kmax)
               TG%hn(i,j,TG%kmax)=TG%zin(i,j)-zz
            end if
         end if
      end do
   end do
   call cpu_time(Lstop)
   write(57,*) Lstop-Lstart
   end associate TGrid

   !! if gga, zin and H are updated in halo zones - extend to all domain
   !! what about mask
   UGrid: associate( UG => self%U )
   do j=UG%l(2),UG%u(2)
      do i=UG%l(1),UG%u(1)-1
         if (UG%mask(i,j) > 0) then
            r=cord_relax/dt*UG%H(i,j)/self%Dmax
            zz=-UG%H(i,j)
            HH=UG%zin(i,j)+UG%H(i,j)
            do k=1,UG%kmax-1
               UG%ho(i,j,k)=UG%hn(i,j,k)
               UG%hn(i,j,k)=(UG%ho(i,j,k)*r+HH*0.5*(gga(i  ,j,k)-gga(i  ,j,k-1) &
                                                   +gga(i+1,j,k)-gga(i+1,j,k-1)))/(r+1._real64)
               zz=zz+UG%hn(i,j,k)
            end do
            UG%ho(i,j,UG%kmax)=UG%hn(i,j,UG%kmax)
            UG%hn(i,j,UG%kmax)=UG%zin(i,j)-zz
         end if
      end do
   end do
   end associate UGrid

   !! if gga, zin and H are updated in halo zones - extend to all domain
!KB   TG%ho=TG%hn
   VGrid: associate( VG => self%V )
   do j=VG%l(2),VG%u(2)-1
      do i=VG%l(1),VG%u(1)
         if (VG%mask(i,j) > 0) then
            r=cord_relax/dt*VG%H(i,j)/self%Dmax
            zz=-VG%H(i,j)
            HH=VG%zin(i,j)+VG%H(i,j)
            do k=1,VG%kmax-1
               VG%ho(i,j,k)=VG%hn(i,j,k)
               VG%hn(i,j,k)=(VG%ho(i,j,k)*r+HH*0.5*(gga(i,j  ,k)-gga(i,j  ,k-1) &
                                                   +gga(i,j+1,k)-gga(i,j+1,k-1)))/(r+1._real64)
            HH=max(VG%zio(i,j)+VG%H(i,j),self%Dmin)
            VG%ho(i,j,k)=HH*0.5*(gga(i,j,k)-gga(i,j,k-1)+gga(i,j+1,k)-gga(i,j+1,k-1))
            VG%hn(i,j,k)=VG%ho(i,j,k)
               zz=zz+VG%hn(i,j,k)
            end do
            VG%ho(i,j,VG%kmax)=VG%hn(i,j,VG%kmax)
            VG%hn(i,j,VG%kmax)=VG%zin(i,j)-zz
         end if
      end do
   end do
   end associate VGrid
END SUBROUTINE do_gvc

!---------------------------------------------------------------------------

END SUBMODULE vertical_gvc_smod
