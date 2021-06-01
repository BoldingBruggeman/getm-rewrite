! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

SUBMODULE (getm_domain : vertical_coordinates_smod) vertical_sigma_smod

!  Module types and variables
   real(real64), dimension(:), allocatable  :: dga

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

MODULE SUBROUTINE init_sigma(self)
   !! A wrapper for vertical coordinate calculations

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_domain), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: stat
   real(real64), dimension(:), allocatable  :: ga
   integer :: k
!-----------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('init_sigma()',level=3)

   allocate(dga(1:self%T%kmax),stat=stat)
   if (self%ddl <= 0._real64 .and. self%ddu <= 0._real64) then
      ! Equidistant sigma coordinates
      dga(:) = 1._real64/self%T%kmax
   else
      ! Non-equidistant sigma coordinates
      if (self%ddu < 0._real64) self%ddu=0._real64
      if (self%ddl < 0._real64) self%ddl=0._real64
      allocate(ga(0:self%T%kmax),stat=stat)
      ga(0)= -1._real64
      do k=1,self%T%kmax
         ! This zooming routine is from Antoine Garapon, ICCH, DK
         ga(k)=tanh((self%ddl+self%ddu)*k/float(self%T%kmax)-self%ddl)+tanh(self%ddl)
         ga(k)=ga(k)/(tanh(self%ddl)+tanh(self%ddu)) - 1._real64
         dga(k)=ga(k)-ga(k-1)
      end do
   end if
END SUBROUTINE init_sigma

!---------------------------------------------------------------------------

MODULE SUBROUTINE do_sigma(self)
   !! A wrapper for vertical coordinate calculations

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_domain), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: i,j,k
   integer :: stat
!-----------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('do_sigma()',level=3)

   !! why not ho=hn as zio=zin
   TGrid: associate( TG => self%T )
!KB
   TG%ho=TG%hn
   do j=TG%l(2),TG%u(2)
      do i=TG%l(1),TG%u(1)
         if (TG%mask(i,j) > 0) then
!KB            TG%ho(i,j,:)=(TG%zio(i,j)+TG%H(i,j))*dga(:)
            TG%hn(i,j,:)=(TG%zin(i,j)+TG%H(i,j))*dga(:)
         end if
      end do
   end do
   end associate TGrid

   !! why not ho=hn as zio=zin
   !! if zin and H are updated in halo zones - extend to all domain
   !! what about mask
   UGrid: associate( UG => self%U )
!KB
   UG%ho=UG%hn
   do j=UG%l(2),UG%u(2)
      do i=UG%l(1),UG%u(1) ! requires zin is HALO-updated (-1)
         if (UG%mask(i,j) > 0) then
!KB            UG%ho(i,j,:)=(UG%zio(i,j)+UG%H(i,j))*dga(:)
            UG%hn(i,j,:)=(UG%zin(i,j)+UG%H(i,j))*dga(:)
         end if
      end do
   end do
   end associate UGrid

   !! if zin and H are updated in halo zones - extend to all domain
   VGrid: associate( VG => self%V )
!KB
   VG%ho=VG%hn
   do j=VG%l(2),VG%u(2)  ! requires zin is HALO-updated (-1)
      do i=VG%l(1),VG%u(1)
         if (VG%mask(i,j) > 0) then
!KB            VG%ho(i,j,:)=(VG%zio(i,j)+VG%H(i,j))*dga(:)
            VG%hn(i,j,:)=(VG%zin(i,j)+VG%H(i,j))*dga(:)
         end if
      end do
   end do
   end associate VGrid
END SUBROUTINE do_sigma

!---------------------------------------------------------------------------

END SUBMODULE vertical_sigma_smod
