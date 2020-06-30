! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!SUBMODULE (getm_domain) vertical_coordinates_smod
SUBMODULE (getm_domain : vertical_coordinates_smod) vertical_sigma_smod

!  Module types and variables
   real(real64), dimension(:), allocatable  :: dga

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

module SUBROUTINE init_sigma(self)
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
   call self%logs%info('init_sigma()',level=3)
   allocate(dga(1:self%T%kmax),stat=stat)
   if (self%ddl .le. 0._real64 .and. self%ddu .le. 0._real64) then
#if 0
      ga(0) = -1._real64
      do k=1,self%T%kmax
         ga(k) = ga(k-1) + 1._real64/self%T%kmax
      end do
      ga(self%T%kmax) = 0._real64
#endif
      dga(:) = 1._real64/self%T%kmax
      write(*,*) 'aaaaa', dga
   else
      allocate(ga(0:self%T%kmax),stat=stat)
      ! Non-equidistant sigma coordinates
      ! This zooming routine is from Antoine Garapon, ICCH, DK
      if (self%ddu .lt. 0._real64) self%ddu=0._real64
      if (self%ddl .lt. 0._real64) self%ddl=0._real64
      allocate(dga(self%T%kmax),stat=stat)
      if (stat /= 0) STOP 'coordinates: Error allocating (dga)'
      ga(0)= -1._real64
      do k=1,self%T%kmax
         ga(k)=tanh((self%ddl+self%ddu)*k/float(self%T%kmax)-self%ddl)+tanh(self%ddl)
         ga(k)=ga(k)/(tanh(self%ddl)+tanh(self%ddu)) - 1._real64
         dga(k)=ga(k)-ga(k-1)
      end do
      deallocate(ga)
   end if
   return
END SUBROUTINE init_sigma

!---------------------------------------------------------------------------

module SUBROUTINE do_sigma(self)
   !! A wrapper for vertical coordinate calculations

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_domain), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: i,j,k
   integer :: stat
!-----------------------------------------------------------------------------
   call self%logs%info('do_sigma()',level=3)
   !! why not ho=hn as sseo=ssen
   do j=self%T%l(2),self%T%u(2)
      do i=self%T%l(1),self%T%u(1)
         if (self%T%mask(i,j) > 0) then
            self%T%ho(i,j,:)=(self%T%sseo(i,j)+self%T%H(i,j))*dga(:)
            self%T%hn(i,j,:)=(self%T%ssen(i,j)+self%T%H(i,j))*dga(:)
!            self%T%ho(i,j,:)=(self%T%H(i,j))*dga(:)
!            self%T%hn(i,j,:)=(self%T%H(i,j))*dga(:)
         end if
      end do
   end do

   !! why not ho=hn as sseo=ssen
   !! if ssen and H are updated in halo zones - extend to all domain
   !! what about mask
   do j=self%U%l(2),self%U%u(2)
      do i=self%U%l(1),self%U%u(1)-1
         if (self%U%mask(i,j) > 0) then
            self%U%ho(i,j,:)=(self%U%sseo(i,j)+self%U%H(i,j))*dga(:)
            self%U%hn(i,j,:)=(self%U%ssen(i,j)+self%U%H(i,j))*dga(:)
            self%U%ho(i,j,:)=(self%U%H(i,j))*dga(:) ! KB
            self%U%hn(i,j,:)=(self%U%H(i,j))*dga(:) ! KB
         end if
      end do
   end do

   !! if ssen and H are updated in halo zones - extend to all domain
   do j=self%V%l(2),self%V%u(2)-1
      do i=self%V%l(1),self%V%u(1)
         if (self%V%mask(i,j) > 0) then
            self%V%ho(i,j,:)=(self%U%sseo(i,j)+self%V%H(i,j))*dga(:)
            self%V%hn(i,j,:)=(self%U%ssen(i,j)+self%V%H(i,j))*dga(:)
            self%V%ho(i,j,:)=(self%V%H(i,j))*dga(:) ! KB
            self%V%hn(i,j,:)=(self%V%H(i,j))*dga(:) ! KB
         end if
      end do
   end do
   return
END SUBROUTINE do_sigma

!---------------------------------------------------------------------------

END SUBMODULE vertical_sigma_smod
