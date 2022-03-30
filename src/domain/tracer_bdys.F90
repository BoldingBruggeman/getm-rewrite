! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!! Update open boundaries - either via externally provided
!! boundary data - or as 0-gradient.

SUBMODULE (getm_domain) tracer_bdys_smod

!---------------------------------------------------------------------------

CONTAINS

!---------------------------------------------------------------------------

!---------------------------------------------------------------------------

module SUBROUTINE tracer_bdy_3d(self,g,f,bdytype,bdy)

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_domain), intent(inout) :: self
   class(type_getm_grid), intent(in) :: g
   real(real64), intent(inout) :: f(g%l(1):,g%l(2):,g%l(3):)
   integer, intent(in) :: bdytype
   real(real64), intent(in), optional :: bdy(g%l(3):,:)

!  Local constants

!  Local variables
   real(real64), parameter :: sp(4)=(/1._real64, 0.5625_real64, &
                                      0.25_real64, 0.0625_real64/)
   integer :: i,j,k,l,n
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('tracer_bdy_3d()',level=2)
   l=0
   do n=1,self%nwb
      l=l+1
      k=self%bdy_index(l)
      i=self%wi(n)
      do j=self%wfj(n),self%wlj(n)
         select case (bdytype)
            case (zero_gradient)
               f(i,j,:)=f(i+1,j,:)
            case (clamped)
               f(i,j,:)=bdy(k,:)
            case (sponge)
!KB - update
               f(i,j,:)=bdy(k,:)
         end select
         k= k+1
      end do
   end do

   do n=1,self%nnb
      l=l+1
      k=self%bdy_index(l)
      j=self%nj(n)
      do i=self%nfi(n),self%nli(n)
         select case (bdytype)
            case (zero_gradient)
               f(i,j,:)=f(i,j-1,:)
            case (clamped)
               f(i,j,:)=bdy(k,:)
            case (sponge)
!KB - update
               f(i,j,:)=bdy(k,:)
         end select
         k=k+1
      end do
   end do

   do n=1,self%neb
      l=l+1
      k=self%bdy_index(l)
      i=self%ei(n)
      do j=self%efj(n),self%elj(n)
         select case (bdytype)
            case (zero_gradient)
               f(i,j,:)=f(i-1,j,:)
            case (clamped)
               f(i,j,:)=bdy(k,:)
            case (sponge)
!KB - update
               f(i,j,:)=bdy(k,:)
         end select
         k=k+1
      end do
   end do

   do n=1,self%nsb
      l=l+1
      k=self%bdy_index(l)
      j=self%sj(n)
      do i=self%sfi(n),self%sli(n)
         select case (bdytype)
            case (zero_gradient)
               f(i,j,:)=f(i,j+1,:)
            case (clamped)
               f(i,j,:)=bdy(k,:)
            case (sponge)
!KB - update
               f(i,j,:)=bdy(k,:)
         end select
         k=k+1
      end do
   end do
END SUBROUTINE tracer_bdy_3d

!---------------------------------------------------------------------------

END SUBMODULE tracer_bdys_smod
