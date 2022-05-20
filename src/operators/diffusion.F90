! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!!{!./code/operators/diffusion.md!}

! https://www.tek-tips.com/viewthread.cfm?qid=1726162
! https://pdfs.semanticscholar.org/0798/fa452cda22b0b501cf1388a021931efe1686.pdf

!> @note
!> ckeck dimension order of auxo and auxn
!> ckeck dimension order of a1, a2, a3, a4
!> @endnote

SUBMODULE (getm_operators) diffusion_smod

!-----------------------------------------------------------------------------

   logical :: is_initialized=.false.

CONTAINS

#define _NORMAL_ORDER_
!-----------------------------------------------------------------------------

MODULE SUBROUTINE vertical_diffusion_initialize_grid(self,grid)

   !! Initialize the vertical diffusion operator from a grid object

   IMPLICIT NONE

   ! Subroutine arguments
   class(type_vertical_diffusion), intent(inout) :: self
   type(type_getm_grid), intent(in) :: grid
      !! grid dimensions in case of dynamic memory allocation

!  Local constants

!  Local variables
!---------------------------------------------------------------------------
#ifndef _STATIC_
   self%halo=grid%halo
   self%imin = grid%imin; self%imax = grid%imax
   self%jmin = grid%jmin; self%jmax = grid%jmax
   self%kmin = grid%kmin; self%kmax = grid%kmax
   allocate(self%auxo(grid%l(1):grid%u(1),grid%l(2):grid%u(2),grid%l(3):grid%u(3)))
   allocate(self%auxn(grid%l(1):grid%u(1),grid%l(2):grid%u(2),grid%l(3):grid%u(3)))
#ifdef _NORMAL_ORDER_
#define ORDER i,j,k
   allocate(self%a1(grid%l(1):grid%u(1),grid%l(2):grid%u(2),grid%l(3):grid%u(3)))
   allocate(self%a2(grid%l(1):grid%u(1),grid%l(2):grid%u(2),grid%l(3):grid%u(3)))
   allocate(self%a3(grid%l(1):grid%u(1),grid%l(2):grid%u(2),grid%l(3):grid%u(3)))
   allocate(self%a4(grid%l(1):grid%u(1),grid%l(2):grid%u(2),grid%l(3):grid%u(3)))
#else
#define ORDER k,i,j
   allocate(self%a1(grid%l(3):grid%u(3),grid%l(1):grid%u(1),grid%l(2):grid%u(2)))
   allocate(self%a2(grid%l(3):grid%u(3),grid%l(1):grid%u(1),grid%l(2):grid%u(2)))
   allocate(self%a3(grid%l(3):grid%u(3),grid%l(1):grid%u(1),grid%l(2):grid%u(2)))
   allocate(self%a4(grid%l(3):grid%u(3),grid%l(1):grid%u(1),grid%l(2):grid%u(2)))
#endif
#endif
   self%matrix_time = 0._real64
   self%tridiag_time = 0._real64
   is_initialized=.true.
END SUBROUTINE vertical_diffusion_initialize_grid

!---------------------------------------------------------------------------

MODULE SUBROUTINE vertical_diffusion_initialize_field(self,f,halo)

   !! Initialize the vertical diffusion operator from a 3D field

   IMPLICIT NONE

   ! Subroutine arguments
   class(type_vertical_diffusion), intent(inout) :: self
   real(real64), dimension(:,:,:), intent(in) :: f
      !! grid dimensions in case of dynamic memory allocation
   integer, dimension(3), intent(in), optional :: halo

!  Local constants

!  Local variables
!---------------------------------------------------------------------------
   if (present(halo)) self%halo=halo
#ifndef _STATIC_
   self%imin = lbound(f,1); self%imax = ubound(f,1)
   self%jmin = lbound(f,2); self%jmax = ubound(f,2)
   self%kmin = lbound(f,3); self%kmax = ubound(f,3)
   allocate(self%auxo(self%imin:self%imax,self%jmin:self%jmax,self%kmin:self%kmax))
   allocate(self%auxn(self%imin:self%imax,self%jmin:self%jmax,self%kmin:self%kmax))
#ifdef _NORMAL_ORDER_
#define ORDER i,j,k
   allocate(self%a1(self%imin:self%imax,self%jmin:self%jmax,self%kmin:self%kmax))
   allocate(self%a2(self%imin:self%imax,self%jmin:self%jmax,self%kmin:self%kmax))
   allocate(self%a3(self%imin:self%imax,self%jmin:self%jmax,self%kmin:self%kmax))
   allocate(self%a4(self%imin:self%imax,self%jmin:self%jmax,self%kmin:self%kmax))
#else
#define ORDER k,i,j
   allocate(self%a1(self%kmin:self%kmax,self%imin:self%imax,self%kmin:self%kmax))
   allocate(self%a2(self%kmin:self%kmax,self%imin:self%imax,self%kmin:self%kmax))
   allocate(self%a3(self%kmin:self%kmax,self%imin:self%imax,self%kmin:self%kmax))
   allocate(self%a4(self%kmin:self%kmax,self%imin:self%imax,self%kmin:self%kmax))
#endif
#endif
   self%matrix_time = 0._real64
   self%tridiag_time = 0._real64
   is_initialized=.true.
END SUBROUTINE vertical_diffusion_initialize_field

!---------------------------------------------------------------------------

MODULE SUBROUTINE vertical_diffusion_calculate(self,dt,cnpar,mask,dzo,dzn,molecular,nuh,var,ea2,ea4)

   !! Vertical diffusion

   IMPLICIT NONE

   ! Subroutine arguments
   class(type_vertical_diffusion), intent(inout) :: self
   real(real64), intent(in) :: dt
   real(real64), intent(in) :: cnpar
#define _T2_ self%imin-self%halo(1):,self%jmin-self%halo(2):
   integer, intent(in) :: mask(_T2_)
#undef _T2_
#define _T3_ self%imin-self%halo(1):,self%jmin-self%halo(2):,self%kmin:
   real(real64), intent(in) :: dzo(_T3_)
   real(real64), intent(in) :: dzn(_T3_)
   real(real64), intent(in) :: molecular
   real(real64), intent(in) :: nuh(_T3_)   ! Note: this is diffusivity defined at the interior interfaces, from kmin to kmax-1. surface/bottom are absent
   real(real64), intent(inout) :: var(_T3_)
   real(real64), intent(in), optional :: ea2(_T3_)
   real(real64), intent(in), optional :: ea4(_T3_)
#undef _T3_

!  Local constants

!  Local variables
   integer :: i,j,k
!---------------------------------------------------------------------------
   if (.not. is_initialized) stop 'vertical_diffusion is not initialized'
   if (self%kmax == 1) return

   ! Setting up the matrix
   matrix: block
   real(real64) :: x
   real(real64) :: matrix_start, matrix_end
   call cpu_time(matrix_start)
   ! Auxilury terms, old and new time level
   do k=self%kmin,self%kmax-1
      do j=self%jmin,self%jmax
         do i=self%imin,self%imax
            if (mask(i,j) ==  1) then
               x = 2._real64*dt*(nuh(i,j,k)+molecular)
               self%auxo(i,j,k)=(1-cnpar)*x/(dzo(i,j,k+1)+dzo(i,j,k))
               self%auxn(i,j,k)=   cnpar *x/(dzn(i,j,k+1)+dzn(i,j,k))
            end if
         end do
      end do
   end do

   ! Matrix elements for surface layer
   k=self%kmax
   do j=self%jmin,self%jmax
      do i=self%imin,self%imax
         if (mask(i,j) ==  1) then
            self%a1(ORDER)=-self%auxn(i,j,k-1)
            self%a2(ORDER)=dzn(i,j,k)+self%auxn(i,j,k-1)
            self%a4(ORDER)=var(i,j,k)*(dzo(i,j,k)-self%auxo(i,j,k-1))+var(i,j,k-1)*self%auxo(i,j,k-1)
         end if
      end do
   end do

   ! Matrix elements for inner layers
   ! do k=kmin,kmax-1
   do k=self%kmin+1,self%kmax-1
      do j=self%jmin,self%jmax
         do i=self%imin,self%imax
            if (mask(i,j) ==  1) then
               self%a3(ORDER)=-self%auxn(i,j,k  )
               self%a1(ORDER)=-self%auxn(i,j,k-1)
               self%a2(ORDER)=dzn(i,j,k)+self%auxn(i,j,k)+self%auxn(i,j,k-1)
               self%a4(ORDER)=var(i,j,k+1)*self%auxo(i,j,k) &
                             +var(i,j,k  )*(dzo(i,j,k)-self%auxo(i,j,k)-self%auxo(i,j,k-1)) &
                             +var(i,j,k-1)*self%auxo(i,j,k-1)
            end if
         end do
      end do
   end do

   ! Matrix elements for bottom layer
   k=self%kmin
   do j=self%jmin,self%jmax
      do i=self%imin,self%imax
         if (mask(i,j) ==  1) then
            self%a3(ORDER)=-self%auxn(i,j,k)
            self%a2(ORDER)=dzn(i,j,k)+self%auxn(i,j,k)
            self%a4(ORDER)=var(i,j,k+1)*self%auxo(i,j,k)+var(i,j,k)*(dzo(i,j,k)-self%auxo(i,j,k))
         end if
      end do
   end do
   if (present(ea2)) self%a2=self%a2-ea2
   if (present(ea4)) self%a4=self%a4+ea4
   call cpu_time(matrix_end)
   self%matrix_time = self%matrix_time + matrix_end - matrix_start
   end block matrix

   tridiagonal: block
   ! Solving the tridiagonal system
   real(real64) :: ru(self%kmax), qu(self%kmax)
   real(real64) :: tridiag_start, tridiag_end
   call cpu_time(tridiag_start)
   do j=self%jmin,self%jmax
      do i=self%imin,self%imax
         if (mask(i,j) ==  1) then
#ifdef _NORMAL_ORDER_
            ru(self%kmax)=self%a1(i,j,self%kmax)/self%a2(i,j,self%kmax)
            qu(self%kmax)=self%a4(i,j,self%kmax)/self%a2(i,j,self%kmax)
            do k=self%kmax-1,self%kmin+1,-1
               ru(k)=self%a1(i,j,k)/(self%a2(i,j,k)-self%a3(i,j,k)*ru(k+1))
               qu(k)=(self%a4(i,j,k)-self%a3(i,j,k)*qu(k+1))/(self%a2(i,j,k)-self%a3(i,j,k)*ru(k+1))
            end do
            qu(self%kmin)=(self%a4(i,j,self%kmin)-self%a3(i,j,self%kmin)*qu(self%kmin+1)) &
                         /(self%a2(i,j,self%kmin)-self%a3(i,j,self%kmin)*ru(self%kmin+1))
#else
            ru(self%kmax)=self%a1(self%kmax,i,j)/self%a2(self%kmax,i,j)
            qu(self%kmax)=self%a4(self%kmax,i,j)/self%a2(self%kmax,i,j)

            do k=self%kmax-1,self%kmin+1,-1
               ru(k)=self%a1(k,i,j)/(self%a2(k,i,j)-self%a3(k,i,j)*ru(k+1))
               qu(k)=(self%a4(k,i,j)-self%a3(k,i,j)*qu(k+1))/(self%a2(k,i,j)-self%a3(k,i,j)*ru(k+1))
            end do
            qu(self%kmin)=(self%a4(self%kmin,i,j)-self%a3(self%kmin,i,j)*qu(self%kmin+1)) &
                         /(self%a2(self%kmin,i,j)-self%a3(self%kmin,i,j)*ru(self%kmin+1))
#endif
            var(i,j,self%kmin)=qu(self%kmin)
            do k=self%kmin+1,self%kmax
               var(i,j,k)=qu(k)-ru(k)*var(i,j,k-1)
            end do
         end if
      end do
   end do
   call cpu_time(tridiag_end)
   self%tridiag_time = self%tridiag_time + tridiag_end - tridiag_start
   end block tridiagonal
END SUBROUTINE vertical_diffusion_calculate

#undef _NORMAL_ORDER_

!---------------------------------------------------------------------------

END SUBMODULE diffusion_smod
