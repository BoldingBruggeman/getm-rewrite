! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

! https://www.tek-tips.com/viewthread.cfm?qid=1726162
! https://pdfs.semanticscholar.org/0798/fa452cda22b0b501cf1388a021931efe1686.pdf

!> @note
!> ckeck dimension order of auxo and auxn
!> ckeck dimension order of a1, a2, a3, a4
!> @endnote

SUBMODULE (getm_operators) diffusion_smod

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

#define _NORMAL_ORDER_

module SUBROUTINE vertical_diffusion_initialize(self,var)

   !! Initialize the salinity field

   IMPLICIT NONE

   ! Subroutine arguments
   class(type_vertical_diffusion), intent(inout) :: self
   real(real64), dimension(:,:,:), intent(in) :: var
      !! grid dimensions in case of dynamic memory allocation

!  Local constants

!  Local variables
!---------------------------------------------------------------------------
#ifndef _STATIC_
   self%imin = lbound(var,1); self%imax = ubound(var,1)
   self%jmin = lbound(var,2); self%jmax = ubound(var,2)
   self%kmin = lbound(var,3); self%kmax = ubound(var,3)
   allocate(self%auxn(self%imin:self%imax,self%jmin:self%jmax,self%kmin:self%kmax))
   allocate(self%auxo(self%imin:self%imax,self%jmin:self%jmax,self%kmin:self%kmax))
#ifdef _NORMAL_ORDER_
#define ORDER i,j,k
   allocate(self%a1(self%imin:self%imax,self%jmin:self%jmax,self%kmin:self%jmax))
   allocate(self%a2(self%imin:self%imax,self%jmin:self%jmax,self%kmin:self%jmax))
   allocate(self%a3(self%imin:self%imax,self%jmin:self%jmax,self%kmin:self%jmax))
   allocate(self%a4(self%imin:self%imax,self%jmin:self%jmax,self%kmin:self%jmax))
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
   return
END SUBROUTINE vertical_diffusion_initialize

!---------------------------------------------------------------------------

module SUBROUTINE vertical_diffusion_calculate(self,mask,dz,dt,cnpar,avmol,nuh,var)

   !! Vertical diffusion

   IMPLICIT NONE

   ! Subroutine arguments
   class(type_vertical_diffusion), intent(inout) :: self
   integer, dimension(:,:), intent(in) :: mask
   real(real64), dimension(:,:,:), intent(in) :: dz
   real(real64), intent(in) :: dt
   real(real64), intent(in) :: cnpar
   real(real64), intent(in) :: avmol
   real(real64), dimension(:,:,:), intent(in) :: nuh
   real(real64), dimension(:,:,:), intent(inout) :: var

!  Local constants

!  Local variables
   integer :: i,j,k
!---------------------------------------------------------------------------

   if (self%kmax == 1) return

matrix: block
   ! Setting up the matrix
   real(real64) :: x
   real(real64) :: matrix_start, matrix_end
   call cpu_time(matrix_start)
   ! Auxilury terms, old and new time level
   do k=self%kmin,self%kmax-1
      do j=self%jmin,self%jmax
         do i=self%imin,self%imax
            if (mask(i,j) ==  1) then
               x = dt*(nuh(i,j,k)+avmol)/(dz(i,j,k+1)+dz(i,j,k))
               self%auxo(i,j,k)=2._real64*(1-cnpar)*x
               self%auxn(i,j,k)=2._real64*   cnpar *x
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
            self%a2(ORDER)=dz(i,j,k)+self%auxn(i,j,k-1)
            self%a4(ORDER)=var(i,j,k  )*(dz(i,j,k)-self%auxo(i,j,k-1)) &
                          +var(i,j,k-1)*self%auxo(i,j,k-1)
         end if
      end do
   end do

   ! Matrix elements for inner layers
   ! do k=kmin,kmax-1
   do k=2,self%kmax-1
      do j=self%jmin,self%jmax
         do i=self%imin,self%imax
            if (mask(i,j) ==  1) then
               self%a3(ORDER)=-self%auxn(i,j,k  )
               self%a1(ORDER)=-self%auxn(i,j,k-1)
               self%a2(ORDER)=dz(i,j,k)+self%auxn(i,j,k)+self%auxn(i,j,k-1)
               self%a4(ORDER)=var(i,j,k+1)*self%auxo(i,j,k) &
                             +var(i,j,k  )*(dz(i,j,k)-self%auxo(i,j,k)-self%auxo(i,j,k-1)) &
                             +var(i,j,k-1)*self%auxo(i,j,k-1)
            end if
         end do
      end do
   end do

   ! Matrix elements for bottom layer
   k=1
   do j=self%jmin,self%jmax
      do i=self%imin,self%imax
         if (mask(i,j) ==  1) then
            self%a3(ORDER)=-self%auxn(i,j,k  )
            self%a2(ORDER)=dz(i,j,k)+self%auxn(i,j,k)
            self%a4(ORDER)=var(i,j,k+1)*self%auxo(i,j,k)+var(i,j,k  )*(dz(i,j,k)-self%auxo(i,j,k))
         end if
      end do
   end do
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
   return
END SUBROUTINE vertical_diffusion_calculate

#undef _NORMAL_ORDER_

!---------------------------------------------------------------------------

END SUBMODULE diffusion_smod
