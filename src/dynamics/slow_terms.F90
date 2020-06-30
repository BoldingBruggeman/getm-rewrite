! Copyright (C) 2020 Bolding & Bruggeman
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
   real(real64) :: vertsum,ip_fac
!---------------------------------------------------------------------------
   call self%logs%info('slow_terms()',level=2)

#define UG self%domain%U
   do j=UG%jmin,UG%jmax
      do i=UG%imin,UG%imax
         if (UG%mask(i,j) .ge. 1) then
            self%SlUx(i,j)=-self%UEx(i,j)+SUM(self%uuEx(i,j,1:))-SUM(idpdx(i,j,1:))
         end if
      end do
   end do

   do j=UG%jmin,UG%jmax
      do i=UG%imin,UG%imax
         if (UG%mask(i,j) .ge. 1) then
            !k=kumin(i,j)
            k=1
            if (UG%kmax .gt. 1) then
               self%Slru(i,j)=-self%ru(i,j)*self%Uint(i,j) &
                              /(0.5_real64*(UG%sseo(i,j)+UG%ssen(i,j))+UG%H(i,j)) &
                              +self%rru(i,j)*self%uu(i,j,k) &
                              /(0.5_real64*(UG%ho(i,j,k)+UG%hn(i,j,k)))
#ifdef STRUCTURE_FRICTION
               do k=1,kmax
                  Slru(i,j)=Slru(i,j)+uu(i,j,k)*_HALF_*(sf(i,j,k)+sf(i+1,j,k))
               end do
#endif
            else
               self%Slru(i,j)= 0._real64
            end if
         end if
      end do
   end do
#undef UG

#define VG self%domain%V
   do j=VG%jmin,VG%jmax
      do i=VG%imin,VG%imax
         if (VG%mask(i,j) .ge. 1) then
            self%SlVx(i,j)=-self%VEx(i,j)+SUM(self%vvEx(i,j,1:))-SUM(idpdy(i,j,1:))
         end if
      end do
   end do

   do j=VG%jmin,VG%jmax
      do i=VG%imin,VG%imax
         if (VG%mask(i,j) .ge. 1) then
            !k=kumin(i,j)
            k=1
            if (VG%kmax .gt. 1) then
               self%Slrv(i,j)=-self%rv(i,j)*self%Vint(i,j) &
                              /(0.5_real64*(VG%sseo(i,j)+VG%ssen(i,j))+VG%H(i,j)) &
                              +self%rrv(i,j)*self%vv(i,j,k) &
                              /(0.5_real64*(VG%ho(i,j,k)+VG%hn(i,j,k)))
#ifdef STRUCTURE_FRICTION
               do k=1,kmax
                  Slru(i,j)=Slru(i,j)+uu(i,j,k)*_HALF_*(sf(i,j,k)+sf(i+1,j,k))
               end do
#endif
            else
               self%Slru(i,j)= 0._real64
            end if
         end if
      end do
   end do
#undef VG
   return
END SUBROUTINE slow_terms

!---------------------------------------------------------------------------

END SUBMODULE slow_terms_smod
