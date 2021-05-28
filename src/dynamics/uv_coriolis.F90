! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!KB!!{!./pages/momentum_3d.md!}

!> @bug
!> check indices in work2d array calculations
!> @endbug

SUBMODULE (getm_momentum) coriolis_smod

CONTAINS

!---------------------------------------------------------------------------

MODULE SUBROUTINE coriolis_fu(self)

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
      !! GETM momentum type

!  Local constants

!  Local variables
   real(real64) :: cord_curv=0._real64
   integer :: i,j
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('coriolis_fu()',level=3)
   ! Semi-implicit treatment of Coriolis force for V-momentum eq.
   TGrid: associate( TG => self%domain%T )
   UGrid: associate( UG => self%domain%U )
   VGrid: associate( VG => self%domain%V )
   XGrid: associate( XG => self%domain%X )
   select case (self%coriolis_scheme)
      case (0)
         return
      case (1)
         do j=VG%jmin,VG%jmax
            do i=VG%imin,VG%imax
               self%work2d(i,j)=0._real64
               if(VG%mask(i,j) > 0) then
                  self%work2d(i,j)=0.25_real64*(self%U(i-1,j)+self%U(i,j)+self%U(i-1,j+1)+self%U(i,j+1))
               end if
            end do
         end do
      case (2) ! Espelid et al. [2000], IJNME 49, 1521-1545
         do j=VG%jmin,VG%jmax
            do i=VG%imin,VG%imax
               self%work2d(i,j)=0._real64
               if(VG%mask(i,j) > 0) then
                  self%work2d(i,j)=0.25_real64*sqrt(VG%D(i,j)) &
                                  *(self%U(i-1,j  )/sqrt(UG%D(i-1,j  ))+self%U(i,j  )/sqrt(UG%D(i,j  ))  &
                                   +self%U(i-1,j+1)/sqrt(UG%D(i-1,j+1))+self%U(i,j+1)/sqrt(UG%D(i,j+1)))
               end if
            end do
         end do
   end select

   do j=VG%jmin,VG%jmax
      do i=VG%imin,VG%imax
         if(VG%mask(i,j) > 0) then
            if (self%domain%domain_type /= 1) then
               cord_curv=(self%V(i,j)*(XG%dy(i,j)-XG%dy(i-1,j)) &
                         -self%work2d(i,j)*(TG%dx(i,j+1)-TG%dx(i,j)))/VG%D(i,j)*VG%iarea(i,j)
            end if
            self%fU(i,j)=(cord_curv+VG%cor(i,j))*self%work2d(i,j)
         else
            self%fU(i,j)=0._real64
         end if
      end do
   end do
   end associate XGrid
   end associate VGrid
   end associate UGrid
   end associate TGrid
END SUBROUTINE coriolis_fu

!---------------------------------------------------------------------------

MODULE SUBROUTINE coriolis_fv(self)

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
      !! GETM momentum type

!  Local constants

!  Local variables
   integer :: i,j
   real(real64) :: cord_curv
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('coriolis_fv()',level=3)
   ! Semi-implicit treatment of Coriolis force for U-momentum eq.
   TGrid: associate( TG => self%domain%T )
   UGrid: associate( UG => self%domain%U )
   VGrid: associate( VG => self%domain%V )
   XGrid: associate( XG => self%domain%X )
   select case (self%coriolis_scheme)
      case (0)
         return
      case (1)
         do j=UG%jmin,UG%jmax
            do i=UG%imin,UG%imax
               self%work2d(i,j)=0._real64
               if(UG%mask(i,j) .ge. 1) then
                  self%work2d(i,j)=0.25_real64*( self%V(i,j-1)+self%V(i+1,j-1)+self%V(i,j)+self%V(i+1,j))
               end if
            end do
         end do
      case (2) ! Espelid et al. [2000], IJNME 49, 1521-1545
         do j=UG%jmin,UG%jmax
            do i=UG%imin,UG%imax
               self%work2d(i,j)=0._real64
               if(UG%mask(i,j) > 0) then
                  self%work2d(i,j)=0.25_real64*sqrt(UG%D(i,j)) &
                                  *(self%V(i,j-1)/sqrt(VG%D(i,j-1))+ self%V(i+1,j-1)/sqrt(VG%D(i+1,j-1))  &
                                   +self%V(i,j  )/sqrt(VG%D(i,j  ))+ self%V(i+1,j  )/sqrt(VG%D(i+1,j  )))
               end if
            end do
         end do
   end select

   do j=UG%jmin,UG%jmax
      do i=UG%imin,UG%imax
         if(VG%mask(i,j) > 0) then
            if (self%domain%domain_type /= 1) then
               cord_curv=(self%work2d(i,j)*(TG%dy(i+1,j)-TG%dy(i,j))) &
                         +self%U(i,j)*(XG%dx(i,j)-XG%dx(i,j-1))/UG%D(i,j)*UG%iarea(i,j)
            end if
            self%fV(i,j)=(cord_curv+UG%cor(i,j))*self%work2d(i,j)
         else
            self%fV(i,j)=0._real64
         end if
      end do
   end do
   end associate XGrid
   end associate VGrid
   end associate UGrid
   end associate TGrid
END subroutine coriolis_fv


!---------------------------------------------------------------------------

MODULE SUBROUTINE coriolis_fpk(self)

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
      !! GETM momentum type

!  Local constants

!  Local variables
   real(real64) :: cord_curv=0._real64
   integer :: i,j,k
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('coriolis_fpk()',level=3)
   TGrid: associate( TG => self%domain%T )
   UGrid: associate( UG => self%domain%U )
   VGrid: associate( VG => self%domain%V )
   XGrid: associate( XG => self%domain%X )
   do k=VG%kmin,VG%kmax
      select case (self%coriolis_scheme)
         case (0)
            return
         case (1)
            do j=VG%jmin,VG%jmax
               do i=VG%imin,VG%imax
                  self%work2d(i,j)=0._real64
                  if(VG%mask(i,j) > 0) then
                     self%work2d(i,j)=0.25_real64*(self%pk(i,j,k)+self%pk(i+1,j,k)+self%pk(i,j-1,k)+self%pk(i+1,j-1,k))
                  end if
               end do
            end do
         case (2) ! Espelid et al. [2000], IJNME 49, 1521-1545
            do j=UG%jmin,UG%jmax
               do i=UG%imin,UG%imax
                  self%work2d(i,j)=0._real64
                  if (UG%mask(i,j) == 1 .or. UG%mask(i,j) == 2) then
                     self%work2d(i,j)=0.25_real64*sqrt(UG%ho(i,j,k))*( &
                                      +self%pk(i  ,j  ,k)/sqrt(UG%ho(i  ,j  ,k)) &
                                      +self%pk(i+1,j  ,k)/sqrt(UG%ho(i+1,j  ,k)) &
                                      +self%pk(i  ,j-1,k)/sqrt(UG%ho(i  ,j-1,k)) &
                                      +self%pk(i+1,j-1,k)/sqrt(UG%ho(i+1,j-1,k)))
                     end if
                  end do
               end do
      end select

      do j=VG%jmin,VG%jmax
         do i=VG%imin,VG%imax
            if(VG%mask(i,j) > 0) then
               if (self%domain%domain_type /= 1) then
                  cord_curv=(self%qk(i,j,k)*(XG%dy(i,j)-XG%dy(i-1,j)) &
                            -self%work2d(i,j)*(TG%dx(i,j+1)-TG%dx(i,j))) &
                            /VG%ho(i,j,k)*VG%iarea(i,j)
               end if
               self%fpk(i,j,k)=(cord_curv-VG%cor(i,j))*self%work2d(i,j)
            end if
         end do
      end do
   end do
   end associate XGrid
   end associate VGrid
   end associate UGrid
   end associate TGrid
END SUBROUTINE coriolis_fpk

!---------------------------------------------------------------------------

MODULE SUBROUTINE coriolis_fqk(self)

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_momentum), intent(inout) :: self
      !! GETM momentum type

!  Local constants

!  Local variables
   real(real64) :: cord_curv=0._real64
   integer :: i,j,k
!---------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('coriolis_fqk()',level=3)
   TGrid: associate( TG => self%domain%T )
   UGrid: associate( UG => self%domain%U )
   VGrid: associate( VG => self%domain%V )
   XGrid: associate( XG => self%domain%X )
   do k=UG%kmin,UG%kmax
      select case (self%coriolis_scheme)
         case (0)
            return
         case (1)
            do j=UG%jmin,UG%jmax
               do i=UG%imin,UG%imax
                  self%work2d(i,j)=0._real64
                  if(UG%mask(i,j) > 0) then
                     self%work2d(i,j)=0.25_real64*(self%qk(i,j,k)+self%qk(i+1,j,k)+self%qk(i,j-1,k)+self%qk(i+1,j-1,k))
                  end if
               end do
            end do
         case (2) ! Espelid et al. [2000], IJNME 49, 1521-1545
            do j=UG%jmin,UG%jmax
               do i=UG%imin,UG%imax
                  self%work2d(i,j)=0._real64
                  if (UG%mask(i,j) == 1 .or. UG%mask(i,j) == 2) then
                     self%work2d(i,j)=0.25_real64*sqrt(UG%ho(i,j,k))*( &
                                      +self%qk(i  ,j  ,k)/sqrt(VG%ho(i  ,j  ,k)) &
                                      +self%qk(i+1,j  ,k)/sqrt(VG%ho(i+1,j  ,k)) &
                                      +self%qk(i  ,j-1,k)/sqrt(VG%ho(i  ,j-1,k)) &
                                      +self%qk(i+1,j-1,k)/sqrt(VG%ho(i+1,j-1,k)))
                     end if
                  end do
               end do
      end select

      do j=UG%jmin,UG%jmax
         do i=UG%imin,UG%imax
            if(UG%mask(i,j) > 0) then
               if (self%domain%domain_type /= 1) then
                  cord_curv=(self%work2d(i,j)*(TG%dy(i+1,j)-TG%dy(i,j)) &
                            -self%pk(i,j,k)*(XG%dx(i,j)-XG%dx(i,j-1))) &
                            /UG%ho(i,j,k)*UG%iarea(i,j)
               end if
               self%fqk(i,j,k)=(cord_curv+UG%cor(i,j))*self%work2d(i,j)
            end if
         end do
      end do
   end do
   end associate XGrid
   end associate VGrid
   end associate UGrid
   end associate TGrid
END SUBROUTINE coriolis_fqk

!---------------------------------------------------------------------------

END SUBMODULE coriolis_smod

#if 0
domain,u,v,fu - height
   select case (self%coriolis_scheme)
      case (0)
         return
      case (1)
         do j=VG%jmin,VG%jmax
            do i=VG%imin,VG%imax
               self%work2d(i,j)=0._real64
               if(VG%mask(i,j) > 0) then
                  self%work2d(i,j)=0.25_real64*(self%U(i-1,j)+self%U(i,j)+self%U(i-1,j+1)+self%U(i,j+1))
               end if
            end do
         end do
      case (2) ! Espelid et al. [2000], IJNME 49, 1521-1545
         do j=VG%jmin,VG%jmax
            do i=VG%imin,VG%imax
               self%work2d(i,j)=0._real64
               if(VG%mask(i,j) > 0) then
                  self%work2d(i,j)=0.25_real64*sqrt(VG%D(i,j)) &
                                  *(self%U(i,j  )/sqrt(UG%D(i,j  ))+ self%U(i-1,j  )/sqrt(UG%D(i-1,j  ))  &
                                   +self%U(i,j+1)/sqrt(UG%D(i,j+1))+ self%U(i-1,j+1)/sqrt(UG%D(i-1,j+1)))
               end if
            end do
         end do
   end select

   do j=VG%jmin,VG%jmax
      do i=VG%imin,VG%imax
         if(VG%mask(i,j) > 0) then
            if (self%domain%domain_type /= 1) then
               cord_curv=(self%V(i,j)*(XG%dy(i,j)-XG%dy(i-1,j)) &
                         -self%work2d(i,j)*(TG%dx(i,j+1)-TG%dx(i,j)))/VG%D(i,j)*VG%iarea(i,j)
            end if
            self%fU(i,j)=(cord_curv+VG%cor(i,j))*self%work2d(i,j)
         end if
      end do
   end do
#endif
