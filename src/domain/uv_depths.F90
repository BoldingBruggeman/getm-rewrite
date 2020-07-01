! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

!> In this routine which is called once during the model initialisation,
!> the bathymetry value in the U- and the V-points are calculated from the
!> bathymetry values in the T-points. The interpolation depends on the value
!> which is given to {\tt vel\_depth\_method}:
!>
!> \begin{equation}
!> H^u_{i,j} = \left\{
!> \begin{array}{ll}
!> \displaystyle
!> \frac12 \left(H_{i,j}+H_{i+1,j}\right), &
!> \displaystyle
!> \mbox{ for vel_depth_method } =0, \\ \\
!> \displaystyle
!> \min\left\{H_{i,j}+H_{i+1,j}\right\}, &
!> \displaystyle
!> \mbox{ for vel_depth_method } =1, \\ \\
!> \displaystyle
!> \min\left\{H_{i,j}+H_{i+1,j}\right\}, &
!> \displaystyle
!> \mbox{ for vel_depth_method } =2 \mbox{ and } \min\{H_{i,j}i,H_{i+1,j}\}<D_{crit} \\ \\
!> \displaystyle
!> \frac12 \left(H_{i,j}+H_{i+1,j}\right), &
!> \displaystyle
!> \mbox{ for vel_depth_method } =2 \mbox{ and } \min\{H_{i,j},H_{i+1,j}\}\geq D_{crit} \\ \\
!> \end{array}
!> \right.
!> \end{equation}
!>
!> The calculation of  \(H^v_{i,j}\) is done accordingly.
!>
!> The options 1 and 2 for **vel\_depth\_method** may help to stabilise
!> calculations when drying and flooding is involved.
!>


SUBMODULE (getm_domain) uv_depths_smod

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

module SUBROUTINE uv_depths(self)

  IMPLICIT NONE

! Subroutine arguments
   class(type_getm_domain), intent(inout) :: self

! Local constants

! Local variables
   integer :: i,j
!-----------------------------------------------------------------------------
   call self%logs%info('uv_depths()',level=2)
   ! U-mask
   do j=self%U%jmin,self%U%jmax
      do i=self%U%imin,self%U%imax-1
         if (self%T%mask(i,j) == 1 .and. self%T%mask(i+1,j) == 1) then
            self%U%mask(i,j)=1
         end if
         if ((self%T%mask(i,j) == 1 .and. self%T%mask(i+1,j) == 2) .or. &
             (self%T%mask(i,j) == 2 .and. self%T%mask(i+1,j) == 1)) then
            self%U%mask(i,j)=2
         end if
         if (self%T%mask(i,j) == 2 .and. self%T%mask(i+1,j) == 2) then
            self%U%mask(i,j)=3
         end if
      end do
   end do
   ! U-depths
   do j=self%U%jmin,self%U%jmax
      do i=self%U%imin,self%U%imax-1
         if (self%U%mask(i,j) > 0) then
            select case (vel_depth_method)
               case (0)
                  self%U%H(i,j) = 0.5_real64 * ( self%T%H(i,j) + self%T%H(i+1,j) )
               case (1)
                  self%U%H(i,j)=min(self%T%H(i,j),self%T%H(i+1,j))
               case (2)
                  if (self%T%H(i,j) .lt. Dcrit .or. self%T%H(i+1,j) .lt. Dcrit) then
                     self%U%H(i,j)=min(self%T%H(i,j),self%T%H(i+1,j))
                  else
                     self%U%H(i,j) = 0.5_real64 * ( self%T%H(i,j) + self%T%H(i+1,j) )
                  end if
            end select
         end if
      end do
   end do

   ! V-mask
   do j=self%V%jmin,self%V%jmax
      do i=self%V%imin,self%V%imax-1
         if (self%T%mask(i,j) == 1 .and. self%T%mask(i,j+1) == 1) then
            self%V%mask(i,j)=1
         end if
         if ((self%T%mask(i,j) == 1 .and. self%T%mask(i,j+1) == 2) .or. &
             (self%T%mask(i,j) == 2 .and. self%T%mask(i,j+1) == 1)) then
            self%V%mask(i,j)=2
         end if
         if (self%T%mask(i,j) == 2 .and. self%T%mask(i,j+1) == 2) then
            self%V%mask(i,j)=3
         end if
      end do
   end do
   ! V-depths
   do j=self%V%jmin,self%V%u(2)-1
      do i=self%V%imin,self%V%imax
         if (self%V%mask(i,j) > 0) then
            select case (vel_depth_method)
               case (0)
                  self%V%H(i,j) = 0.5_real64 * ( self%T%H(i,j) + self%T%H(i,j+1) )
               case (1)
                  self%V%H(i,j)=min(self%T%H(i,j),self%T%H(i,j+1))
               case (2)
                  if (self%T%H(i,j) .lt. Dcrit .or. self%T%H(i,j+1) .lt. Dcrit) then
                     self%V%H(i,j)=min(self%T%H(i,j),self%T%H(i,j+1))
                  else
                     self%V%H(i,j) = 0.5_real64 * ( self%T%H(i,j) + self%T%H(i,j+1) )
                  end if
            end select
         end if
      end do
   end do
   ! Initialize time varying depths in U and V points
   self%T%D = self%T%H
   self%U%D = self%U%H
   self%V%D = self%V%H
   call self%logs%info('done',level=2)
   return
END SUBROUTINE uv_depths

!---------------------------------------------------------------------------

END SUBMODULE uv_depths_smod
