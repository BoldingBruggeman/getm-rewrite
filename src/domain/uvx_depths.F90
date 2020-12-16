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


SUBMODULE (getm_domain) uvx_depths_smod

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

module SUBROUTINE uvx_depths(self)

  IMPLICIT NONE

! Subroutine arguments
   class(type_getm_domain), intent(inout) :: self

! Local constants

! Local variables
   integer :: i,j
!-----------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('uvx_depths()',level=2)
   TGrid: associate( TG => self%T )
   ! Initialize time varying depths at T-points - U, V and X depths done in loops
   TG%D = TG%H

   ! U-mask
   UGrid: associate( UG => self%U )
   do j=UG%jmin,UG%jmax
      do i=UG%imin,UG%imax
         if (TG%mask(i,j) == 1 .and. TG%mask(i+1,j) == 1) then
            UG%mask(i,j)=1
         end if
         if ((TG%mask(i,j) == 1 .and. TG%mask(i+1,j) == 2) .or. &
             (TG%mask(i,j) == 2 .and. TG%mask(i+1,j) == 1)) then
            UG%mask(i,j)=2
         end if
         if (TG%mask(i,j) == 2 .and. TG%mask(i+1,j) == 2) then
            UG%mask(i,j)=3
         end if
      end do
   end do
   ! U-depths
   do j=UG%jmin,UG%jmax
      do i=UG%imin,UG%imax
         if (UG%mask(i,j) > 0) then
            select case (vel_depth_method)
               case (0)
                  UG%H(i,j) = 0.5_real64 * ( TG%H(i,j) + TG%H(i+1,j) )
               case (1)
                  UG%H(i,j)=min(TG%H(i,j),TG%H(i+1,j))
               case (2)
                  if (TG%H(i,j) .lt. self%Dcrit .or. TG%H(i+1,j) .lt. self%Dcrit) then
                     UG%H(i,j)=min(TG%H(i,j),TG%H(i+1,j))
                  else
                     UG%H(i,j) = 0.5_real64 * ( TG%H(i,j) + TG%H(i+1,j) )
                  end if
            end select
            UG%D(i,j) = UG%H(i,j)
         end if
      end do
   end do
   end associate UGrid

   ! V-mask
   VGrid: associate( VG => self%V )
   do j=VG%jmin,VG%jmax
      do i=VG%imin,VG%imax
         if (TG%mask(i,j) == 1 .and. TG%mask(i,j+1) == 1) then
            VG%mask(i,j)=1
         end if
         if ((TG%mask(i,j) == 1 .and. TG%mask(i,j+1) == 2) .or. &
             (TG%mask(i,j) == 2 .and. TG%mask(i,j+1) == 1)) then
            VG%mask(i,j)=2
         end if
         if (TG%mask(i,j) == 2 .and. TG%mask(i,j+1) == 2) then
            VG%mask(i,j)=3
         end if
      end do
   end do
   ! V-depths
   do j=VG%jmin,VG%jmax
      do i=VG%imin,VG%imax
         if (VG%mask(i,j) > 0) then
            select case (vel_depth_method)
               case (0)
                  VG%H(i,j) = 0.5_real64 * ( TG%H(i,j) + TG%H(i,j+1) )
               case (1)
                  VG%H(i,j)=min(TG%H(i,j),TG%H(i,j+1))
               case (2)
                  if (TG%H(i,j) .lt. self%Dcrit .or. TG%H(i,j+1) .lt. self%Dcrit) then
                     VG%H(i,j)=min(TG%H(i,j),TG%H(i,j+1))
                  else
                     VG%H(i,j) = 0.5_real64 * ( TG%H(i,j) + TG%H(i,j+1) )
                  end if
            end select
            VG%D(i,j) = VG%H(i,j)
         end if
      end do
   end do
   end associate VGrid

   XGrid: associate( XG => self%X )
   ! X-mask
   do j=XG%jmin+1,XG%jmax
      do i=XG%imin+1,XG%imax
         if (TG%mask(i  ,j  ) > 0 .and. TG%mask(i+1,j) > 0 .and. &
             TG%mask(i+1,j+1) > 0 .and. TG%mask(i,j+1) > 0) then
            XG%mask(i,j)=1
         end if
      end do
   end do
   ! X-depths
   do j=XG%jmin,XG%jmax
      do i=XG%imin,XG%imax
         if (XG%mask(i,j) > 0) then
            XG%H(i,j) = 0.25_real64*(TG%H(i,j)+TG%H(i+1,j)+TG%H(i+1,j+1)+TG%H(i,j+1))
            XG%D(i,j) = XG%H(i,j)
         end if
      end do
   end do
   end associate XGrid
   end associate TGrid
END SUBROUTINE uvx_depths

!---------------------------------------------------------------------------

END SUBMODULE uvx_depths_smod
