! Copyright (C) 2020 Bolding & Bruggeman and Hans Burchard

SUBMODULE (getm_domain) metrics_smod

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

module SUBROUTINE metrics(self)
   !! Allocate all domain related variables

   IMPLICIT NONE

!  Subroutine arguments
   class(type_getm_domain), intent(inout) :: self

!  Local constants

!  Local variables
   integer :: i,j,k,n
   real(real64) :: dx,dy
   real(real64) :: dlon,dlat
!-----------------------------------------------------------------------------
   if (associated(self%logs)) call self%logs%info('metrics()',level=2)

!  might need to be changed if we go for a different C-grid layout
   if (self%domain_type == cartesian .or. self%domain_type .eq. spherical) then
      dx = (self%T%c1(self%T%imax) - self%T%c1(self%T%imin))/(self%T%imax-self%T%imin)
      dy = (self%T%c2(self%T%jmax) - self%T%c2(self%T%jmin))/(self%T%jmax-self%T%jmin)
      dlon = dx; dlat = dy
      do n=1,size(self%T%c1)
         self%U%c1(n) = self%T%c1(n)-dx/2._real64
         self%V%c1(n) = self%T%c1(n)
         self%X%c1(n) = self%T%c1(n)-dx/2._real64
      end do
      do n=1,size(self%T%c2)
         self%U%c2(n) = self%T%c2(n)
         self%V%c2(n) = self%T%c2(n)-dx/2._real64
         self%X%c2(n) = self%T%c2(n)-dx/2._real64
      end do
   end if

!  calculate missing fields - will depend on the domain_type
!  and what is read from bathymetry.nc
   select case (self%domain_type)
      case(cartesian)
         self%T%dx = dx; self%T%dy = dy
         self%U%dx = dx; self%U%dy = dy
         self%V%dx = dx; self%V%dy = dy
         self%X%dx = dx; self%X%dy = dy
         if (allocated(self%T%lon) .and. allocated (self%T%lat)) then
            do j=self%T%jmin,self%T%jmax
               do i=self%T%imin,self%T%imax-1
                  self%U%lat(i,j) = 0.5_real64 * ( self%T%lat(i,j) + self%T%lat(i+1,j) )
               end do
            end do
            do j=self%T%jmin,self%T%jmax-1
               do i=self%T%imin,self%T%imax
                  self%V%lat(i,j) = 0.5_real64 * ( self%T%lat(i,j) + self%T%lat(i,j+1) )
               end do
            end do
         end if
      case(spherical)
         do j=self%T%jmin,self%T%jmax
            do i=self%T%imin,self%T%imax
               self%T%lon(i,j) = self%T%c1(i)
               self%T%lat(i,j) = self%T%c2(j)
               self%T%dx(i,j) = deg2rad*dlon*rearth*cos(deg2rad*self%T%lat(i,j))
!               self%V%dx(i,j) = deg2rad*dlon*rearth*cos(deg2rad*self%X%lat(i,j))
            end do
         end do
         self%T%dy = deg2rad*dlat*rearth

         self%U%dx = self%T%dx
         self%U%dy = self%T%dy

         self%V%dy = self%T%dy

         self%X%dx = self%V%dx
         self%X%dy = self%T%dy
#if 0
      case(curvilinear)
#endif
      case default
         stop 'emil'
   end select

!  calculate various metrics- will depend on the domain_type
   select case (self%domain_type)
      case(cartesian)
         self%T%area = self%T%dx*self%T%dy
         self%U%area = self%U%dx*self%U%dy
         self%V%area = self%V%dx*self%V%dy
         self%X%area = self%X%dx*self%X%dy
      case(spherical)
         self%T%area = self%T%dx*self%T%dy
         self%U%area = self%U%dx*self%U%dy
         self%V%area = self%V%dx*self%V%dy
         self%X%area = self%X%dx*self%X%dy
#if 0
      case(curvilinear)
#endif
      case default
         stop 'egon'
   end select

   self%T%inv_area = 1._real64/self%T%area
   self%U%inv_area = 1._real64/self%U%area
   self%V%inv_area = 1._real64/self%V%area
   self%X%inv_area = 1._real64/self%X%area
   if (associated(self%logs)) call self%logs%info('done',level=2)
   return
END SUBROUTINE metrics

!---------------------------------------------------------------------------

END SUBMODULE metrics_smod
