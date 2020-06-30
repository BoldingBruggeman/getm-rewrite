!>
!> Interpolate (latx,lonx) and (xx,yx) to the u, v, and T-points.
!>

SUBROUTINE x2uvc()

   !! This routine interpolates (latx,lonx), (xx,yx), and convx to the
   !! u-points, v-points, and the central T-points. The data at the T-points
   !! are only updated from values of the X-points if the logical flags {\tt updateXYC},
   !! {\tt updateXYC}, and {\tt updateXYC} are {\tt .true.}. This is not necessary
   !! if data at the T-points have been read from the topo input file.
   !!
   !! Original author(s): Lars Umlauf

   IMPLICIT NONE

!  Subroutine arguments

!  Local constants

!  Local variables
   integer :: i,j,n
   real(rk) :: x

!------------------------------------------------------------------------

#ifdef _CARRTESIAN_
   if ( have_lonlat ) then

      do j=jll,jhl
         do i=ill,ihl-1
            latu(i,j) = _HALF_ * ( latc(i,j) + latc(i+1,j) )
         end do
      end do

      do j=jll,jhl-1
         do i=ill,ihl
            latv(i,j) = _HALF_ * ( latc(i,j) + latc(i,j+1) )
         end do
      end do

! this is just a check and can be deleted if nobody experiences problems
! Note (KK): there should only be problems for periodic domains (mask=1)
!            (the mask=2 condition is always wrong)
#if 1
      if ( joff+jhl .eq. jextr ) then ! most northern subdomain
         do i=ill,ihl
            if ( av(i,jhl) .eq. 1 .or. av(i,jhl) .eq. 2 ) then
               latv(i,jhl) = 1000
               LEVEL0 'x2uvc() - warning: latv is set to illegal value'
               LEVEL0 'please report the problem on getm-users'
               stop
            end if
         end do
      end if

      if ( ioff+ihl .eq. iextr ) then ! most eastern subdomain
         do j=jll,jhl
            if ( au(ihl,j) .eq. 1 .or. au(ihl,j) .eq. 2 ) then
               latu(ihl,j) = 1000
               LEVEL0 'x2uvc() - warning: latu is set to illegal value'
               LEVEL0 'please report the problem on getm-users'
               stop
            end if
         end do
      end if
#endif

   end if
#endif

#ifdef _SPHERICAL_

!        we need latx to calculate dxv - utilize equidistance
   latx(ill:ihl,jll-1) = latc(ill:ihl,jll) - dlat/2.
   n=1
   do j=jll,jhl
      latx(ill:ihl,j) = latx(ill:ihl,jll-1) + n*dlat
      n=n+1
   end do

   latu = latc
   latv(ill:ihl,jll:jhl) = latx(ill:ihl,jll:jhl)
#endif

#ifdef _CURVILINEAR_
   do j=jll,jhl
      do i=ill,ihl
         xu(i,j)   = ( xx(i,j) +   xx(i,j-1) ) / 2
         yu(i,j)   = ( yx(i,j) +   yx(i,j-1) ) / 2

         xv(i,j)   = ( xx(i,j) +   xx(i-1,j) ) / 2
         yv(i,j)   = ( yx(i,j) +   yx(i-1,j) ) / 2
      end do
   end do

   do j=jll,jhl
      do i=ill+1,ihl
         xc(i,j)   = ( xu(i,j) +  xu(i-1,j) ) / 2
      end do
   end do

   do j=jll+1,jhl
      do i=ill,ihl
         yc(i,j)   = ( yv(i,j) +   yv(i,j-1) ) / 2
      end do
   end do


   if ( have_lonlat ) then

      do j=jll,jhl
         do i=ill,ihl

            latu(i,j)  = ( latx(i,j) + latx(i,j-1) ) / 2

            latv(i,j)  = ( latx(i,j) + latx(i-1,j) ) / 2

            lonc(i,j)  = ( lonx(i-1,j-1) + lonx(i-1,j) &
                         + lonx(i  ,j-1) + lonx(i,j  ) ) / 4

            latc(i,j)  = ( latx(i-1,j-1) + latx(i-1,j) &
                         + latx(i  ,j-1) + latx(i,j  ) ) / 4

            convc(i,j) = ( convx(i-1,j-1) + convx(i-1,j) &
                         + convx(i  ,j-1) + convx(i,j  ) ) / 4
         end do
      end do
   end if
#endif

   return
END SUBROUTINE x2uvc
