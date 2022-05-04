module mod_longwave_radiation

   use mod_airsea_variables, only: emiss,bolz,rk

   implicit none

   private

   public clark, hastenrath, bignami, berliand, josey1, josey2

contains

   real(rk) elemental function cloud_correction_factor(dlat)
      real(rk), intent(in)  :: dlat

      real(rk), parameter, dimension(91)  :: ccf = (/ &
        0.497202,     0.501885,     0.506568,     0.511250,     0.515933, &
        0.520616,     0.525299,     0.529982,     0.534665,     0.539348, &
        0.544031,     0.548714,     0.553397,     0.558080,     0.562763, &
        0.567446,     0.572129,     0.576812,     0.581495,     0.586178, &
        0.590861,     0.595544,     0.600227,     0.604910,     0.609593, &
        0.614276,     0.618959,     0.623641,     0.628324,     0.633007, &
        0.637690,     0.642373,     0.647056,     0.651739,     0.656422, &
        0.661105,     0.665788,     0.670471,     0.675154,     0.679837, &
        0.684520,     0.689203,     0.693886,     0.698569,     0.703252, &
        0.707935,     0.712618,     0.717301,     0.721984,     0.726667, &
        0.731350,     0.736032,     0.740715,     0.745398,     0.750081, &
        0.754764,     0.759447,     0.764130,     0.768813,     0.773496, &
        0.778179,     0.782862,     0.787545,     0.792228,     0.796911, &
        0.801594,     0.806277,     0.810960,     0.815643,     0.820326, &
        0.825009,     0.829692,     0.834375,     0.839058,     0.843741, &
        0.848423,     0.853106,     0.857789,     0.862472,     0.867155, &
        0.871838,     0.876521,     0.881204,     0.885887,     0.890570, &
        0.895253,     0.899936,     0.904619,     0.909302,     0.913985, &
        0.918668 /)

      ! Calculate cloud correction factor,fortran counts from 1 !
      ! Of course abs(latitude) should always lie between 0 and 90.
      ! But if the caller includes points with a missing/masked value,
      ! latitude may lie outside the natural range. That should not trigger
      ! out-of-bound array access.
      if (dlat >= -90._rk .and. dlat <= 90._rk) then
         cloud_correction_factor = ccf(nint(abs(dlat))+1)
      else
         cloud_correction_factor = 1._rk
      end if
   end function

   ! Clark et al. (1974) formula.
   elemental subroutine clark(dlat,tw,ta,cloud,ea,ql)
      real(rk), intent(in)  :: dlat,tw,ta,cloud,ea
      real(rk), intent(out) :: ql
      real(rk)              :: x1,x2,x3,ccf

      ccf = cloud_correction_factor(dlat)

      ! unit of ea is Pascal, must hPa
      ! Black body defect term, clouds, water vapor correction
      x1=(1.0-ccf*cloud*cloud)*(tw**4)
      x2=(0.39-0.05*sqrt(ea*0.01))
      ! temperature jump term
      x3=4.0*(tw**3)*(tw-ta)
      ql=-emiss*bolz*(x1*x2+x3)
   end subroutine

   ! Hastenrath and Lamb (1978) formula.
   elemental subroutine hastenrath(dlat,tw,ta,cloud,qa,ql)
      real(rk), intent(in)  :: dlat,tw,ta,cloud,qa
      real(rk), intent(out) :: ql
      real(rk)              :: x1,x2,x3,ccf

      ccf = cloud_correction_factor(dlat)

      ! qa in g(water)/kg(wet air)
      x1=(1.0-ccf*cloud*cloud)*(tw**4)
      x2=(0.39-0.056*sqrt(1000.0*qa))
      x3=4.0*(tw**3)*(tw-ta)
      ql=-emiss*bolz*(x1*x2+x3)
   end subroutine

   ! Bignami et al. (1995) formula (Med Sea).
   elemental subroutine bignami(dlat,tw,ta,cloud,ea,ql)
      real(rk), intent(in)  :: dlat,tw,ta,cloud,ea
      real(rk), intent(out) :: ql
      real(rk)              :: x1,x2,x3,ccf

      ccf = cloud_correction_factor(dlat)

      ! unit of ea is Pascal, must hPa
      ccf=0.1762
      x1=(1.0+ccf*cloud*cloud)*ta**4
      x2=(0.653+0.00535*(ea*0.01))
      x3= emiss*(tw**4)
      ql=-bolz*(-x1*x2+x3)
   end subroutine

   ! Berliand & Berliand (1952) formula (ROMS).
   elemental subroutine berliand(tw,ta,cloud,ea,ql)
      real(rk), intent(in)  :: tw,ta,cloud,ea
      real(rk), intent(out) :: ql
      real(rk)              :: x1,x2,x3

      x1=(1.0-0.6823*cloud*cloud)*ta**4
      x2=(0.39-0.05*sqrt(0.01*ea))
      x3=4.0*ta**3*(tw-ta)
      ql=-emiss*bolz*(x1*x2+x3)
   end subroutine

   ! Josey et.al. 2003 - (J1,9)
   elemental subroutine josey1(tw,ta,cloud,ql)
      real(rk), intent(in)  :: tw,ta,cloud
      real(rk), intent(out) :: ql
      real(rk)              :: x1,x2,x3

      x1=emiss*tw**4
      x2=(10.77*cloud+2.34)*cloud-18.44
      x3=0.955*(ta+x2)**4
      ql=-bolz*(x1-x3)
   end subroutine

   ! Josey et.al. 2003 - (J2,14)
   elemental subroutine josey2(tw,ta,cloud,ea,ql)
      real(rk), intent(in)  :: tw,ta,cloud,ea
      real(rk), intent(out) :: ql
      real(rk)              :: x1,x2,x3

      x1=emiss*tw**4
      ! AS avoid zero trap, limit to about 1% rel. humidity ~ 10Pa
      x2=34.07+4157.0/(log(2.1718e10/max(ea, 10._rk)))
      x2=(10.77*cloud+2.34)*cloud-18.44+0.84*(x2-ta+4.01)
      x3=0.955*(ta+x2)**4
      ql=-bolz*(x1-x3)
   end subroutine

end module
!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
