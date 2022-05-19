module mod_longwave_radiation

   use mod_airsea_variables, only: emiss,bolz,rk

   implicit none

   private

   public clark, hastenrath, bignami, berliand, josey1, josey2

contains

   real(rk) elemental function cloud_correction_factor(dlat)
      real(rk), intent(in)  :: dlat

      real(rk), parameter, dimension(91)  :: ccf = (/ &
        0.497202_rk,     0.501885_rk,     0.506568_rk,     0.511250_rk,     0.515933_rk, &
        0.520616_rk,     0.525299_rk,     0.529982_rk,     0.534665_rk,     0.539348_rk, &
        0.544031_rk,     0.548714_rk,     0.553397_rk,     0.558080_rk,     0.562763_rk, &
        0.567446_rk,     0.572129_rk,     0.576812_rk,     0.581495_rk,     0.586178_rk, &
        0.590861_rk,     0.595544_rk,     0.600227_rk,     0.604910_rk,     0.609593_rk, &
        0.614276_rk,     0.618959_rk,     0.623641_rk,     0.628324_rk,     0.633007_rk, &
        0.637690_rk,     0.642373_rk,     0.647056_rk,     0.651739_rk,     0.656422_rk, &
        0.661105_rk,     0.665788_rk,     0.670471_rk,     0.675154_rk,     0.679837_rk, &
        0.684520_rk,     0.689203_rk,     0.693886_rk,     0.698569_rk,     0.703252_rk, &
        0.707935_rk,     0.712618_rk,     0.717301_rk,     0.721984_rk,     0.726667_rk, &
        0.731350_rk,     0.736032_rk,     0.740715_rk,     0.745398_rk,     0.750081_rk, &
        0.754764_rk,     0.759447_rk,     0.764130_rk,     0.768813_rk,     0.773496_rk, &
        0.778179_rk,     0.782862_rk,     0.787545_rk,     0.792228_rk,     0.796911_rk, &
        0.801594_rk,     0.806277_rk,     0.810960_rk,     0.815643_rk,     0.820326_rk, &
        0.825009_rk,     0.829692_rk,     0.834375_rk,     0.839058_rk,     0.843741_rk, &
        0.848423_rk,     0.853106_rk,     0.857789_rk,     0.862472_rk,     0.867155_rk, &
        0.871838_rk,     0.876521_rk,     0.881204_rk,     0.885887_rk,     0.890570_rk, &
        0.895253_rk,     0.899936_rk,     0.904619_rk,     0.909302_rk,     0.913985_rk, &
        0.918668_rk /)

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
