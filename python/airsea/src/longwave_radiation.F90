module mod_longwave_radiation

   use mod_airsea_variables, only: emiss,bolz,rk

   implicit none

   private

   public clark, hastenrath, bignami, berliand, josey1, josey2

   contains

   ! According to Clark et al. (1974), this table is from:
   ! JOHNSON, J. H.. G. A. FLITTNER, and M. W. CLINE (1965)
   ! Automatic data processing program for marine synoptic radio weather reports.
   ! U.S. Fish Wildl. Serv., Spec. Sci. Rep. Fish. 503, iv + 74 p.
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
         cloud_correction_factor = ccf(nint(abs(dlat)) + 1)
      else
         cloud_correction_factor = 1._rk
      end if
   end function

   ! Clark et al. (1974), https://swfsc-publications.fisheries.noaa.gov/publications/CR/1974/7406.PDF
   elemental subroutine clark(dlat, tw, ta, cloud, ea, ql)
      real(rk), intent(in)  :: dlat, tw, ta, cloud, ea
      real(rk), intent(out) :: ql
      real(rk)              :: x1, x2, x3, ccf, tw3

      ccf = cloud_correction_factor(dlat)

      ! unit of vapor pressure ea is Pascal, must hPa (= mbar)
      ! Black body defect term, clouds, water vapor correction
      tw3 = tw**3
      x1 = (1.0_rk - ccf * cloud * cloud) * tw3 * tw
      x2 = 0.39_rk - 0.05_rk * sqrt(ea * 0.01_rk)
      ! temperature jump term
      x3 = 4.0_rk * tw3 * (tw - ta)
      ql = -emiss * bolz * (x1 * x2 + x3)
   end subroutine

   ! Hastenrath and Lamb (1978)
   ! Heat Budget Atlas of the Tropical Atlantic and Eastern Pacific Oceans
   ! 90 pp., University of Wisconsin Press, Madison
   elemental subroutine hastenrath(dlat, tw, ta, cloud, qa, ql)
      real(rk), intent(in)  :: dlat, tw, ta, cloud, qa
      real(rk), intent(out) :: ql
      real(rk)              :: x1, x2, x3, ccf, tw3

      ccf = cloud_correction_factor(dlat)

      ! unit of specific humidity qa is kg/kg, must be g(water)/kg(wet air)
      tw3 = tw**3
      x1 = (1.0_rk - ccf * cloud * cloud) * tw3 * tw
      x2 = 0.39_rk - 0.056_rk * sqrt(1000.0_rk * qa)
      x3 = 4.0_rk * tw3 * (tw - ta)
      ql = -emiss * bolz * (x1 * x2 + x3)
   end subroutine

   ! Bignami et al. (1995) (Med Sea), https://doi.org/10.1029/94JC02496
   elemental subroutine bignami(dlat, tw, ta, cloud, ea, ql)
      real(rk), intent(in)  :: dlat, tw, ta, cloud, ea
      real(rk), intent(out) :: ql
      real(rk)              :: x1, x2, x3, ccf

      ! unit of vapor pressure ea is Pascal, must hPa (= mbar)
      ccf = 0.1762_rk
      x1 = (1.0_rk + ccf * cloud * cloud) * ta**4
      x2 = 0.653_rk + 0.00535_rk * (ea * 0.01_rk)
      x3 = emiss * tw**4
      ql = -bolz * (-x1 * x2 + x3)
   end subroutine

   ! Berliand & Berliand (1952) (ROMS)
   ! Measurement of the effective radiation of the Earth with varying cloud amounts
   ! Izv. Akad. Nauk. SSSR, Ser. Geofiz., 1, 1952.
   elemental subroutine berliand(tw, ta, cloud, ea, ql)
      real(rk), intent(in)  :: tw, ta, cloud, ea
      real(rk), intent(out) :: ql
      real(rk)              :: x1, x2, x3, ta3

      ! Note: the B&B equation given in Table 3 of Bignami et al. (https://doi.org/10.1029/94JC02496)
      ! has x1 linear in cloud cover instead of quadratic below
      ! unit of vapor pressure ea is Pascal, must hPa (= mbar)
      ta3 = ta**3
      x1 = (1.0_rk - 0.6823_rk * cloud * cloud) * ta3 * ta
      x2 = 0.39_rk - 0.05_rk * sqrt(0.01_rk * ea)
      x3 = 4.0_rk * ta3 * (tw - ta)
      ql = -emiss * bolz * (x1 * x2 + x3)
   end subroutine

   ! Josey et.al. 2003 - (J1,9), https://doi.org/10.1029/2002jc001418
   elemental subroutine josey1(tw, ta, cloud, ql)
      real(rk), intent(in)  :: tw, ta, cloud
      real(rk), intent(out) :: ql
      real(rk)              :: x1, x2, x3

      x1 = emiss * tw**4
      x2 = (10.77_rk * cloud + 2.34_rk) * cloud - 18.44_rk  ! Josey et al Eq 8
      x3 = 0.955_rk * (ta + x2)**4
      ql = -bolz * (x1 - x3)
   end subroutine

   ! Josey et.al. 2003 - (J2,14), https://doi.org/10.1029/2002jc001418
   elemental subroutine josey2(tw, ta, cloud, ea, ql)
      real(rk), intent(in)  :: tw, ta, cloud, ea
      real(rk), intent(out) :: ql
      real(rk)              :: x1,x2,x3

      x1 = emiss * tw**4
      ! AS avoid zero trap, limit to about 1% rel. humidity ~ 10Pa
      x2 = 34.07_rk + 4157.0_rk / (log(2.1718e10_rk / max(ea, 10._rk))) ! Dew point temperature (Josey et al Eq 10)
      x2 = (10.77_rk * cloud + 2.34_rk) * cloud - 18.44_rk + 0.84_rk * (x2 - ta + 4.01_rk)  ! Josey et al Eq 13
      x3 = 0.955_rk * (ta + x2)**4
      ql = -bolz * (x1 - x3)
   end subroutine

end module
!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
