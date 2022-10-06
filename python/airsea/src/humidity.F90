module mod_humidity

   use mod_airsea_variables, only: kelvin, mw_per_ma, Ra, Rw, rk

   implicit none

   public

contains

   ! Saturation vapor pressure (Pa) at given temperature (Celsius)
   ! Lowe (1977), https://doi.org/10.1175/1520-0450(1977)016<0100:AAPFTC>2.0.CO;2
   elemental function saturation_vapor_pressure(t) result(qa)
      real(rk), intent(in) ::  t
      real(rk) :: qa

      ! Note shift of indices for coefficients compared to Lowe (1977)
      real(rk), parameter :: a1=6.107799961_rk
      real(rk), parameter :: a2=4.436518521e-1_rk
      real(rk), parameter :: a3=1.428945805e-2_rk
      real(rk), parameter :: a4=2.650648471e-4_rk
      real(rk), parameter :: a5=3.031240396e-6_rk
      real(rk), parameter :: a6=2.034080948e-8_rk
      real(rk), parameter :: a7=6.136820929e-11_rk

      qa = a1 + t * (a2 + t * (a3 + t * (a4 + t * (a5 + t * (a6 + t * a7)))))
      qa = qa * 100.0_rk ! Conversion millibar --> Pascal
   end function

   ! Specific humidity (kg water/kg moist air) at given vapor pressure (Pa) and air pressure (Pa)
   elemental function specific_humidity(ea, airp) result(qa)
      real(rk), intent(in) ::  ea, airp
      real(rk) :: qa

      ! Mixing ratio r (kg water/kg dry air) = mw_per_ma * ea / (airp - ea)
      ! This relates to specific humidity (kg water/kg moist air) qa = r / (1 + r)
      ! After rearranging we obtain:
      qa = mw_per_ma * ea / (airp - (1.0_rk - mw_per_ma) * ea)
   end function

   ! vapor pressure (Pa) at given specific humidity (kg water/kg moist air) and air pressure (Pa)
   elemental function vapor_pressure(qa, airp) result(ea)
      real(rk), intent(in) ::  qa, airp
      real(rk) :: ea

      ! Expression below obtained by isolating ea from expression for specific_humidity in previous routine
      ea = (qa * airp) / (mw_per_ma + (1.0_rk - mw_per_ma) * qa)
   end function

   ! Air density at given air pressure (Pa), temperature (Celsius) and specific humidity (kg water/kg moist air)
   elemental function air_density(airp, ta, qa) result(rhoa)
      real(rk), intent(in) :: airp, ta, qa
      real(rk) :: rhoa

      rhoa = airp / (((1.0_rk - qa) * Ra + qa * Rw) * (ta + kelvin))
   end function

end module
!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
