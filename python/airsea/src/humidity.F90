module mod_humidity

   use mod_airsea_variables, only: kelvin,const06,rgas,rk

   implicit none

   public

   contains

   ! Saturation vapor pressure (Pa) at given temperature (Celsius)
   elemental function saturation_vapor_pressure(t) result(qa)
      real(rk), intent(in) ::  t
      real(rk) :: qa

      ! Note shift of indices for coefficients compared to Lowe (1977, J. Appl. Met.)
      real(rk), parameter       :: a1=6.107799961_rk
      real(rk), parameter       :: a2=4.436518521e-1_rk
      real(rk), parameter       :: a3=1.428945805e-2_rk
      real(rk), parameter       :: a4=2.650648471e-4_rk
      real(rk), parameter       :: a5=3.031240396e-6_rk
      real(rk), parameter       :: a6=2.034080948e-8_rk
      real(rk), parameter       :: a7=6.136820929e-11_rk

      qa = a1 +t*(a2+t*(a3+t*(a4+t*(a5+t*(a6+t*a7)))))
      qa = qa * 100.0_rk ! Conversion millibar --> Pascal
   end function

   ! Specific humidity (kg/kg) at given vapor pressure (Pa) and air pressure
   elemental function specific_humidity(ea, airp) result(qa)
      real(rk), intent(in) ::  ea, airp
      real(rk) :: qa
      qa = const06 * ea / (airp - 0.377_rk * ea)
   end function

   ! Air density at given air pressure, temperature (Celsius) and specific humidity (kg/kg)
   elemental function air_density(airp, ta, qa) result(rhoa)
      real(rk), intent(in) :: airp, ta, qa
      real(rk) :: rhoa
      rhoa = airp / (rgas * (ta + kelvin) * (1.0_rk + const06 * qa))
   end function

end module
!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
