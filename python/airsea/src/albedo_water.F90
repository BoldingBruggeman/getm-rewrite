module mod_albedo_water

   use mod_airsea_variables, only: rk

   implicit none

   private

contains

   subroutine albedo_water_2d(nx, ny, imin, imax, jmin, jmax, method, zen, yday, albedo) bind(c)
      integer, intent(in), value                 :: nx, ny, imin, imax, jmin, jmax, method, yday
      real(rk), intent(in),    dimension(nx, ny) :: zen
      real(rk), intent(inout), dimension(nx, ny) :: albedo

      select case (method)
         case (1)
            albedo(imin:imax, jmin:jmax) = albedo_payne(zen(imin:imax, jmin:jmax))
         case (2)
            albedo(imin:imax, jmin:jmax) = albedo_cogley(zen(imin:imax, jmin:jmax), yday)
         case default
            stop 'albedo_2d()'
      end select
   end subroutine

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Albedo over water asccording to Payne(1972)
!
! !INTERFACE:
   real(rk) elemental function albedo_payne(zen)
!
! !DESCRIPTION:
!  The albedo monthly values from \cite{Payne72} are given  here
!  as means of the values between at 30$^{\circ}$ N and 40$^{\circ}$ N
!  for the Atlantic Ocean (hence the same latitudinal band of the
!  Mediterranean Sea).
!  Linear interpolation is done in latitude - time variation is included
!  via the solar zenith angle.
!
! !INPUT PARAMETERS:
   real(rk), intent(in)      :: zen
!
! !DEFINED PARAMETERS:
!  albedo values for 0.7 transmittance
   real(rk), parameter       :: alb1(20) = &
                 (/.719_rk,.656_rk,.603_rk,.480_rk,.385_rk,.300_rk,.250_rk,.193_rk,.164_rk,.131_rk, &
                   .103_rk,.084_rk,.071_rk,.061_rk,.054_rk,.043_rk,.039_rk,.036_rk,.034_rk,.034_rk /)
   real(rk), parameter       :: za(20) = &
                 (/90._rk,88._rk,86._rk,84._rk,82._rk,80._rk,78._rk,76._rk,74._rk,70._rk,  &
                   66._rk,62._rk,58._rk,54._rk,50._rk,40._rk,30._rk,20._rk,10._rk,0.0_rk /)
   real(rk), parameter       :: dza(19) = &
                 (/2.0_rk, 2.0_rk, 2.0_rk, 2.0_rk, 2.0_rk, 2.0_rk, 2.0_rk, 2.0_rk, 4.0_rk, 4.0_rk, &
                   4.0_rk, 4.0_rk, 4.0_rk, 4.0_rk, 10.0_rk, 10.0_rk, 10.0_rk, 10.0_rk, 10.0_rk/)
!
! !LOCAL VARIABLES:
   integer                   :: jab
   real(rk)                  :: dzen
!EOP
!-----------------------------------------------------------------------
!BOC
   if (.not. (zen >= 0. .and. zen <= 90._rk)) then
      ! catches NaN
      jab = 1
   else if(zen .ge. 74.)then
      jab=.5*(90.-zen)+1.
   else if (zen .ge. 50.) then
      jab=.23*(74.-zen)+9.
   else
      jab=.10*(50.-zen)+15.
   endif
   if (jab .eq. 20) then
      dzen = 0._rk
      albedo_payne=alb1(jab)
   else
      dzen=(za(jab)-zen)/dza(jab)
      albedo_payne=alb1(jab)+dzen*(alb1(jab+1)-alb1(jab))
   end if

   end function albedo_payne
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Albedo over water asccording to Cogley()
!
! !INTERFACE:
   real(rk) elemental function albedo_cogley(zen,yd)
!
! !DESCRIPTION:
!  The albedo monthly values are from  Graham Cogley (1979)
!  applied bilinear interpolation in latitude and time.
!  Likely the search of the indices used to look up values in the
!  different tables could be optimized. This is a reference implementation
!  based on code provided by Adolf Stips, JRC.
!
! !INPUT PARAMETERS:
   real(rk), intent(in)      :: zen
   integer, intent(in)       :: yd
!
! !DEFINED PARAMETERS:
   real(rk), parameter       :: a1(12,10) = reshape( &
!  Maybe 0.293 should be changed to 0.239 to give a smooth curve - KB 2015-05-06
!  From Table 5 in Cogley 1979
      (/1._rk   ,1._rk   ,0.301_rk,0.293_rk,0.171_rk,0.148_rk,0.160_rk,0.246_rk,0.342_rk,1._rk   ,1._rk   ,1._rk   ,  &
        1._rk   ,0.301_rk,0.319_rk,0.225_rk,0.160_rk,0.131_rk,0.145_rk,0.206_rk,0.294_rk,0.305_rk,1._rk   ,1._rk   ,  &
        0.301_rk,0.338_rk,0.229_rk,0.148_rk,0.116_rk,0.112_rk,0.114_rk,0.134_rk,0.202_rk,0.313_rk,0.301_rk,1._rk   ,  &
        0.339_rk,0.240_rk,0.155_rk,0.105_rk,0.088_rk,0.084_rk,0.086_rk,0.098_rk,0.136_rk,0.216_rk,0.321_rk,0.355_rk,  &
        0.220_rk,0.161_rk,0.108_rk,0.084_rk,0.075_rk,0.073_rk,0.074_rk,0.080_rk,0.099_rk,0.144_rk,0.210_rk,0.241_rk,  &
        0.145_rk,0.111_rk,0.085_rk,0.073_rk,0.068_rk,0.067_rk,0.068_rk,0.071_rk,0.080_rk,0.103_rk,0.138_rk,0.161_rk,  &
        0.103_rk,0.086_rk,0.073_rk,0.067_rk,0.065_rk,0.064_rk,0.064_rk,0.066_rk,0.071_rk,0.082_rk,0.100_rk,0.111_rk,  &
        0.083_rk,0.074_rk,0.067_rk,0.064_rk,0.063_rk,0.063_rk,0.063_rk,0.064_rk,0.066_rk,0.072_rk,0.081_rk,0.087_rk,  &
        0.072_rk,0.067_rk,0.064_rk,0.063_rk,0.064_rk,0.064_rk,0.064_rk,0.063_rk,0.063_rk,0.066_rk,0.071_rk,0.074_rk,  &
        0.066_rk,0.064_rk,0.063_rk,0.064_rk,0.066_rk,0.068_rk,0.067_rk,0.064_rk,0.063_rk,0.064_rk,0.066_rk,0.068_rk/) &
       , shape(a1))
   real(rk), parameter       :: za(10) = (/90._rk,80._rk,70._rk,60._rk,50._rk,40._rk,30._rk,20._rk,10._rk,0.0_rk/)
   real(rk), parameter       :: dz =  10.0_rk
#if 0
! This is the original provided by AS - as Fortran code
   real(rk), parameter       :: tim(12)= (/ &
                                           15.21_rk,45.62_rk,76.03_rk,106.44_rk,136.85_rk,167.26_rk, &
                                          197.67_rk,228.08_rk,258.49_rk,288.90_rk, 319.31_rk, 349.72_rk &
                                         /)
#else
!  The corresponds to the R test code  - also from AS
   real(rk), parameter       :: tim(12)= (/ &
                                 1._rk,32._rk,60._rk,91._rk,121._rk,152._rk,182._rk,213._rk,244._rk,274._rk,305._rk,335._rk &
                                         /)
#endif
   real(rk), parameter       :: dt= 365.25_rk/12._rk
!
! !LOCAL VARIABLES:
   integer                   :: jab,tab,jab1,tab1
   real(rk)                  :: dzen,dti,r1,r2,x
!EOP
!-----------------------------------------------------------------------
!BOC
   
   if (zen >= 0. .and. zen <= 90._rk) then
      jab=min(max(((90._rk-zen)/dz+1._rk),1._rk),10._rk)
   else
      ! catches NaN
      jab = 1
   end if
   tab=min(max(floor(yd/dt)+1._rk,1._rk),12._rk)

   dzen=(za(jab)-zen)/dz
   dti=(yd-tim(tab))/dt

!  interpolate the two latitudes
   jab1=min(jab+1,10)
   tab1 = tab+1
   if (tab1 .gt. 12) tab1 = 1
   r1=a1(tab,jab) +dzen*(a1(tab,jab1) -a1(tab,jab))
   r2=a1(tab1,jab)+dzen*(a1(tab1,jab1)-a1(tab1,jab))

!  interpolate the time
   x=r1+dti*(r2-r1)
   albedo_cogley = max(min(x, 1.0_rk), 0.0_rk)

   end function albedo_cogley
!EOC
end module
!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
