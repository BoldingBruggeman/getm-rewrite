module mod_t10m
   ! Re-written by Per Berg from the COARE 3.0 algorithm
   !     https://www7320.nrlssc.navy.mil/nasec/Module_coare.F
   ! with the purpose of obtaining t,q @10m when the t,q @2m are available.
   !
   ! Please note, "pure" fortran subroutines were developed in order to ease
   ! code optimization (inlining, etc). The PURE prefix-spec is not standard
   ! Fortran 90 but was introduced from standard Fortran 95.


   !- modules: -----------------------------------------------------------------
   ! none used here

   !- implicit directives: -----------------------------------------------------
   use mod_airsea_variables, only: rk

   implicit none

   private

   !- Public routines: ---------------------------------------------------------
   public  :: t10m, t10m_init

   !- Private constants: -------------------------------------------------------
   real(rk), parameter :: qra = 0.61
   real(rk), parameter :: g=9.81, vkarm=0.41
   real(rk), parameter :: qrt = 0.25
   real(rk), parameter :: pi=3.14159265358979323846
   real(rk), parameter :: zero = 0.0, one = 1.0, half = 0.5
   real(rk), parameter :: two = 2.0, three = 3.0
   real(rk), parameter :: four = 4.0, ten = 10.0
   real(rk), parameter :: onethird=1.0/3.0, twothird=2.0/3.0

   integer,  parameter :: nitmin = 1, nitmax = 3

   real(rk), parameter :: Beta = 1.2, von = vkarm, fdg = 1.0, Rgas = 287.1
   real(rk), parameter :: tdk = 273.16, grav = g
   real(rk), parameter :: tsle = -2370.0, Le0 = 2501000.0, wetc0 = 0.622
   real(rk), parameter :: vsa0 = 1.326e-5, vsa1 = 6.542e-3
   real(rk), parameter :: vsa2 = 8.301e-6, vsa3 = -4.84e-9
   real(rk), parameter :: b4 = 0.004
   real(rk), parameter :: charnl = 0.011, charnh = 0.018, utup = 18.0
   real(rk), parameter :: zetuhigh = 50.0
   real(rk), parameter :: neutraldt = -0.0098, dter0 = 0.3, usr0 = 0.035
   real(rk), parameter :: zouu = 0.011, zovu = 0.11, Ch10 = 0.00115
   real(rk), parameter :: xlpow2 = 0.333
   real(rk), parameter :: ug0 = 0.2
   real(rk), parameter :: zoqmax = 0.000115, zoqfac = 0.000055, zoqpow = 0.6
   real(rk), parameter :: cmax = 50.0, cslope = 0.35
   real(rk), parameter :: c0 = 14.28, c1 = 8.525
   real(rk), parameter :: p1  = 10.15, p2  = 1.5, p3 = 34.15
   real(rk), parameter :: fifteen = 15.0, oneph = 1.5

   !- Private variables: -------------------------------------------------------
   real(rk), save :: sq3, a1, a2, Ribcu

   ! zu      wind speed measurement height (m)
   ! zt      air T measurement height (m)
   ! zq      air q measurement height (m)
   ! zi      PBL depth (m), usually 600m
   real(rk), save :: zu, zt, zq, zi

contains

   !----------------------------------------------------------------------------

   pure subroutine t10m (u, ts, t, qs, q, tzs, qzs)
      ! Get temperature @ std ht 10m.
      ! Based on the original cor30a.
      !
      !
      ! input:
      ! u       wind speed (m/s) at height zu (m)
      ! ts      bulk water temperature (C)
      ! t       bulk air temperature (C) at height zt
      ! qs      bulk water spec hum (kg/kg)
      ! q       bulk air spec hum (kg/kg) at height zq
      !
      !
      ! output:
      ! tzs     bulk air temperature (C) at standard height zs=10m
      ! qzs     humidity at standard height


      implicit none

      !- Arguments -------------------------------------------------------------
      real(rk), intent(in)  :: u, ts, t, Qs, Q
      real(rk), intent(out) :: tzs, qzs

      !- Local vars ------------------------------------------------------------
      real(rk)    :: Le, wetc
      real(rk)    :: dt, dq
      real(rk)    :: bf, CC, charn, Ct
      real(rk)    :: L10, Ribu, rr, ta, ut, zet, zetu, zo10, zot10
      real(rk)    :: zo, zot, zoq, L, usr, tsr, qsr
      real(rk)    :: Cd, ug 
      real(rk)    :: xx, visa, psiuo, psiq, psit, psik, psic, z, f
      real(rk)    :: lzu10, lzt10, lzq10, dtjcl
      integer  :: i, nits

      ! air parameters ---------------------------------------------------------
      Le   = Le0 + tsle*ts
      visa = vsa0*(one + (vsa1 + (vsa2 + vsa3*t)*t)*t) 

      ! cool skin parameters ---------------------------------------------------
      wetc = wetc0*Le*qs/(Rgas*(ts+tdk)**2) 
     
      !***************   Begin bulk loop ***************************************
      ! first guess ------------------------------------------------------------
      dt    = ts - t + neutraldt*zt 
      dq    = qs - q 
      ta    = t + tdk 
      ut    = sqrt(u*u + qrt) 
      usr   = usr0*ut 
      zo10  = zouu*usr*usr/grav + zovu*visa/usr 
      xx    = von/log(ten/zo10)
      zot10 = ten*exp(-(von/Ch10)*xx) 
      lzu10 = log(zu/zo10)
      lzt10 = log(zt/zot10)
      lzq10 = log(zq/zot10)
      Cd    = (von/lzu10)**2 
      Ct    = von/lzt10 
      CC    = von*Ct/Cd 
      dtjcl = dter0
      Ribu  = -grav*zu*((dt-dtjcl)/ta + qra*dq)/ut**2 
      if (Ribu < zero) then 
         zetu = CC*Ribu/(one + Ribu/Ribcu) 
      else 
         zetu =    Ribu*(CC  + three*Ribu)
      endif 
      L10 = zu/zetu 
      if (zetu > zetuhigh) then 
         nits = nitmin 
      else
         nits = nitmax
      endif 

      ! wave parameters --------------------------------------------------------
      if (ut > utup) then
         charn = charnh
      elseif (ut > ten) then
         charn = charnl + (ut-ten)/(utup-ten)*(charnh-charnl)
      else
         charn = charnl
      endif

      ! first guesses for usr, sr, qst:
      if (L10 <= zero) then
         ! wind
         z     = zu/L10
         xx    = sqrt(sqrt(one-fifteen*z))
         psik  = log((one+xx)*(one+xx)*(one+xx*xx)) - atan(xx) + a1
         xx    = (one - p1*z)**onethird
         psic  = p2*log(one+xx+xx*xx) - sq3*atan((one+two*xx)/sq3) + a2
         f     = z*z/(one+z*z)
         psiuo = (one-f)*psik + f*psic
         ! temperature
         z     = zt/L10
         xx    = sqrt(one-fifteen*z)
         psik  = two*log(one+xx) - log(four)
         xx    = (one - p3*z)**onethird
         psic  = oneph*log(one+xx+xx*xx) - sq3*atan((one+two*xx)/sq3) + a2
         f     = z*z/(one+z*z)
         psit  = (one-f)*psik + f*psic
         ! humidity
         z     = zq/L10
         xx    = sqrt(one-fifteen*z)
         psik  = two*log(one+xx) - log(four)
         xx    = (one - p3*z)**onethird
         psic  = oneph*log(one+xx+xx*xx) - sq3*atan((one+two*xx)/sq3) + a2
         f     = z*z/(one+z*z)
         psiq  = (one-f)*psik + f*psic
      else
         ! wind
         z     = zu/L10
         psiuo = -((one+z) + twothird*(z-c0)*exp(-min(cmax, cslope*z)) + c1)
         ! temperature
         z     = zt/L10
         xx    = (one+twothird*z)
         psit  = -(xx*sqrt(xx) + twothird*(z-c0)*exp(-min(cmax, cslope*z)) + c1)
         ! humidity
         z     = zq/L10
         xx    = (one+twothird*z)
         psiq  = -(xx*sqrt(xx) + twothird*(z-c0)*exp(-min(cmax, cslope*z)) + c1)
      endif
      usr = ut*von/(lzu10 - psiuo)
      tsr = -(dt-dtjcl)*von*fdg/(lzt10 - psit) 
      qsr = -(dq-wetc*dtjcl)*von*fdg/(lzq10 - psiq) 
        
      ! bulk loop --------------------------------------------------------------
      do i=1,nits-1
         zet = von*grav*zu/ta*(tsr*(one+qra*q)+qra*ta*qsr)/(usr*usr)/(one+qra*q) 
         zo  = charn*usr*usr/grav+zovu*visa/usr
         rr  = zo*usr/visa 
         L   = zu/zet 
         if (L <= zero) then
            ! wind
            z     = zu/L
            xx    = sqrt(sqrt(one-fifteen*z))
            psik  = log((one+xx)*(one+xx)*(one+xx*xx)) - atan(xx) + a1
            xx    = (one - p1*z)**onethird
            psic  = p2*log(one+xx+xx*xx) - sq3*atan((one+two*xx)/sq3) + a2
            f     = z*z/(one+z*z)
            psiuo = (one-f)*psik + f*psic
            ! temperature
            z     = zt/L
            xx    = sqrt(one-fifteen*z)
            psik  = two*log(one+xx) - log(four)
            xx    = (one - p3*z)**onethird
            psic  = oneph*log(one+xx+xx*xx) - sq3*atan((one+two*xx)/sq3) + a2
            f     = z*z/(one+z*z)
            psit  = (one-f)*psik + f*psic
            ! humidity
            z     = zq/L
            xx    = sqrt(one-fifteen*z)
            psik  = two*log(one+xx) - log(four)
            xx    = (one - p3*z)**onethird
            psic  = oneph*log(one+xx+xx*xx) - sq3*atan((one+two*xx)/sq3) + a2
            f     = z*z/(one+z*z)
            psiq  = (one-f)*psik + f*psic
         else
            ! wind
            z     = zu/L
            psiuo = -((one+z) + twothird*(z-c0)*exp(-min(cmax, cslope*z)) + c1)
            ! temperature
            z    =zt/L
            xx   =(one+twothird*z)
            psit =-(xx*sqrt(xx) + twothird*(z-c0)*exp(-min(cmax,cslope*z)) + c1)
            ! humidity
            z    =zq/L
            xx   =(one+twothird*z)
            psiq =-(xx*sqrt(xx) + twothird*(z-c0)*exp(-min(cmax,cslope*z)) + c1)
         endif

         usr = ut*von/(lzu10 - psiuo) 
         qsr = -(dq-wetc*dtjcl)*von*fdg/(lzq10 - psiq) 
         tsr = -(dt-dtjcl)*von*fdg/(lzt10 - psit) 

         if (i == 1) then
            Bf = -grav*usr*(tsr + qra*ta*qsr)/ta
            if (Bf > zero) then
               ug = Beta*(Bf*zi)**xlpow2 
            else
               ug = ug0
            endif
            ut = sqrt(u*u+ug*ug) 
         endif
      enddo !bulk iter loop
      zet = von*grav*zu/ta*(tsr*(one+qra*q)+qra*ta*qsr)/(usr*usr)/(one+qra*q) 
      zo  = charn*usr*usr/grav+zovu*visa/usr
      rr  = zo*usr/visa 
      L   = zu/zet 
      zoq = min(zoqmax, zoqfac/rr**zoqpow) 
      zot = zoq 
      if (L <= zero) then
         ! temperature
         z    = zt/L
         xx   = sqrt(one-fifteen*z)
         psik = two*log(one+xx) - log(four)
         xx   = (one - p3*z)**onethird
         psic = oneph*log(one+xx+xx*xx) - sq3*atan((one+two*xx)/sq3) + a2
         f    = z*z/(one+z*z)
         psit = (one-f)*psik + f*psic
         ! humidity
         z    = zq/L
         xx   = sqrt(one-fifteen*z)
         psik = two*log(one+xx) - log(four)
         xx   = (one - p3*z)**onethird
         psic = oneph*log(one+xx+xx*xx) - sq3*atan((one+two*xx)/sq3) + a2
         f    = z*z/(one+z*z)
         psiq = (one-f)*psik + f*psic
      else
         ! temperature
         z    = zt/L
         xx   = (one+twothird*z)
         psit = -(xx*sqrt(xx) + twothird*(z-c0)*exp(-min(cmax,cslope*z)) + c1)
         ! humidity
         z    = zq/L
         xx   = (one+twothird*z)
         psiq = -(xx*sqrt(xx) + twothird*(z-c0)*exp(-min(cmax,cslope*z)) + c1)
      endif
      !******************* End bulk loop ***************************************
    
      ! Calculate and return the temperature value at std ht (10m) -------------
      tzs = ts -      dtjcl                                                    &
               -      (dt-dtjcl)*fdg*(log(ten/zot)-psit)/(lzt10-psit)
      qzs = qs - wetc*dtjcl                                                    &
               - (dq-wetc*dtjcl)*fdg*(log(ten/zoq)-psiq)/(lzq10-psiq)

   end subroutine t10m

   !----------------------------------------------------------------------------

   subroutine t10m_init()

      implicit none

      sq3 = sqrt(three)
      a1  = pi*half - three*log(two)
      a2  = pi/sq3 - log(three*sq3)

      zu = ten
      zq = two
      zt = two
      zi = 600.0

      Ribcu = -zu/(zi*b4*Beta**3) 

   end subroutine t10m_init

   !----------------------------------------------------------------------------

end module mod_t10m
