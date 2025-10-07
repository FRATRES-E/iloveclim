#include "choixcomposantes.h"

!=======================================================================
      MODULE carbonate_speciation_mod
!=======================================================================

! Adapted and extended from routines taken from the SolveSAPHE library
! version 1.0.3.
! For details about the adopted algorithms and the original source code,
! please refer to
!  - Munhoven, G.: Mathematics of the total alkalinity-pH equation --
!    pathway to robust and universal solution algorithms: the SolveSAPHE
!    package v1.0.1, Geosci. Model Dev., 6, 1367-1388,
!    doi:10.5194/gmd-6-1367-2013, 2013.
!  - Munhoven, Guy. (2020, April 15). SolveSAPHE (Solver Suite for
!    Alkalinity-PH Equations) (Version 1.0.3). Geosci. Model Dev.
!    Zenodo. doi:10.5281/zenodo.3752633
!
! SolveSAPHE is free software and distributed under the GNU LGPL3 license.

! Author: Guy Munhoven
! Last modified: 18th May 2020


      PRIVATE

      PUBLIC :: incche
      PUBLIC :: calc_co3sat_arag, calc_co3sat_calc


                                    ! Set the internally used type of
                                    ! real floating point variables
      INTEGER, PARAMETER :: wp = KIND(1.0D+00)


! --------------------------------------
! Parameters for usage within the module
! --------------------------------------

! Gas constant
! ------------

      REAL(KIND=wp), PARAMETER :: gasconst_bar_cm3_o_mol_k = 83.14472_wp ! Handbook (2007)

! 0 degrees centigrade in Kelvin
! ------------------------------

      REAL(KIND=wp), PARAMETER :: t_k_zerodegc = 273.15_wp ! Handbook (2007)


! ----------------------------------------------------------------
! Chemical constants' products: for usage by members of the module
! ----------------------------------------------------------------

      REAL(KIND=wp) :: api1_dic, api2_dic
      REAL(KIND=wp) :: api1_bor
      REAL(KIND=wp) :: api1_wat



      CONTAINS

!=======================================================================
      FUNCTION A_BTOT_SALIN(s)
!=======================================================================

      !!Function returns total borate concentration in mol/kg-SW
      ! given the salinity of a sample

      ! References: Uppström (1974), cited by  Dickson et al. (2007, chapter 5, p 10)
      !             Millero (1982) cited in Millero (1995)
      ! pH scale  : N/A


      IMPLICIT NONE

      REAL(KIND=wp) :: A_BTOT_SALIN


      ! ------------------
      ! Argument variables
      ! ------------------

      REAL(KIND=wp), INTENT(IN) :: s


      A_BTOT_SALIN = 0.000416_wp*(s/35._wp)

      RETURN

!=======================================================================
      END FUNCTION A_BTOT_SALIN
!=======================================================================



!=======================================================================
      FUNCTION A_CATOT_SALIN(s)
!=======================================================================

      !!Function returns total calcium concentration in mol/kg-SW
      ! given the salinity of a sample

      ! References: Culkin and Cox (1966)
      ! pH scale  : N/A


      IMPLICIT NONE

      REAL(KIND=wp) :: A_CATOT_SALIN


      ! ------------------
      ! Argument variables
      ! ------------------

      REAL(KIND=wp), INTENT(IN) :: s

      ! (g Ca/kg)/Cl_permil values
      ! Culkin (1967):                         0.0213
      ! Culkin and Cox (DSR 13 1966):          0.02126 +/- 0.00004 (stdev)
      ! Riley and Tongudai (Chem Geol 2 1967): 0.02128 +/- 0.00006 (stdev)
      ! Handbook (2007):                       0.02127
      !    with reference to Riley and Tongudai (1967) (???)

      ! Here:
      ! (g Ca/kg)/Cl_permil = 0.02127
      ! (g Ca)/(mol Ca)     = 40.078
      !  Cl_permil          = S/1.80655
      ! mol Ca/kg = (0.02127/40.078) * (35/1.80655)
      !A_CATOT_SALIN = 0.010282_wp*(s/35._wp)
      A_CATOT_SALIN = (0.02127_wp/40.078_wp) * (s/1.80655_wp)


      RETURN

!=======================================================================
      END FUNCTION A_CATOT_SALIN
!=======================================================================



!=======================================================================
      FUNCTION A_FTOT_SALIN(s)
!=======================================================================

      !!Function returns total calcium concentration in mol/kg-SW
      ! given the salinity of a sample

      ! References: Culkin (1965) (???)
      ! pH scale  : N/A


      IMPLICIT NONE

      REAL(KIND=wp) :: A_FTOT_SALIN


      ! ------------------
      ! Argument variables
      ! ------------------

      REAL(KIND=wp), INTENT(IN) :: s


      A_FTOT_SALIN = 0.000068_wp*(s/35._wp)


      RETURN

!=======================================================================
      END FUNCTION A_FTOT_SALIN
!=======================================================================



!=======================================================================
      FUNCTION A_SO4TOT_SALIN(s)
!=======================================================================

      !!Function returns total sulfate concentration in mol/kg-SW
      ! given the salinity of a sample

      ! References: Morris, A.W. and Riley, J.P. (1966) quoted in Handbook (2007)
      ! pH scale  : N/A


      IMPLICIT NONE

      REAL(KIND=wp) :: A_SO4TOT_SALIN


      ! ------------------
      ! Argument variables
      ! ------------------

      REAL(KIND=wp), INTENT(IN) :: s


      !A_SO4TOT_SALIN = 0.028234_wp*(s/35._wp) ! in libthdyct and Thesis
      !A_SO4TOT_SALIN = 0.02824_wp*(s/35._wp)                ! Handbook (2007, chap 6, p 10, tab 2, col 3)
      A_SO4TOT_SALIN = (0.1400_wp/96.062_wp)*(s/1.80655_wp)  ! Handbook (2007, chap 6, p 10)

      RETURN

!=======================================================================
      END FUNCTION A_SO4TOT_SALIN
!=======================================================================




      !=======================================================================
       FUNCTION AK_CARB_0_WEIS74(t_k, s)
      !=======================================================================

      !!Function calculates K0 in (mol/kg-SW)/atmosphere

      ! References: Weiss (1979) [(mol/kg-SW)/atm]
      ! pH scale  : N/A
      ! Note      : currently no pressure correction


      IMPLICIT NONE

      REAL(KIND=wp) :: AK_CARB_0_WEIS74


      ! ------------------
      ! Argument variables
      ! ------------------

      !     s      : salinity
      !     t_k    : temperature in K

      REAL(KIND=wp), INTENT(IN) :: t_k
      REAL(KIND=wp), INTENT(IN) :: s


      ! ---------------
      ! Local variables
      ! ---------------

      !     zt_k_o_100   : zt_k/100

      REAL(KIND=wp) :: zt_k_o_100



      zt_k_o_100 = t_k/100._wp

      AK_CARB_0_WEIS74                                                       &
                = EXP( -60.2409_wp + 93.4517_wp/zt_k_o_100                   &
                      + 23.3585_wp*LOG(zt_k_o_100)                           &
                      + (   0.023517_wp - 0.023656_wp*zt_k_o_100             &
                         + 0.0047036_wp*zt_k_o_100*zt_k_o_100)*s )


      RETURN

      !=======================================================================
       END FUNCTION AK_CARB_0_WEIS74
      !=======================================================================



      !=======================================================================
       FUNCTION AK_CARB_1_MILL95(t_k, s, p_bar)
      !=======================================================================
      !!Function calculates first dissociation constant of carbonic acid
      ! in mol/kg-SW on the SWS pH-scale.

      ! References: Millero (1995, eq 50 -- ln K1(COM))
      !             Millero (1982) pressure correction
      ! pH scale:   SWS


      IMPLICIT NONE

      REAL(KIND=wp) :: AK_CARB_1_MILL95


      ! ------------------
      ! Argument variables
      ! ------------------

      !     t_k    : temperature in Kelvin
      !     s      : salinity
      !     p_bar  : applied pressure in bar

      REAL(KIND=wp), INTENT(IN) :: t_k
      REAL(KIND=wp), INTENT(IN) :: s
      REAL(KIND=wp), INTENT(IN) :: p_bar


      ! ---------------
      ! Local variables
      ! ---------------

      !     zrt            : R*t_k, R in bar*cm3/(mol*K)
      !     zt_degc        : temperature in degrees Celsius
      !     zdvi           : volume change for ionization
      !     zdki           : compressibility change for ionization
      !     zsqrts         : square root of salinity
      !     zds            : salinity-34.8
      !     zln_kc1_p0     : ln(K_C1) at p_bar = 0
      !     zln_kc1_pp     : pressure correction for p_bar /= 0

      REAL(KIND=wp) :: zrt, zt_degc, zdvi, zdki, zds, zsqrts
      REAL(KIND=wp) :: zln_kc1_p0, zln_kc1_pp


      ! ln(K_C1) value at p_bar = 0

      zsqrts     = SQRT(s)

      zln_kc1_p0 =      2.18867_wp - 2275.0360_wp/t_k - 1.468591_wp*LOG(t_k) &
                   + (-0.138681_wp -   9.33291_wp/t_k)*zsqrts                &
                   +  0.0726483_wp*s                                         &
                   - 0.00574938_wp*s*zsqrts


      ! Pressure correction

      zt_degc    = t_k - t_k_zerodegc
      zds        = s - 34.8_wp
      zrt        = gasconst_bar_cm3_o_mol_k * t_k

      zdvi       =  -25.50_wp - 0.151_wp*zds + 0.1271_wp*zt_degc
      zdki       = ( -3.08_wp - 0.578_wp*zds + 0.0877_wp*zt_degc)*1.0E-03_wp

      zln_kc1_pp = (-zdvi + zdki*p_bar/2._wp)*p_bar/zrt


      ! Final K_C1 value

      AK_CARB_1_MILL95 = EXP( zln_kc1_p0 + zln_kc1_pp )

      RETURN

      !=======================================================================
       END FUNCTION AK_CARB_1_MILL95
      !=======================================================================



      !=======================================================================
       FUNCTION AK_CARB_2_MILL95(t_k, s, p_bar)
      !=======================================================================
      !!Function calculates second dissociation constant K1
      ! in mol/kg-SW on the SWS pH-scale.

      ! References: Millero (1995, eq 51 -- ln K2(COM))
      !             Millero (1979) pressure correction
      ! pH scale:   SWS


      IMPLICIT NONE

      REAL(KIND=wp) :: AK_CARB_2_MILL95


      ! Argument variables
      ! ------------------

      !     t_k    : temperature in Kelvin
      !     s      : salinity
      !     p_bar  : applied pressure in bar

      REAL(KIND=wp), INTENT(IN) :: t_k
      REAL(KIND=wp), INTENT(IN) :: s
      REAL(KIND=wp), INTENT(IN) :: p_bar


      ! Local variables
      ! ---------------

      !     zrt            : R*t_k, R in bar*cm3/(mol*K)
      !     zt_degc        : temperature in degrees Celsius
      !     zdvi           : volume change for ionization
      !     zdki           : compressibility change for ionization
      !     zsqrts         : square root of salinity
      !     zds            : salinity-34.8
      !     zln_kc2_p0     : ln(K_C2) at p_bar = 0
      !     zln_kc2_pp     : pressure correction for p_bar /= 0

      REAL(KIND=wp) :: zrt, zt_degc, zdvi, zdki, zds, zsqrts
      REAL(KIND=wp) :: zln_kc2_p0, zln_kc2_pp


      ! ln(K_C2) value at p_bar = 0

      zsqrts     = SQRT(s)

      zln_kc2_p0 =     -0.84226_wp - 3741.1288_wp/t_k - 1.437139_wp*LOG(t_k) &
                   + (-0.128417_wp -  24.41239_wp/t_k)*zsqrts                &
                   +  0.1195308_wp*s                                         &
                   - 0.00912840_wp*s*zsqrts


      ! Pressure correction

      zt_degc    = t_k - t_k_zerodegc
      zds        = s - 34.8_wp
      zrt        = gasconst_bar_cm3_o_mol_k * t_k

      zdvi       =  -15.82_wp + 0.321_wp*zds - 0.0219_wp*zt_degc
      zdki       =  ( 1.13_wp - 0.314_wp*zds - 0.1475_wp*zt_degc)*1.0E-03_wp

      zln_kc2_pp =  (-zdvi + zdki*p_bar/2._wp)*p_bar/zrt


      ! Final K_C2 value

      AK_CARB_2_MILL95  = EXP( zln_kc2_p0 + zln_kc2_pp )

      RETURN

      !=======================================================================
       END FUNCTION AK_CARB_2_MILL95
      !=======================================================================



      !=======================================================================
       FUNCTION AK_BORA_DICK90(t_k, s, p_bar)
      !=======================================================================
      !!Function calculates boric acid dissociation constant KB
      ! in mol/kg-SW on the total pH-scale.

      ! References: Dickson (1990, eq. 23) -- also Handbook (2007, eq. 37)
      !             Millero (1979) pressure correction
      ! pH scale  : total


      IMPLICIT NONE

      REAL(KIND=wp) :: AK_BORA_DICK90


      ! ------------------
      ! Argument variables
      ! ------------------

      !     t_k    : temperature in Kelvin
      !     s      : salinity
      !     p_bar  : applied pressure in bar

      REAL(KIND=wp), INTENT(IN) :: t_k
      REAL(KIND=wp), INTENT(IN) :: s
      REAL(KIND=wp), INTENT(IN) :: p_bar


      ! ---------------
      ! Local variables
      ! ---------------

      !     zrt          : R*t_k, R in bar*cm3/(mol*K)
      !     zt_degc      : temperature in degrees Celsius
      !     zdvi         : volume change for ionization
      !     zdki         : compressibility change for ionization
      !     zsqrts       : square root of salinity
      !     zds          : salinity-34.8
      !     zln_kb_p0    : K_b at p_bar = 0
      !     zln_kb_pp    : pressure correction for p_bar /= 0

      REAL(KIND=wp) :: zrt, zt_degc, zdvi, zdki, zds, zsqrts
      REAL(KIND=wp) :: zln_kb_p0, zln_kb_pp


      ! ln(K_B) value at p_bar = 0

      zsqrts     = SQRT(s)

      zln_kb_p0  = ( -8966.90_wp                                                     &
                                + zsqrts*( -2890.53_wp                               &
                                + zsqrts*(  -77.942_wp                               &
                                + zsqrts*(    1.728_wp - 0.0996_wp*zsqrts)))) / t_k  &
                     +  148.0248_wp + zsqrts*(137.1942_wp + zsqrts*1.62142_wp)       &
                     + (-24.4344_wp + zsqrts*(-25.085_wp - zsqrts*0.2474_wp)) * LOG(t_k) &
                     + 0.053105_wp*zsqrts*t_k


      ! Pressure correction

      zt_degc   = t_k - t_k_zerodegc
      zds       = s - 34.8_wp
      zrt       = gasconst_bar_cm3_o_mol_k * t_k

      zdvi      = -29.48_wp + 0.295_wp*zds + 0.1622_wp*zt_degc - 0.002608_wp*zt_degc*zt_degc
      zdki      = (-2.84_wp + 0.354_wp*zds)*1.0E-03_wp

      zln_kb_pp =  (-zdvi + zdki*p_bar/2._wp)*p_bar/zrt


      ! Final K_B value

      AK_BORA_DICK90   = EXP( zln_kb_p0 + zln_kb_pp )

      !=======================================================================
       END FUNCTION AK_BORA_DICK90
      !=======================================================================



      !=======================================================================
       FUNCTION AK_W_MILL95(t_k, s, p_bar)
      !=======================================================================

      !!Function calculates water dissociation constant Kw in (mol/kg-SW)^2

      ! References: Millero (1995) for value at p_bar = 0
      !             Millero (pers. comm. 1996) for pressure correction
      ! pH scale  : SWS


      IMPLICIT NONE

      REAL(KIND=wp) :: AK_W_MILL95


      ! ------------------
      ! Argument variables
      ! ------------------

      !     t_k    : temperature in K
      !     s      : salinity
      !     p_bar  : applied pressure in bar

      REAL(KIND=wp), INTENT(IN) ::  t_k, s, p_bar


      ! ---------------
      ! Local variables
      ! ---------------

      !     zrt        : R*t_k
      !     zt_degc    : temperature in degrees Celsius
      !     zdvi       : volume change for ionization
      !     zdki       : compressibility change for ionization
      !     zln_kw_p0  : ln(K_w) at p_bar = 0
      !     zln_kw_pp  : pressure correction for p_bar /= 0

      REAL(KIND=wp) :: zrt, zt_degc, zdvi, zdki
      REAL(KIND=wp) :: zln_kw_p0, zln_kw_pp


      ! ln(K_w) value at p_bar = 0

      zln_kw_p0 =    148.9802_wp                                                     &
                   - 13847.26_wp/t_k                                                 &
                   -  23.6521_wp*LOG(t_k)                                            &
                   + ( -5.977_wp + 118.67_wp/t_k + 1.0495_wp*LOG(t_k))*SQRT(s)       &
                   - 0.01615_wp*s


      ! Pressure correction

      zt_degc = t_k - t_k_zerodegc
      zrt     = gasconst_bar_cm3_o_mol_k * t_k

      zdvi    =  -20.02_wp + 0.1119_wp*zt_degc - 0.1409E-02_wp*zt_degc*zt_degc
      zdki    = ( -5.13_wp + 0.0794_wp*zt_degc)*1.0E-03_wp

      zln_kw_pp =  (-zdvi + zdki*p_bar/2._wp)*p_bar/zrt


      ! Final K_w value

      AK_W_MILL95 = EXP( zln_kw_p0 + zln_kw_pp )


      RETURN

      !=======================================================================
       END FUNCTION AK_W_MILL95
      !=======================================================================



      !=======================================================================
       FUNCTION ASP_CALC_MUCC83(t_k, s, p_bar)
      !=======================================================================

      !!Function returns stoechiometric solubility product
      ! of calcite in seawater

      ! References: Mucci (1983)
      !             Millero (1995) for pressure correction
      ! pH scale  : N/A
      ! Units     : (mol/kg-SW)^2


      IMPLICIT NONE

      REAL(KIND=wp) :: ASP_CALC_MUCC83

      ! ------------------
      ! Argument variables
      ! ------------------

      !     t_k    : temperature in K
      !     s      : salinity
      !     p_bar  : applied pressure in bar

      REAL(KIND=wp), INTENT(IN) ::  t_k, s, p_bar

      ! ---------------
      ! Local variables
      ! ---------------

      !     zrt          : R*t_k, R in bar*cm3/(mol*K)
      !     zsqrts       : square root of salinity
      !     zt_degc      : temperature in degrees Celsius
      !     zdvi         : volume change for ionization
      !     zdki         : compressibility change for ionization
      !     zln_kp1_p0   : ln(K_p1) at p_bar = 0
      !     zln_kp1_pp   : pressure correction for p_bar /= 0

      REAL(KIND=wp) :: zrt, zsqrts, zt_degc, zdvi, zdki
      REAL(KIND=wp) :: zlog10_kspcalc_p0, zln_kspcalc_pp


      zsqrts    = SQRT(s)

      ! log10(Ksp_Calc) for p_bar = 0
      zlog10_kspcalc_p0 = &
                    -171.9065_wp -   0.077993_wp*t_k                         &
                  +  2839.319_wp/t_k + 71.595_wp*LOG10(t_k)                  &
                  + (-0.77712_wp +  0.0028426_wp*t_k + 178.34_wp/t_k)*zsqrts &
                  -   0.07711_wp*s                                           &
                  + 0.0041249_wp*s*zsqrts


      ! Pressure correction
      zt_degc   = t_k - t_k_zerodegc
      zrt       = gasconst_bar_cm3_o_mol_k * t_k

      zdvi      =  -48.76_wp + 0.5304_wp*zt_degc
      zdki      = (-11.76_wp + 0.3692_wp*zt_degc)*1.0E-03_wp

      zln_kspcalc_pp = (-zdvi + zdki*p_bar/2._wp)*p_bar/zrt


      ! Final value of Ksp_Calc

      ASP_CALC_MUCC83 = 10._wp**(zlog10_kspcalc_p0) * EXP(zln_kspcalc_pp)

      RETURN

      !=======================================================================
       END FUNCTION ASP_CALC_MUCC83
      !=======================================================================



      !=======================================================================
       FUNCTION ASP_ARAG_MUCC83(t_k, s, p_bar)
      !=======================================================================

      !!Function returns stoechiometric solubility product
      ! of aragonite in seawater

      ! References: Mucci (1983)
      !             Millero (1979) for pressure correction
      ! pH scale  : N/A
      ! Units     : (mol/kg-SW)^2


      IMPLICIT NONE

      REAL(KIND=wp) :: ASP_ARAG_MUCC83

      ! ------------------
      ! Argument variables
      ! ------------------

      !     t_k    : temperature in K
      !     s      : salinity
      !     p_bar  : applied pressure in bar

      REAL(KIND=wp), INTENT(IN) ::  t_k, s, p_bar

      ! ---------------
      ! Local variables
      ! ---------------

      !     zrt          : R*t_k, R in bar*cm3/(mol*K)
      !     zsqrts       : square root of salinity
      !     zt_degc      : temperature in degrees Celsius
      !     zdvi         : volume change for ionization
      !     zdki         : compressibility change for ionization
      !     zln_kp1_p0   : ln(K_p1) at p_bar = 0
      !     zln_kp1_pp   : pressure correction for p_bar /= 0

      REAL(KIND=wp) :: zrt, zsqrts, zt_degc, zdvi, zdki
      REAL(KIND=wp) :: zlog10_ksparag_p0, zln_ksparag_pp


      zsqrts    = SQRT(s)

      ! log10(Ksp_Arag) for p_bar = 0
      zlog10_ksparag_p0 = &
                      -171.945_wp   - 0.077993_wp*t_k                        &
                  +   2903.293_wp/t_k + 71.595_wp*LOG10(t_k)                 &
                  + (-0.068393_wp +  0.0017276_wp*t_k + 88.135_wp/t_k)*zsqrts &
                  -    0.10018_wp*s                                          &
                  +  0.0059415_wp*s*zsqrts


      ! Pressure correction
      zt_degc   = t_k - t_k_zerodegc
      zrt       = gasconst_bar_cm3_o_mol_k * t_k

      zdvi      =  -48.76_wp + 0.5304_wp*zt_degc  + 2.8_wp
      zdki      = (-11.76_wp + 0.3692_wp*zt_degc)*1.0E-03_wp

      zln_ksparag_pp = (-zdvi + zdki*p_bar/2._wp)*p_bar/zrt


      ! Final value of Ksp_Arag

      ASP_ARAG_MUCC83 = 10._wp**(zlog10_ksparag_p0) * EXP(zln_ksparag_pp)

      RETURN

      !=======================================================================
       END FUNCTION ASP_ARAG_MUCC83
      !=======================================================================



      !=======================================================================
       FUNCTION ABETA_HF_DIRI79(t_k, s, p_bar)
      !=======================================================================

      !!Function calculates association constant \beta_{HF} [(mol/kg-SW)^{-1}]
      ! in (mol/kg-SW)^{-1}, where
      !   \beta_{HF} = \frac{ [HF] }{ [H^{+}] [F^{-}] }

      ! References: Dickson and Riley (1979)
      !             Millero (1995) for pressure correction
      ! pH scale  : free
      ! Note      : converted here from mol/kg-H2O to mol/kg-SW

      IMPLICIT NONE

      REAL(KIND=wp) :: ABETA_HF_DIRI79


      ! ------------------
      ! Argument variables
      ! ------------------

      !     t_k    : temperature in K
      !     s      : salinity
      !     p_bar  : applied pressure in bar

      REAL(KIND=wp), INTENT(IN) ::  t_k, s, p_bar


      ! ---------------
      ! Local variables
      ! ---------------

      !     zrt            : R*t_k, R in bar*cm3/(mol*K)
      !     zt_degc        : temperature in degrees Celsius
      !     zdvi           : volume change for ionization
      !     zdki           : compressibility change for ionization
      !     zionst         : ionic strength [mol/kg-H2O]
      !     zcvt_to_kgsw   : mass of pure water in 1kg of seawater as a fct. of salinity
      !     zln_bhf_p0     : \beta_HF at p_bar = 0
      !     zln_khf_pp     : pressure correction for k_HF = 1/\beta_HF at p_bar /= 0


      REAL(KIND=wp) :: zrt, zt_degc, zdvi, zdki
      REAL(KIND=wp) :: zionst, zcvt_to_kgsw
      REAL(KIND=wp) :: zln_bhf_p0, zln_khf_pp


      ! \beta_HF at p_bar = 0
      ! ---------------------

      zcvt_to_kgsw    = ACVT_KGH2O_O_KGSW(s)
      zionst          = A_IONSTRENGTH_SALIN(s)/zcvt_to_kgsw

      zln_bhf_p0      = -1590.2_wp/t_k + 12.641_wp - 1.525_wp*SQRT(zionst)


      ! Pressure correction
      ! -------------------

      zt_degc      = t_k - t_k_zerodegc
      zrt          = gasconst_bar_cm3_o_mol_k * t_k

      zdvi         =   -9.78_wp + zt_degc*(-0.0090_wp - zt_degc*0.942E-03_wp)
      zdki         = ( -3.91_wp + zt_degc*0.054_wp)*1.0E-03_wp

      zln_khf_pp   = (-zdvi + zdki*p_bar/2._wp)*p_bar/zrt


      ! Final \beta_HF value
      ! --------------------
      !  notice that  ln(k_HF(P)) = ln(k_HF(0)) + zln_khf_pp
      !         <=>  -ln(\beta_HF(P)) = -ln(\beta_HF(0)) + zln_khf_pp
      !         <=>   ln(\beta_HF(P)) =  ln(\beta_HF(0)) - zln_khf_pp

      ABETA_HF_DIRI79 = EXP(zln_bhf_p0 - zln_khf_pp ) / zcvt_to_kgsw

      RETURN

      !=======================================================================
       END FUNCTION ABETA_HF_DIRI79
      !=======================================================================



      !=======================================================================
       FUNCTION AK_HSO4_DICK90(t_k, s, p_bar)
      !=======================================================================

      !!Function returns the dissociation constant of hydrogen sulfate (bisulfate)

      ! References: Dickson (1990) -- also Handbook (2007)
      !             Millero (1995) for pressure correction
      ! pH scale  : free
      ! Note      : converted here from mol/kg-H2O to mol/kg-SW

      IMPLICIT NONE

      REAL(KIND=wp) :: AK_HSO4_DICK90


      ! ------------------
      ! Argument variables
      ! ------------------

      !     t_k    : temperature in K
      !     s      : salinity
      !     p_bar  : applied pressure in bar

      REAL(KIND=wp), INTENT(IN) ::  t_k, s, p_bar


      ! ---------------
      ! Local variables
      ! ---------------

      !     zrt            : R*t_k, R in bar*cm3/(mol*K)
      !     zt_degc        : temperature in degrees Celsius
      !     zdvi           : volume change for ionization
      !     zdki           : compressibility change for ionization
      !     zionst         : ionic strength in mol/-kg-H2O
      !     zsqrti         : square root og ion strength
      !     zcvt_to_kgsw   : mass of pure water in 1kg of seawater as a fct. of salinity
      !     zln_khso4_p0   : K_HSO4 at p_bar = 0
      !     zln_khso4_pp   : pressure correction for p_bar /= 0

      REAL(KIND=wp) :: zrt, zt_degc, zdvi, zdki
      REAL(KIND=wp) :: zcvt_to_kgsw, zionst, zsqrti
      REAL(KIND=wp) :: zln_khso4_p0, zln_khso4_pp


      ! ln(K_HSO4) at p_bar = 0

      zcvt_to_kgsw = ACVT_KGH2O_O_KGSW(s)
      zionst       = A_IONSTRENGTH_SALIN(s)/zcvt_to_kgsw
      zsqrti       = SQRT(zionst)

      zln_khso4_p0 =    -4276.1_wp/t_k + 141.328_wp -  23.093_wp*LOG(t_k)           &
                     + (-13856._wp/t_k +  324.57_wp -  47.986_wp*LOG(t_k)) * zsqrti &
                     + ( 35474._wp/t_k -  771.54_wp + 114.723_wp*LOG(t_k)) * zionst &
                     - (  2698._wp/t_k)*zsqrti * zionst                             &
                     + (  1776._wp/t_k)*zionst*zionst


      ! Pressure correction

      zt_degc      = t_k - t_k_zerodegc
      zrt          = gasconst_bar_cm3_o_mol_k * t_k

      zdvi         =  -18.03_wp + zt_degc*(0.0466_wp + zt_degc*0.316E-03_wp)
      zdki         = ( -4.53_wp + zt_degc*0.0900_wp)*1.0E-03_wp

      zln_khso4_pp = (-zdvi + zdki*p_bar/2._wp)*p_bar/zrt


      ! ln(K_HSO4) at p_bar = 0

      AK_HSO4_DICK90 = zcvt_to_kgsw * EXP( zln_khso4_p0 + zln_khso4_pp )

      RETURN

      !=======================================================================
       END FUNCTION AK_HSO4_DICK90
      !=======================================================================



      !=======================================================================
       FUNCTION ACVT_HSWS_O_HTOT(t_k, s, p_bar)
      !=======================================================================

      !!Function returns the ratio H_SWS/H_Tot as a function of salinity s

      ! Reference:  Munhoven
      ! pH scale:   all


      IMPLICIT NONE

      REAL(KIND=wp) :: ACVT_HSWS_O_HTOT


      ! ------------------
      ! Argument variables
      ! ------------------

      !     t_k    : temperature in K
      !     s      : salinity
      !     p_bar  : applied pressure in bar

      REAL(KIND=wp), INTENT(IN) ::  t_k
      REAL(KIND=wp), INTENT(IN) ::  s
      REAL(KIND=wp), INTENT(IN) ::  p_bar


      ! ---------------
      ! Local variables
      ! ---------------

      !     zso4_tot: total sulfate concentration in mol/kg-SW
      !     zf_tot  : total fluoride concentration in mol/kg-SW

      REAL(KIND=wp) :: zso4_tot, zf_tot

      !-----------------------------------------------------------------------


      zso4_tot = A_SO4TOT_SALIN(s)
      zf_tot   = A_FTOT_SALIN(s)


      ACVT_HSWS_O_HTOT = 1._wp +  (zf_tot*ABETA_HF_DIRI79(t_k, s, p_bar)) &
                                 /(1._wp + zso4_tot/AK_HSO4_DICK90(t_k,s, p_bar))

      RETURN


      !=======================================================================
       END FUNCTION ACVT_HSWS_O_HTOT
      !=======================================================================



      !=======================================================================
       FUNCTION ACVT_KGH2O_O_KGSW(s)
      !=======================================================================

      !!Function returns the mass of pure water in one kg of seawater
      ! of salinity s

      ! References: "libthdyct" -- derived by Munhoven (1997) from data by Millero (1982)
      !             "Handbook (2007)" -- Handbook (2007)
      ! pH scale:   N/A


      IMPLICIT NONE

      REAL(KIND=wp) :: ACVT_KGH2O_O_KGSW

      REAL(KIND=wp), INTENT(IN) :: s

      !ACVT_KGH2O_O_KGSW = 1._wp - 0.0010049_wp*s ! libthdyct
      ACVT_KGH2O_O_KGSW = 1._wp - 0.001005_wp*s ! Handbook (2007)

      RETURN


      !=======================================================================
       END FUNCTION ACVT_KGH2O_O_KGSW
      !=======================================================================



      !=======================================================================
       FUNCTION A_IONSTRENGTH_SALIN(s)
      !=======================================================================

      !!Function calculates ionic strength in mol/kg-SW, for given salinity.

      ! References: "libthdyct" -- derived by Munhoven (1997) from data by Millero (1982)
      !             "Handbook (2007)" -- Handbook (2007)
      ! pH scale:   N/A


      IMPLICIT NONE

      REAL(KIND=wp) :: A_IONSTRENGTH_SALIN


      ! ------------------
      ! Argument variables
      ! ------------------

      REAL(KIND=wp), INTENT(IN) :: s



      !A_IONSTRENGTH_SALIN = (0.019920D+00*s) ! libthdyct
      A_IONSTRENGTH_SALIN = (0.019924D+00*s) ! Handbook (2007)

      RETURN

      !=======================================================================
       END FUNCTION A_IONSTRENGTH_SALIN
      !=======================================================================



      !=======================================================================
       SUBROUTINE SETUP_API4PHSWS(t_k, s, p_bar)
      !=======================================================================

      ! ------------------
      ! Argument variables
      ! ------------------

      !     t_k    : temperature in Kelvin
      !     s      : salinity
      !     p_bar  : applied pressure in bar

      REAL(KIND=wp), INTENT(IN) :: t_k
      REAL(KIND=wp), INTENT(IN) :: s
      REAL(KIND=wp), INTENT(IN) :: p_bar

      ! ---------------
      ! Local variables
      ! ---------------

      REAL(KIND=wp) :: zcvt_hsws_o_htot


      zcvt_hsws_o_htot  = ACVT_HSWS_O_HTOT(t_k, s, p_bar)


      api1_dic =            AK_CARB_1_MILL95(t_k, s, p_bar)
      api2_dic = api1_dic * AK_CARB_2_MILL95(t_k, s, p_bar)

      api1_bor =            AK_BORA_DICK90(t_k, s, p_bar) * zcvt_hsws_o_htot

      api1_wat =            AK_W_MILL95(t_k, s, p_bar)


      RETURN

      !=======================================================================
       END SUBROUTINE SETUP_API4PHSWS
      !=======================================================================




      !=======================================================================
       FUNCTION SOLVE_ACBW_POLY(p_alkcbw, p_dictot, p_bortot, p_hini, p_val)
      !=======================================================================

      !!Function
      !  - determines the positive root of the DIC - B_T - OH-H - A_CBW
      !    equation for [H+].
      !  - returns -1 if divergent (may only happen if the maximum number
      !    of iterations is exceeded.

      ! The solution is based upon the resolution of the quintic polynomial
      ! equation. Safe bracketting is implemented, and the resolution is
      ! carried out in the pH-Alk space, similarly to SOLVE_ACBW_GENERAL.

      ! Author: Guy Munhoven
      ! Last modified: 18th May 2020


      IMPLICIT NONE


      !------------------------------------!
      ! General parameters of the function !
      !------------------------------------!

      ! Type the function itself
      REAL(KIND=wp)            :: SOLVE_ACBW_POLY

      REAL(KIND=wp), PARAMETER :: pp_rdel_ah_target     = 1.E-6_wp
      INTEGER, PARAMETER       :: jp_maxniter_acbw_poly = 50
      INTEGER                  :: niter_acbw_poly       = 0


      !--------------------!
      ! Argument variables !
      !--------------------!

      REAL(KIND=wp), INTENT(IN)            :: p_alkcbw
      REAL(KIND=wp), INTENT(IN)            :: p_dictot
      REAL(KIND=wp), INTENT(IN)            :: p_bortot
      REAL(KIND=wp), INTENT(IN),  OPTIONAL :: p_hini
      REAL(KIND=wp), INTENT(OUT), OPTIONAL :: p_val


      !-----------------!
      ! Local variables !
      !-----------------!

      REAL(KIND=wp)            ::  za4, za3, za2, za1, za0
      REAL(KIND=wp)            ::  zh_ini, zh_min, zh_max, zh, zh_prev
      REAL(KIND=wp)            ::  zeqn, zdeqndh, zh_delta, zh_lnfactor, zeqn_absmin

      LOGICAL                  :: l_exitnow

                                              ! Threshold value for switching from
                                              ! pH-space to [H^+]-space iterations.
      REAL(KIND=wp), PARAMETER :: pz_exp_threshold = 1.0_wp


      !==============================================================================


      niter_acbw_poly = 0


      IF(PRESENT(p_hini)) THEN                ! If an initial value is provided, use it, ...

         zh_ini = p_hini

      ELSE                                    ! if not, get it from AHINI_FOR_ACBW

         CALL AHINI_FOR_ACBW(p_alkcbw, p_dictot, p_bortot, zh_ini)

      ENDIF

                                              ! Get the brackets for the root
       CALL ACBW_HINFSUP(p_alkcbw, p_dictot, p_bortot, zh_min, zh_max)



      zh   = MAX(MIN(zh_max, zh_ini), zh_min) ! Bracket the initial value

                                              ! Coefficients of the quintic polynomial
                                              ! Use a leading coefficient of -1, so that
                                              ! polynomial is oriented the same way as
                                              ! the rational function equation for H>>
                                              ! (decreasing towards -\infty
      za4 = -(p_alkcbw + api1_dic + api1_bor)
      za3 =   api1_dic*(p_dictot - p_alkcbw - api1_bor) - api2_dic &
            + api1_bor*(p_bortot - p_alkcbw) + api1_wat
      za2 =   api1_dic*api1_bor*(p_dictot + p_bortot - p_alkcbw) &
            + api2_dic*(p_dictot+p_dictot - p_alkcbw - api1_bor) &
            + api1_wat*(api1_dic + api1_bor)
      za1 =   api2_dic*(api1_bor*(p_dictot+p_dictot + p_bortot - p_alkcbw) + api1_wat) &
            + api1_dic*api1_bor*api1_wat
      za0 =   api2_dic*api1_bor*api1_wat



      zeqn_absmin = HUGE(1._wp)               ! Preset the current best minimum value
                                              ! found to the highest possible value

      DO

         IF(niter_acbw_poly >= jp_maxniter_acbw_poly) THEN

            zh = -1._wp
            EXIT

         ENDIF

         zh_prev = zh

         zeqn = za0 + zh*(za1 + zh*(za2 + zh*(za3 + zh*(za4 - zh))))

         IF(zeqn > 0._wp) THEN
            zh_min = zh_prev
         ELSEIF(zeqn < 0._wp) THEN
            zh_max = zh_prev
         ELSE
            EXIT                              ! zh is the root; unlikely but, one never knows
         ENDIF


         ! Now start to calculate the next iterate zh

         niter_acbw_poly = niter_acbw_poly + 1


         zdeqndh = za1 + zh*(2._wp*za2 + zh*(3._wp*za3 + zh *(4._wp*za4 - 5._wp*zh)))


         IF(ABS(zeqn) >= 0.5_wp*zeqn_absmin) THEN

            zh = SQRT(zh_max * zh_min)

                                              ! zh_lnfactor required to test convergence below
            zh_lnfactor = (zh - zh_prev)/zh_prev

         ELSE

            zh_lnfactor = -zeqn/(zdeqndh*zh_prev)

            IF(ABS(zh_lnfactor) > pz_exp_threshold) THEN
               zh       = zh_prev*EXP(zh_lnfactor)
            ELSE
               zh_delta = zh_lnfactor*zh_prev
               zh       = zh_prev + zh_delta
            ENDIF


            IF( zh < zh_min ) THEN

               zh          = SQRT(zh_min * zh_max)

                                               ! zh_lnfactor required to test convergence below
               zh_lnfactor = (zh - zh_prev)/zh_prev

            ENDIF


            IF( zh > zh_max ) THEN

               zh          = SQRT(zh_min * zh_max)

                                               ! zh_lnfactor required to test convergence below
               zh_lnfactor = (zh - zh_prev)/zh_prev

            ENDIF


         ENDIF

         zeqn_absmin = MIN( ABS(zeqn), zeqn_absmin)

         l_exitnow = (ABS(zh_lnfactor) < pp_rdel_ah_target)

         IF(l_exitnow) EXIT

      ENDDO


      SOLVE_ACBW_POLY = zh


      IF(zh > -1._wp) THEN

         IF(PRESENT(p_val)) p_val = za0 + zh*(za1 + zh*(za2 + zh*(za3 + zh*(za4 - zh))))

      ELSE

         IF(PRESENT(p_val)) p_val = HUGE(1._wp)

      ENDIF


      RETURN


      CONTAINS


      !-----------------------------------------------------------------------
       SUBROUTINE AHINI_FOR_ACBW(p_alkcb, p_dictot, p_bortot, p_hini)
      !-----------------------------------------------------------------------

      ! Subroutine returns the root for the 2nd order approximation of the
      ! DIC -- B_T -- A_CB equation for [H+] (reformulated as a cubic polynomial)
      ! around the local minimum, if it exists.

      ! Returns * 1E-03_wp if p_alkcb <= 0
      !         * 1E-10_wp if p_alkcb >= 2*p_dictot + p_bortot
      !         * 1E-07_wp if 0 < p_alkcb < 2*p_dictot + p_bortot
      !                    and the 2nd order approximation does not have a solution


      IMPLICIT NONE


      !--------------------!
      ! Argument variables !
      !--------------------!

      REAL(KIND=wp), INTENT(IN)             ::  p_alkcb, p_dictot, p_bortot
      REAL(KIND=wp), INTENT(OUT)            ::  p_hini


      !-----------------!
      ! Local variables !
      !-----------------!

      REAL(KIND=wp)  ::  zca, zba
      REAL(KIND=wp)  ::  zd, zsqrtd, zhmin
      REAL(KIND=wp)  ::  za2, za1, za0


      !-----------------------------------------------------------------------


      IF (p_alkcb <= 0._wp) THEN
        p_hini = 1.e-3_wp
      ELSEIF (p_alkcb >= (2._wp*p_dictot + p_bortot)) THEN
        p_hini = 1.e-10_wp
      ELSE
        zca = p_dictot/p_alkcb
        zba = p_bortot/p_alkcb

        ! Coefficients of the cubic polynomial
        za2 = api1_bor*(1._wp - zba) + api1_dic*(1._wp-zca)
        za1 = api1_dic*api1_bor*(1._wp - zba - zca) + api2_dic*(1._wp - (zca+zca))
        za0 = api2_dic*api1_bor*(1._wp - zba - (zca+zca))


                                              ! Taylor expansion around the minimum

        zd = za2*za2 - 3._wp*za1              ! Discriminant of the quadratic equation
                                              ! for the minimum close to the root

        IF(zd > 0._wp) THEN                   ! If the discriminant is positive

          zsqrtd = SQRT(zd)

          IF(za2 < 0) THEN
            zhmin = (-za2 + zsqrtd)/3._wp
          ELSE
            zhmin = -za1/(za2 + zsqrtd)
          ENDIF

          p_hini = zhmin + SQRT(-(za0 + zhmin*(za1 + zhmin*(za2 + zhmin)))/zsqrtd)

        ELSE

          p_hini = 1.e-7_wp

        ENDIF

      ENDIF


      RETURN

      !-----------------------------------------------------------------------
       END SUBROUTINE AHINI_FOR_ACBW
      !-----------------------------------------------------------------------



      !-----------------------------------------------------------------------
       SUBROUTINE ACBW_HINFSUP(p_alkcbw, p_dictot, p_bortot, p_hinf, p_hsup)
      !-----------------------------------------------------------------------

      ! Subroutine returns the lower and upper brackets of the root


      IMPLICIT NONE


      !--------------------!
      ! Argument variables !
      !--------------------!

      REAL(KIND=wp), INTENT(IN)  :: p_alkcbw
      REAL(KIND=wp), INTENT(IN)  :: p_dictot
      REAL(KIND=wp), INTENT(IN)  :: p_bortot
      REAL(KIND=wp), INTENT(OUT) :: p_hinf
      REAL(KIND=wp), INTENT(OUT) :: p_hsup


      !-----------------!
      ! Local variables !
      !-----------------!

      REAL(KIND=wp) :: zalknw_inf, zalknw_sup
      REAL(KIND=wp) :: zsqrtdelta


      !-----------------------------------------------------------------------


      ! Infimum and supremum for the ALK_ACBW not related
      ! to water self-ionization.

      zalknw_inf = 0._wp
      zalknw_sup = p_dictot+p_dictot + p_bortot


      ! Lower bound for the root

      zsqrtdelta = SQRT((p_alkcbw-zalknw_inf)**2 + 4._wp*api1_wat)

      IF(p_alkcbw >= zalknw_inf) THEN
         p_hinf = 2._wp*api1_wat /( p_alkcbw-zalknw_inf + zsqrtdelta )
      ELSE
         p_hinf = (-(p_alkcbw-zalknw_inf) + zsqrtdelta ) / 2._wp
      ENDIF


      ! Upper bound for the root

      zsqrtdelta = SQRT((p_alkcbw-zalknw_sup)**2 + 4._wp*api1_wat)

      IF(p_alkcbw <= zalknw_sup) THEN
         p_hsup = (-(p_alkcbw-zalknw_sup) + zsqrtdelta ) / 2._wp
      ELSE
         p_hsup = 2._wp*api1_wat /( p_alkcbw-zalknw_sup + zsqrtdelta )
      ENDIF


      RETURN

      !-----------------------------------------------------------------------
       END SUBROUTINE ACBW_HINFSUP
      !-----------------------------------------------------------------------


      !=======================================================================
       END FUNCTION SOLVE_ACBW_POLY
      !=======================================================================



      !=======================================================================
       SUBROUTINE SPECIATION_DIC(p_dictot, p_h, p_co2, p_hco3, p_co3)
      !=======================================================================

      ! Subroutine returns the speciation of the carbonate system


      !------------------------------!
      ! Chemical constants' products !
      !------------------------------!
      ! - api1_dic = K_1
      ! - api2_dic = K_1*K_2


      IMPLICIT NONE


      !--------------------!
      ! Argument variables !
      !--------------------!

      REAL(KIND=wp), INTENT(IN)  :: p_dictot
      REAL(KIND=wp), INTENT(IN)  :: p_h
      REAL(KIND=wp), INTENT(OUT) :: p_co2
      REAL(KIND=wp), INTENT(OUT) :: p_hco3
      REAL(KIND=wp), INTENT(OUT) :: p_co3


      !-----------------!
      ! Local variables !
      !-----------------!

      REAL(KIND=wp) :: z_dictot_over_denom


      !=======================================================================


      z_dictot_over_denom = p_dictot / ( api2_dic + p_h*( api1_dic + p_h) )


      IF (p_h < api1_dic)  THEN              ! CO_2 is dominant

        p_co3  = api2_dic     * z_dictot_over_denom
        p_hco3 = api1_dic*p_h * z_dictot_over_denom

        p_co2  = p_dictot - (p_hco3 + p_co3)

      ELSEIF(api1_dic*p_h < api2_dic) THEN   ! HCO_3^- is dominant

        p_co3  = api2_dic     * z_dictot_over_denom
        p_co2  =      p_h*p_h * z_dictot_over_denom

        p_hco3 = p_dictot - (p_co2 + p_co3)

      ELSE                                   ! CO_3^2- is dominant

        p_hco3 = api1_dic*p_h * z_dictot_over_denom
        p_co2  =      p_h*p_h * z_dictot_over_denom

        p_co3  = p_dictot - (p_co2 + p_hco3)

      ENDIF


      !=======================================================================
       END SUBROUTINE SPECIATION_DIC
      !=======================================================================



      !=======================================================================
       FUNCTION CALC_CO3SAT_CALC(t_k, s, p_bar)
      !=======================================================================


      !!Function returns the CO3 concentration at saturation with respect
      ! to calcite, in [mol/kg-SW], assuming that the total Ca
      ! concentration is conservative (i.e., can be derived from salinity)

      IMPLICIT NONE

      ! Type the function itself
      REAL(KIND=wp)            :: CALC_CO3SAT_CALC

      !--------------------!
      ! Argument variables !
      !--------------------!

      REAL, INTENT(IN)  :: t_k             ! temperature [K]
      REAL, INTENT(IN)  :: s               ! salinity [-]
      REAL, INTENT(IN)  :: p_bar           ! applied pressure [p_bar]


      !-----------------!
      ! Local variables !
      !-----------------!

      REAL(KIND=wp) :: z_t_k, z_s, z_p_bar


      z_t_k   = t_k
      z_s     = s
      z_p_bar = p_bar


      CALC_CO3SAT_CALC = ASP_CALC_MUCC83(z_t_k, z_s, z_p_bar) / A_CATOT_SALIN(z_s)


      RETURN

      !=======================================================================
       END FUNCTION CALC_CO3SAT_CALC
      !=======================================================================



      !=======================================================================
       FUNCTION CALC_CO3SAT_ARAG(t_k, s, p_bar)
      !=======================================================================

      !!Function returns the CO3 concentration at saturation with respect
      ! to aragonite, in [mol/kg-SW], assuming that the total Ca
      ! concentration is conservative (i.e., can be derived from salinity)


      IMPLICIT NONE

      ! Type the function itself
      REAL(KIND=wp)            :: CALC_CO3SAT_ARAG

      !--------------------!
      ! Argument variables !
      !--------------------!

      REAL, INTENT(IN)  :: t_k             ! temperature [K]
      REAL, INTENT(IN)  :: s               ! salinity [-]
      REAL, INTENT(IN)  :: p_bar           ! applied pressure [p_bar]


      !-----------------!
      ! Local variables !
      !-----------------!

      REAL(KIND=wp) :: z_t_k, z_s, z_p_bar


      z_t_k   = t_k
      z_s     = s
      z_p_bar = p_bar


      CALC_CO3SAT_ARAG = ASP_ARAG_MUCC83(z_t_k, z_s, z_p_bar) / A_CATOT_SALIN(z_s)


      RETURN

      !=======================================================================
       END FUNCTION CALC_CO3SAT_ARAG
      !=======================================================================



      !=======================================================================
       SUBROUTINE INCCHE(t_k, s, p_bar, xDIC,xALK, sCO2,xpCO2,xCO2,xHCO3,xCO3, z_h)
      !=======================================================================

      ! pH solver and carbonate system speciation calculation for iLOVECLIM-OCCYC
      ! Author: G. Munhoven
      ! Last modified: 18th May 2020


      IMPLICIT NONE


      !--------------------!
      ! Argument variables !
      !--------------------!

      REAL, INTENT(IN)  :: t_k             ! temperature [K]
      REAL, INTENT(IN)  :: s               ! salinity [-]
      REAL, INTENT(IN)  :: p_bar           ! applied pressure [p_bar]

      REAL, INTENT(IN)  :: xDIC            ! Total Dissolved Inorganic Carbon [mol/kg-SW]
      REAL, INTENT(IN)  :: xALK            ! Total alkalinity [eq/kg-SW, mol/kg-SW]
      REAL, INTENT(OUT) :: xpCO2           ! CO2 partial pressure [atm]
      REAL, INTENT(OUT) :: sCO2            ! Henry's constant for CO2 [(mol/kg-SW)/atm]
      REAL, INTENT(OUT) :: xCO2            ! (CO_2)_(aq) concentration [mol/kg-SW]
      REAL, INTENT(OUT) :: xHCO3           ! HCO_3 concentration [mol/kg-SW]
      REAL, INTENT(OUT) :: xCO3            ! CO_3 concentration [mol/kg-SW]
      REAL, INTENT(OUT) :: z_h             ! [H+] (for PH)


      !-----------------!
      ! Local variables !
      !-----------------!

      REAL(KIND=wp) :: z_t_k, z_s, z_p_bar
      REAL(KIND=wp) :: z_dictot, z_alkcbw, z_bortot
      !REAL(KIND=wp) :: z_h
      REAL(KIND=wp) :: z_co2, z_hco3, z_co3
      REAL(KIND=wp) :: z_k0



      z_t_k    = t_k
      z_s      = s
      z_p_bar  = p_bar

      z_dictot = xDIC
      z_alkcbw = xALK

                                          ! Prepare the constants' products
                                          ! for the given T, S, P conditions
      CALL SETUP_API4PHSWS(z_t_k, z_s, z_p_bar)

      z_bortot = A_BTOT_SALIN(z_s)        ! Total boron concentration

                                          ! Solve for pH (here as [H^+])
      z_h = SOLVE_ACBW_POLY(z_alkcbw, z_dictot, z_bortot)

                                          ! Calculate CO2-HCO3-CO3 speciation
      CALL SPECIATION_DIC(z_dictot, z_h, z_co2, z_hco3, z_co3)

                                          ! Henry's constant for CO2
      z_k0 = AK_CARB_0_WEIS74(z_t_k, z_s)


      sCO2  = z_k0                        ! Finalize rsults and transfer
      xpCO2 = z_co2/z_k0                  ! to dummy variables (possibly type change)
      xCO2  = z_co2
      xHCO3 = z_hco3
      xCO3  = z_co3


      RETURN

      !=======================================================================
       END SUBROUTINE incche
      !=======================================================================

      !=======================================================================
       END MODULE carbonate_speciation_mod
      !=======================================================================

