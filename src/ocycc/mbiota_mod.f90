!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:21 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:21 CET 2009

! Modif pour le Pa/Th, dmr&lim
! Updated to be cleaner and with an eye on performance - 2020-06-12

!nb TO CHECK:
!delete sigma_m?


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
module mbiota_mod
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

use declars_mod, only: lt, jt, noc_cbr, jx
use global_constants_mod, only: dblp=>dp, ip
use mod_sync_time, only: tstoc, tday

use para0_mod, only: NISOO2
use iso_dioxygen_mod, only: iair

#if ( OOISO ==1 )
use iso_dioxygen_mod, only: iair16, iair17, iair18
use iso_dioxygen_mod, only: Ototal_min, r18smow, r17smow
use iso_dioxygen_mod, only: compute_ISOO2_mbiodyn
use iso_dioxygen_mod, only: GPPO2_func, Ray_respO2
use iso_dioxygen_mod, only: OO2_saturation
#endif

implicit none

REAL(kind=dblp) :: n0_m, PI_m, PAR_m, lef_m ,pmin_m, dp_m, er_doc, a_m, b_m, c_m, p0_m, g0_m, ex_doc, dz_m, zmin_m, l0_m, d0_m, &
                   kd_m, rdoc_ms, O2min1, O2min2

REAL(kind=dblp), dimension(LT,JT,NOC_CBR) :: PHYTO_M, ZOO_M
REAL(kind=dblp), dimension(LT,JT,NOC_CBR) :: PHYTO_M13, ZOO_M13

REAL(kind=dblp), dimension(LT,JT,NOC_CBR,NISOO2) :: prod_O2, resp_O2

REAL(kind=dblp), dimension(LT,NOC_CBR) :: TPP_m = 0.0_dblp, TPP_D13C = 0.0_dblp, caco3_d13C = 0.0_dblp
REAL(kind=dblp), dimension(LT,NOC_CBR) :: caco3_m = 0.0_dblp
REAL(kind=dblp), dimension(LT,NOC_CBR) :: caco3_m_b_sh = 0.0_dblp

REAL(kind=dblp), dimension(JX) :: SUE_MCA !nb tbd?, fhypso

!REFACTORING DONE: OrgCFlxAttFactor replaces SUE_M
REAL(kind=dblp), dimension(LT, JX, NOC_CBR) :: OrgCFlxAttFactor


#if ( REMIN == 1 )
REAL(kind=dblp), dimension(LT,JT,NOC_CBR) :: SUE_3D !3D remineralisation
!REAL(kind=dblp), dimension(JT) :: R !remineralisation rate
!REAL(kind=dblp) :: kremin !remineralisation rate
REAL(kind=dblp), dimension(LT,JT,NOC_CBR) :: kremin !remineralisation rate
REAL(kind=dblp) :: Ea, Rgaz, betaPOM
#endif
#if ( REMIN == 1 || REMIN_CACO3 == 1 )
REAL (kind=dblp) :: dt !day
REAL(kind=dblp), dimension(JT) :: w_sink !sinking particle velocity
#endif
#if ( REMIN_CACO3 == 1 )
REAL(kind=dblp), dimension(LT,JT,NOC_CBR) :: kremin_ca !remineralisation rate
REAL(kind=dblp), dimension(LT,JT,NOC_CBR) :: kremin_ar !remineralisation rate
REAL(kind=dblp), dimension(LT,JT,NOC_CBR) :: SUE_ca_3D !3D remineralisation
REAL(kind=dblp), dimension(LT,JT,NOC_CBR) :: SUE_ar_3D !3D remineralisation
REAL(kind=dblp) :: k_diss
#endif
REAL(kind=dblp), dimension(JT) :: SWR_FRAC !dmr --- transferred to do it once
!cnb  caco3 pour exp.f
REAL(kind=dblp), dimension(LT,JT,NOC_CBR) :: caco3

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

REAL(kind=dblp) :: prod_sum = 0.0_dblp, graz_sum = 0.0_dblp, remin_sum = 0.0_dblp, exc_sum = 0.0_dblp, exu_sum = 0.0_dblp,      &
                   barem_sum = 0.0_dblp, tpp_sum = 0.0_dblp, ort_sum = 0.0_dblp, pel_sum = 0.0_dblp &
                  ,total_car, total_phos  &
                  ,total_adv, SCALE_M,SCALE_B,C14RA,RR,SCANU, caco3_sum = 0.0_dblp, fPOC_top = 0.0_dblp, fPOC_1000 = 0.0_dblp,  &
                   fPOC_bot = 0.0_dblp, fCAL_top = 0.0_dblp, fCAL_2000 = 0.0_dblp, fCAL_bot = 0.0_dblp

REAL(kind=dblp), dimension(LT,JT,NOC_CBR) :: phyto_day_sum = 0.0_dblp, prodO2_sum = 0.0_dblp, respO2_sum = 0.0_dblp,            &
                                             reminO2 = 0.0_dblp, NCP_sum = 0.0_dblp
#if ( ARAG == 1 )
REAL(kind=dblp) :: RR_ar, Kmax
REAL(kind=dblp), dimension(LT,NOC_CBR) :: caco3_d13C_ar = 0.0_dblp
REAL(kind=dblp), dimension(JX) ::SUE_MAR
#if ( CORAL == 0 )
REAL(kind=8), dimension(LT,JT,NOC_CBR), target :: omega_arag3D
#endif
#endif
#if ( OXNITREUX == 1 )
REAL(kind=dblp) :: err_ON2O
#endif

REAL(kind=dblp) :: total_tpp

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr [NOTA] POINTERS for MEDUSA coupling ... is this the right place to put them?

REAL(kind=dblp), dimension(:,:,:), pointer :: temp_ma, salt_ma, odic_ma, oalk_ma, ooxy_ma
REAL(kind=dblp), dimension(:,:,:), pointer :: ono3_ma, opo4_ma, oc13_ma, oc14_ma

!dmr [NOTA] TPP_ma & caco3_ma are likely incorrect when not using MEDUSA since
!                 the zero reset is done in MEDUSA only(!!)
REAL(kind=dblp), dimension(LT,JT,NOC_CBR) :: clay_ma, tracer01_ma, TPP_ma = 0.0_dblp, caco3_ma = 0.0_dblp
REAL(kind=dblp), dimension(LT,NOC_CBR)    :: caco3_mabot = 0.0_dblp

REAL(kind=dblp), dimension(LT,JT,NOC_CBR):: TPPC13_ma, caco3C13_ma, TPPC14_ma, caco3C14_ma, caco3O18_ma

REAL(kind=dblp) :: total_caco3

REAL(kind=dblp) :: b_sh_fr ! dmr == big shells fraction


REAL(kind=dblp), dimension(LT,JT,NOC_CBR) :: TPP_mave, caco3_mave, temp_mave, salt_mave, odic_mave, oalk_mave, clay_mave   &
               , tracer01_mave, ooxy_mave, ono3_mave, opo4_mave, TPPC13_mave,caco3C13_mave, TPPC14_mave,caco3C14_mave      &
               , oc13_mave, oc14_mave, oo17_mave, oo18_mave, oohd_mave, caco3O18_mave

REAL(kind=dblp), dimension(LT,NOC_CBR) :: TPP_mafond, caco3_mafond, temp_mafond, salt_mafond, ooxy_mafond, tracer01_mafond &
               , ono3_mafond, opo4_mafond, ohco3_mafond, oco3_mafond, oco2_mafond, caco3O18_mafond, rdf_c_mafond           &
               , rdf_n_mafond, rdf_p_mafond, rdf_ro2_mafond, clay_mafond, TPPC13_mafond, caco3C13_mafond, TPPC14_mafond    &
               , caco3C14_mafond, oc13_mafond, oc14_mafond, oo17_mafond, oo18_mafond,oohd_mafond

REAL(kind=dblp) :: summary_flux_O2S_calc, summary_flux_O2S_orgm, summary_flux_S2O_dic ! in Tmol.yr-1

!tbd REAL(kind=dblp) :: riverine_input_dic = 0.1833D+2                     &               ! in Tmols.yr-1
!tbd                 , riverine_input_alk = 0.1202D+2                     &               ! in Tmols.yr-1
!tbd                 , riverine_input_PO4 = 0.597674 !????

REAL(kind=dblp) :: riverine_dic_input_tot         &
                 , riverine_alk_input_tot         &
                 , riverine_o2_input_tot          &
                 , riverine_NO3_input_tot         &
                 , riverine_PO4_input_tot
#ifdef WITH_C13
REAL(kind=dblp) :: riverine_dic13_input_tot != 0.0
#endif
#ifdef WITH_C14
REAL(kind=dblp) :: riverine_dic14_input_tot != 0.0
#endif


!! variables to store fluxes from sediment to ocean !mohr

REAL(kind=dblp), dimension(LT,NOC_CBR) :: ODOCS_sed2oc, ODIC_sed2oc, OALK_sed2oc, OO2_sed2oc, ONO3_sed2oc, OPO4_sed2oc     &
               , OO17_sed2oc, OO18_sed2oc, OOHD_sed2oc, ODIC13_sed2oc, ODIC14_sed2oc

!! organic matter loss through sediment bottom for reinjection by rivers
!! flowing from continent ("loop-back")

REAL(kind=dblp), dimension(LT,NOC_CBR) :: orgm_dic_loopback, orgm_no3_loopback, orgm_po4_loopback, orgm_oxyg_loopback      &
               , orgm_alk_loopback, calc_loopback, orgm_dic13_loopback, orgm_dic14_loopback, calc13_loopback               &
               , calc14_loopback


#if ( PATH >= 1 )
REAL(kind=dblp), dimension(LT,JT,NOC_CBR) :: TPP_path, caco3_path
#endif

LOGICAL, DIMENSION(LT,JT,NOC_CBR) :: oc_bottom_cell

!REFACTORING: too many hardcoded constants here about time control
REAL(kind=dblp),  PARAMETER :: TIME_STM = TDAY/10._dblp
INTEGER(kind=ip), PARAMETER :: nb_step_bio=NINT(86400.0_dblp/TIME_STM+0.001_dblp) !TSTOC = 86400 s / TIME_STM=T_DAY/10


contains


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
  subroutine mbiodyn(im, nm, RAD_M, watTempArray, watTempSurf, OXCO2, oc13bio, OC13, ODIC, &
            OPO4, ODOC, ODOCS, ODOC13, ODOCS13, PHYTO_PROD, OC14, ONO3, OO2, OALK, DVOL, ZZ, ZX, &
            FRICE, MGT, O2_sat_thistime)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!REFACTORING: the following need to be moved from marine_bio_mod to this
!             module as they are only being used here;
!             also transfer their initialisation, if not done at declaration
!             into an initialisation subroutine herein.
    use marine_bio_mod, ONLY: eher, zinges, ecan, sigma_md

    use marine_bio_mod, ONLY: OetaC_POMoxid, OetaN_POMoxid, OetaO2_POMoxid,          &
                              OetaC_DOMoxid_1D, OetaN_DOMoxid_1D, OetaO2_DOMoxid_1D, &
                              sigma_m

! [???] TDAY is already used above in the parameter declaration of TIME_STM
    use mod_sync_time, only: tday


    use marine_bio_mod, only: JPROD

#if ( OOISO == 1 )
    use para0_mod, ONLY: NISOO2     ! required to set the dimensions of the dummy argument OO2
#endif

#if (REMIN_CACO3 == 1 )
    use omega_mod, only: calc_omega_ca, omega_calc3D
#endif

#if ( ARAG == 1 )
    use omega_mod, only: calc_omega_ar, omega_arag3D
#endif


#if ( IRON_LIMITATION == 1 )
       use loveclim_transfer_mod, only: IRON_LIM
#endif

    implicit none


    ! Dummy arguments
    integer(kind=ip), intent(in) :: IM, NM
    real(kind=dblp),  intent(in) :: RAD_M, watTempSurf
    real(kind=dblp), dimension(JT), intent(in) :: watTempArray
    real, dimension(LT,NOC_CBR), intent(in) :: OXCO2
    real, intent(out) :: oc13bio
    real, dimension(LT,JT,NOC_CBR), intent(inout) :: OC13
    real(kind=8), dimension(LT,JT,NOC_CBR), intent(inout) :: ODIC
    real, dimension(LT,JT,NOC_CBR), intent(inout) :: OPO4
    real, dimension(LT,JT,NOC_CBR), intent(inout) :: ODOC
    real, dimension(LT,JT,NOC_CBR), intent(inout) :: ODOCS
    real, dimension(LT,JT,NOC_CBR), intent(inout) :: ODOC13
    real, dimension(LT,JT,NOC_CBR), intent(inout) :: ODOCS13
    double precision, dimension(LT,JT,NOC_CBR), intent(out) :: phyto_prod
    real, dimension(LT,JT,NOC_CBR), intent(inout) :: OC14
    real, dimension(LT,JT,NOC_CBR), intent(inout) :: ONO3

    real(kind=8), dimension(LT,JT,NOC_CBR,NISOO2), intent(inout) :: OO2

    real(kind=8), dimension(LT,JT,NOC_CBR), intent(inout) :: OALK
    real(kind=dblp), dimension(LT,JT,NOC_CBR), intent(in) :: DVOL
    real(kind=dblp), dimension(JT), intent(in) :: ZZ
    real(kind=dblp), dimension(JX), intent(in) :: ZX
    real(kind=dblp), dimension(LT,NOC_CBR), intent(in) :: FRICE
    integer(kind=ip), dimension(LT,JT,NOC_CBR), intent(in) :: MGT
    real(kind=dblp), dimension(lt,jprod,noc_cbr), intent(in) :: O2_sat_thistime


    integer(kind=ip):: timeloop, j
    real(kind=dblp) :: FUNT_M, FUNL_M, RATE_M, G_M, RDOC_M, O2_sat, watTemp
    real(kind=dblp) :: avanut     ! avanut: AVAilable NUTrients

    ! Local variables
    real(kind=dblp) :: deltaC13, DeltaC14
    real(kind=dblp) :: phyto_growth, phyto_senesc, phyto_DOCexudat
    real(kind=dblp) :: grazing, grazing_max
    real(kind=dblp) :: phyto_zoousage, phyto_egest, phyto_remin
    real(kind=dblp) :: zoo_growth, zoo_mortal, zoo_DOCexcret
    real(kind=dblp) :: zoo_egest, zoo_remin
    real(kind=dblp) :: doc_productot
    real(kind=dblp) :: DOC_product, DOC_decay, DOC13_decay
    real(kind=dblp) :: DOCS_product, DOCS_decay, DOCS13_decay
    real(kind=dblp) :: DIC_dif_fromPOC, DIC_dif_fromDOCtot, DIC_dif_fromPZD
    real(kind=dblp) :: NPP_O2, NCP, OO2_before
    real(kind=dblp) :: cvt_mumolC_kg_to_GtC
    real(kind=dblp) :: residual_O2

! temporary var
    REAL(kind=dblp) :: PHYTO_dif, ZOO_dif, TPP_dif, ODOC_dif, ODOCS_dif, &
                       ODOC13_dif, ODOCS13_dif, calred, caco3_dif, vdic_POC, vdic_DOC
#if ( ARAG == 1 )
!    REAL(kind=dblp) :: caco3_m_ar
      REAL(kind=dblp) :: caco3_dif_ar
      REAL(kind=dblp), dimension(LT,NOC_CBR) :: caco3_m_ar = 0.0_dblp
#endif

#if ( OOISO == 1 )
    REAL(kind=dblp), dimension(NISOO2) :: OO2_netflux
    REAL(kind=dblp), dimension(NISOO2) :: OO2_netflux_prod, OO2_netflux_consum
    REAL(kind=dblp), dimension(NISOO2) :: OO2_sati
#endif
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


    ! N-P-Z-D Model following Six and Maier-Reimer (1996,
    ! Paleoceanography, DOI:10.1029/96GB02561).
    !
    ! Structure:
    ! We start with the P-Z-D parts and derive the N parts afterwards.

    ! A summary after the P-Z-D parts are done shows how the different
    ! exchange terms must add up to have a consistent system
    ! (search for "Summary of the P-Z-D model terms").
    ! Relevant equations are numbered in comments in the code
    ! (search for "--> Summary Equation ($n$)" with $n$ replaced
    ! by the actual number to find them from the summary itself).


    do j = 1, JPROD
!REFACTORING: the following "if (MGT ...) then" will soon be removed;
!             the code indentation applied already takes this into account.
      if (MGT(im,j,nm).eq.1) then   ! if wet cell

#if ( ARAG == 1 )
      call calc_omega_ar(im,j,nm)
#endif
#if ( REMIN_CACO3 == 1 )
      call calc_omega_ca(im,j,nm)
#endif

      watTemp = watTempArray(j)

                                    ! Maximum phytoplanktion rate
                                    ! ---------------------------
                                    ! function r(T,L) [S&MR96, Eq. (1)] -> rate_m
                                    !  - component f(T) [S&MR96, Eq. (2)]
      funt_m = a_m*b_m**(c_m*watTemp)

      funl_m = RAD_M * SWR_FRAC(j)  !  - component f(L) [S&MR96, Eq. (3)]

                                    ! [???] Not sure why the must be safe-guarded:
                                    !       all factors are positive
      rate_m = max((funt_m*funl_m/sqrt(funt_m*funt_m+funl_m*funl_m)/TDAY), 0.0_dblp)

                                    ! Potential grazing
                                    ! -----------------
                                    ! function g(T): g_m

                                    ! g0_m = 1.0 / TDAY [S&MR96, Eq. (4)]
                                    ! c_m = 1.0 [S&MR96, Eq. (4)] same as for f(T) above

      g_m = g0_m * (1.01_dblp)**(c_m*watTemp)

                                    ! Carbonate flux
                                    ! --------------
                                    ! Limitation by local temperature
                                    ! [???] reference
      calred = exp(watTempSurf-3.5_dblp)/(5._dblp + exp(watTempSurf-3.5_dblp))

                                    ! **_day_sum reset to 0.0 umolC/day
                                    ! (start of the day)
      phyto_day_sum(im,j,nm) = 0.0
      NCP_sum(im,j,nm) = 0.0
      prodO2_sum(im,j,nm) = 0.0
      respO2_sum(im,j,nm) = 0.0

! REFACTORING -- it would be referable not to have any time loop in this
! subroutine, and certainly not with a somehow implicit time step
! (OK it is parametrized above, but it is hardcoded).
! If deemed indispensable: move the timestep length and the number of
! timesteps to the dummy arguments of this subroutine

      do timeloop = 1, nb_step_bio

                                    ! Phyto_day_sum reset to 0.0 umolC/day
                                    ! (start of the day)

        avanut = max(0.0_dblp, OPO4(im,j,nm))

!REFACTORING: better move the oc13bio calculation to a function
!REFACTORING: Notice: the types of all constants have already been
!             changed to real(dblp), although the variables oxCO2 and oc13bio
!             are still plain real
                                    ! Dependence of C13 fractionation from
                                    ! Rau et al. (1991, Paleoceanography,
                                    ! doi:10.1029/91PA00321 -- à verifier)
                                    ! oxCO2 = xCO2 dans ocn_bio.f = [CO2aq]
                                    ! unites requises pour [CO2]: umol/L
                                    ! oxCO2 = xCO2 from ocn_bio.f = [CO2aq]
        if (oxCO2(im,nm) > 1.0e-10_dblp) then
          oc13bio = -9866._dblp/(273._dblp+watTemp) + 24.12_dblp &
                    - 17._dblp*log10(oxCO2(im,nm)*1000._dblp*1000._dblp) &
                    + 3.4_dblp + OC13(im,J,nm)/ODIC(im,J,nm)
        else
          oc13bio = -25._dblp
        endif

                                    ! Convenience variables for the delta^13C
                                    ! Delta^14C of the water DIC - these are
                                    ! used in numerous places.
        deltaC13 = OC13(im,J,nm)/ODIC(im,J,nm)
        DeltaC14 = OC14(im,J,nm)/ODIC(im,J,nm)


        ! ==== N-*P*-Z-D: LIFE AND DEATH OF PHYTOPLANKTON ====
        !
        ! Notice: no phytoplankton underneath sea ice,
        ! hence application of (1 - frice) factor in phyto_prod

        phyto_prod(im,j,nm) = rate_m * avanut/(avanut + n0_m) &
                                     * phyto_m(im,j,nm)  &
                                     * (1._dblp-FRICE(im,nm)) * TIME_STM

                                    ! Units of phyto_prod: mumol/kg
                                    ! [???] The limiting term below was said to be
                                    !       essential in setting the prod_sum annually.
                                    !       --> WHY?
                                    ! [???] Was 0.1 instead of 1.0 in 2004
                                    !       --> WHY 1 now?

!REFACTORING: hardcoded not obvious absolute constant
!             --> parametrize, or use fractio of something
        if (phyto_prod(im,j,nm) > 1._dblp) phyto_prod(im,j,nm) = 1._dblp


                                    ! [???] WHY this limitation now?
                                    !       Dimensionally strange (depends - again -
                                    !       on implicit time information ???)
        if (phyto_prod(im,j,nm) >= avanut*OetaC_POMoxid(im,j,nm)) then
          phyto_prod(im,j,nm) = avanut * OetaC_POMoxid(im,j,nm)
        endif

#if ( IRON_LIMITATION == 1)
                                    ! Iron limitation (LFe) from PISCES model
                                    ! Pre-industriel period
        phyto_prod(im,j,nm) = phyto_prod(im,j,nm) * IRON_LIM(im,j,nm)
#endif


!REFACTORING: This explanation better had to be elsewhere  (or removed)
                                    ! TIME_STM is length of time step per day (1/10. probably for us ...)
!REFACTORING -- the comment says it: "probably".
!             IDEALLY: No time loop in this subroutine!

!REFACTORING: phyto_growth should be used throughout once phyto_prod is not global anymore
        phyto_growth = phyto_prod(im,j,nm)

                                    ! phyto_senesc: death due to age
                                    ! dp_m = 0.008 / TDAY
                                    ! pmin_m = minimum phytoplankton concentration = 0.01
        phyto_senesc = dp_m * (PHYTO_M(im,j,nm) - pmin_m) * TIME_STM

                                    ! phyto_DOCexudat: Phyto C loss by exudation of DOC
                                    ! er_doc = 0.03 / TDAY (was parameter gamma_P in S&MR96)
        phyto_DOCexudat = er_doc * (PHYTO_M(im,j,nm) - pmin_m) * TIME_STM

                                    ! grazing --> loss of Phyto C as a result of the
                                    ! activity of Zooplankton (direct and indirect)
                                    ! p0_m = half-saturation constant for grazing
                                    !        (value: 3.5 mumol/kg - was 4 mumolC/L in S&MR96)

        grazing = g_m * (PHYTO_M(im,j,nm) - pmin_m)  &
                            / (PHYTO_M(im,j,nm) + p0_m)  &
                            * ZOO_M(im,j,nm) * TIME_STM

                                    ! [???] Maximum amount of phytoplankton available for grazing
                                    !  - dimensionally strange (adds concentrations and variations)
                                    !  - grazing_max could theoretically even be negative ...
        grazing_max = PHYTO_M(im,j,nm) - pmin_m  &
                      + phyto_growth - phyto_senesc - phyto_DOCexudat

                                    ! If grazing exceeds available Phyto, reduce grazing
        if (grazing > grazing_max)  grazing = grazing_max


                                    ! Phytoplankton dynamics
                                    ! ----------------------

        PHYTO_dif = phyto_growth - grazing - phyto_senesc - phyto_DOCexudat
                                    ! --> Summary Equation (1)

        PHYTO_M(im,j,nm) = PHYTO_M(im,j,nm) + PHYTO_dif

                                    ! [???] With the grazing limitation above,
                                    !       this should never happen
        if (PHYTO_M(im,j,nm) < pmin_m)  then
          PHYTO_dif = pmin_m - PHYTO_M(im,j,nm)
          PHYTO_M(im,j,nm) = pmin_m
        endif

                                    ! [???] We assume C-13 equilibration
        PHYTO_M13(im,j,nm) = PHYTO_M(im,j,nm) * oc13bio


        ! ==== N-P-*Z*-D: LIFE AND DEATH OF ZOOPLANKTON ====

                                    ! Partitioning of the grazing flux:
                                    !  + phyto_egest: part of the grazing that
                                    !                 zooplankton egests as fecal
                                    !                 pellets of digested phytoplankton
                                    !  + phyto_zoousage: remainder of the grazing,
                                    !                 either used for organic growth
                                    !                 or remineralized/respired
                                    !                 by zooplankton
                                    ! Formulation:
                                    !  - phyto_zoousage = eher * grazing
                                    !    with 0.8 < eher < 1.0  (here 0.9)
                                    !    (eher was eps_her in S&MR96)
                                    !  - phyto_egest = (1 - eher) * grazing
        phyto_zoousage = eher * grazing
        phyto_egest = grazing - phyto_zoousage      ! = (1 - eher) * grazing
                                    ! --> Summary Equation (6)


                                    ! zoo_growth: organic growth of zooplankton,
                                    ! derived from the partitioning of
                                    ! phyto_zoousage into
                                    !  + zoo_growth = zinges * phyto_zoousage
                                    !    where zinges = assimilation efficiency
                                    !    zinges = 6./10
                                    !  + phyto_remin, the remainder of phyto_zoousage
                                    !    that is remineralized/respired

        zoo_growth = zinges * phyto_zoousage
        phyto_remin = phyto_zoousage - zoo_growth ! = (1 - zinges) * phyto_zoousage
                                    ! --> Summary Equation (7)

                                    ! zoo_mortal: disparition par deces
        zoo_mortal = dz_m * (ZOO_M(im,j,nm) - zmin_m) * TIME_STM

                                    ! For later usage, zoo_mortal is
                                    ! partitioned onto
                                    !  - zoo_remin = ecan * zoo_mortal,
                                    !    supposed to be remineralized as a
                                    !    zooplankton ingestion by other zooplankton
                                    !    0 < ecan < 0.05 (was eps_can in S&MR96)
                                    !  - zoo_egest, the remainder of zoo_mortal,
                                    !    egested by the same carnivore zooplankton
        zoo_remin = ecan * zoo_mortal
        zoo_egest = zoo_mortal - zoo_remin ! = (1 - ecan) * zoo_mortal
                                    ! --> Summary Equation (8)

                                    ! zoo_DOCexcret: Zooplankton C loss by DOC excretion
                                    ! ex_doc (was parameter gamma_Z in S&MR96)
        zoo_DOCexcret = ex_doc * (ZOO_M(im,j,nm) - zmin_m) * TIME_STM


                                    ! Zooplankton dynamics
                                    ! --------------------

        ZOO_dif = zoo_growth - zoo_mortal - zoo_DOCexcret
                                    ! --> Summary Equation (2)

        ZOO_M(im,j,nm) = ZOO_M(im,j,nm) + ZOO_dif

                                    ! [???] HERE WE ARE REALLY ASKING FOR REAL TROUBLE:
                                    !       when the limitation becomes effective, we
                                    !       do not know which term to apply it to - for
                                    !       phytoplankton at least, this was "grazing.
                                    !       Beware of possible conservation issues that
                                    !       could result from this limitation.
        if (ZOO_M(im,j,nm) < zmin_m)  then
          ZOO_dif = zmin_m - ZOO_M(im,j,nm)
          ZOO_M(im,j,nm) = zmin_m
        endif

                                    ! [???] We assume C-13 equilibration
                                    !       between seawter and Phyto- plus
                                    !       Zooplankton. Justification?
        ZOO_M13(im,j,nm) = ZOO_M(im,j,nm) * oc13bio


                                    ! TPP = Total Particle Production in this layer.
                                    !       Units claimed to be mumolC/kg
                                    ! [???] should rather be mumolC/kg / sub-time-step
        TPP_dif = phyto_egest + phyto_senesc + zoo_egest
                                    ! --> Summary Equation (3)


                                    ! Accumulate TPP_dif
                                    !  - over the sub-timesteps
                                    !  - over the photic zone
                                    ! to derive the export production,
                                    ! for the whole photic zone over the
                                    ! complete time step for usage in maphot
                                    ! afterwards

                                    ! SCALE_B is the conversion factor from mumol/kg
                                    ! to Tmol/m^3 (10^3 [kg-> m^3] * 10^-18 [mumol -> Tmol])
                                    ! Hence TPP_dif*SCAL_B is in Tmol.m^-3 and TPP_m is in
                                    !PROBLEM: includes hidden seawater density of 1000 kg/m^3
        TPP_m(im, nm) = TPP_m(im, nm) &
                          + TPP_dif * SCALE_B * DVOL(im, J, nm)

        TPP_D13C(im, nm) = TPP_D13C(im, nm) &
                             + TPP_dif*SCALE_B * DVOL(im, J, nm) * oc13bio

!nb CaCO3
        caco3_m(im,nm)=RR*calred*TPP_m(im,nm)*(1-b_sh_fr)/(1-sigma_m)
        caco3_m_b_sh(im,nm)=RR*calred*TPP_m(im,nm)*b_sh_fr/(1-sigma_m)
#if( ARAG == 1 )
        caco3_m_ar(im,nm)=-RR_ar*calred*TPP_m(im,nm)/(1-sigma_m)
#endif


        ! ==== N-P-Z-*D*: EVOLUTION OF DOC and DOCS ====

                                    ! Total DOC production, from
                                    ! phytoplankton exudation and
                                    ! zooplankton excretion
        doc_productot = phyto_DOCexudat + zoo_DOCexcret
                                    ! --> Summary Equation (9)


        ! fast DOC dynamics
        ! -----------------

        rdoc_m = d0_m * avanut/(avanut + kd_m)

        ! [???] Not sure what is going on here: if rdoc_m is small,
        !       reset it to a higher value --> recipe for disaster!!!
        ! [???] Should possibly be >= 0.1*d0_m ???
        !       Was, however, like this in the original
        if (rdoc_m <= (0.1_dblp*d0_m)) then
          rdoc_m = 0.1_dblp * d0_m
        endif

                                  ! production of fast DOC is only a fraction
                                  ! (1 - sigma_md) of the total DOC production
        DOC_product = (1._dblp - sigma_md) * doc_productot

                                  ! decay of fast DOC
        DOC_decay = rdoc_m * ODOC(im,j,nm) * TIME_STM

        ODOC_dif = DOC_product - DOC_decay
                                  ! --> Summary Equation (4)


        DOC13_decay = rdoc_m * ODOC13(im,j,nm) * TIME_STM

        ODOC13_dif = DOC_product * oc13bio - DOC13_decay


        ODOC(im, j, nm) = ODOC(im, j, nm) + ODOC_dif
        ODOC13(im, j, nm) = ODOC13(im, j, nm) + ODOC13_dif


        ! slow DOC dynamics
        ! -----------------

        DOCS_product = doc_productot - DOC_product
                                  ! --> Summary Equation (10)
        DOCS_decay = rdoc_ms * ODOCS(im,j,nm) * TIME_STM

        ODOCS_dif = DOCS_product - DOCS_decay
                                  ! --> Summary Equation (5)

        DOCS13_decay = rdoc_ms * ODOCS13(im,j,nm) * TIME_STM
        ODOCS13_dif= DOCS_product * oc13bio - DOCS13_decay

        ODOCS(im,j,nm) = ODOCS(im,j,nm) + ODOCS_dif
        ODOCS13(im,j,nm) = ODOCS13(im,j,nm) + ODOCS13_dif



        ! Summary of the P-Z-D model terms
        ! ================================

        !  (1) PHYTO_dif = phyto_growth - grazing - phyto_senesc - phyto_DOCexudat
        !  (2) ZOO_dif = zoo_growth - zoo_mortal - zoo_DOCexcret
        !  (3) TPP_dif = phyto_egest + phyto_senesc + zoo_egest
        !  (4) ODOC_dif = DOC_product - DOC_decay
        !  (5) ODOCS_dif = DOCS_product - DOCS_decay

        !  (6) phyto_egest = grazing - phyto_zoousage
        !  (7) phyto_remin = phyto_zoousage - zoo_growth
        !  (8) zoo_egest = zoo_mortal - zoo_remin
        !  (9) doc_productot = phyto_DOCexudat + zoo_DOCexcret
        ! (10) DOCS_product = doc_productot - DOC_product

        ! Notice: search for "--> Summary Equation ($n$)" with $n$
        ! replaced by the actual number to find the equations back
        ! in the code.

        ! Global conservation check:

        ! Add (6) + (7):
        !  (A) phyto_egest + phyto_remin = grazing - zoo_growth

        ! PHYTO_dif + ZOO_dif + ODOC_dif + ODOCS_dif + TPP_dif
        !   = phyto_growth - grazing - phyto_senesc - phyto_DOCexudat   | expand terms
        !     + zoo_growth - zoo_mortal - zoo_DOCexcret
        !     + DOC_product - DOC_decay
        !     + DOCS_product - DOCS_decay
        !     + phyto_egest + phyto_senesc + zoo_egest
        !   = phyto_growth - (grazing - zoo_growth) - phyto_senesc      | collect siblings
        !     - (phyto_DOCexudat + zoo_DOCexcret) - zoo_mortal
        !     + (DOC_product + DOCS_product) - (DOC_decay + DOCS_decay)
        !     + phyto_egest + phyto_senesc + zoo_egest
        !   = phyto_growth - (phyto_egest + phyto_remin) - phyto_senesc | introduce (A)
        !     - doc_productot - zoo_mortal                              | introduce (9)
        !     + doc_productot - (DOC_decay + DOCS_decay)                | introduce (10) and simplify
        !     + phyto_egest + phyto_senesc + zoo_egest
        !   = phyto_growth - phyto_remin - zoo_mortal                   | simplify all
        !     - (DOC_decay + DOCS_decay)
        !     + zoo_mortal - zoo_remin                                  | introduce (8)
        !   = phyto_growth - phyto_remin                                | simplify all
        !     - (DOC_decay + DOCS_decay) - zoo_remin

        ! Hence, by C conservation, we get
        !
        ! DIC_dif_fromNPZ
        !   = -(PHYTO_dif + ZOO_dif + ODOC_dif + ODOCS_dif + TPP_dif)
        !   = -(phyto_growth - phyto_remin - (DOC_decay + DOCS_decay) - zoo_remin)
        !   = phyto_remin + zoo_remin - phyto_growth - (DOC_decay + DOCS_decay)
        !
        ! to be separated in POC and DOCtot contributions to be able to
        ! consistently derive concomitant nutrient and O2 fluxes (done below)



        ! ==== N-P-Z-D PIGGYBACKING: CARBONATE PRDUCTION ====

!REFACTORING: use a positive value for caco3_m and chose a more meaningful name
                                    ! CaCO3 formation following a Rain Ratio "RR"
!nb        caco3_m = -RR * calred * TPP_dif
!nb        caco3_d13C(im, nm) = caco3_d13C(im, nm) - caco3_m * deltaC13 * DVOL(im,J,nm)*SCALE_B
!dmr --- formation suivant le Rain Ratio "RR"
         caco3_dif=-RR*calred*TPP_dif/(1-sigma_m)
         caco3_d13C(im,nm)=caco3_d13C(im,nm)-caco3_dif*(OC13(im,J,nm)/ODIC(im,J,nm))*(DVOL(im,J,nm)*SCALE_B)

                                    ! [???] Not sure the following is correct
#if ( ARAG == 1 )
! "juste une partie de CaCO3:"
! [???] Hmmm, je ne sais pas quelle est l'intention ici: "juste une partie de CaCO3"
!       alors, ce serait mieux d'utiliser caco3_m_ar = frac_ar * caco3_m, et de
!       supprimer en gros toutes les caco3_m_ar qui suivent, car caco3_m englobe alors caco3_m_ar
!       ou de n'utiliser caco3 que pour les carbonates toutes minérlogies confondues, et
!       faire clairemenbt la distinction entre calicte et aragonite.
!nb        caco3_m_ar = -RR_ar * calred * TPP_dif
!nb        caco3_d13C_ar(im,nm) = caco3_d13C_ar(im,nm)-caco3_m_ar*deltaC13*(DVOL(im,J,nm)*SCALE_B)
         caco3_dif_ar=-RR_ar*calred*TPP_dif/(1-sigma_m)
         caco3_d13C_ar(im,nm)=caco3_d13C_ar(im,nm)-caco3_dif_ar*(OC13(im,J,nm)/ODIC(im,J,nm))*(DVOL(im,J,nm)*SCALE_B)
! ou Formation d aragonite principalement des pteropodes, et depend de omega, cf Gangsto
! et al 2008
!         caco3_dif_ar=-RR_ar*TPP_dif*(omega_arag3D(im,j,nm)-1)/(Kmax+(omega_arag3D(im,j,nm)-1))
!         caco3_d13C_ar(im,nm)=caco3_d13C_ar(im,nm)-caco3_dif_ar*(OC13(im,J,nm)/ODIC(im,J,nm))*(DVOL(im,J,nm)*SCALE_B)
#endif


        ! ==== *N*-P-Z-D: THE FATE OF NUTRIENTS (DIC, PO4, NO3) ====

        ! DIC dynamics
        ! ------------

        ! Notice: the follwoing are deliberately named DIC_dif... and
        ! not ODIC_dif..., because their units are mumol/kg (same as OC)
        ! and not mol/kg (as ODIC).
        DIC_dif_fromPOC = phyto_remin + zoo_remin - phyto_growth
        DIC_dif_fromDOCtot = doc_decay + docs_decay
        DIC_dif_fromPZD = DIC_dif_fromPOC + DIC_dif_fromDOCtot


#if ( ARAG == 0 )
                                    ! DIC_dif_fromPZD > 0 is a net source of DIC
                                    ! abs(caco3_m) is a sink of DIC, but since
                                    ! caco3_m < 0 it is net source
!nb        ODIC(im,J,nm) = ODIC(im,J,nm) + (DIC_dif_fromPZD + caco3_m)*SCANU
        ODIC(im,J,nm) = ODIC(im,J,nm) + (DIC_dif_fromPZD + caco3_dif)*SCANU

                                    ! [???] I doubt the following gives good results
                                    !       as it neglects phytoplankton, zooplankton
                                    !       stock changes.
        OC13(im,J,nm) = OC13(im,J,nm) &
!nb                          + ( DIC_dif_fromPOC*oc13bio + caco3_m * deltaC13 &
                          + ( DIC_dif_fromPOC*oc13bio + caco3_dif * deltaC13 &
                                + (DOC13_decay + DOCS13_decay) ) &
                            * SCANU

                                    ! [???] C14 decay missing ???
                                    !       What about DOC14 and DOCS14 ???
!nb        OC14(im,J,nm) = OC14(im,J,nm) + (DIC_dif_fromPZD + caco3_m) * DeltaC14 * SCANU
        OC14(im,J,nm) = OC14(im,J,nm) + (DIC_dif_fromPZD + caco3_dif) * DeltaC14 * SCANU
#else
        ODIC(im,J,nm) = ODIC(im,J,nm) &
!nb                          + (DIC_dif_fromPZD + caco3_m + caco3_m_ar) * SCANU
                          + (DIC_dif_fromPZD + caco3_dif + caco3_dif_ar) * SCANU

        OC13(im,J,nm) = OC13(im,J,nm) &
!nb                          + ( DIC_dif_fromPOC*oc13bio + (caco3_m+caco3_m_ar)*deltaC13 &
                          + ( DIC_dif_fromPOC*oc13bio + (caco3_dif+caco3_dif_ar)*deltaC13 &
                                + (DOC13_decay + DOCS13_decay) ) &
                            * SCANU

        OC14(im,J,nm) = OC14(im,J,nm) &
!nb                          + (DIC_dif_fromPZD + caco3_m + caco3_m_ar) * DeltaC14 * SCANU
                          + (DIC_dif_fromPZD + caco3_dif + caco3_dif_ar) * DeltaC14 * SCANU
#endif


        ! PO4 dynamics
        ! ------------

        OPO4(im,j,nm) = &
          OPO4(im,j,nm) + DIC_dif_fromPOC/OetaC_POMoxid(im,j,nm) &
                        + DIC_dif_fromDOCtot/OetaC_DOMoxid_1D(j)


        ! NO3 dynamics
        ! ------------

        ONO3(im,j,nm) = &
          ONO3(im,j,nm) + DIC_dif_fromPOC/OetaC_POMoxid(im,j,nm)*OetaN_POMoxid(im,j,nm)                                     &
                        + DIC_dif_fromDOCtot/OetaC_DOMoxid_1D(j)*OetaN_DOMoxid_1D(j)




        ! O2 dynamics
        ! ------------
                                    ! Diagnostic variable for OOISO
        NPP_O2 = phyto_growth*ABS(OetaO2_POMoxid(im,j,nm)/OetaC_POMoxid(im,j,nm))

                                    ! Oxygen sat at this place
        O2_sat = O2_sat_thistime(im,j,nm)

                                    ! Oxygen at time t-1
        OO2_before = OO2(im,j,nm,1)


                                    ! Guy Calculation : Justification: when DIC_dif_fromPOC > 0
                                    ! (POC is a C source), O2 gets consumed and since OetaO2... < 0
                                    ! we have to add DIC_dif_fromPOC/OetaC... * OetaO2...

        !NCP = DIC_dif_fromPOC/OetaC_POMoxid(im,j,nm)*OetaO2_POMoxid(im,j,nm)  &
        !      + DIC_dif_fromDOCtot/OetaC_DOMoxid_1D(j)*OetaO2_DOMoxid_1D(j)


                                    ! NCP =  Net Community Production (NCP), or in other words
                                    ! : the oxygen balance. Termes of production : phyto_growth.
                                    ! Termes of consumption : phyto_remin, zoo_remin, DOC_decay and
                                    ! DOCS_decay -> Converted by Oeta factors (With OetaO2 < 0).

        NCP = (phyto_growth*ABS(OetaO2_POMoxid(im,j,nm)/OetaC_POMoxid(im,j,nm)))  &
              + ((phyto_remin +  zoo_remin)/OetaC_POMoxid(im,j,nm)*OetaO2_POMoxid(im,j,nm))  &
              + (DIC_dif_fromDOCtot/OetaC_DOMoxid_1D(j)*OetaO2_DOMoxid_1D(j))

        OO2(im,j,nm,iair) = OO2(im,j,nm,iair) + NCP

#if ( OOISO == 1 )

#if ( RAYLEIGH == 0 )
        call compute_ISOO2_mbiodyn(NPP_O2,NCP,OO2(im,j,nm,:),OO2_netflux(:)) !watTemp optional
        OO2(im,j,nm,iair16:NISOO2) = OO2(im,j,nm,iair16:NISOO2)+OO2_netflux(iair16:NISOO2)
#else
        ! 1) Update OO2 after photosynthesis
        OO2_netflux_prod(:) = GPPO2_func(NPP_O2)
        OO2(im,j,nm,iair16:NISOO2) = OO2(im,j,nm,iair16:NISOO2)+OO2_netflux_prod(iair16:NISOO2)

        ! 2) Update OO2 after respiration according to Rayleigh Equation
        OO2_netflux_consum(:) = Ray_respO2(OO2_netflux_prod(iair),NCP,OO2_before,OO2(im,j,nm,:)) !watTemp optional
        OO2(im,j,nm,iair16:NISOO2) = OO2(im,j,nm,iair16:NISOO2)-OO2_netflux_consum(iair16:NISOO2)
#endif

! Residual --> Verification de la conservation: OO2(im,j,nm,iair) = OO2(im,j,nm,iair16) + OO2(im,j,nm,iair17) + OO2(im,j,nm,iair18)
        residual_O2 = OO2(im,j,nm,iair) - sum(OO2(im,j,nm,iair16:NISOO2))

! Saturation in the surface in Oxygene Isotopes :
!        OO2_sati(:) = OO2_saturation(OO2(im,1,nm,:),O2_sat)
!        OO2(im,1,nm,iair16:NISOO2) = OO2_sati(iair16:NISOO2)
#endif

                                    ! -- ecl : The saturation need to be applied after
                                    ! the isotopes calculation and on the surface :
                                    ! only if the OO2 > O2_sat.

!       OO2(im,1,nm,iair) = min(O2_sat, OO2(im,1,nm,iair))

                                    ! -- ecl : Variables pour sorties NEWGEN
       phyto_day_sum(im,j,nm) = phyto_day_sum(im,j,nm) + phyto_growth
       prodO2_sum(im,j,nm) = prodO2_sum(im,j,nm) + (2*NPP_O2)
       respO2_sum(im,j,nm) = respO2_sum(im,j,nm) + ( (2*NPP_O2) - NCP )
       NCP_sum(im,j,nm) = NCP_sum(im,j,nm) + NCP

        ! ALK dynamics
        ! ------------
#if ( ARAG == 0 )
        OALK(im,j,nm) =  &
!nb          OALK(im,j,nm) + ( 2._dblp*caco3_m  &
          OALK(im,j,nm) + ( 2._dblp*caco3_dif  &
                           - DIC_dif_fromPOC/OetaC_POMoxid(im,j,nm)*(OetaN_POMoxid(im,j,nm) + 1._dblp)  &
                           - DIC_dif_fromDOCtot/OetaC_DOMoxid_1D(j)*(OetaN_DOMoxid_1D(j) + 1._dblp))  &
                          * SCANU
#else
        OALK(im,j,nm) =  &
!nb          OALK(im,j,nm) + ( 2._dblp*(caco3_m + caco3_m_ar)  &
          OALK(im,j,nm) + ( 2._dblp*(caco3_dif + caco3_dif_ar)  &
                           - DIC_dif_fromPOC/OetaC_POMoxid(im,j,nm)*(OetaN_POMoxid(im,j,nm) + 1._dblp)  &
                           - DIC_dif_fromDOCtot/OetaC_DOMoxid_1D(j)*(OetaN_DOMoxid_1D(j) + 1._dblp))
                          * SCANU
#endif

!REFACTORING: On all the following lines, the conversion of the carbon fluxes is truly made from
!  mumol/kg to GtC, i.e., an explicit water density of 1.028 kg/L is taken into account.
!  Elsewhere SCALE_B is used, which implies a water density of 1 kg/L = 1000 kg/m3.

        cvt_mumolC_kg_to_GtC = DVOL(im,J,nm) * SCALE_M * 12._dblp * 1.028_dblp
        !GtC     =  mumol/kg * m3            * 1e18    * gC/mol * kg/L = GtC (1e18= 1e15 *1e3 L/m3)

        ! All the following are expressed in GtC
        prod_sum = prod_sum + phyto_growth * cvt_mumolC_kg_to_GtC
        graz_sum = graz_sum + grazing * cvt_mumolC_kg_to_GtC
        tpp_sum = tpp_sum + TPP_dif * cvt_mumolC_kg_to_GtC
        exu_sum = exu_sum + phyto_DOCexudat * cvt_mumolC_kg_to_GtC
        exc_sum = exc_sum + zoo_DOCexcret * cvt_mumolC_kg_to_GtC
        ort_sum = ort_sum + phyto_senesc * cvt_mumolC_kg_to_GtC
        pel_sum = pel_sum + (phyto_egest + zoo_egest) * cvt_mumolC_kg_to_GtC
        remin_sum = remin_sum + (phyto_remin + zoo_remin) * cvt_mumolC_kg_to_GtC
        barem_sum = barem_sum + doc_decay * cvt_mumolC_kg_to_GtC
        caco3_sum = caco3_sum + caco3_dif * cvt_mumolC_kg_to_GtC


      enddo ! timeloop


      endif ! on MGT == 1
    enddo ! on j, JPROD

    return

  end subroutine mbiodyn


!==================================
!REFACTORING: NOTHING CHANGED BELOW
!==================================


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
#if ( REMIN == 1 )
      SUBROUTINE compute_remin(i,j,n)
! Computes the remineralisation 3D values for POC (TPP) depending on temperature
! Dependance to temperature as in Crichton 202

      use loveclim_transfer_mod, only: TM, zz
      use marine_bio_mod, only: jprod, OO2

      integer i,j,n
      real deltaO2

      !kremin(i,j,n)=betaPOM*exp(-Ea/(Rgaz*(TM(i,j,n)+273.15))) ! From Crichton 2021
!      deltaO2=min(1.,max(0.,0.4*(O2min1-OO2(i,j,n))/(O2min2+OO2(i,j,n))))
      deltaO2=min(1.,max(0.,0.8*(O2min1-OO2(i,j,n))/(O2min2+OO2(i,j,n))))
      kremin(i,j,n)=betaPOM*exp(-Ea/(Rgaz*(TM(i,j,n)+273.15)))*(1-0.45*deltaO2)
      !kremin(i,j,n)=betaPOM*exp(-Ea/(Rgaz*(TM(i,j,n)+273.15)))*(1-0.9*deltaO2)
      !write(*,*) 'kremin', kremin(i,j,n)

      if (j.le.jprod+1) then
        SUE_3D(i,j,n)=1.
      else
        !SUE_3D(i,j,n)= dt/zz(j)*w_sink(j)*kremin(i,j,n)*TDAY !in /day (in /s-> in /day)
        SUE_3D(i,j,n)= SUE_3D(i,j-1,n)/(1+kremin(i,j-1,n)*zz(j)/w_sink(j)) !in /day
      endif

      END SUBROUTINE compute_remin

#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
#if (REMIN_CACO3 == 1 )
      SUBROUTINE compute_remin_ca(i,j,n,omega, kremin_carb, kremin_carbp)
! Computes the remineralisation 3D values for PIC (CacO3) depending on omega

      use loveclim_transfer_mod, only: zz
      use marine_bio_mod, only: jprod

      integer i,j,n
      real omega, kremin_carb, kremin_carbp

      if (omega .le. 1) then
         kremin_carb=k_diss*(1-omega)!**n_reac
      else
         kremin_carb=0
      endif
      !write(*,*) 'kremin', kremin_carb

      if (j.le.jprod+1) then
        SUE_ca_3D(i,j,n)=1.
      else
        !SUE_ca_3D(i,j,n)= dt/zz(j)*w_sink(j)*kremin_caco3*TDAY !in /day (in /s-> in /day)
        SUE_ca_3D(i,j,n)= SUE_ca_3D(i,j-1,n)/(1+kremin_carbp*zz(j)/w_sink(j)) !in /day kremin_carbp pour j-1
      endif

      return
      END SUBROUTINE compute_remin_ca
#endif
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
end module mbiota_mod
