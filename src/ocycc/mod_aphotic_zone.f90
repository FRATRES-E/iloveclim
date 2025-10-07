!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!   Copyright 2024 Didier M. Roche (a.k.a. dmr)

!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at

!       http://www.apache.org/licenses/LICENSE-2.0

!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:50 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:50 CET 2009


      MODULE MOD_APHOTIC_ZONE


       CONTAINS

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       SUBROUTINE MAPHOT(I,J,N, oc_bottom_cell, TM, OO2, OPO4, ONO3, OALK, ODIC, OC13, OC14, ODOC, ODOC13, ODOCS, ODOCS13, &
                OetaC_POMoxid, OetaN_POMoxid, OetaO2_POMoxid, TPP_ma, caco3_ma, DVOL_col)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       use global_constants_mod, only: dblp=>dp

       use declars_mod, ONLY: JT ! , LT, NOC_CBR

       use loveclim_transfer_mod, ONLY: SQRO2 ! , DVOL
       use mod_sync_time, ONLY:  TSTOC

       use mbiota_mod, ONLY: TPP_m, caco3_d13C, TPP_D13C, SUE_MCA, SCANU, RR, rdoc_ms, barem_sum   &
                           , d0_m, kd_m, SCALE_B, SCALE_M, fPOC_top, fPOC_1000, fPOC_bot, fCAL_top,&
                             fCAL_2000, fCAL_bot, b_sh_fr &
                           , caco3_mabot, O2min1, O2min2,OrgCFlxAttFactor, caco3_m, caco3_m_b_sh,  &
                           reminO2
#if ( RAYLEIGH == 1 )
       use iso_dioxygen_mod, ONLY: Ray_reminO2
#endif

#if ( OOISO == 1 )
       use iso_dioxygen_mod, ONLY: compute_ISOO2_maphot
#endif

#if ( REMIN == 1 )
       use mbiota_mod, ONLY: SUE_3D
#endif

#if ( REMIN_CACO3 == 1 )
       use mbiota_mod, ONLY: SUE_ca_3D, SUE_ar_3D
#endif

#if ( ARAG == 1 )
       use mbiota_mod, ONLY: caco3_d13C_ar, RR_ar, SUE_MAR, caco3_m_ar
#endif

       use marine_bio_mod, ONLY: pr, OetaC_DOMoxid_1D, OetaN_DOMoxid_1D, OetaO2_DOMoxid_1D, sigma_m

       use marine_bio_mod, only: JPROD

#if ( OXNITREUX == 1 )
       use marine_bio_mod, ONLY: ON2O
#endif

       use para0_mod, ONLY: NISOO2
       use iso_dioxygen_mod, ONLY: iair, Ototal_min

#if ( OOISO == 1 )
       use iso_dioxygen_mod, ONLY: iair16, iair17, iair18
       use iso_dioxygen_mod, ONLY: r18smow, r17smow
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      IMPLICIT NONE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


      INTEGER, INTENT(in) :: I, J, N

      REAL(kind=dblp), dimension(JT, NISOO2), INTENT(inout) :: OO2

      LOGICAL,         DIMENSION(JT), INTENT(in)    :: oc_bottom_cell
      REAL(kind=dblp), dimension(JT), INTENT(in)    :: TM
      REAL(kind=dblp), dimension(JT), INTENT(in)    :: DVOL_col

      REAL(kind=dblp), dimension(JT), INTENT(inout) :: OPO4
      REAL(kind=dblp), dimension(JT), INTENT(inout) :: ONO3
      REAL(kind=dblp), dimension(JT), INTENT(inout) :: OALK
      REAL(kind=dblp), dimension(JT), INTENT(inout) :: ODIC
      REAL(kind=dblp), dimension(JT), INTENT(inout) :: OC13
      REAL(kind=dblp), dimension(JT), INTENT(inout) :: OC14
      REAL(kind=dblp), dimension(JT), INTENT(inout) :: ODOC
      REAL(kind=dblp), dimension(JT), INTENT(inout) :: ODOC13
      REAL(kind=dblp), dimension(JT), INTENT(inout) :: ODOCS
      REAL(kind=dblp), dimension(JT), INTENT(inout) :: ODOCS13

      REAL(pr), dimension(JT), INTENT(in)           :: OetaC_POMoxid
      REAL(pr), dimension(JT), INTENT(in)           :: OetaN_POMoxid
      REAL(pr), dimension(JT), INTENT(in)           :: OetaO2_POMoxid

      REAL(kind=dblp), dimension(JT), INTENT(inout) :: TPP_ma
      REAL(kind=dblp), dimension(JT), INTENT(inout) :: caco3_ma


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      REAL(kind=dblp) :: OPOC_dif    ! POC loss from flux though the water column (OPOC_dif > 0: POC loss!!!)
      REAL(kind=dblp) :: ODOC_dif    ! [DOC] variation (ODOC_dif < 0: DOC loss), added to ODOC(i,j,n)
      REAL(kind=dblp) :: ODOCS_dif   ! [DOCS] variation (ODOCS_dif < 0: DOC loss), added to ODOCS(i,j,n)
      REAL(kind=dblp) :: OPO4_dif    ! [PO4] variation (OPO4_dif > 0: PO4 gain), added to OPO4(i,j,n)
      REAL(kind=dblp) :: ONO3_dif    ! [NO3] variation (ONO3_dif > 0: NO3 gain), added to ONO3(i,j,n)
      REAL(kind=dblp) :: OO2_dif     ! [O2] variation (OO2_dif < 0: O2 loss), added to OO2(i,j,n)
      REAL(kind=dblp) :: OPOC13_dif
      REAL(kind=dblp) :: ODOCS13_dif
      REAL(kind=dblp) :: ODOC13_dif
      REAL(kind=dblp) :: caco3_dif    ! caco3 loss from flux though the water column (caco3_dif > 0: caco3 loss!!!)
#if ( OXNITREUX == 1 )
      REAL(kind=dblp) :: ON2O_dif
#endif

      REAL(kind=dblp) :: caco3_13, TPP_compute, caco3_compute ! calred unused a priori
#if ( ARAG == 1 )
      REAL(kind=dblp) :: caco3_13_ar, caco3_compute_ar
#endif
      REAL(kind=dblp) :: RDOC_M
      REAL(kind=dblp) :: VDIC        ! [DIC] variation due to oxidation of DOC, DOCS and POC (VDIC > 0: DIC gain),
                                     !   added to ODIC(i,j,n) after multiplication by SCANU
      REAL(kind=dblp) :: lo2_m       ! Replaces OPOC_dif from line 155 (notice that Six & Maier-Reimer (1996)
                                     !   denote the remineralization of POC by $l(O_2)$); OPOC_dif not changed afterwards
      REAL(kind=dblp) :: deltaO2

      REAL(kind=dblp)                    :: OO2_before_remin
#if ( OOISO == 1 )
      REAL :: residual_O2
      REAL(kind=dblp), dimension(NISOO2) :: OO2_flux_isoremin
#endif
!nb #define bottom_remin 0 ! to be removed
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!--- Init of local variables
        OPOC_dif=0
        ODOC_dif=0
        ODOCS_dif=0
        caco3_dif=0
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8-|
! dmr
!              [NOTE] Guy & Didier :
!                                   [TPP_m] = Tmols
!                                   [TPP_m*OrgCFlxAttFactor] = Tmols
!                                   [TPP_ma] = Tmols.m-2
!                     2021-10-14 : changed TPP_ma so it is the flux through the
!                                  bottom interface of the cell
!                                  ... unit is per timestep of maphot = day-1
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8-|
#if ( REMIN == 1 )
       OrgCFlxAttFactor(i, j, n) = SUE_3D(i, j, n)
       OrgCFlxAttFactor(i, j+1, n) = SUE_3D(i, j+1, n)
#endif

        TPP_compute = OrgCFlxAttFactor(i, j+1, n)*TPP_m(i,n)                      ! [FIXEDINPUT] : OrgCFlxAttFactor
        TPP_ma(j) = TPP_ma(j) + TPP_compute/SQRO2(i,n)                            ! [INOUTPUT]   : tpp_ma(i,j,n)
                                                                                  ! [FIXEDINPUT] : SQRO2(i,n)

!RFC: The following looks like a diagnostic. Correct? Comment would be helpful.
        if (j.eq.(JPROD+1)) then
          fPOC_top = fPOC_top + TPP_compute
        elseif ((j.eq.13).or.(j.eq.14)) then                 ! index 8 on CLIO is 850 m, index 7 is 1200 meters, so I take the mean of the two (on Carbon => 21-8 = 13 and 21-7 = 14)
          fPOC_1000 = fPOC_1000 + TPP_compute/2.
!nb test        elseif (oc_bottom_cell(j)) then
        elseif (oc_bottom_cell(j).and.(j.gt.13)) then   ! j=13 limite de 1000m
          fPOC_bot  = fPOC_bot + TPP_compute  ! here fPOC is in Tmols
        endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! POC dynamics
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!if not medusa, no sediments, all is remineralised/dissolved in last layer
!nb #if ( bottom_remin == 1 )
#if ( MEDUSA == 0 )
          if (oc_bottom_cell(j)) then
!REFACTORING NOTE (IMPORTANT):
! since the OrgCFlxAttFactor(i,j+1,n) when (i,j+1,n) is a bottom cell
! is already set to zero, or not if MEDUSA is coupled,
! there is not need to process the bottom cell separately.
! When MEDUSA is coupled, it is actually wrong to process
! it separately, as the OrgCFlxAttFactor(i,j+1,n) is not zero in this case

! dmr [OPOC_dif] == mumol/kg
!PROBLEM: for TPP_m in Tmol, [OPOC_dif] is actually in mumol/L: the conversion factor SCALE_B
!         used in the calculation implicitly includes a water density of 1 kg/L or 1000 kg/m3,
!         This is in contradiction with a density of 1.028 kg/L used elsewhere
!         (e.g., mbiota_mod.f90:570ff)
            OPOC_dif=TPP_m(i,n)*OrgCFlxAttFactor(i, j, n)/(DVOL_col(J)*SCALE_B)      ! [FIXEDINPUT] : OrgCFlxAttFactor, DVOL, SCALE_B
            OPOC13_dif=TPP_D13C(i,n)*OrgCFlxAttFactor(i, j, n)/(DVOL_col(J)*SCALE_B) ! [INPUT]      : TPP_D13C(i,n)
            !write*,* 'maphot, TPP_D13C' TPP_D13C(i,n), OrgCFlxAttFactor(i, j, n)

          else ! .not. bottom_cell
#endif

            OPOC_dif=TPP_m(i,n)*(OrgCFlxAttFactor(i, J, n)-OrgCFlxAttFactor(i, J+1, n))/(DVOL_col(j)*SCALE_B)
            OPOC13_dif=TPP_D13C(i,n)*(OrgCFlxAttFactor(i, J, n)-OrgCFlxAttFactor(i, J+1, n))/(DVOL_col(j)*SCALE_B)

!nb #if ( bottom_remin == 1 )
#if ( MEDUSA == 0 )
          endif ! .is. bottom_cell
#endif
! missing oxygen concentration -ox_m
! changes - POC decay is constant (limited by oxygen concentration)

         lo2_m=OPOC_dif

!nb only if under euphotic zone
         if (j.ge.JPROD+1) then

!dmr --- Par construction
!dmr --- 0.1*d0_m*TSTOC < rdoc_m < d0_m*TSTOC
!dmr --- rdoc : remineralisation / dissolution par rapport au PO4

!        rdoc_m=d0_m*max(0.0,OPO4(i,j,n))/(max(0.0,OPO4(i,j,n))+kd_m)*TSTOC      ! [INPUT]   : OPO4(i,j,n)

        deltaO2=min(1.,max(0.,0.4*(O2min1-OO2(j,iair))/(O2min2+OO2(j,iair))))
        rdoc_m=d0_m*max(0.0,OPO4(j))/(max(0.0,OPO4(j))+kd_m)*TSTOC*(1-deltaO2)

        if (rdoc_m.le.0.1*d0_m*TSTOC) then
          rdoc_m=0.1*d0_m*TSTOC
        endif

        barem_sum=barem_sum+rdoc_m*ODOC(j)*DVOL_col(j)*SCALE_M*12           ! [INPUT]   : ODOC(j)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! fast DOC dynamics
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        ODOC_dif=-rdoc_m*ODOC(j)
        ODOC13_dif=-rdoc_m*ODOC13(j)

        ODOC(j)=ODOC(j)+ODOC_dif                                       ! [INOUTPUT]   : ODOC(j)
        ODOC13(j)=ODOC13(j)+ODOC13_dif                                 ! [INOUTPUT]   : ODOC13(j)

!~ !dmr --- gestion de l "underflow"
!~         if( ODOC(j).lt.0)  then
!~ !nb        ODOC(j)=0
!~           err_ODOC=err_ODOC-ODOC(j)/OetaC_DOMoxid_1D(j)*DVOL(j)*SCALE_M
!~         endif

! slow DOC dynamics

        ODOCS_dif=-rdoc_ms*TSTOC*ODOCS(j)                                  ! [FIXEDINPUT] : TSTOC, rdoc_ms
        ODOCS13_dif=-rdoc_ms*TSTOC*ODOCS13(j)

        ODOCS(j)=ODOCS(j)+ODOCS_dif                                    ! [INOUTPUT]   : ODOCS(j)
        ODOCS13(j)=ODOCS13(j)+ODOCS13_dif                              ! [INOUTPUT]   : ODOCS13(j)

!~ !dmr --- gestion de l "underflow"
!~         if( ODOCS(j).lt.0)  then
!~ !nb        ODOCS(j)=0
!~          err_ODOCS=err_ODOCS-ODOCS(j)/OetaC_DOMoxid_1D(j)*DVOL_col(j)*SCALE_M
!~         endif

!nb
      endif ! j.ge.JT+1

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! CaCO3 - caco3_m -> calcite, caco3_m_ar -> aragonite
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if (REMIN_CACO3 == 1 )
         SUE_MCA(j)=SUE_ca_3D(i,j,n)
         SUE_MCA(j+1)=SUE_ca_3D(i,j+1,n)
#if (ARAG == 1 )
         SUE_MAR(j)=SUE_ar_3D(i,j,n)
         SUE_MAR(j+1)=SUE_ar_3D(i,j+1,n)
#endif
#endif

!nb #if ( bottom_remin == 1 )
#if ( MEDUSA == 0 )
         if (oc_bottom_cell(j)) then

!nb remove           caco3_m=(RR*calred*SUE_MCA(j)*TPP_m(i,n)/(1-sigma_m)    &!/(DVOL_col(j)*SCALE_B) ! [FIXEDINPUT] : SUE_MCA
           caco3_dif=(SUE_MCA(j)*caco3_m(i,n)    &!/(DVOL(j)*SCALE_B) ! [FIXEDINPUT] : SUE_MCA
!                   + b_sh_fr*RR*calred*TPP_m/(1-sigma_m))     &!nb test
                   + caco3_m_b_sh(i,n))     &!nb test
                   /(DVOL_col(j)*SCALE_B) ! [FIXEDINPUT] : SUE_MCA

           caco3_13=SUE_MCA(j)*caco3_d13C(i,n)/(DVOL_col(j)*SCALE_B)
           caco3_d13C(i,n)=0

#if ( ARAG == 1 )
           caco3_dif_ar=SUE_MAR(j)*caco3_m_ar/(DVOL_col(j)*SCALE_B) ! [FIXEDINPUT] : SUE_MCA
           caco3_13_ar=SUE_MAR(j)*caco3_d13C_ar/(DVOL_col(j)*SCALE_B)
           caco3_d13C_ar=0
#endif

         else ! .not. bottom_cell
#endif

!nb remove            caco3_m=RR*calred*(SUE_MCA(J)-SUE_MCA(J+1))*TPP_m(i,n)/(1-sigma_m)/(DVOL_col(j)*SCALE_B)
            caco3_dif=(SUE_MCA(J)-SUE_MCA(J+1))*caco3_m(i,n)/(DVOL_col(j)*SCALE_B)
            caco3_13=(SUE_MCA(J)-SUE_MCA(J+1))*caco3_d13C(i,n)/(DVOL_col(j)*SCALE_B)

#if ( ARAG == 1 )
            caco3_m_ar=(SUE_MAR(J)-SUE_MCA(J+1))*caco3_m_ar/(DVOL_col(j)*SCALE_B)
            caco3_13_ar=(SUE_MAR(J)-SUE_MAR(J+1))*caco3_d13C_ar/(DVOL_col(j)*SCALE_B)
#endif

! #if ( bottom_remin == 1 )
#if ( MEDUSA == 0 )
         endif ! .is. bottom_cell
#endif


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8-|
! dmr   (Same procedure as for TPP)
!              [NOTE] Guy & Didier :
!                                   [TPP_m] = Tmols
!                                   [TPP_m*OrgCFlxAttFactor] = Tmols
!                                   [CaCO3_ma] = Tmols.m-2
!                     2021-10-23 : changed caco3_ma so it is the flux through the
!                                  bottom interface of the cell
!                                  ... unit is per timestep of maphot = day-1
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8-|

!nb remove        caco3_compute = RR*calred*SUE_MCA(j+1)*TPP_m(i,n)/(1-sigma_m)
        caco3_compute = SUE_MCA(j+1)*caco3_m(i,n)
#if ( ARAG == 1 )
!nb remove        caco3_compute_ar = RR_ar*calred*SUE_MAR(j+1)*TPP_m(i,n)/(1-sigma_m)
        caco3_compute_ar = SUE_MAR(j+1)*caco3_m_ar(i,n)
#endif
        if ( oc_bottom_cell(j) ) then ! Need to add the big_shell_fraction that is sedimented directly
!nb remove          caco3_compute = caco3_compute + b_sh_fr*RR*calred*TPP_m/(1-sigma_m)
          caco3_compute = caco3_compute + caco3_m_b_sh(i,n)
        endif
#if ( ARAG == 0 )
        caco3_ma(j) = caco3_ma(j) + caco3_compute/SQRO2(i,n)                 ! [INOUTPUT]   : caco3_ma(j)
#else
!a voir : une variable caco3_ma_ar?
        caco3_ma(j) = caco3_ma(j) + (caco3_compute+caco3_compute_ar)/SQRO2(i,n)                 ! [INOUTPUT]   : caco3_ma(j)
#endif
                                                                                 ! [FIXEDINPUT] : SQRO2(i,n)

        if (j.eq.(JPROD+1)) then
          fCAL_top = fCAL_top + caco3_compute
        elseif ((j.eq.16).or.(j.eq.15)) then                   ! index 6 on CLIO is 1700 m, index 5 is 2300 meters, so I take the mean of the two (on Carbon => 21-5 = 16 and 21-6 = 15)
          fCAL_2000 = fCAL_2000 + caco3_compute/2.
!nb        elseif (oc_bottom_cell(i,j,n)) then
        elseif (oc_bottom_cell(j).and.(j.gt.13)) then
          fCAL_bot  = fCAL_bot + caco3_compute  ! here fCAL is in Tmols
          caco3_mabot(i,n) = caco3_mabot(i,n) + caco3_compute/SQRO2(i,n) ! in Tmols.m-2.day-1
        endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! DIC dynamics
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!dmr --- ODOC  = Fast Organic Carbon Pool
!dmr --- ODOCS = Slow Organic Carbon Pool
!gm  --- lo2_m = OPOC_dif (name switch at some place)

!PROBLEM: on all the following lines, the conversion of the carbon fluxes is not made from
!         mumol/kg to GtC, i.e., water density is taken into account.
!         elsewhere, where SCANU is used, water density is implicitly set to 1 kg/L = 1000 kg/m3

        vdic=-(ODOC_dif+ODOCS_dif)+lo2_m
#if ( ARAG == 0 )
        OC13(j)=OC13(j)+(OPOC13_dif+caco3_13-ODOC13_dif-ODOCS13_dif)*SCANU
        !write(*,*) 'maphot OC13', OC13(j), OPOC13_dif, caco3_13, -ODOC13_dif, -ODOCS13_dif
        ODIC(j)=ODIC(j)+(vdic+caco3_dif)*SCANU
        OC14(j)=OC14(j)+(vdic+caco3_dif)*OC14(j)/ODIC(j)*SCANU
#else
        OC13(j)=OC13(j)+(OPOC13_dif+(caco3_13+caco3_13_ar)-ODOC13_dif-ODOCS13_dif)*SCANU
        ODIC(j)=ODIC(j)+(vdic+caco3_dif+caco3_dif_ar)*SCANU !to be modified -> nb TODO ???
        OC14(j)=OC14(j)+(vdic+caco3_dif+caco3_dif_ar)*OC14(j)/ODIC(j)*SCANU
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! PO4 dynamics
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

         if (oc_bottom_cell(j)) then
           OPO4_dif = -(ODOC_dif+ODOCS_dif)/OetaC_DOMoxid_1D(j) + lo2_m/OetaC_POMoxid(j)
           OO2_dif = (-(ODOC_dif+ODOCS_dif)/OetaC_DOMoxid_1D(j))*OetaO2_DOMoxid_1D(j) + &       ! note that OetaO2_DOMoxid_1D < 0
                     lo2_m/OetaC_POMoxid(j)*OetaO2_POMoxid(j)                           ! and OetaO2_POMoxid < 0

#if ( OXNITREUX == 1 )
           ON2O_dif=(-(ODOC_dif+ODOCS_dif)/OetaC_DOMoxid_1D(j))*OetaO2_DOMoxid_1D(j)
#endif
           TPP_m(i,n)=0
           TPP_D13C(i,n)=0
           caco3_d13C(i,n)=0 !nb sediments
#if (ARAG == 1 )
           caco3_d13C_ar(i,n)=0 !nb sediments
#endif
         else ! .not. bottom_cell

           OPO4_dif = -(ODOC_dif+ODOCS_dif)/OetaC_DOMoxid_1D(j) + lo2_m/OetaC_POMoxid(j)
           OO2_dif = (-(ODOC_dif+ODOCS_dif)/OetaC_DOMoxid_1D(j)*OetaO2_DOMoxid_1D(j) + &
                     lo2_m/OetaC_POMoxid(j)*OetaO2_POMoxid(j))
#if ( OXNITREUX == 1 )
           ON2O_dif=(-(ODOC_dif+ODOCS_dif)/OetaC_DOMoxid_1D(j))*OetaO2_DOMoxid_1D(j)
#endif

         endif ! .is. bottom_cell

         OPO4(j)=OPO4(j)+OPO4_dif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! NO3 dynamics
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

         ONO3_dif = -(ODOC_dif+ODOCS_dif)/OetaC_DOMoxid_1D(j)*OetaN_DOMoxid_1D(j) + &
                    lo2_m/OetaC_POMoxid(j)*OetaN_POMoxid(j)

         ONO3(j) = ONO3(j) + ONO3_dif


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! Three isotopes of dissolved O2 (activated if OOISO == 1)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! O2 dynamics
!jwy --- total O2 consumption = OO2_dif
!jwy --- 18O16O/16O16O ratio of OO2_dif is alpha(resp) * current R * OO2_dif

!ecl ---  [TO BE CHECKED] Normalement: pas de min mais une remin positive a-t-elle un sens ???
!ecl --- OO2_diff est en majorité négatif tout le temps
         OO2_before_remin = OO2(j,iair)
         OO2_dif = min(0.0, OO2_dif)

!ecl -- Remin limitation to manage ZMO and upwelling zone:
         if ( (OO2(j,iair) + OO2_dif) < Ototal_min ) then
            OO2_dif = 0.0
         endif

         OO2(j,iair)=OO2(j,iair)+OO2_dif

#if ( OOISO == 1 )

#if ( RAYLEIGH == 0 )
         call compute_ISOO2_maphot(OO2_dif,OO2(j,:),OO2_flux_isoremin(:)) !TM(j) optional
         OO2(j,iair16:NISOO2) = max(0.0, OO2(j,iair16:NISOO2)+OO2_flux_isoremin(iair16:NISOO2))
#else
         OO2_flux_isoremin(:) = Ray_reminO2(OO2_dif,OO2_before_remin, OO2(j,:)) !Tm(j) optional argument
         OO2(j,iair16:NISOO2) = max(0.0, OO2(j,iair16:NISOO2) + OO2_flux_isoremin(iair16:NISOO2))
#endif

! Residual --> Verification de la conservation des isotopes:
         residual_O2 = OO2(j,iair) - sum(OO2(j,iair16:NISOO2))
#endif

! Newgen output
        reminO2(i,j,n) = OO2_dif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! N2O dynamics
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!vm D apres la parametrisation utilisee dans le modele PISCES

#if ( OXNITREUX == 1 )

         ON2O(j)=ON2O(j)+ON2O_dif*(1.e-4+4.e-3*exp(-0.1*OO2(j,1)*1.e6-1.))
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ALK dynamics
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!REFACTORING DONE: effect of PO4 variations on Alkalinity changes added;
!                  effect of NO3 variations now explicit (not derived from PO4)
#if ( ARAG == 0 )
         OALK(j)=OALK(j)+SCANU*(2.*caco3_dif - OPO4_dif - ONO3_dif)

#else
         OALK(j)=OALK(j)+SCANU*(2.*(caco3_dif+caco3_dif_ar) - OPO4_dif - ONO3_dif) !nb to be modified
#endif
        return

        END SUBROUTINE MAPHOT

      END MODULE MOD_APHOTIC_ZONE
