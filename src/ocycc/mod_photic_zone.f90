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


      MODULE MOD_PHOTIC_ZONE


       CONTAINS


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      SUBROUTINE MBIOTA(tm, frice, mgt, dvol, oc_bottom_cell, phyto_m, zoo_m, phyto_m13, zoo_m13, phyto_prod,          &
                        oo2, opo4, ono3, oalk, odic, oc13, oc14, odoc, odoc13, odocs, odocs13, oxCO2,                  &
                        OetaC_POMoxid, OetaN_POMoxid, OetaO2_POMoxid,  OetaC_DOMoxid_1D, OetaN_DOMoxid_1D,             &
                        tpp_ma, caco3_ma)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!*******************MARINE BIOTA ************************************
!
!  Code written by: V.Brovkin
!  First version: 03.07.98
!  Last modification: 16.07.99
!  Last modification: 08.11.07, D.M. Roche & N. Bouttes ; 11.07.19, idem ; 02.12.21, dmr, 2024-11-08
!
!  Purpose: simulation of nutrients/phytoplankton/zooplankton/DOC/POC
!           dynamics in the ocean
!
!  Description:  K.Six and E.Maier-Reimer, GBC, 10: 559-583, 1996
!********************************************************************
       use global_constants_mod, only: dblp=>dp, ip

       use declars_mod, only: lt, jt, noc_cbr

       use loveclim_transfer_mod, only: ovol, SABST_O, TM_surface, ZZ, ZX
       use mod_sync_time, only: tstoc, kendy, NYR
       use para0_mod, ONLY: NISOO2
#if ( CORAL == 1 )
       use mod_sync_time, only: TYER, TDAY, KWEEK, NMONTH, NDAY
!#elif ( REMIN == 1 || REMIN_CACO3 == 1 )
!     >                          , zz, TYER, TDAY
#endif
!~        use mbiota_mod, only: time_stm, tpp_mas, caco3_mas, total_car
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       use mbiota_mod, only: time_stm, total_car, total_phos, total_tpp, scale_m, scale_b, total_caco3, scanu, prod_sum     &
                     , graz_sum, remin_sum, exc_sum, exu_sum, barem_sum, tpp_sum, ort_sum, pel_sum, mbiodyn, caco3_sum      &
                     , fPOC_top, fPOC_1000, fPOC_bot, fCAL_top, fCAL_2000, fCAL_bot, TPP_m, TPP_D13C

!pb #if ( PATH == 1 )
!pb       use mbiota_mod, only: particlesFields
!pb #endif
#if ( REMIN == 1 )
!       use mbiota_mod, only: SUE_3D, kremin, Ea, Rgaz, betaPOM
       use mbiota_mod, only: compute_remin
#endif

!#if ( REMIN == 1 || REMIN_CACO3 == 1 )
!       use mbiota_mod, only: dt, w_sink
!#endif
#if ( REMIN_CACO3 == 1 )
       use omega_mod, only : calc_omega_ca, omega_calc3D
       use mbiota_mod, only: kremin_ca, kremin_ar, compute_remin_ca
#endif
       use marine_bio_mod, only: jprod, pr, oalk_ini, opo4_ini, ono3_ini, oc13bio


       use C_res_mod, only: c_ocean_fich, fPOC_fich, fCAL_fich

#if ( CORAL == 1 )
       use coral_mod, only : total_area_coral_an, total_prod_coral_an,
     >              total_mass_coral_an, coral_area_out, coral_area,
     >              coral_prod_out, coral_prod, coral_co2, c_riv,
     >              a_riv, c_car_a, out_coral_global,
     >              mbiota_corals, calc_monthly_temp,
     >              calc_DHW, tau_bleach, tau_bleach_out,
     >              DHW_nb, DHW_out, window_MMM,
     >              temp_too_low_all, calc_temp_variability,
     >              k_mbiota_rand, calc_nino3, indice_hs, nb_hs,
     >              coral_mass_out, coral_cum_mass
       use omega_mod, only : calc_omega_ar
#endif

#if ( ARAG == 1 )
       use omega_mod, only : calc_omega_ar, omega_arag3D
#endif
       use O2SAT_mod, only : O2_sat_thistime

       use mod_aphotic_zone, only: maphot

      IMPLICIT NONE


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=dblp), dimension(LT,JT,NOC_CBR), intent(in)   :: TM        ! from loveclim_transfer_mod
      real(kind=dblp), dimension(LT,NOC_CBR), intent(in)      :: FRICE        ! from loveclim_transfer_mod
      integer(kind=ip), dimension(LT,JT,NOC_CBR), intent(in)  :: MGT      ! from loveclim_transfer_mod
      real(kind=dblp),  dimension(LT,JT,NOC_CBR), intent(in)  :: DVOL     ! from loveclim_transfer_mod
      logical, dimension(LT,JT,NOC_CBR), intent(in)           :: oc_bottom_cell    ! from mbiota_mod
      real(kind=dblp), dimension(LT,JT,NOC_CBR), intent(in)   :: PHYTO_M   ! from mbiota_mod
      real(kind=dblp), dimension(LT,JT,NOC_CBR), intent(in)   :: ZOO_M     ! from mbiota_mod
      real(kind=dblp), dimension(LT,JT,NOC_CBR), intent(in)   :: PHYTO_M13 ! from mbiota_mod
      real(kind=dblp), dimension(LT,JT,NOC_CBR), intent(in)   :: ZOO_M13   ! from mbiota_mod
      double precision, dimension(LT,JT,NOC_CBR), intent(out) ::  phyto_prod ! from marine_cio_mod
      real, dimension(LT,JT,NOC_CBR), intent(inout)           :: ono3              ! from marine_bio_mod
      real, dimension(LT,JT,NOC_CBR), intent(inout)           :: opo4              ! from marine_bio_mod
      real, dimension(LT,JT,NOC_CBR), intent(inout)           :: oalk              ! from marine_bio_mod
      real, dimension(LT,JT,NOC_CBR), intent(inout)           :: oc13              ! from marine_bio_mod
      real, dimension(LT,JT,NOC_CBR), intent(inout)           :: oc14              ! from marine_bio_mod
      real, dimension(LT,JT,NOC_CBR), intent(inout)           :: odoc              ! from marine_bio_mod
      real, dimension(LT,JT,NOC_CBR), intent(inout)           :: odoc13            ! from marine_bio_mod
      real, dimension(LT,JT,NOC_CBR), intent(inout)           :: odocs             ! from marine_bio_mod
      real, dimension(LT,JT,NOC_CBR), intent(inout)           :: odocs13           ! from marine_bio_mod
      real, dimension(LT,NOC_CBR), intent(in)                 :: oxCO2                   ! from marine_bio_mod
      real(pr), dimension(LT,JT,NOC_CBR), intent(in)          :: OetaC_POMoxid    ! from marine_bio_mod
      real(pr), dimension(LT,JT,NOC_CBR), intent(in)          :: OetaN_POMoxid    ! from marine_bio_mod
      real(pr), dimension(LT,JT,NOC_CBR), intent(in)          :: OetaO2_POMoxid   ! from marine_bio_mod
      real(pr), dimension(JT), intent(in)                     :: OetaC_DOMoxid_1D            ! from marine_bio_mod
      real(pr), dimension(JT), intent(in)                     :: OetaN_DOMoxid_1D            ! from marine_bio_mod
      real(kind=dblp), dimension(LT,JT,NOC_CBR), intent(out)  :: TPP_ma   ! from mbiota_mod
      real(kind=dblp), dimension(LT,JT,NOC_CBR), intent(out)  :: caco3_ma ! from mbiota_mod

      real(kind=8), dimension(LT,JT,NOC_CBR,NISOO2), intent(inout), target :: oo2 ! from marine_bio_mod
      real(kind=8), dimension(LT,JT,NOC_CBR), intent(inout), target        :: odic ! from marine_bio_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      integer(kind=ip):: N, I, J

      real(kind=dblp) :: TOTAL_PHYTO, TOTAL_ZOO, TOTAL_DOC, TOTAL_SDOC, TOTAL_POC, TOTAL_ALK, TOTAL_NIT, TOTAL_D13     &
                       , TEMP_DIFA, TEMP_DIFP, TEMP_DIFN

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      real(kind=dblp) :: total_phos_1,total_phos_2, tot_err, tot_err2, vtmp_POC, vtmp_DOC,error_phos, vol_col,         &
                         total_alk_1, total_alk_2, error_alk, total_nit_1, total_nit_2, error_nit, total_car_1,        &
                         total_car_2, error_car, total_car13_1, total_car13_2, error_car13

      real(kind=dblp) :: cvt_mumolC_kg_to_GtC

! input data:
!                ocean surface temperature - TM(i,1,n)
!               total incoming radiation to the surface - SABST(i,1,n (W/m2))
!                oxygen concentration - Ox_m (not used)
!               simulation of advection and diffusion of
!               nutrients, DOC, POC - in file ocn.f
!                area free from ice - (1.-FRICE(i,n))
!
! output data:         OPO4(i,j,n)- phosphorus concentration in the ocean
!                PHYTO_M(i,j,n)-phytoplankton concentr. in the ocean
!                ZOO_M(i,j,n)-zooplankton concentration in the ocean
!                ODOC(i,j,n)-DOC concentration in the ocean
!                OPOC(i,j,n)-POC concentration in the ocean-> this is
!                zero and is removed !nb
!nb ODIC in  [mol/kg]
!nb PHYTO_M (and others) in [mumol/kg]
!
! COMMON block described in mbiota.inc
!
! Details:
!                time step of MBIOTA is 1/10 day;
!                POC distribution is calculated every time step.
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( CORAL == 1 )
! for temperature variability
      k_mbiota_rand = 0
! for bleaching
      if (indice_hs.lt.nb_hs) then ! nb_hs=84
        indice_hs=(NMONTH-1)*30+NDAY
        !write(*,*) "indice_hs", indice_hs
      endif
#endif

!       call for MBIOTA model
!...    SPATIAL AND TEMPORAL LOOP

!nb TSTOC=86400s, TIME_STM=T_DAY/10, T_DAY= 86400s
!nb -> itimemax=10+0.001
!~       itimemax=TSTOC/TIME_STM+0.001 !TSTOC = 86400 s / TIME_STM=T_DAY/10

!nb currently mbiota called once per day by emic (through ocyc)
!nb so mbiodyn called 10 times /day
!nb maphot called once per day
!#if ( CORAL ==1 )
!       total_area_coral_an=0
!#endif


#if ( REMIN == 1 )
!recompute remineralisation profile dependant on temperature
!before call to maphot
      do n=1,NOC_CBR
       do j=1,JT
        do i=1,LT
          call compute_remin(i,j,n)
        enddo
       enddo
      enddo
#endif

#if ( REMIN_CACO3 == 1 )
!recompute remineralisation profile for CaCO3 dependant on omega
      do n=1,NOC_CBR
       do j=1,JT
        do i=1,LT
           call calc_omega_ca(i,j,n)
           call compute_remin_ca(i,j,n,omega_calc3D(i,j,n), kremin_ca(i,j,n), kremin_ca(i,j-1,n))
#if ( ARAG == 1 )
           call calc_omega_ar(i,j,n)
           call compute_remin_ca(i,j,n,omega_arag3D(i,j,n), kremin_ar(i,j,n), kremin_ar(i,j-1,n))
#endif
        enddo
       enddo
      enddo
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       do n=1,NOC_CBR
        do i=1,LT

        if (MGT(i,1,n).eq.1) then

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
                call MBIODYN(i,n, SABST_O(i,n),TM(i,:,n),TM_surface(i,n), OXCO2, oc13bio, OC13, ODIC,OPO4,   &
                             ODOC, ODOCS, ODOC13, ODOCS13,PHYTO_PROD, OC14, ONO3, OO2, OALK, DVOL ,ZZ, ZX,   &
                             FRICE, MGT, O2_sat_thistime)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!dmr --- [WARNING] CORALS ARE BROKEN -> need to internalize the computing in MBIODYN in a next stage
#if ( CORAL == 1 )
           !in the euphotic zone
           do j=1,JPROD
             if (MGT(i,j,n).eq.1) then
              !nb for bleaching
              ! add variability to temperature for bleaching
              call calc_temp_variability(i,j,n)

              !replace end of matrix with monthly mean
              !temperature and compute max of montly climatology
              if (NYR.le.window_MMM) then ! if fixed reference for
                                          !bleaching (si commente
                                          !rolling mean)
              call calc_monthly_temp(i,j,n)
              endif
             !nb once the reference has been computed we can start using it for
             !bleaching
              if (NYR.gt.window_MMM) then
                call calc_DHW(i,j,n) ! computes degree heating weeks DHW
              endif
!              if (NYR.eq.1) then !for test with fixed omega to have temperature effect
                call calc_omega_ar(i,j,n)
!              endif
                call mbiota_corals(i,j,n)
             endif !MGT
           enddo !j
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!dmr&nb [NOTA] Partie de MAPHOT à faire sur la partie photique aussi
!dmr&nb          do j=JPROD+1,JT
          do j=1,JT
            if (MGT(i,j,n).eq.1) then
#if ( CORAL == 1 || ARAG == 1 )
              call calc_omega_ar(i,j,n)
#endif
              call MAPHOT(i,j,n, oc_bottom_cell(i,:,n), TM(i,:,n), OO2(i,:,n,1:NISOO2), OPO4(i,:,n), ONO3(i,:,n), OALK(i,:,n), &
                          ODIC(i,:,n), OC13(i,:,n), OC14(i,:,n), ODOC(i,:,n), ODOC13(i,:,n), ODOCS(i,:,n), ODOCS13(i,:,n),     &
                          OetaC_POMoxid(i,:,n), OetaN_POMoxid(i,:,n), OetaO2_POMoxid(i,:,n), TPP_ma(i,:,n), caco3_ma(i,:,n),   &
                          DVOL(i,:,n))
            endif
          enddo ! j

           endif ! MGT(i,1,n)
          enddo !i
        enddo !n


!#if ( CORAL == 1 )
!      if (KENDY.eq.1) then !if last day of year
!      print*, ' '
!      print*, 'corals: mass_carb= ', mass_carb
!      print*, 'corals: total_area_coral_an (km2)= ', total_area_coral_an
!      total_area_coral_an=0
!      print*, 'corals: total_prod_coral_an= ', total_prod_coral_an
!      !print*, 'corals: total_mass_coral_an= ', total_mass_coral_an
!      endif
!#endif

#if ( CORAL == 1 )
      call calc_nino3
#endif

#if ( PATH >= 1 )
!        call particlesFields(TPP_mas,caco3_mas,JPROD,MGT)
#endif




!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr --- en fin d annee
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!         calculation of total inventory of tracers every year
!         calculation of total variables and calculational errors


      if (KENDY.eq.1) then

!dmr --- on reinitialise les variables
        total_phyto=0
        total_zoo=0
!dmr --- DOC total
        total_doc=0
!dmr --- DOC total "lent" (slow carbon pool)
        total_sdoc=0
        total_poc=0
!dmr --- Carbone total oceanique
        total_car=0
        total_phos=0
        total_alk=0
        total_nit=0
        total_d13=0

       do n=1,NOC_CBR
        do i=1,LT
          do j=1,JT
          if (MGT(i,j,n).eq.1) then

!nb  carbone organique
!PROBLEM: Not sure what is done here: Oeta(j,4) = POC:POP and Oeta(j,5) = DOC:DOP, Why?
!         Inconsistent mixture of POM and DOM C:P and O2:P ratios
!RFC         vtmp=((PHYTO_M(i,J,n)+ZOO_M(i,J,n)+
!     <          ODOC(i,J,n)+ODOCS(i,J,n))*Oeta(j,4)/Oeta(j,5)
!     <       +OPOC(i,j,n))*DVOL(i,J,n)/OVOL
!nb we think this is wrong as mixing poc and doc
!         vtmp=((PHYTO_M(i,J,n)+ZOO_M(i,J,n)+
!     <          ODOC(i,J,n)+ODOCS(i,J,n))*
!     <          OetaC_POMoxid(i,J,n)/OetaC_DOMoxid_1D(j)
!     <       +OPOC(i,j,n))*DVOL(i,J,n)/OVOL

         vtmp_POC=(PHYTO_M(i,J,n)+ZOO_M(i,J,n))*DVOL(i,J,n)/OVOL

         vtmp_DOC=(ODOC(i,J,n)+ODOCS(i,J,n))*DVOL(i,J,n)/OVOL

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!PROBLEM: 1.028 density conversion factor included; leads to consistent units,
!         but inconsistent with other conversions
        total_phyto=total_phyto+PHYTO_M(i,J,n)*DVOL(i,J,n)*SCALE_M*12*1.028 !nb en GtC
        total_zoo=total_zoo+ZOO_M(i,J,n) * DVOL(i,J,n)*SCALE_M*12*1.028 !nb en GtC
        total_doc=total_doc+ODOC(i,J,n) * DVOL(i,J,n)*SCALE_M*12*1.028 !nb en GtC
        total_sdoc=total_sdoc+ODOCS(i,J,n)*DVOL(i,J,n)*SCALE_M*12*1.028 !nb en GtC
!        total_poc=total_poc+OPOC(i,J,n)*
!     <  DVOL(i,J,n)*SCALE_M*12*1.028 !nb en GtC
!dmr --- Carbone total = DIC partout dans l ocean
        total_car=total_car+ODIC(i,J,n)*DVOL(i,J,n)*SCALE_M*12*1.028*1e6 !nb en GtC

! modified - 09.06.07
!-----|--1--------2---------3---------4---------5---------6---------7-|
!dmr --- Comment: since we are in a 3D ocean, fhypt is 1.0 everywhere  (!!!)
!-----|--1--------2---------3---------4---------5---------6---------7-|
        total_alk=total_alk+OALK(i,J,n)*DVOL(i,J,n)/OVOL*1000

! modified - 09.06.07
!PROBLEM: vtmp is a mixture of DIC from POC and DOC, so a unique
!  conversion factor is inadequate
!        total_phos=total_phos+vtmp/OetaC_POMoxid(i,j,n)+OPO4(i,J,n)*
!     >  DVOL(i,J,n)/OVOL

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!nb try to separate POC and DOC
        total_phos=total_phos+vtmp_POC/OetaC_POMoxid(i,j,n)+vtmp_DOC/OetaC_DOMoxid_1D(j)+OPO4(i,J,n)*DVOL(i,J,n)/OVOL

! modified - 09.06.07
!PROBLEM: vtmp is a mixture of DIC from POC and DOC, so a unique
!  conversion factor is inadequate
!        total_nit=total_nit+vtmp/OetaC_POMoxid(i,j,n)*
!     >  OetaN_POMoxid(i,j,n)+ONO3(i,J,n) * DVOL(i,J,n)/OVOL

!nb try to separate POC and DOC
        total_nit=total_nit+vtmp_POC/OetaC_POMoxid(i,j,n)*OetaN_POMoxid(i,j,n)+vtmp_DOC/OetaC_DOMoxid_1D(j)*OetaN_DOMoxid_1D(j) &
                           +ONO3(i,J,n) * DVOL(i,J,n)/OVOL

        total_d13=total_d13+OC13(i,J,n)*DVOL(i,J,n)*1000/OVOL
          endif ! if MGT ...
!dmr --- Boucles spatiales
          enddo
        enddo
       enddo

!        print*,'total_phos, OPO4_ini', total_phos, OPO4_ini
!        print*,'total_nit, ONO3_ini', total_nit, ONO3_ini
!        print*,'total_alk, OALK_ini', total_alk/1000, OALK_ini

!         correction of calculation errors: ALC, PO4, NO3
! ....
!dmr --- si deuxieme annee du siecle ... !!!
        if(NYR-(NYR/100)*100.eq.2) then

        print*,'total_alk, OALK_ini', total_alk/1000, OALK_ini

          temp_difa=total_alk/1000.-OALK_ini

          tot_err=tot_err-temp_difa

!nb !!!OPO4_ini au lieu de 2.08 !!
!        temp_difp=total_phos-(OPO4_ini+total_sdoc/OetaC_DOMoxid_1D(1))
        temp_difp=total_phos-OPO4_ini
!        temp_difp=total_phos-(2.08 +4.434/OetaC_DOMoxid_1D(1))

!PROBLEM: inconsistent mixture of POM and DOM C:P and O2:P ratios
!  If total_sdoc is from DOM, replace OetaN_POMoxid(:,1,:)) by OetaN_DOMoxid_1D(1)
        tot_err2=tot_err2-temp_difp
        !temp_difn=total_nit-(ONO3_ini+total_sdoc/OetaC_DOMoxid_1D(1)*OetaN_POMoxid(:,1,:))
        temp_difn=total_nit-ONO3_ini
!        temp_difn=total_nit-(33.6+4.434/OetaC_DOMoxid_1D(1)*OetaN_POMoxid(:,1,:))

        print*, 'temp_difa, temp_difp, temp_difn', temp_difa, temp_difp,temp_difn

! homogeneous correction
!#if ( MEDUSA == 0 )

!nb test       do n=1,NOC_CBR
!nb test        do i=1,LT
!nb test          do j=1,JT
!nb test          if (MGT(i,j,n).eq.1) then


! dissolution of removed alkalinity in the surface water


!nb test        if(OALK(i,j,n).gt.temp_difa) OALK(i,j,n)=OALK(i,j,n)-temp_difa
!nb test        if(OPO4(i,j,n).gt.temp_difp) OPO4(i,j,n)=OPO4(i,j,n)-temp_difp
!nb test        if(ONO3(i,j,n).gt.temp_difn) ONO3(i,j,n)=ONO3(i,j,n)-temp_difn

!nb test          endif
!nb test          enddo
!nb test        enddo
!nb test       enddo

!#endif
        endif ! 2eme annee du siecle

!nb - sediments
        !dmr --- [DELETE] -> suppressed the whole lot since caco3 and
        !                    tpp have changed unit and are in Tmol.m-2 now

       do n=1,NOC_CBR
        do i=1,LT
          do j=1,JT
          if (MGT(i,j,n).eq.1) then

!dmr ###         total_tpp=total_tpp+TPP_ma(i,j,n)*12*1.E-3
!dmr ###     <   *(SQRO2(i,n)*1.E-8)

!         if (j.LT.JT) then
!           if (MGT(i,j+1,n).LT.1) then


!           total_tpp=total_tpp+TPP_ma(i,j,n)*SCALE_M/SCALE_B*12
!           total_tpp=total_tpp+TPP_ma(i,j,n)*DVOL(i,j,n)*SCALE_M*12
           total_tpp=total_tpp+TPP_ma(i,j,n)*DVOL(i,j,n)*SCALE_M*12*1.028
           total_caco3=total_caco3+caco3_ma(i,j,n)*(DVOL(i,j,n)*SCALE_M*12) ! remplacer 12 par mass molaire CacO3 ?

!           endif

!         else if (j.eq.JT) then
           !total_tpp=total_tpp+TPP_ma(i,j,n)*SCALE_M/SCALE_B*12
!           total_tpp=total_tpp+TPP_ma(i,j,n)*DVOL(i,j,n)*SCALE_M*12
!           total_caco3=total_caco3+caco3_ma(i,j,n)*(DVOL(i,j,n)*SCALE_M*12)

!         endif ! on j==JT

         endif ! on MGT ...

            enddo
          enddo
        enddo

!nb       dans C_ocean.txt
       write (c_ocean_fich%id,'(5F12.6,1F20.3,5F12.6,9F8.3,3F20.3)') &
!       1           2         3         4
        total_phyto,total_zoo,total_doc,total_sdoc,                  &
!       5         6         7         8         9
        total_poc,total_car,total_alk,total_nit,total_phos,          &
!       10
        total_d13/total_car,                                         &
!       11
        total_d13/(total_car+(total_doc+total_sdoc)*SCANU),          &
!       12        13        14
        prod_sum, graz_sum, remin_sum,                               &
!       15      16       17         18       19       20
        exc_sum,exu_sum, barem_sum, tpp_sum, ort_sum, pel_sum,       &
!       21         22           23
        total_tpp, total_caco3, caco3_sum

        prod_sum=0
        graz_sum=0
        remin_sum=0
        exc_sum=0
        exu_sum=0
        barem_sum=0
        tpp_sum=0
        ort_sum=0
        pel_sum=0
        caco3_sum=0

        total_tpp = 0.0
!        FORALL (n=1:NOC_CBR, i=1:LT, j=1:JT)
!           TPP_ma(i,j,n)=0.0
!        END FORALL

        total_caco3 = 0.0


       write (fPOC_fich%id, '(5F8.3,1F20.3,5F8.3,9F8.3,3F20.3)')    &
!           1         2         3
        fPOC_top, fPOC_1000, fPOC_bot

        fPOC_top  = 0.0
        fPOC_1000 = 0.0
        fPOC_bot  = 0.0


       write (fCAL_fich%id, '(5F8.3,1F20.3,5F8.3,9F8.3,3F20.3)')    &
!           1         2         3
        fCAL_top, fCAL_2000, fCAL_bot

        fCAL_top  = 0.0
        fCAL_2000 = 0.0
        fCAL_bot  = 0.0

#if ( CORAL == 1 )
!nb corals
      print*, ' '
      print*, 'CORAL at last day of the year'
      print*, 'total_area_coral_an in 1e3 km2 = ',  total_area_coral_an*1e-3
      print*, 'total prod_coral_an in Pmol/year = ' , total_prod_coral_an
      print*, 'total mass_coral_an in Pg/year = ', total_mass_coral_an
      print*, 'coral_CO2=', coral_CO2, 'gC/an'


      temp_too_low_all(:,:)=0.0

#if ( 1 )
!nb To fix input from rivers
      !C_riv= 2* C_sed
!      C_riv = 2* total_prod_coral_an/(TYER/TDAY) ! in Pmol/day
      A_riv = C_riv
      write(*,*) 'weathering computed as output: ', C_riv, total_prod_coral_an
      C_car_a=C_riv/2.
#endif

! write annual global mean in file coral_output.txt
      CALL out_coral_global

!nb remise a 0 a la fin de l annee
! and save 3d output
      total_area_coral_an=0.0
      total_prod_coral_an=0.0
      total_mass_coral_an=0.0
      coral_CO2=0
      coral_area_out(:,:,:)=coral_area(:,:,:)
      coral_prod_out(:,:,:)=coral_prod(:,:,:)
      coral_mass_out(:,:,:)=coral_cum_mass(:,:,:)
      tau_bleach_out(:,:,:)=tau_bleach(:,:,:)
      DHW_out(:,:,:)=DHW_nb(:,:,:)
      ! do n=1,NOC_CBR
      !  do i=1,LT
      !    do j=1,JT
      !    if (MGT(i,j,n).eq.1) then
      !      if (coral_prod_out(i,j,n).ne.0) then
      !  write(*,*) 'coral_prod_out mbiota ', coral_prod_out(i,j,n)
      !      endif
      !    endif
      !    enddo
      !  enddo
      !enddo
      coral_area(:,:,:)=0.0
      coral_prod(:,:,:)=0.0
      DHW_nb(:,:,:)=0.0
#endif
        endif ! de KENDY == 1

        return
        END SUBROUTINE MBIOTA

      END MODULE MOD_PHOTIC_ZONE
