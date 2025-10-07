!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:21 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:21 CET 2009
c********************************************************************
                    SUBROUTINE INIT_OCC
c********************************************************************
c              OCEAN CARBON CYCLE MODEL INITIALISATION
c
c  Purpose: preparation of initial conditions and parameters
c           for ocean carbon cycle model
c
c  By V.Brovkin
c  Last modification: 04.04.2001
c********************************************************************
        use global_constants_mod, only: dblp=>dp, ip
        use declars_mod, only: LT, JT, NOC_CBR
        use loveclim_transfer_mod, only: KLSR, MGT, DVOL, OVOL
        use marine_bio_mod, only: FOPO4, FONO3, FOSI, FOO2, FOALK, FODIC
     >                    , FODOC, FODOCS, FOC13, FOC14, FODOC13
     >                    , FODOCS13, FOC14, OALK, ODOC, ODOCS
     >                    , OPO4, ONO3, OALK_INI, OPO4_INI
     >                    , ONO3_INI, ODIC, OC13, OC14
     >                    , OetaC_POMoxid, OetaC_DOMoxid_1D
     >                    , OetaN_POMoxid, OetaN_DOMoxid_1D
#if ( BATHY >= 1 )
     >                    , global_dic, global_alk, global_po4
     >                    , global_no3, global_oc13, global_oc14
     >                    , tot_alk_prev,tot_phos_prev, tot_nit_prev
        use loveclim_transfer_mod, only: OVOL_prev
#endif
        use mbiota_mod, only: PHYTO_M, ZOO_M
        use para0_mod, ONLY: NISOO2
        use C_res_mod, ONLY: alk_oc_rest

#if ( CORAL == 1 )
       use coral_mod, only: ini_coral
#endif

      IMPLICIT NONE

      INTEGER(kind=ip) :: i,j,n, km
      REAL(kind=dblp)  :: tot_alk, tot_phos, vtmp, vtmp_POC, vtmp_DOC
      REAL(kind=dblp)  :: tot_nit, vol_cor
      INTEGER(kind=ip) :: bgc_hyps_id

#if ( BATHY >= 1 )      
      REAL(kind=dblp)  :: global_dic_time, diff_dic 
      REAL(kind=dblp)  :: global_po4_time, diff_po4 
      REAL(kind=dblp)  :: global_no3_time, diff_no3 
      REAL(kind=dblp)  :: global_alk_time, diff_alk 
      REAL(kind=dblp)  :: global_oc13_time, diff_oc13 
      REAL(kind=dblp)  :: global_oc14_time, diff_oc14 
#endif
c********************************************************************


c     =================
       call INITMBIOPAR
c     =================

#if ( CORAL == 1 )
! Corals
      call ini_coral
#endif

      if(KLSR.eq.0) then
        call INIT_MB

        FOPO4(:,:) = 0.0
        FONO3(:,:) = 0.0
        FOSI(:,:)  = 0.0

        FOO2(:,:,:)= 0.0

        FOALK(:,:) = 0.0
        FODIC(:,:) = 0.0
        FODOC(:,:) = 0.0
        FODOCS(:,:) = 0.0
        FOC13(:,:) = 0.0
#if ( KC14 == 1 )
        FOC14(:,:) = 0.0
#endif
        FODOC13(:,:) = 0.0
        FODOCS13(:,:) = 0.0
#if ( OXNITREUX == 1 )
        FON2O(:,:) = 0.0
#endif

      else
        call restart_mb(1)
      endif

#if ( BATHY >=1 )
cnb volume change from last state

      global_dic_time=0.
      global_po4_time=0.
      global_no3_time=0.
      global_alk_time=0.
      global_oc13_time=0.
      global_oc14_time=0.
      do j=1,JT
       do i=1,LT
        do n=1,NOC_CBR
          if (MGT(i,j,n).eq.1) then
          global_dic_time=global_dic_time+ODIC(i,j,n)*DVOL(i,j,n)
          global_po4_time=global_po4_time+OPO4(i,j,n)*DVOL(i,j,n)
          global_no3_time=global_no3_time+ONO3(i,j,n)*DVOL(i,j,n)
          global_alk_time=global_alk_time+OALK(i,j,n)*DVOL(i,j,n)
          global_oc13_time=global_oc13_time+OC13(i,j,n)*DVOL(i,j,n)
          global_oc14_time=global_oc14_time+OC14(i,j,n)*DVOL(i,j,n)
           endif
          enddo
        enddo
      enddo

      write(*,*) 'global dic prev, now ', global_dic, global_dic_time

      diff_dic=global_dic_time-global_dic
      diff_po4=global_po4_time-global_po4
      diff_no3=global_no3_time-global_no3
      diff_alk=global_alk_time-global_alk
      diff_oc13=global_oc13_time-global_oc13
      diff_oc14=global_oc14_time-global_oc14
      write(*,*) 'dic correction diff ', diff_dic

      do j=1,JT
       do i=1,LT
        do n=1,NOC_CBR
          if (MGT(i,j,n).eq.1)  then
           ODIC(i,j,n)=ODIC(i,j,n)-diff_dic/OVOL !*DVOL(i,j,n)/DVOL(i,j,n)
           OPO4(i,j,n)=OPO4(i,j,n)-diff_po4/OVOL !*DVOL(i,j,n)/DVOL(i,j,n)
           ONO3(i,j,n)=ONO3(i,j,n)-diff_no3/OVOL !*DVOL(i,j,n)/DVOL(i,j,n)
           OALK(i,j,n)=OALK(i,j,n)-diff_alk/OVOL !*DVOL(i,j,n)/DVOL(i,j,n)
           OC13(i,j,n)=OC13(i,j,n)-diff_oc13/OVOL !*DVOL(i,j,n)/DVOL(i,j,n)
           OC14(i,j,n)=OC14(i,j,n)-diff_oc14/OVOL !*DVOL(i,j,n)/DVOL(i,j,n)
          endif
        enddo
       enddo
      enddo

!verif
      global_dic_time=0.
      global_po4_time=0.
      global_no3_time=0.
      global_alk_time=0.
      global_oc13_time=0.
      global_oc14_time=0.
      do j=1,JT
       do i=1,LT
        do n=1,NOC_CBR
          if (MGT(i,j,n).eq.1)  then
          global_dic_time=global_dic_time+ODIC(i,j,n)*DVOL(i,j,n)
          global_po4_time=global_po4_time+OPO4(i,j,n)*DVOL(i,j,n)
          global_no3_time=global_no3_time+ONO3(i,j,n)*DVOL(i,j,n)
          global_alk_time=global_alk_time+OALK(i,j,n)*DVOL(i,j,n)
          global_oc13_time=global_oc13_time+OC13(i,j,n)*DVOL(i,j,n)
          global_oc14_time=global_oc14_time+OC14(i,j,n)*DVOL(i,j,n)
          endif
          enddo
        enddo
      enddo

      write(*,*) 'verification global_dic_time ', global_dic_time



cnb correction to account for volume change due to bathymetry change
!       do n=1,NOC_CBR
!         do i=1,LT
!           do j=1,JT
!           if (MGT(i,j,n).eq.1) then 
!             ODIC(i,j,n)=ODIC(i,j,n)*OVOL_chge ! manque * vol grid cell /
!                                               ! volume total
!             OPO4(i,j,n)=OPO4(i,j,n)*OVOL_chge
!             ONO3(i,j,n)=ONO3(i,j,n)*OVOL_chge
!             OALK(i,j,n)=OALK(i,j,n)*OVOL_chge
!             ODOC(i,j,n)=ODOC(i,j,n)*OVOL_chge
!             ODOCs(i,j,n)=ODOCs(i,j,n)*OVOL_chge
!             OC13(i,j,n)=OC13(i,j,n)*OVOL_chge
!             OC14(i,j,n)=OC14(i,j,n)*OVOL_chge
!             ODOC13(i,j,n)=ODOC13(i,j,n)*OVOL_chge
!             ODOCs13(i,j,n)=ODOCs13(i,j,n)*OVOL_chge
!           endif
!           enddo
!         enddo
!       enddo
#endif

!-----|--1--------2---------3---------4---------5---------6---------7-|
c  ---   Initial alkalinity from restart file
!-----|--1--------2---------3---------4---------5---------6---------7-|
       tot_alk=0.
       tot_phos=0.
       tot_nit=0.
       do n=1,NOC_CBR
        do i=1,LT
          do j=1,JT
          if (MGT(i,j,n).eq.1) then
        tot_alk=tot_alk+OALK(i,J,n)*
     >  DVOL(i,J,n)*1000

!nb vtmp separated into POC and DOC
         vtmp_POC=(PHYTO_M(i,J,n)+ZOO_M(i,J,n))*DVOL(i,J,n)/OVOL

         vtmp_DOC=(ODOC(i,J,n)+ODOCS(i,J,n))*DVOL(i,J,n)/OVOL

!nb try to separate POC and DOC
        tot_phos=tot_phos+vtmp_POC/OetaC_POMoxid(i,j,n)
     >  +vtmp_DOC/OetaC_DOMoxid_1D(j)+OPO4(i,J,n)*
     >  DVOL(i,J,n)/OVOL

!nb try to separate POC and DOC
        tot_nit=tot_nit+vtmp_POC/OetaC_POMoxid(i,j,n)*
     >  OetaN_POMoxid(i,j,n)+vtmp_DOC/OetaC_DOMoxid_1D(j)*
     >  OetaN_DOMoxid_1D(j)+ONO3(i,J,n)*DVOL(i,J,n)/OVOL


          endif
         enddo
        enddo
       enddo

        OALK_ini=tot_alk/1000/OVOL
        !OALK_ini=alk_oc_rest ! now read in restart file directly
        !OALK_ini=2.365136591371019E-003 ! test
        OPO4_ini=tot_phos/OVOL
        ONO3_ini=tot_nit/OVOL
        print *,'OALK_ini ori',OALK_ini
        print *,'OPO4_ini ori',OPO4_ini
        print *,'ONO3_ini ori',ONO3_ini

#if ( BATHY >=1 )
        OALK_ini=tot_alk_prev/1000./OVOL
        !OALK_ini=alk_oc_rest*OVOL_prev/OVOL
        OPO4_ini=tot_phos_prev/OVOL
        ONO3_ini=tot_nit_prev/OVOL
        print *,'OALK_ini modif',OALK_ini
        print *,'OPO4_ini modif',OPO4_ini
        print *,'ONO3_ini modif',ONO3_ini
#endif

c weathering flux on the ocean grid mask

c      do i=1,LT
c        l=(i+3)/4
c       do n=1,NOC_CBR
c       if (MASKT(i,n).eq.1) then

c        RUNOC_WEAV_BIC_oc(i,n)=0.25*RUNOC_WEAV_BIC(l,n)
c        RUNOC_WEAV_13BIC_oc(i,n)=0.25*RUNOC_WEAV_13BIC(l,n)

c       endif
c       enddo
c      enddo

      return
      END SUBROUTINE INIT_OCC
