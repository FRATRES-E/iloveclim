!dmr -- Ajout du choix optionnel des composantes - Wed Aug 18 16:16:38 CEST 2010
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Wed Aug 18 16:16:38 CEST 2010
!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Cette routine restart_mb sert a ecrire / lire le redemarrage
!       des valeurs carbones dans un fichier
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche 
!      Date   : 18 Aout 2010
!      Derniere modification : 01 Decembre 2025 (nbouttes & ecl)
!-----|--1--------2---------3---------4---------5---------6---------7-|
       SUBROUTINE restart_mb(choix)

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree  : 
!         Choix : ecriture ou lecture du restart ... 
!       Variables de sortie : 
!         Neant
!-----|--1--------2---------3---------4---------5---------6---------7-|

!   INSERER ICI LES EVENTUELS "USE MODULE"
       USE declars_mod !cnb

#if ( BATHY >= 1) 
       USE loveclim_transfer_mod, ONLY:MGT, MGT_prev, DVOL, OVOL        &
                            , OVOL_prev, DVOL_prev
#else
       USE loveclim_transfer_mod, ONLY:MGT
#endif

       USE marine_bio_mod, ONLY:OPO4, ONO3, OSI, OO2, OALK,  ODIC, ODOC,&
                     OPOC, OC13, OC14, ODOCS, ODOC13, ODOCS13,          &
                     FOPO4, FONO3, FOSI, FOO2, FOALK, FODIC, FODOC,     &

#if ( OXNITREUX == 1 )
                     ON2O, FON2O,                                       &
#endif

#if ( BATHY >= 1 )
                     global_dic,                                        &
                     global_po4,                                        &
                     global_no3,                                        &
                     global_alk,                                        &
                     global_oc13,                                       &
                     global_oc14,                                       &
                     tot_alk_prev,                                      &
                     tot_phos_prev,                                     &
                     tot_nit_prev,                                      &
!nb                     Oeta,                                              &
                     OetaC_POMoxid, OetaC_DOMoxid_1D,                   &
                     OetaN_POMoxid, OetaN_DOMoxid_1D
#endif

                     FODOCS, FOC13, FODOC13, FODOCS13 

       USE mbiota_mod, ONLY: PHYTO_M, ZOO_M, PHYTO_M13, ZOO_M13

!cnb       USE veget_iso, ONLY: C13ATM

       USE C_res_mod, ONLY: C13ATM

       USE carbone_co2, ONLY: PA0_C, PA_C, C14ATM0, C14ATM

!nb & fl in the future faire avec 3D
!       use update_clio_bathy_tools, only: mean_neighbours_with_mask
#if ( KC14 == 1 )
       USE C_res_mod, ONLY: cav_oc14_b, cav_oc_b
#endif

#if ( BATHY >= 1 )
       USE update_clio_bathy_tools, only: mean_neighbours_with_mask_CC
       USE para0_mod, ONLY: NISOO2
#endif

#if ( OOISO == 1 )
       USE para0_mod, ONLY: NISOO2
       USE iso_dioxygen_mod, ONLY: iair16, iair17, iair18, iair
#endif

       IMPLICIT NONE


       INTEGER :: choix

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Variables locales
!-----|--1--------2---------3---------4---------5---------6---------7-|


       INTEGER :: fich_num, nrecl
       CHARACTER*11, PARAMETER :: fich_res_name="rest_mb.dat"
       CHARACTER*21, PARAMETER ::                                       &
                              fich_res_name_old="startdata/rest_mb.dat"
       LOGICAL :: existe
       INTEGER :: n,i,j

!nb & fl
       logical :: logic_elmt

       integer :: km

#if ( BATHY >= 1 )
!nb    DOUBLE PRECISION vtmp
       DOUBLE PRECISION vtmp_POC, vtmp_DOC
#endif

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Determination d'un numero de fichier libre ...
!-----|--1--------2---------3---------4---------5---------6---------7-|

       fich_num=297
       print*, "init_mb_fich", fich_num
       existe=.TRUE.

       DO WHILE (existe)
        INQUIRE(fich_num,OPENED=existe)
        print*, "init_mb_fich", fich_num, existe
        fich_num=fich_num+1 
       END DO 

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Cas ou l'on ecrit le restart
!-----|--1--------2---------3---------4---------5---------6---------7-|
       IF (choix.EQ.0) THEN

!       Determination de la taille de l'ecriture
!       ----------------------------------------
!       nrecl = 13*SIZE(OPO4)+11*SIZE(FOPO4)+SIZE(PHYTO_M)+SIZE(ZOO_M)+  &
       nrecl = 13*SIZE(OPO4)+10*SIZE(FOPO4)+SIZE(PHYTO_M)+SIZE(ZOO_M)   &
               +3 ! PA0_C, PA_C, C13ATM

#if ( OOISO == 1 )
       nrecl =3*SIZE(OPO4)+SIZE(OO2)+6*SIZE(OALK)                       &
                  +SIZE(PHYTO_M)+SIZE(ZOO_M)+3*SIZE(ODOCS)              &
                  +3+3*SIZE(FOPO4)+SIZE(FOO2)+6*SIZE(FOALK)
#endif

#if ( OXNITREUX == 1 )
       nrecl = nrecl+SIZE(ON2O)+SIZE(FON2O)
#endif

#if ( KC14 == 1 )
       nrecl = nrecl+3  ! SIZE(C14ATM)+SIZE(cav_oc14_b)+SIZE(cav_oc_b)
#endif

       nrecl=nrecl+SIZE(FODOCS13)

       nrecl = nrecl * KIND(PA0_C)

!      write restart
!      -------------
        OPEN(UNIT=fich_num, FILE=fich_res_name, STATUS='unknown',       &
              ACCESS='direct', RECL=nrecl, ACTION='write')

        WRITE(UNIT=fich_num, REC=1)                                     &
             OPO4,  ONO3,  OSI,                                         & 

#if ( OOISO == 1 )
             OO2(:,:,:,1),                                              &
#else
             OO2,                                                       &
#endif

             OALK,  ODIC, ODOC,                                         &
             OPOC, OC13, OC14, PHYTO_M, ZOO_M, ODOCS, ODOC13,           &
             ODOCS13, PA0_C, PA_C, C13ATM, FOPO4, FONO3, FOSI,          &

#if ( OOISO == 1 )
             FOO2(:,:,1),                                               &
#else
             FOO2,                                                      &
#endif

             FOALK, FODIC, FODOC, FODOCS, FOC13, FODOC13,               &

#if ( OXNITREUX == 1 )
                     ON2O ,FON2O,                                       &
#endif

#if ( KC14 == 1 )
                     C14ATM, cav_oc14_b, cav_oc_b,                      &
#endif

#if ( OOISO == 1 )
                     FODOCS13, OO2(:,:,:,2:NISOO2), FOO2(:,:,2:NISOO2)
#else
                     FODOCS13
#endif

        CLOSE(UNIT=fich_num)


!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Cas ou l'on lit le restart
!-----|--1--------2---------3---------4---------5---------6---------7-|
       ELSE IF (choix.EQ.1) THEN

!       Determination de la taille de l'ecriture
!       ----------------------------------------
!       nrecl = 13*SIZE(OPO4)+11*SIZE(FOPO4)+SIZE(PHYTO_M)+SIZE(ZOO_M)+  &
       nrecl = 13*SIZE(OPO4)+10*SIZE(FOPO4)+SIZE(PHYTO_M)+SIZE(ZOO_M)   &
               +3 ! PA0_C, PA_C, C13ATM

#if ( OOISO == 1 )
       nrecl =3*SIZE(OPO4)+SIZE(OO2)+6*SIZE(OALK)                       &
                  +SIZE(PHYTO_M)+SIZE(ZOO_M)+3*SIZE(ODOCS)              &
                  +3+3*SIZE(FOPO4)+SIZE(FOO2)+6*SIZE(FOALK)
#endif

#if ( OXNITREUX == 1 )
       nrecl = nrecl+SIZE(ON2O)+SIZE(FON2O)
#endif

#if ( KC14 == 1 )
       nrecl = nrecl+3  ! SIZE(C14ATM)+SIZE(cav_oc14_b)+SIZE(cav_oc_b)
#endif

       nrecl=nrecl+SIZE(FODOCS13)

       nrecl = nrecl * KIND(PA0_C)

!     read restart
!     ------------
      WRITE(*,*) "Lecture du redemarrage carbone a partir du fichier : "
      WRITE(*,*) fich_res_name_old

        OPEN(UNIT=fich_num, FILE=fich_res_name_old, STATUS='old',       &
              ACCESS='direct', RECL=nrecl, ACTION='read')

        READ(UNIT=fich_num,REC=1)                                       &
             OPO4,  ONO3,  OSI,                                         &

#if ( OOISO == 1 )
             OO2(:,:,:,1),                                              &
#else
             OO2,                                                       &
#endif

             OALK,  ODIC, ODOC,                                         &
             OPOC, OC13, OC14, PHYTO_M, ZOO_M, ODOCS, ODOC13,           &
             ODOCS13, PA0_C, PA_C, C13ATM, FOPO4, FONO3, FOSI,          &

#if ( OOISO == 1 )
             FOO2(:,:,1),                                               &
#else
             FOO2,                                                      &
#endif

             FOALK, FODIC, FODOC, FODOCS, FOC13, FODOC13,               &

#if ( OXNITREUX == 1 )
                     ON2O ,FON2O,                                       &
#endif

#if ( KC14 == 1 )
                     C14ATM0, cav_oc14_b, cav_oc_b,                     &
#endif

#if ( OOISO == 1 )
                     FODOCS13, OO2(:,:,:,2:NISOO2), FOO2(:,:,2:NISOO2)
#else
                     FODOCS13
#endif
                     
        CLOSE(UNIT=fich_num)

!nb later : add in restart
      PHYTO_M13(:,:,:)=PHYTO_M(:,:,:)*OC13(:,:,:)/ODIC(:,:,:)
      ZOO_M13(:,:,:)=ZOO_M(:,:,:)*OC13(:,:,:)/ODIC(:,:,:)

#if ( BATHY >=1 )
!nb computes previous global values
      global_dic=0.
      global_po4=0.
      global_no3=0.
      global_alk=0.
      global_oc13=0.
      global_oc14=0.
      tot_alk_prev=0.
      tot_phos_prev=0.
      tot_nit_prev=0.
      do j=1,JT
       do i=1,LT
        do n=1,NOC_CBR
          if (MGT_prev(i,j,n).eq.1) then
          global_dic=global_dic+ODIC(i,j,n)*DVOL_prev(i,j,n)
          global_po4=global_po4+OPO4(i,j,n)*DVOL_prev(i,j,n)
          global_no3=global_no3+ONO3(i,j,n)*DVOL_prev(i,j,n)
          global_alk=global_alk+OALK(i,j,n)*DVOL_prev(i,j,n)
          global_oc13=global_oc13+                                      &
                      OC13(i,j,n)*DVOL_prev(i,j,n)
          global_oc14=global_oc14+                                      &
                      OC14(i,j,n)*DVOL_prev(i,j,n)

          tot_alk_prev=tot_alk_prev+OALK(i,j,n)                         &
                      *DVOL_prev(i,j,n)*1000
!nb          vtmp=((PHYTO_M(i,j,n)+ZOO_M(i,j,n)+ ODOC(i,j,n)+ODOCS(i,j,n)) &
!nb                      *OetaC_POMoxid_1D(j)/OetaC_DOMoxid_1D(j))         &
!nb                      *DVOL_prev(i,j,n)
         vtmp_POC=(PHYTO_M(i,j,n)+ZOO_M(i,j,n))*DVOL_prev(i,j,n)/OVOL

         vtmp_DOC=(ODOC(i,J,n)+ODOCS(i,J,n))*DVOL(i,J,n)/OVOL

!nb          tot_phos_prev=tot_phos_prev+vtmp/Oeta(j,4)+OPO4(i,j,n)        &
!nb                      *DVOL_prev(i,j,n)
          tot_phos_prev=tot_phos_prev+vtmp_POC/OetaC_POMoxid(i,j,n)     &
          +vtmp_DOC/OetaC_DOMoxid_1D(j)+OPO4(i,j,n)*                    &
                      *DVOL_prev(i,j,n)

!nb          tot_nit_prev=tot_nit_prev+vtmp/Oeta(j,4)*Oeta(j,1)+ONO3(i,j,n)&
!nb                      *DVOL_prev(i,j,n)

          tot_nit_prev=tot_nit_prev+vtmp_POC/OetaC_POMoxid(i,j,n)*      &
          OetaN_POMoxid(i,j,n)+vtmp_DOC/OetaC_DOMoxid_1D(j)*            &
          OetaN_DOMoxid_1D(j)+ONO3(i,j,n) * DVOL_prev(i,j,n)

           endif
          enddo
        enddo
      enddo

      write(*,*) 'global dic in restart_mb ', global_dic

     write(*,*) 'tot_alk_prev dans restart_mb', tot_alk_prev
     write(*,*) 'OALK_ini avec OVOL_prev=', tot_alk_prev/OVOL_prev/1000.
     write(*,*) 'OALK_ini avec OVOL=', tot_alk_prev/OVOL/1000.



      write(*,*) 'UPDATE in restart_mb '

!nb initialise values for new grid cells
       logic_elmt= mean_neighbours_with_mask_CC(OPO4(:,:,:),1) ! last argument 1 = tms2D
       logic_elmt= mean_neighbours_with_mask_CC(ONO3(:,:,:),1) ! last argument 1 = tms2D
       logic_elmt= mean_neighbours_with_mask_CC(OSI(:,:,:),1) ! last argument 1 = tms2D

#if ( ISOO2 == 0 )
       logic_elmt= mean_neighbours_with_mask_CC(OO2(:,:,:),1) ! last argument 1 = tms2D
#else
       do km=1,NISOO2
         logic_elmt= mean_neighbours_with_mask_CC(OO2(:,:,:,km),1) ! last argument 1 = tms2D
       enddo
#endif

       logic_elmt= mean_neighbours_with_mask_CC(OALK(:,:,:),1) ! last argument 1 = tms2D
       logic_elmt= mean_neighbours_with_mask_CC(ODIC(:,:,:),1) ! last argument 1 = tms2D
       logic_elmt= mean_neighbours_with_mask_CC(ODOC(:,:,:),1) ! last argument 1 = tms2D
!       logic_elmt= mean_neighbours_with_mask_CC(OPOC(:,:,:),1) ! last argument 1 = tms2D
       logic_elmt= mean_neighbours_with_mask_CC(OC13(:,:,:),1) ! last argument 1 = tms2D
       logic_elmt= mean_neighbours_with_mask_CC(OC14(:,:,:),1) ! last argument 1 = tms2D
       logic_elmt= mean_neighbours_with_mask_CC(PHYTO_M(:,:,:),1) ! last argument 1 = tms2D
       logic_elmt= mean_neighbours_with_mask_CC(ZOO_M(:,:,:),1) ! last argument 1 = tms2D
       logic_elmt= mean_neighbours_with_mask_CC(ODOCs(:,:,:),1) ! last argument 1 = tms2D
       logic_elmt= mean_neighbours_with_mask_CC(ODOC13(:,:,:),1) ! last argument 1 = tms2D
       logic_elmt= mean_neighbours_with_mask_CC(ODOCs13(:,:,:),1) ! last argument 1 = tms2D


#endif
!--- !cnb modif valeur oc13 globale (a supprimer ensuite)
!       print*, 'Modification OC13 ini'
!       do n=1,NOC_CBR
!         do i=1,LT
!           do j=1,JT
!           if (MGT(i,j,n).eq.1) then 
!---             OC13(i,j,n)=OC13(i,j,n)+0.65*ODIC(i,j,n) !0.35
!             ODIC(i,j,n)=ODIC(i,j,n)*1.023 !1.000
!             OC13(i,j,n)=OC13(i,j,n)+1.0*ODIC(i,j,n) !
!              OPO4(i,j,n)=OPO4(i,j,n)-0.1
! re initialisation des variables
!             OO2(i,j,n,1) = 250
!             OPO4(i,j,n) = 2.08
!             ONO3(i,j,n) = 33.6
!             OALK(i,j,n) = 2373.0e-6
!             ODOCs(i,j,n) = 4.434
!           endif
!           enddo
!         enddo
!       enddo

!!--- dmr modif valeur odic globale (a supprimer ensuite)
!       print*, 'Modification ODIC ini'
!       do n=1,NOC_CBR
!         do i=1,LT
!           do j=1,JT
!           if (MGT(i,j,n).eq.1) then 
!             ODIC(i,j,n)=ODIC(i,j,n)*0.998 !1.000 0.0995
!           endif
!           enddo
!         enddo
!       enddo

!cnb moved to restart_veccarb       PA0_C = PA_C

!cnb modify values for LGM : increase of 3.3percent due to sea level
!change.
!#if ( LGMSWITCH == 1 )
!       print*, 'Modification OPO4, ONO3, ALK for LGM in restart_mb'
!       do n=1,NOC_CBR
!         do i=1,LT
!           do j=1,JT
!           if (MGT(i,j,n).eq.1) then 
!             OPO4(i,j,n)=OPO4(i,j,n)*1.033
!             ONO3(i,j,n)=ONO3(i,j,n)*1.033
!             ALK(i,j,n)=ALK(i,j,n)*1.033
!           endif
!           enddo
!         enddo
!       enddo
!#endif

       ENDIF

       END SUBROUTINE restart_mb
