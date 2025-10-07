!dmr -- Ajout du choix optionnel des composantes - Wed Aug 18 16:16:38 CEST 2010
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Wed Aug 18 16:16:38 CEST 2010
!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Cette routine restart_veccarb sert a ecrire / lire le redemarrage
!       des valeurs globales de carbone dans un fichier
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : N Bouttes
!      Date   : 26 February 2019
!      Derniere modification :
!-----|--1--------2---------3---------4---------5---------6---------7-|

       SUBROUTINE restart_cc(choix)

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree  :
!         Choix : ecriture ou lecture du restart ...
!       Variables de sortie :
!         Neant
!-----|--1--------2---------3---------4---------5---------6---------7-|

!   INSERER ICI LES EVENTUELS "USE MODULE"
!       USE declars_mod !cnb
       USE C_res_mod, ONLY: ca_oc_rest, ca_la_rest, cav_oc, cav_la      &
                         ,ca13_oc_rest,ca13_la_rest,cav_oc13,cav_la13   &
                         ,c13atm,c13atm_rest, alk_oc_rest
#if ( KC14 == 1 )
       use C_res_mod, only: ca14_oc_rest,ca14_la_rest,cav_oc14,cav_la14
#endif

       USE carbone_co2, ONLY: PA0_C, PA_C, C14ATM, C14ATM_rest

       USE marine_bio_mod, only: OALK_ini
#if ( BATHY >= 1 )
       USE loveclim_transfer_mod, ONLY: OVOL, OVOL_prev
#endif

       IMPLICIT NONE


       INTEGER :: choix

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Variables locales
!-----|--1--------2---------3---------4---------5---------6---------7-|


       INTEGER :: fich_num, nrecl
       CHARACTER*11, PARAMETER :: fich_res_name="rest_cc.dat"
       CHARACTER*21, PARAMETER ::                                       &
                         fich_res_name_old="startdata/rest_cc.dat"
       LOGICAL :: existe
       INTEGER :: n,i,j
!nb & fl
       logical :: logic_elmt

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Determination d'un numero de fichier libre ...
!-----|--1--------2---------3---------4---------5---------6---------7-|

       fich_num=298
       print*, "init_CC_fich", fich_num
       existe=.TRUE.

       DO WHILE (existe)
        INQUIRE(fich_num,OPENED=existe)
        print*, "init_CC_fich", fich_num, existe
        fich_num=fich_num+1
       END DO

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Determination de la taille de l'ecriture ...
!-----|--1--------2---------3---------4---------5---------6---------7-|
       nrecl = 10 !nb of variables written

       nrecl = nrecl * KIND(PA0_C)
!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Cas ou l'on ecrit le restart
!-----|--1--------2---------3---------4---------5---------6---------7-|
       IF (choix.EQ.0) THEN
        ca_oc_rest=cav_oc
        ca_la_rest=cav_la
        ca13_oc_rest=cav_oc13
        ca13_la_rest=cav_la13
        c13atm_rest=c13atm
        alk_oc_rest=OALK_ini
#if ( KC14 == 1 )
        ca14_oc_rest=cav_oc14
        ca14_la_rest=cav_la14
        c14atm_rest=c14atm
#endif

        OPEN(UNIT=fich_num, FILE=fich_res_name, STATUS='unknown',       &
              ACCESS='direct', RECL=nrecl, ACTION='write')

#if( KC14 == 1 )
        WRITE(UNIT=fich_num, REC=1)                                     &
                     ca_oc_rest, ca_la_rest, PA_C                       &
                     ,ca13_oc_rest, ca13_la_rest, c13atm_rest           &
                     ,alk_oc_rest                                       &
                     ,ca14_oc_rest, ca14_la_rest, c14atm_rest           
#else
        WRITE(UNIT=fich_num, REC=1)                                     &
                     ca_oc_rest, ca_la_rest, PA_C                       &
                     ,ca13_oc_rest, ca13_la_rest, c13atm_rest           &
                     ,alk_oc_rest
#endif

        CLOSE(UNIT=fich_num)

        write(*,*), "write carbon values in rest_cc.dat"
        write(*,*), ca_oc_rest,ca_la_rest,ca13_oc_rest,ca13_la_rest
        write(*,*), PA_C, c13atm_rest, alk_oc_rest

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Cas ou l'on lit le restart
!-----|--1--------2---------3---------4---------5---------6---------7-|
       ELSE IF (choix.EQ.1) THEN

!       nrecl = 9 !9 variables written
!       nrecl = nrecl * KIND(PA0_C)

      WRITE(*,*) "Lecture du redemarrage carbone a partir du fichier : "
      WRITE(*,*) fich_res_name_old

        OPEN(UNIT=fich_num, FILE=fich_res_name_old, STATUS='old',       &
              ACCESS='direct', RECL=nrecl, ACTION='read')

#if ( KC14 == 1 )
        READ(UNIT=fich_num,REC=1)                                       &
                     ca_oc_rest, ca_la_rest, PA_C                       &
                     ,ca13_oc_rest, ca13_la_rest, c13atm_rest           &
                     ,alk_oc_rest                                       &
                     ,ca14_oc_rest, ca14_la_rest, c14atm_rest
#else
        READ(UNIT=fich_num,REC=1)                                       &
                     ca_oc_rest, ca_la_rest, PA_C                       &
                     ,ca13_oc_rest, ca13_la_rest, c13atm_rest           &
                     ,alk_oc_rest                                       
#endif

        CLOSE(UNIT=fich_num)

       PA0_C = PA_C
       c13atm=c13atm_rest
#if ( KC14 == 1 )
       c14atm=c14atm_rest
#endif
       write(*,*), "initialisation of atmospheric carbon"
       write(*,*), PA0_C, c13atm, c14atm

!alkalinity in ocean
       OALK_ini=alk_oc_rest
!nb test
#if ( BATHY >= 1 )
       !write(*,*) 'OVOL_prev, OVOL ', OVOL_prev, OVOL
       OALK_ini=alk_oc_rest*OVOL_prev/OVOL
#endif
       !write(*,*) 'OALK_ini from restart: ', OALK_ini

       ENDIF

       END SUBROUTINE restart_cc
