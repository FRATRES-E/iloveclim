#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:21 CET 2009
!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Cette routine sync_lcm_cc sert a synchroniser les variables
!       interne du cycle du carbone ocycc et les variables CLIO
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche 
!      Date   : 10 fevrier 2010
!      Derniere modification : 12 juillet 2010, Didier M. Roche
!-----|--1--------2---------3---------4---------5---------6---------7-|

       SUBROUTINE sync_lcm_ocycc(sens_sync)

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree  :  les variables ocycc a advecter
!       Variables de sortie :  le champ scal(:,:,:,3:nsmax)
!-----|--1--------2---------3---------4---------5---------6---------7-|


!! START_OF_USE_SECTION

       USE loveclim_transfer_mod, ONLY: TM, SM
       USE marine_bio_mod, ONLY: ODOC, ODOCS, ODIC, OPO4, ONO3, OALK
     & ,OO2, OC13, ODOC13, ODOCS13, OC14, oxpCO2
       USE iso_dioxygen_mod, ONLY: iair

       USE bloc0_mod, ONLY: scal, oxpco2_clio, TPPma_clio, CACO3ma_clio
     & , fPOC_flx_clio, fCAL_flx_clio
       USE para0_mod, only: imax, jmax, kmax, nsmax

       USE declars_mod, only: LT, JT, NOC_CBR 
       USE loveclim_transfer_mod, ONLY:MGT


#if ( OXNITREUX == 1 )
       USE marine_bio_mod, ONLY: ON2O
#endif
c~        USE mbiota_mod, ONLY: TPP_mas, caco3_mas

#if ( PATH >= 1 )
       USE path_mod, only: Aactivity, scalstart_path, scalend_path
       USE path_mod, only: particles_fluxes_field, ncaco3, npoc, nopal
#endif

#if ( NEOD > 0 )
!       use neodymium_mod, only: neodymium,
!      & scalstart_neodymium, scalend_neodymium
       USE neodymium_mod, only: neodymium
       USE neodymium_mod, only: scalstart_neodymium
       USE neodymium_mod, only: scalend_neodymium
#endif

#if ( OCYCC == 1 )
       USE mbiota_mod, only: TPP_ma, caco3_ma
#endif

#if ( CORAL == 1 )
       USE coral_mod, only: coral_area_out, coral_prod_out,
     >                     tau_bleach_out,
     >                     DHW_out,
     >                     coral_mass_out
      use  omega_mod, only: PH, omega_arag3D
      use marine_bio_mod, only: OCO3
#endif

#if ( REMIN == 1 )
      use mbiota_mod, only: kremin 
#endif

#if ( OOISO == 1 )
      USE para0_mod, only: oo2iso16, oo2iso17, oo2iso18
      USE iso_dioxygen_mod, ONLY: iair16, iair17, iair18
      USE mbiota_mod, ONLY: phyto_day_sum, prodO2_sum, respO2_sum
      USE mbiota_mod, ONLY : reminO2, NCP_sum
      USE bloc0_mod, ONLY: phyto_clio, prodO2_clio, respO2_clio
      USE bloc0_mod, ONLY: reminO2_clio, NCP_clio
#endif

!! END_OF_USE_SECTION
      IMPLICIT NONE

!! START_OF_INCLUDE_SECTION

!! #include "type.com"
! [SCRPTCOM] #include "para.com"
! [SCRPTCOM] #include "bloc.com"

!! END_OF_INCLUDE_SECTION

       INTEGER :: sens_sync

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Variables locales
!-----|--1--------2---------3---------4---------5---------6---------7-|

       REAL, PARAMETER :: tzero = 273.15 
       INTEGER :: k, n, i, j

!#if ( CORAL == 1 )
!       INTEGER :: i,j,n
!#endif


!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Description de la variable "scal" : 
!
!       Scal(1) => Temperature
!       Scal(2) => Salinite
!      Ces deux-ci sont recopiees initialement dans "fait_pointer" puis
!       MaJ ici (suppression de fait pointer in fine ?)

!      Ajoutees dans le cadre de OCYCC
!
!        Scal(3) => ODOC
!        Scal(4) => ODOCS
!        Scal(5) => ODIC
!        Scal(6) => OPO4
!        Scal(7) => ONO3
!        Scal(8) => OALK
!        Scal(9) => OO2
!        Scal(10) => OC13
!        Scal(11) => ODOC13
!        Scal(12) => ODOCS13
!        Scal(13) => OC14
!nb si oxygen isotopes
!        Scal(oo2iso16) => OO2(:,:,2)
!        Scal(oo2iso17) => OO2(:,:,3)
!        Scal(oo2iso18) => OO2(:,:,4)

!-----|--1--------2---------3---------4---------5---------6---------7-|

!      Les dimensions de scal sont : 
!       Scal(Lon, Lat, Levels, nb_tracers) == (imax, jmax, kmax, nsmax)
!       ... et les variables climber sont :
!       odoc(Lat, Levels, Lon) == (LT, JT, NOC_CBR)

!      Pour passer de l'un a l'autre a traceur et profondeur constants
!       il faut donc transposer lat,lon. 
!      Note a Benets : en plus les profondeurs sont inversees avec le 
!       premier niveau OCYCC en surface et au fond en Scal

!      Pour passer de l'un a l'autre il s'en suit donc aisement : 

!       FORALL (k=1:kmax)
!         odoc(:,(kmax+1-k),:) = TRANSPOSE(scal(:,:,k,3))
!       END FORALL

!nb     on ne veut pas les premiere et derniere colonnes de scal car ce
!nb     sont des copies des colonnes avant derniere et 2 (halo points)
!nb     du coup avec NOC_CBR = 120 et plus 122 :

!       FORALL (k=1:kmax)
!        FORALL (n=1:NOC_CBR)
!         odoc(:,(kmax+1-k),n) = scal(n+1,:,k,3)
!       END FORALL

!-----|--1--------2---------3---------4---------5---------6---------7-|

       if (sens_sync.EQ.1) then ! LCM -> OCYCC

!        FORALL (k=1:kmax)
!         FORALL (n=1:NOC_CBR)
        
      do n=1,NOC_CBR
       do k=1,kmax
        do i=1,LT

         if (MGT(i,(kmax+1-k),n).eq.1) then


!tbd         TM(:,(kmax+1-k),:) = TRANSPOSE(scal(:,:,k,1)) - tzero
!         TM(i,(kmax+1-k),n) = scal(n+1,i,k,1) - tzero ! dans
!         fait_pointer
!         SM(i,(kmax+1-k),n) = scal(n+1,i,k,2)
         ODOC(i,(kmax+1-k),n) = scal(n+1,i,k,3)
         ODOCS(i,(kmax+1-k),n) = scal(n+1,i,k,4)
         ODIC(i,(kmax+1-k),n) = scal(n+1,i,k,5)
         OPO4(i,(kmax+1-k),n) = scal(n+1,i,k,6)
!vm: on enleve le traceur ONO3
#if ( OXNITREUX == 1 )
         ON2O(i,(kmax+1-k),n) = scal(n+1,i,k,7)
#else
         ONO3(i,(kmax+1-k),n) = scal(n+1,i,k,7)
#endif
!vm
         OALK(i,(kmax+1-k),n) = scal(n+1,i,k,8)
         OC13(i,(kmax+1-k),n) = scal(n+1,i,k,10)
         ODOC13(i,(kmax+1-k),n) = scal(n+1,i,k,11)
         ODOCS13(i,(kmax+1-k),n) = scal(n+1,i,k,12)
         OC14(i,(kmax+1-k),n) = scal(n+1,i,k,13)
         OO2(i,(kmax+1-k),n,iair) = scal(n+1,i,k,9)

#if ( OOISO >= 1 )
         OO2(i,(kmax+1-k),n,iair16) = scal(n+1,i,k,oo2iso16)
         OO2(i,(kmax+1-k),n,iair17) = scal(n+1,i,k,oo2iso17)
         OO2(i,(kmax+1-k),n,iair18) = scal(n+1,i,k,oo2iso18)
#endif


         
            endif

           enddo
          enddo
         enddo
    
!         END FORALL
!        END FORALL

#if ( PATH >= 1 )
        Aactivity(:,:,:,:) = scal(:,:,:,scalstart_path:scalend_path)
#endif

#if ( NEOD > 0 )
        neodymium(:,:,:,:) = scal(:,:,:,scalstart_neodymium:scalend_neodymium)
#endif
#if ( MEDUSA == 0 )
        caco3_ma(:,:,:) = 0.0d0
        TPP_ma(:,:,:) = 0.0d0
#endif

       else ! OCYCC -> LCM

!        FORALL (k=1:kmax)
!         FORALL (n=1:NOC_CBR)

      do n=1,NOC_CBR
       do k=1,kmax
        do i=1,LT

         if (MGT(i,(kmax+1-k),n).eq.1) then

!tbd         scal(:,:,k,3) = TRANSPOSE(ODOC(:,(kmax+1-k),:))
         scal(n+1,i,k,3) = ODOC(i,(kmax+1-k),n)
         scal(n+1,i,k,4) = ODOCS(i,(kmax+1-k),n)
         scal(n+1,i,k,5) = ODIC(i,(kmax+1-k),n)
         scal(n+1,i,k,6) = OPO4(i,(kmax+1-k),n)
!vm: on enleve le traceur ONO3
#if ( OXNITREUX == 1 )
         scal(n+1,i,k,7) = ON2O(i,(kmax+1-k),n)
#else
         scal(n+1,i,k,7) = ONO3(i,(kmax+1-k),n)
#endif
!vm
         scal(n+1,i,k,8) = OALK(i,(kmax+1-k),n)
         scal(n+1,i,k,10) = OC13(i,(kmax+1-k),n)
         scal(n+1,i,k,11) = ODOC13(i,(kmax+1-k),n)
         scal(n+1,i,k,12) = ODOCS13(i,(kmax+1-k),n)
         scal(n+1,i,k,13) = OC14(i,(kmax+1-k),n)
         scal(n+1,i,k,9) = OO2(i,(kmax+1-k),n,iair)

#if ( OOISO >= 1 )
         scal(n+1,i,k,oo2iso16) = OO2(i,(kmax+1-k),n,iair16)
         scal(n+1,i,k,oo2iso17) = OO2(i,(kmax+1-k),n,iair17)
         scal(n+1,i,k,oo2iso18) = OO2(i,(kmax+1-k),n,iair18)
#endif

! ### #if ( MEDUSA == 1 )
! ###          OPOCFlux_clio(:,:,k) = TRANSPOSE(OPOC_flux(:,(kmax+1-k),:))
! ### #endif

#if ( CORAL == 1 )
      coral_area_clio(n+1,i,k) = coral_area_out(i,(kmax+1-k),n)
      coral_prod_clio(n+1,i,k) = coral_prod_out(i,(kmax+1-k),n)
      coral_mass_clio(n+1,i,k) = coral_mass_out(i,(kmax+1-k),n)
      omega_arag_clio(n+1,i,k) = omega_arag3D(i,(kmax+1-k),n)
      tau_bleach_clio(n+1,i,k) = tau_bleach_out(i,(kmax+1-k),n)
      DHW_clio(n+1,i,k)= DHW_out(i,(kmax+1-k),n)
      oco3_clio(n+1,i,k) = oco3(i,(kmax+1-k),n)
      PH_clio(n+1,i,k)= PH(i,(kmax+1-k),n)
#endif

#if ( REMIN == 1 )
      kremin_clio(n+1,i,k) = kremin(i,(kmax+1-k),n)
#endif

#if ( OCYCC == 1 )
!         if (k.lt.14) then
         fPOC_flx_clio(n+1,i,k) = TPP_ma(i,(kmax+1-k),n) 
         fCAL_flx_clio(n+1,i,k) = caco3_ma(i,(kmax+1-k),n)
!         endif
#endif

#if ( OOISO == 1 )
      phyto_clio(n+1,i,k) = phyto_day_sum(i,(kmax+1-k),n)
      prodO2_clio(n+1,i,k) = prodO2_sum(i,(kmax+1-k),n)
      respO2_clio(n+1,i,k) = respO2_sum(i,(kmax+1-k),n)
      reminO2_clio(n+1,i,k) = reminO2(i,(kmax+1-k),n) !OO2_flux_isoremin(1)
      NCP_clio(n+1,i,k) = NCP_sum(i,(kmax+1-k),n) 
#endif      

            endif

           enddo
          enddo
         enddo



!         END FORALL
!        END FORALL

!#if ( CORAL == 1 )
!       write(*,*) 'IN SYNC ', NOC_CBR, LT, JT
!       do n=1,NOC_CBR
!        do i=1,LT
!          do j=1,JT
!          if (MGT(i,j,n).eq.1) then
!           !if (coral_prod_out(i,j,n).ne.0) then 
!           if (coral_prod_clio(n,j,i).ne.0) then 
!              !write(*,*) 'in sync_lcm_ocycc', coral_prod_out(i,j,n)
!              write(*,*) 'in sync_lcm_ocycc 2',
!     >           coral_prod_clio(n,j,i)
!         !write(*,*) 'in sync_lcm_ocycc', coral_area_clio(:,:,:)
!           endif
!          endif
!          enddo
!        enddo
!       enddo
!#endif

!#if ( OCYCC == 1 )
!dmr --- Ajout du transfert de la pression partielle du CO2
!        FORALL (n=2:NOC_CBR+1)

        j=1
        do i=1,LT
         do n=2,NOC_CBR+1
          if (MGT(i,j,n-1).eq.1) then

          oxpco2_clio(n,i) = oxpCO2(i,n-1)
c~ !dmr --- Ajout du transfert de la M.O. exportée sous la zone photique
c~           TPPma_clio(n,i) = TPP_mas(i,n-1)
c~           CACO3ma_clio(n,i) = caco3_mas(i,n-1)
c~ !dmr --- Ajout du transfert de la M.O. exportée sous la zone photique
          endif
         enddo
        enddo

!        END FORALL
!#endif

#if ( PATH >= 1 )
        scal(:,:,:,scalstart_path:scalend_path) = Aactivity(:,:,:,:)
#endif

#if ( NEOD > 0 )
         scal(:,:,:,scalstart_neodymium:scalend_neodymium) = neodymium(:,:,:,:)
#endif


!nb tableaux croises pour scal
!        write(*,*) 'TEST nsmax = ', nsmax
        FORALL (n=3:nsmax)
          scal(1,:,:,n) = scal(NOC_CBR+1,:,:,n)
          scal(NOC_CBR+2,:,:,n) = scal(2,:,:,n)
        END FORALL

#if ( CORAL == 1 )
!nb tableaux croises
      coral_area_clio(1,:,:) = coral_area_clio(NOC_CBR+1,:,:)
      coral_area_clio(NOC_CBR+2,:,:) = coral_area_clio(2,:,:)
      coral_prod_clio(1,:,:) = coral_prod_clio(NOC_CBR+1,:,:)
      coral_prod_clio(NOC_CBR+2,:,:) = coral_prod_clio(2,:,:)
      coral_mass_clio(1,:,:) = coral_mass_clio(NOC_CBR+1,:,:)
      coral_mass_clio(NOC_CBR+2,:,:) = coral_mass_clio(2,:,:)
      omega_arag_clio(1,:,:) = omega_arag_clio(NOC_CBR+1,:,:)
      omega_arag_clio(NOC_CBR+2,:,:) = omega_arag_clio(2,:,:)
      tau_bleach_clio(1,:,:) = tau_bleach_clio(NOC_CBR+1,:,:)
      tau_bleach_clio(NOC_CBR+2,:,:) = tau_bleach_clio(2,:,:)
      DHW_clio(1,:,:) = DHW_clio(NOC_CBR+1,:,:)
      DHW_clio(NOC_CBR+2,:,:) = DHW_clio(2,:,:)
      oco3_clio(1,:,:) = oco3_clio(NOC_CBR+1,:,:)
      oco3_clio(NOC_CBR+2,:,:) = oco3_clio(2,:,:)
      PH_clio(1,:,:) = PH_clio(NOC_CBR+1,:,:)
      PH_clio(NOC_CBR+2,:,:) = PH_clio(2,:,:)
#endif

!nb tableaux croises
      oxpco2_clio(1,:) = oxpco2_clio(NOC_CBR+1,:)
      oxpco2_clio(NOC_CBR+2,:) = oxpco2_clio(2,:)
      TPPma_clio(1,:) = TPPma_clio(NOC_CBR+1,:)
      TPPma_clio(NOC_CBR+2,:) = TPPma_clio(2,:)
      CACO3ma_clio(1,:) = CACO3ma_clio(NOC_CBR+1,:)
      CACO3ma_clio(NOC_CBR+2,:) = CACO3ma_clio(2,:)

#if ( REMIN == 1 )
!nb tableau croise pour remineralisation output
      kremin_clio(1,:,:) = kremin_clio(NOC_CBR+1,:,:)
      kremin_clio(NOC_CBR+2,:,:) = kremin_clio(2,:,:)
#endif


       endif

       END SUBROUTINE sync_lcm_ocycc
