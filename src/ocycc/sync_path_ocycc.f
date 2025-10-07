!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:21 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:21 CET 2009
!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Cette routine sync_path_ocycc sert a synchroniser les variables
!       interne du cycle du carbone ocycc et les variables path
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche & Lise Missiaen
!      Date   : 01 Juin 2017
!      Derniere modification : ~
!-----|--1--------2---------3---------4---------5---------6---------7-|

       SUBROUTINE sync_path_ocycc(sens_sync)

#if ( PATH == 2 )

       use global_constants_mod, only: solar_day ! number of seconds per 24 h

       use path_mod, only: particles_fluxes_field, ncaco3, npoc, nlith,
     &                  nopal, kmax => kmax_loc
       use mbiota_mod, only: TPP_ma, caco3_ma

#endif

       INTEGER :: sens_sync

#if ( PATH == 2 )

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Variables locales
!-----|--1--------2---------3---------4---------5---------6---------7-|

       INTEGER :: k, n

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


!-----|--1--------2---------3---------4---------5---------6---------7-|

       if (sens_sync.EQ.1) then ! LCM -> OCYCC

          ! Not used since PaTh does not affect OCYCC (lim&dmr)

       else ! OCYCC -> PaTh

         
       FORALL (k=1:kmax)
         ! --- lim&dmr caco3_path is in molC.m**-2.s-1
         ! --- particule_fluxes_field is in molC.m**-2.s-1
         particles_fluxes_field(2:NOC_CBR,:,k,ncaco3)
     &               = TRANSPOSE(caco3_ma(:,(kmax+1-k),:))
!          FORALL (n=2:NOC_CBR+1)
!         particles_fluxes_field(n,:,k,ncaco3)
!     &               = TRANSPOSE(caco3_path(:,(kmax+1-k),n-1))
         ! --- lim&dmr TPP_path is in molC.m**-2.s-1
         ! --- particule_fluxes_field is in molC.m**-2.s-1
         particles_fluxes_field(2:NOC_CBR,:,k,npoc)
     &               = TRANSPOSE(TPP_ma(:,(kmax+1-k),:))
!         particles_fluxes_field(n,:,k,npoc)
!     &               = TRANSPOSE(TPP_path(:,(kmax+1-k),n-1))
!          END FORALL
!        END FORALL

!nb tableau croise
       particles_fluxes_field(1,:,k,ncaco3)=
     &    particles_fluxes_field(NOC_CBR+1,:,k,ncaco3)
       particles_fluxes_field(NOC_CBR+2,:,k,ncaco3)=
     &    particles_fluxes_field(2,:,k,ncaco3)
       particles_fluxes_field(1,:,k,npoc)=
     &    particles_fluxes_field(NOC_CBR+1,:,k,npoc)
       particles_fluxes_field(NOC_CBR+2,:,k,npoc)=
     &    particles_fluxes_field(2,:,k,npoc)


         END FORALL
!-----lim opal table to initialize (first step to 0)
!  TODO: get more sophisticated opal table to read from netcdf file !

        particles_fluxes_field(:,:,:,nopal)=0.

! unit conversion from [TmolsC m-2 timestep-1] to [mol m-2 s-1]
        particles_fluxes_field(:,:,:,ncaco3)=
     &      particles_fluxes_field(:,:,:,ncaco3)*1E12/86400
        particles_fluxes_field(:,:,:,npoc)=
     &      particles_fluxes_field(:,:,:,npoc)*1E12/86400

       endif

#endif

       END SUBROUTINE sync_path_ocycc
