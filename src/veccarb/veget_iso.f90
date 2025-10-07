!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:58 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:58 CET 2009
!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Ce module sert à declarer les variables isotopiques du cycle du
!       carbone de vecode (importe de CLIMBER, version 2.4)
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche, Nathaelle Bouttes
!      Date   : 01 decembre 2009
!      Derniere modification : 07 novembre 2012 dmr, vm
!-----|--1--------2---------3---------4---------5---------6---------7-|

       MODULE VEGET_ISO

       use global_constants_mod, only: dblp=>dp, ip        

       IMPLICIT NONE

!-----|--1--------2---------3---------4---------5---------6---------7-|
       INTEGER(kind=ip), PARAMETER :: nnnlat = 32, nnnlon = 64
!-----|--1--------2---------3---------4---------5---------6---------7-|

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Declaration isotopes du carbone
!-----|--1--------2---------3---------4---------5---------6---------7-|
       REAL(kind=dblp), dimension(nnnlat,nnnlon) :: B4T14, B3T14, B4T13,&
         B3T13, B2T14, B1T14, B2T13, B1T13, B4G14, B3G14, B4G13, B3G13, &
         B2G14, B1G14, B2G13, B1G13, BC13=0.0_dblp, BC14, ANUP13       ! [NOTA]: should not BC13 be read in restart somewhere?
       REAL(kind=dblp) :: C13FRAC, C13FRAC4

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Declarations supplementaires pas mise dans "veget.h"
!-----|--1--------2---------3---------4---------5---------6---------7-|
       REAL(kind=dblp), dimension(nnnlat,nnnlon) :: SG4 = 0.0_dblp      &
                      , TATMSMIN
       REAL(kind=dblp) :: G4SHARE_ST
!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Declaration a replace plus tard dans les modules atmosphere
!cnb all moved elsewhere
!-----|--1--------2---------3---------4---------5---------6---------7-|
!       REAL :: C13ATM = 1000.-6.5, C13INIT, C14INIT
!       REAL :: C13ATM, C13INIT, C14INIT
       END MODULE VEGET_ISO
