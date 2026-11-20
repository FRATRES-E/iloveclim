!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Ce module contient les definitions de taille de grille ISM & 
!       atmosphere ECBilt, dans le cadre du developpement du couplage 
!       LOVECLIM --> GRISLI
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche
!      Date   : 21 juillet 2008
!      Derniere modification : 25 mars 2014, Didier M. Roche
!-----|--1--------2---------3---------4---------5---------6---------7-|
!Adding the choice of components through the pre-processing option
#include "choixcomposantes.h"

       MODULE taillesGrilles

       IMPLICIT NONE

!      Definition locale de la taille de la grille ECBilt
       INTEGER, PARAMETER :: iEcb = 32,  jEcb = 64
!      Definition locale du nombre de sous grilles disjointes
       INTEGER, PARAMETER :: sgd = 1

!      Definition locale des sous grilles (nombre de points en x et y de
!      la grille cible
#if ( GRIDS_DOWN == 1)
!Global grid
       INTEGER, PARAMETER,DIMENSION (sgd) :: sgnx = (/720/) 
       INTEGER, PARAMETER,DIMENSION (sgd) :: sgny = (/360/) 
#elif ( GRIDS_DOWN == 2)
!North Hemisphere grid (GRISLI north)
       INTEGER, PARAMETER,DIMENSION (sgd) :: sgnx = (/241/)
       INTEGER, PARAMETER,DIMENSION (sgd) :: sgny = (/241/)
#elif ( GRIDS_DOWN == 3)
!South Hemisphere grid (GRISLI south)
       INTEGER, PARAMETER,DIMENSION (sgd) :: sgnx = (/141/)
       INTEGER, PARAMETER,DIMENSION (sgd) :: sgny = (/141/)
#elif ( GRIDS_DOWN == 4)
!Europe grid       
       INTEGER, PARAMETER,DIMENSION (sgd) :: sgnx = (/256/)
       INTEGER, PARAMETER,DIMENSION (sgd) :: sgny = (/156/) 
!Uncomment to add a specific area
!#elif ( GRIDS_DOWN == 5)       
!       INTEGER, PARAMETER,DIMENSION (sgd) :: sgnx = (/??/)
!       INTEGER, PARAMETER,DIMENSION (sgd) :: sgny = (/??/)

#else
        INTEGER, PARAMETER,DIMENSION (sgd) :: sgnx = (/1/)
        INTEGER, PARAMETER,DIMENSION (sgd) :: sgny = (/1/)


#endif
!      Index max in x and y:
       INTEGER, PARAMETER :: sgnxm = maxval(sgnx)
       INTEGER, PARAMETER :: sgnym = maxval(sgny)
!      Definition locale de la taille de la grille CLIO
       INTEGER, PARAMETER :: CNX = 122, CNY = 65, CNZ = 20
!      Definition locale de la taille de la grille CLIO
       INTEGER, PARAMETER :: CNXI = 144, CNYI = 72
       
! [TODO]  nbmois & nbjours to be DEPRECATED ... cf. global_constants_mod    
   
!      Definition locale du nombre de mois par an
       INTEGER, PARAMETER :: nbmois = 12
!      Definition locale du nombre de jours par an
       INTEGER, PARAMETER :: nbjours = 360

       END MODULE taillesGrilles
