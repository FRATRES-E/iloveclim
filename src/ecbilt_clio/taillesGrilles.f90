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

       MODULE taillesGrilles


       IMPLICIT NONE

!      Definition locale de la taille de la grille ECBilt
       INTEGER, PARAMETER :: iEcb = 32,  jEcb = 64
!      Definition locale du nombre de sous grilles disjointes
       INTEGER, PARAMETER :: sgd = 1
!      Definition locale des sous grilles
       INTEGER, PARAMETER,DIMENSION (sgd) :: sgnx = (/241/) !Hemin40 only
       INTEGER, PARAMETER,DIMENSION (sgd) :: sgny = (/241/) 
       !INTEGER, PARAMETER,DIMENSION (sgd) :: sgnx = (/241 , 141 /) !Hemin40+Ant40
       !INTEGER, PARAMETER,DIMENSION (sgd) :: sgny = (/241 , 141 /)

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
