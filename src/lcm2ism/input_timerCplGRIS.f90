!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Ce module contient les timers en lien avec GRISLI (frequence
!       couplage, frequences sorties  etc.
!       couplage LOVECLIM --> GRISLI
!       (dans l'environnement logiciel LUDUS)
!   
!      Auteur : Didier M. Roche 
!      Date   : 18 Novembre 2008
!      Derniere modification : 18 Novembre 2008
!-----|--1--------2---------3---------4---------5---------6---------7-|

       MODULE input_timerCplGRIS


       IMPLICIT NONE

! Si dessous la frequence de couplage vue par LOVECLIM
       INTEGER, PARAMETER :: timCplGRISday = 360*1
! Frequence de decouplage : GRISLI integre freqdecoup fois le pas
!  de temps de LOVECLIM. En transitoire vrai, freqdecoup = 1
!  En mode decouple, freqdecoup > 1
       INTEGER, PARAMETER :: freqdecoup = 1
! Si dessous la frequence de couplage vue par GRISLI
       INTEGER, PARAMETER :: timCplGRISyr = INT(timCplGRISday/360)      &
                                             *freqdecoup
! La frequence d'ecriture des fichiers de redemarrage GRISLI
!       INTEGER, PARAMETER :: freqsortie = 25*timCplGRISyr
       INTEGER :: freqsortie 

       END MODULE input_timerCplGRIS
