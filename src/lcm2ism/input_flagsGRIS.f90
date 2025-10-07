!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Ce module contient les FLAGS en lien avec GRISLI
!       couplage LOVECLIM --> GRISLI
!       (dans l'environnement logiciel LUDUS)
!   
!      Auteur : Didier M. Roche 
!      Date   : 18 Novembre 2008
!      Derniere modification : 18 Novembre 2008
!-----|--1--------2---------3---------4---------5---------6---------7-|

       MODULE input_flagsGRIS


       IMPLICIT NONE
! dmr On peut faire une construction judicieuse : 0 = pas de calotte, 
!      1 = calotte variable mais forcee, 2 = calotte interactive
       INTEGER, PARAMETER :: nord_GRIS = 2, sud_GRIS = 0 

       END MODULE input_flagsGRIS
