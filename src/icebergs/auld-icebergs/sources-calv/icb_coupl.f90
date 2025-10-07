!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Ce module est une construction pour le couplage des icebergs 
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche, Marianne Buegelmayer
!      Date   : 25 Aout 2011
!      Derniere modification : 
!-----|--1--------2---------3---------4---------5---------6---------7-|


       MODULE icb_coupl


       IMPLICIT NONE
 
       INTEGER, PARAMETER, PRIVATE :: my_CLIOlat = 65, my_CLIOlon = 122
       INTEGER, PARAMETER, PRIVATE :: my_mon = 12
       INTEGER, PARAMETER, PRIVATE :: my_icbmax = 500000


       REAL, DIMENSION(my_CLIOlon,my_CLIOlat) :: startbox
       REAL, DIMENSION(my_CLIOlon,my_CLIOlat) :: calvCLIO, mass_icb
!       REAL, DIMENSION(my_CLIOlon,my_CLIOlat) :: clat, clon
       REAL, DIMENSION(my_CLIOlon,my_CLIOlat) :: mass_rest
       REAL, DIMENSION(my_CLIOlon,my_CLIOlat) :: mass_day,runcalvgris
       REAL, DIMENSION(my_CLIOlon,my_CLIOlat,my_mon) :: MonIcb
       REAL, DIMENSION(my_CLIOlon,my_CLIOlat) :: vol_rest
       REAL :: RO = 910.
       REAL :: testmass_day
       REAL :: testwatcons2

       END MODULE icb_coupl
