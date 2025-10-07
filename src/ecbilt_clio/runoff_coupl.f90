!-----|--1--------2---------3---------4---------5---------6---------7-|
!  Ce module est une construction pour le couplage de GRISLI a ECBILT 
!  This module is used for the coupling between GRISLI and ECBILT 
!  regarding the runoff of GRISLI! (in the LUDUS environment)
!
!      Auteur : Marianne Buegelmayer
!      Date   : 16 Janvier, 2012
!      Derniere modification : 
!-----|--1--------2---------3---------4---------5---------6---------7-|
!mab: runoflECB = runoffland in ECB in m/yr water
!mab: runoECB = runoff oncean in ECB in m/yr ice
!mab: runoECB_vol is used in thersf.f and there made to m3/day (still ice)
!mab: rungris is used in thersf.f and there runoECB is transformed into
!     kg/m2/s so is is then in kg/m2/s

       MODULE runoff_coupl


       IMPLICIT NONE
 
       INTEGER, PARAMETER, PRIVATE :: my_ECBlat = 32, my_ECBlon = 64
       INTEGER, PARAMETER, PRIVATE :: my_CLIOlat = 65, my_CLIOlon = 122

       REAL, DIMENSION(my_ECBlat,my_ECBlon) :: runoflECB
       REAL, DIMENSION(my_CLIOlon,my_CLIOlat) :: runoECB, rungris
       
       END MODULE runoff_coupl
