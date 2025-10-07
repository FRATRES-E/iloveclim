!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Ce module contient les FLAGS en lien avec GRISLI
!       couplage LOVECLIM --> GRISLI
!       (dans l'environnement logiciel LUDUS)
!   
!      Auteur : Didier M. Roche 
!      Date   : 18/Juin/2007 <-- Mask of Precip, Lev Tarasov forcing !
!      Derniere modification : 18 Novembre 2008, Conversion en module
!-----|--1--------2---------3---------4---------5---------6---------7-|

      MODULE input_icemask
! ### dmr ### ci-dessous : sale, mais problemes de compatibilite sinon
      REAL, DIMENSION(32,64), TARGET :: icemask
!###      REAL :: totalprecipcases
! afq -- we duplicate the icemask, for albedo purposes:
!        this is used by ec_co2oc to force a snow albedo is there is
!        an ice sheet over the ocean
      real, dimension(32,64) :: icemaskalb

      END MODULE input_icemask
