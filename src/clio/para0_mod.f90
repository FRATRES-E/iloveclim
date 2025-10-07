!dmr -- Added optional components choice - Tue Jun 14 15:59:05 CEST 2011
#include "choixcomposantes.h"
!dmr -- Added optional components choice - Tue Jun 14 15:59:05 CEST 2011

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Ce module est un portage en FORTRAN 90 du module CLIO initial
!       cree en FORTRAN 77.
!      (dans l'environnement logiciel LUDUS)

!      Auteur : ??, Didier M. Roche
!      Date   : ??
!      Derniere modification : 10 avril 2017
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module para0_mod

      implicit none

      public

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!--parametres lies a la taille du domaine :
! [indispensable pour inclusion des fichiers bloc.com, reper.com, var??.com]
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( HRCLIO == 0 )
      integer, parameter :: imax = 122, jmax = 65,  kmax = 20
#elif ( HRCLIO == 1 )
      integer, parameter :: imax = 242, jmax = 128, kmax = 30
#endif

!nb parameter added for input output
      integer, parameter :: nsmax_TS = 2

#if ( OCYCC == 1 && ISOOCN == 0 && ISOATM == 0 && OOISO== 0)
      integer, parameter :: NISOO2 = 1 !nb

#if ( PATH == 0 && NEOD == 0)
!nb case with only the carbon cycle, no PaTh and no water or O2 isotopes
      integer, parameter :: nsmax = 13

!#elif (PATH == 1 )
#elif (PATH >= 1 && NEOD == 0)
!nb case with carbon cycle and PaTH, but no water or O2 isotopes
      integer, parameter :: nsmax = 17
#elif (PATH >= 1 && NEOD == 1)
      integer, parameter :: nsmax = 19   
#endif

#elif ( OCYCC == 0 && ISOOCN >= 1 && ISOATM >= 2 )
!nb case with no carbon cycle, but with water isotopes, no O2 isotopes (implicit)
!dmr --- Added only three items for the ocean: h217O, h218O, hdo, in that order
!dmr --- Beware for future use: I assume implicitly that these are items 3,4,5 in the arrays
!dmr ---  something that is not a priori straightforward if using the Carbon Cycle in addition
!dmr --- Variables affected are scal & phiss mainly if not only
      integer, parameter :: nsmax = 5 ! [NOTA] if want to have water isotopes after 3,4,5, need to increase that number ...

      integer, parameter :: isoocn_restart=0

      integer, parameter :: owiso = 3, owisostrt=3, owisostop=owisostrt+2
      integer, parameter :: ocnw17=owisostrt, ocnw18=owisostrt+1, ocnw2h=owisostop ! indexation dans les tableaux nsmax

      integer, parameter :: owatert=5

!dmr --- [NOTA] 2018-10-29
!dmr ---        owatert is the number of tracers with water relation T,S,17O,18O,2H
!dmr ---        owatert is *not* the number in nsmax arrays

!dmr ---        nsmax arrays are arranged with anynumber of tracers and the water isotopes are placed
!dmr ---        from owisostrt to owisotop
!dmr ---        indexes to be used are ocnw17, ocnw18, ocnw2h

!dmr ---        conversely, owatert arrays are always 5, since they are mainly water arrays, in link with
!dmr ---        the placement of the water isotopes in the atmosphere
!dmr ---        indexes to be used are the atmospheric ones ... confusing!

#elif ( OCYCC == 1 && ISOOCN >= 1 && ISOATM >= 2 && OOISO == 0 ) /* I assume that if ISOOCN & OCYCC then, PATH is automatically in */
!nb case with carbon cycle and PaTh and water isotopes and O2 isotopes
!dmr --- we are in a case where we have: 17 + 3 = 20 tracers
      integer, parameter :: nsmax = 20

      integer, parameter :: isoocn_restart=0

      integer, parameter :: owiso = 3, owisostrt=18, owisostop=owisostrt+2
      integer, parameter :: ocnw17=owisostrt, ocnw18=owisostrt+1, ocnw2h=owisostop ! indexation dans les tableaux nsmax

      integer, parameter :: owatert=5

      integer, parameter :: NISOO2 = 1 !nb


#elif ( OCYCC == 1 && ISOOCN == 0 && ISOATM == 0 && OOISO == 1 )
!nb with carbon cycle (no PaTH implicit), no water isotopes, with oxygen isotopes in ocean
!nb   2 basic tracers (T and S) + 11 tracers for carbon cycle + 3
!     tracers for oxygen isotopes in ocean = 16
      integer, parameter :: nsmax = 16

      integer, parameter :: isoo2_restart=0 !useful or to be suppressed?

      integer, parameter :: oo2iso = 3, oo2isostrt=14, oo2isostop=oo2isostrt+oo2iso-1
      integer, parameter :: oo2iso16=oo2isostrt, oo2iso17=oo2isostrt+1, oo2iso18=oo2isostop

      integer, parameter ::  NISOO2=4 !test with 1, should be set to 4

#elif ( OCYCC == 1 && PATH == 1 && ISOOCN >= 1 && ISOATM >= 2 && OOISO == 1 ) 
!nb with carbon cycle and PaTh and water isotopes and oxygen isotopes in ocean
!nb   2 basic tracers (T and S) + 11 tracers for carbon cycle + 4 for
!PaTh + 2 for water isotopes + 3 tracers for oxygen isotopes in ocean = 23
      integer, parameter :: nsmax = 23

      integer, parameter :: isoocn_restart=0

      integer, parameter :: owiso = 3, owisostrt=18, owisostop=owisostrt+2
      integer, parameter :: ocnw17=owisostrt, ocnw18=owisostrt+1, ocnw2h=owisostop ! indexation dans les tableaux nsmax

      integer, parameter :: owatert=5

!nb for oxygen isotopes in ocean
      integer, parameter :: isoo2_restart=0 !useful or to be suppressed?

      integer, parameter :: oo2iso = 3, oo2isostrt=21, oo2isostop=oo2isostrt+2
      integer, parameter :: oo2iso16=oo2isostrt, oo2iso17=oo2isostrt+1, oo2iso18=oo2isostop

      integer, parameter :: NISOO2 = 4 !test with 1, should be set to 4

#else

      integer, parameter :: nsmax = 2
#endif

      integer, parameter :: nbpt=imax*jmax
      integer, parameter :: ixjmax = imax*jmax, ijkmax = ixjmax*kmax

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!--Nombre maximum de coins (par niveaux et par type) :
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( HRCLIO == 0 ) /* CL15 */
!      integer, parameter :: ncomax = 100
      integer, parameter :: ncomax = 200
#elif ( HRCLIO == 1 ) /* CL30 */
      integer, parameter :: ncomax = 400
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!--Nombre maximum d'arretes (suivant X, Y, et pour les kmax Niv.)
!      (ordre de grandeur : imax*jmax)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
#if ( HRCLIO == 0 ) /* CL15 */
      integer, parameter :: nlpmax = 6000
#elif ( HRCLIO == 1 ) /* CL30 */
      integer, parameter :: nlpmax = 30000
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!--Nombre maximum de points de grille avec Rap.Inter. (depend de "forc.corr"):
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      integer, parameter :: nrpmax = 10

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!--determine le type de bassin : 0=ferme , 1=cyclique , 2=grilles raccordees
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( HRCLIO == 0 ) /* CL15 */
      integer, parameter :: ltest = 3, jsepar = 50
#elif ( HRCLIO == 1 ) /* CL30 */
      integer, parameter :: ltest = 3, jsepar = 98
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!--parametres lies a la frequence des donnees pour T , S , tau.x et tau.y
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      integer, parameter :: nmois = 12, nseas = 4

#if ( ISOOCN >= 1 )
      contains
      
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!-- Following subroutine returns the atmospheric index of an isotope if the oceanic one is known
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      
      subroutine oc2atwisoindx(ocnwiso, atmwiso)
      
        use iso_param_mod, only: ieau17, ieau18, ieaud
       
        integer, intent(in) :: ocnwiso
        integer, intent(out):: atmwiso
        
        select case (ocnwiso)
           case (ocnw17)
              atmwiso = ieau17
           case (ocnw18)
              atmwiso = ieau18
           case (ocnw2h)
              atmwiso = ieaud
        end select
        
        return
      
      end subroutine oc2atwisoindx
#endif

      end module para0_mod
