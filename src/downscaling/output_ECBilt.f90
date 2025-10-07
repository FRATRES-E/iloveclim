!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:35 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:35 CET 2009
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Ce module contient les variables de l'ISM a destination de
!       l'atmosphere ECBilt, cadre du developpement du couplage
!       LOVECLIM --> GRISLI
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche
!      Date   : 21 juillet 2008
!      Derniere modification : 25 janvier 2016, Didier M. Roche
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       module output_ECBilt

       use taillesGrilles, only: iEcb,jEcb

       implicit none

! pour ECBilt, il faut transposer ...
       double precision, dimension(iEcb,jEcb)        :: topoECB, masqueECB, modiffractn, masqueISM2ECB

#if ( DOWNSTS == 1 )
       integer,          parameter                   :: nb_stat = 3
       double precision, dimension(iEcb,jEcb,nb_stat):: qmount_dif, topdifECB

#endif

#if ( F_PALAEO == 1 )

       integer,     parameter                     :: time_max = 3204

       integer,     dimension(iEcb,jEcb)          :: where_update
       integer,     dimension(iEcb,jEcb,time_max) :: masq_trans
       real(kind=8),dimension(iEcb,jEcb,time_max) :: topo_trans
       real(kind=8),dimension(iEcb,jEcb)          :: topoECB_ti
#endif
       end module output_ECBilt
