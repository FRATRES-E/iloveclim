!dmr -- Adding the choice of components through the pre-processing options
#include "choixcomposantes.h"
!dmr -- Adding the choice of components through the pre-processing options
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!   Copyright 2021 Didier M. Roche (a.k.a. dmr)

!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at

!       http://www.apache.org/licenses/LICENSE-2.0

!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Ce module contient les variables commune pour les isotopes de 
!        l'eau. Module commun a priori aux parties océan et atmosphère
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche 
!      Date   : 07 décembre 2010
!      Derniere modification : 23 juin 2021 (version 21)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       MODULE iso_param_mod



       use global_constants_mod, only: ip, dblp=>dp
#if ( WISOATM == 1 )
       use comatm, only: nlat, nlon, nwisos, iwat18, iwat2h
#endif

       implicit none

       private ip, dblp
#if ( WISOATM == 1 )
       private nlat, nlon, nwisos, iwat18, iwat2h
#endif

!      Définitions de compatibilité ...
!~        INTEGER, PARAMETER :: iso_LAT = 32, iso_LON = 64, iso_TYPS = 3
        INTEGER(kind=ip), PARAMETER :: iso_noc = 1_ip, iso_nld = 3_ip, iso_nse = 2_ip

! STD wisos arrangement has been moved to comatm
       
!~ !      Nombre d'isotopes de l'eau : 4 plus l'eau elle même ... = 5
!~        INTEGER, parameter :: neauiso = 5

!~ !      Arrangement des 6 éléments dans les tableaux concernés
!~ !      ======================================================
!~ !      L'eau totale, conservée par compatibilité avec le modèle de base 
!~        INTEGER, parameter :: ieau   = 1
!~ !      Molécule H216O
!~        INTEGER, parameter :: ieau16 = 2
!~ !      Molécule H217O
!~        INTEGER, parameter :: ieau17 = 3
!~ !      Molécule H218O
!~        INTEGER, parameter :: ieau18 = 4
!~ !      Molécule HD16O
!~        INTEGER, parameter :: ieaud = 5

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Definition of the reference abundances
!      cf.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        REAL(KIND=dblp), PARAMETER :: r18smow = 2005.2D-6
        REAL(KIND=dblp), PARAMETER :: r17smow = 379.9D-6
        REAL(KIND=dblp), PARAMETER :: r2hsmow = 155.76D-6

#if ( WISOATM == 1 )
!--- dmr  In practice in version 21, we will use 1.0 for any one of the isotopes
        REAL(KIND=dblp), PARAMETER :: rsmow(nwisos) = [1.0_dblp,1.0_dblp,1.0_dblp,1.0_dblp,1.0_dblp]

!      Molar masses of the Mendeleiv table

        REAL(KIND=dblp), PARAMETER    :: MH = 1.007825032_dblp, M2H = 2.014101778_dblp                                     &
            , MO16 = 15.994914622_dblp, MO17 = 16.9991315_dblp, MO18 = 17.9991604_dblp 

!      Molar masses of the different types of water

        REAL(KIND=dblp), PARAMETER   :: M18 = 2._dblp * MH + MO18, M17 = 2. * MH + MO17, M16 = 2._dblp * MH + MO16         &
                                      , M19 = MH + M2H + MO16, M32 = MO16 + MO16, M33 = MO16 + MO17, M34 = MO16 + MO18

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      For kinetic fractionation, use diffusivity coefficient as defined
!      in reference document. Please note that the coefficient is in the
!      form [Di/D]^n, the coefficients given here are the Di/D
!      Given the formulation used (cf. ref. document), there is no such
!      coefficient for ieau or ieau16. Thus, they are set to 1.0
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!      Values are coming from Merlivat, 1978 & Barkan and Luz, 2007
       REAL(kind=dblp), DIMENSION(nwisos), PARAMETER :: Diff = [1.0d0, 1.0d0, 0.9855d0, 0.9723d0, 0.9755d0]

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      For the special case of the ocean, we want the formulation to be
!      still in the type of [Di/D], though it is fixed in Merlivat and
!      Jouzel, 1979 that epsilon_k = alpha_k + 1 = 6 per mil for 18O.
!      Thus, n_ka_oc is fixed to ca. 0.214
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!      However, it is also fixed in Merlivat and Jouzel, 1979 that
!      epsilon_k = alpha_k + 1 = 6 *0.88 per mil for D. 
       REAL(kind=dblp), DIMENSION(nwisos), PARAMETER, PRIVATE ::                                                                 &
         n_ka_oc = [1.0_dblp,1.0_dblp,LOG(1.0_dblp-6.0E-3_dblp)/LOG(Diff(iwat18)),  LOG(1.0_dblp-6.0E-3_dblp)/LOG(Diff(iwat18)), &
                                      LOG(1.0_dblp-6.0E-3_dblp)/LOG(Diff(iwat18))]

!      From these considerations, it naturally folows that the kinetic
!      fractionation coefficients are, for the ocean:
       REAL(kind=dblp), DIMENSION(nwisos), PARAMETER :: alpha_diff_oc = Diff**n_ka_oc

!      Kinetic fractionnation coefficient for the land surface
!      evaporation

!      First guess
       REAL(kind=dblp), PARAMETER, PRIVATE           :: n_ka_lnd = 0.57_dblp

       REAL(kind=dblp), DIMENSION(nwisos), PARAMETER :: alpha_diff_lnd = Diff**n_ka_lnd

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Moved here the remnants of isoatm_mod ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       REAL(kind=dblp), DIMENSION(nlat, nlon, nwisos) :: ratio_oceanatm
       REAL(kind=dblp), DIMENSION(nlat, nlon, nwisos) :: ratio_evap_ocean, ratio_evap_land, ratio_evap_snow

       REAL(kind=dblp), PARAMETER                     :: fac_17Oexc = 0.528_dblp

       REAL(kind=dblp), PARAMETER :: d18atmini = -12.E-3_dblp,                                            &
                        datmini(nwisos) = [0._dblp,0._dblp,(d18atmini+1._dblp)**fac_17Oexc-1.0_dblp, d18atmini,d18atmini*8._dblp]

       REAL(kind=dblp), PARAMETER :: d18lbmini = -30.E-3_dblp,                                            &
                        dlbmini(nwisos) = [0.d0,0.d0,(d18lbmini+1.0d0)**0.528d0-1.0, d18lbmini,d18lbmini*8._dblp]

       INTEGER(kind=ip),PARAMETER :: isoatm_restart = WISOATM_RESTART
       
! dmr --- [TODO] isolbm_restart should be moved in an other, more logical place!       
       INTEGER(kind=ip),PARAMETER :: isolbm_restart = WISOLND_RESTART


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

          CONTAINS

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Routine pour le calcul du deuterium-excess
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
        pure elemental function dexcess(d,o) result (h2excess)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Variables d'entree  : 
!        d : ratio deuterium
!        o : ratio 18O
!
!       Variables de sortie : 
!        dexcess : le d-excess
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

          IMPLICIT NONE

          REAL(kind=dblp)             :: h2excess
          REAL(kind=dblp), INTENT(IN) :: d,o

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

          h2excess=(d/rsmow(iwat2h)-1.0_dblp)*1000._dblp-8.0_dblp*(o/rsmow(iwat18)-1.0_dblp)*1000._dblp

        end function dexcess

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Routine pour le calcul du delta
!-----|--1--------2---------3---------4---------5---------6---------7-|
        pure elemental function delta(rwat,riso,n) result (delta_val)
!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree  : 
!        ratio : ratio in moisture
!        n : isotope number
!
!       Variables de sortie : 
!        delta : delta
!-----|--1--------2---------3---------4---------5---------6---------7-|

       IMPLICIT NONE

       REAL(kind=dblp)             :: delta_val
       REAL(kind=dblp), INTENT(IN) :: rwat, riso
       integer(kind=ip), intent(in):: n

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Variables locales
!-----|--1--------2---------3---------4---------5---------6---------7-|

       delta_val = 0._dblp

       if (rwat.gt.epsilon(rwat)) then
         delta_val = (riso/rwat/rsmow(n) - 1._dblp) * 1000._dblp
       endif

       end function delta

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Routine pour le calcul du ratio a partir du delta
!-----|--1--------2---------3---------4---------5---------6---------7-|
        pure elemental function delta_inv(rwat,diso,n) result (riso_val)
!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree  : 
!        rwat : amount of water (moisture)
!        diso : delta of the isotope considered
!        n : isotope number
!
!       Variables de sortie : 
!        riso_val: the ratio of the moisture in the isotope considered
!-----|--1--------2---------3---------4---------5---------6---------7-|

       IMPLICIT NONE

       REAL(kind=dblp)             :: riso_val
       REAL(kind=dblp), INTENT(IN) :: rwat, diso
       integer(kind=ip), intent(in):: n

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Variables locales
!-----|--1--------2---------3---------4---------5---------6---------7-|

       riso_val = rwat * rsmow(n) * (diso/1000._dblp + 1._dblp) 

       end function delta_inv
       
!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Routine pour le calcul du delta a partir du ratio
!-----|--1--------2---------3---------4---------5---------6---------7-|
        pure elemental function deltaR(riso,n) result (delta_val)
!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree  : 
!        ratio : ratio in moisture
!        n : isotope number
!
!       Variables de sortie : 
!        delta : delta
!-----|--1--------2---------3---------4---------5---------6---------7-|

       IMPLICIT NONE

       REAL(kind=dblp)             :: delta_val
       REAL(kind=dblp), INTENT(IN) :: riso
       integer(kind=ip), intent(in):: n

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Variables locales
!-----|--1--------2---------3---------4---------5---------6---------7-|

       delta_val = 0._dblp

!~        if (rwat.gt.epsilon(rwat)) then
         delta_val = (riso/rsmow(n) - 1._dblp) * 1000._dblp
!~        endif

       end function deltaR

#endif

       END MODULE iso_param_mod
