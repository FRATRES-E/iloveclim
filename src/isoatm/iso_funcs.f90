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


      module iso_funcs

#if ( WISOATM == 1 )

       use global_constants_mod, only: ip, dblp=>dp
       use comatm, only: nlat, nlon, nwisos, iwat18, iwat2h

       implicit none


       private nlat, nlon, nwisos, iwat18, iwat2h, ip, dblp


      public

      contains
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Ce fichier contient les routines de calcul des ratios a partir de la quantité d'eau précipitable isotopique, rmoisg
!
!      Auteur : Didier M. Roche
!      Date   : 15 décembre 2010
!      Derniere modification : 17 décembre 2010, 31 mars 2016, 23 juin 2021 (version 21)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


       pure elemental function rayleigh_dist(fract_wat, ratio_in, alpha) result(rayl_dist)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Variables d'entree  :
!         fract_wat : la fraction d'eau initiale restante
!         ratio_in  : le ration isotopique de la vapeur avant distillation
!         alpha     : le coefficient de fractionnement pour l'isotope considéré
!       Variables de sortie :
!         rayl_dist: le ratio isotopique de la vapeur après distillation
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        implicit none

        real(kind=dblp), intent(in) :: fract_wat
        real(kind=dblp), intent(in) :: ratio_in
        real(kind=dblp), intent(in) :: alpha
        
        real(kind=dblp)             :: rayl_dist

        rayl_dist = ratio_in * (fract_wat)**(alpha - 1.0d0)

       end function rayleigh_dist

!~       SUBROUTINE mois2ratios(moisg,R)

!~ !-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!~ !       Variables d'entree  :
!~ !         rmoisg : la quantité d'eau précipitable normale et isotopique
!~ !       Variables de sortie :
!~ !         R : les ratios isotopiques molaires
!~ !-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!~        USE iso_param_mod, ONLY : iso_LAT, iso_LON, neauiso, ieau18,     &
!~                           ieau17, ieaud, ieau16, M16, M17, M18, MD

!~        IMPLICIT NONE

!~        REAL, DIMENSION(iso_LAT,iso_LON,neauiso), INTENT(IN) :: moisg
!~        REAL, DIMENSION(iso_LAT,iso_LON,neauiso), INTENT(OUT) :: R

!~ !-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!~ !      Variables locales
!~ !-----|--1--------2---------3---------4---------5---------6---------7-|

!~        WRITE(*,*)
!~        WRITE(*,*)
!~        WRITE(*,*) "Cette implementation est fausse !! [STOP]"
!~        WRITE(*,*) "=> a reprendre dans sources/iso_funcs.f90"
!~        WRITE(*,*)
!~        WRITE(*,*)
!~        STOP
!~ !-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!~ !       R18
!~ !-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!~        R(:,:,ieau18) = (moisg(:,:,ieau18)/M18)/((moisg(:,:,ieau16)/M16)+&
!~                                                 (moisg(:,:,ieaud)/MD))

!~ !-----|--1--------2---------3---------4---------5---------6---------7-|
!~ !       R17
!~ !-----|--1--------2---------3---------4---------5---------6---------7-|
!~        R(:,:,ieau17) = (moisg(:,:,ieau17)/M17)/((moisg(:,:,ieau16)/M16)+&
!~                                                 (moisg(:,:,ieaud)/MD))

!~ !-----|--1--------2---------3---------4---------5---------6---------7-|
!~ !       Rd
!~ !-----|--1--------2---------3---------4---------5---------6---------7-|
!~        R(:,:,ieaud) = (moisg(:,:,ieaud)/MD)/(                           &
!~                  2.0d0*((moisg(:,:,ieau16)/M16)+(moisg(:,:,ieau17)/M17) &
!~                        +(moisg(:,:,ieau18)/M18)                         &
!~                        )+ (moisg(:,:,ieaud)/MD)                         &
!~                                             )

!~       END SUBROUTINE mois2ratios

!~ !-----|--1--------2---------3---------4---------5---------6---------7-|
!~ !
!~ !-----|--1--------2---------3---------4---------5---------6---------7-|

!~       SUBROUTINE mois2ratiosMJ(moisg,R)

!~ !-----|--1--------2---------3---------4---------5---------6---------7-|
!~ !       Variables d'entree  :
!~ !         rmoisg : la quantité d'eau précipitable normale et isotopique
!~ !       Variables de sortie :
!~ !         R : les ratios isotopiques molaires
!~ !-----|--1--------2---------3---------4---------5---------6---------7-|

!~        USE iso_param_mod, ONLY : iso_LAT, iso_LON, neauiso, ieau18,     &
!~                           ieau17, ieaud, ieau16, ieau

!~        IMPLICIT NONE

!~        REAL, DIMENSION(iso_LAT,iso_LON,neauiso), INTENT(IN) :: moisg
!~        REAL, DIMENSION(iso_LAT,iso_LON,neauiso), INTENT(OUT) :: R

!~ !-----|--1--------2---------3---------4---------5---------6---------7-|
!~ !      Variables locales
!~ !-----|--1--------2---------3---------4---------5---------6---------7-|

!~        WHERE(moisg(:,:,ieau).NE.0.0D0)
!~ !-----|--1--------2---------3---------4---------5---------6---------7-|
!~ !       R18
!~ !-----|--1--------2---------3---------4---------5---------6---------7-|

!~          R(:,:,ieau18) = moisg(:,:,ieau18) / moisg(:,:,ieau)

!~ !-----|--1--------2---------3---------4---------5---------6---------7-|
!~ !       R17
!~ !-----|--1--------2---------3---------4---------5---------6---------7-|

!~          R(:,:,ieau17) = moisg(:,:,ieau17) / moisg(:,:,ieau)

!~ !-----|--1--------2---------3---------4---------5---------6---------7-|
!~ !       Rd
!~ !-----|--1--------2---------3---------4---------5---------6---------7-|

!~          R(:,:,ieaud) = moisg(:,:,ieaud) / moisg(:,:,ieau)

!~        ELSEWHERE
!~          R(:,:,ieau18) = 0.0d0
!~          R(:,:,ieau17) = 0.0d0
!~          R(:,:,ieaud) = 0.0d0
!~        ENDWHERE

!~       END SUBROUTINE mois2ratiosMJ

!~ !-----|--1--------2---------3---------4---------5---------6---------7-|
!~ !
!~ !-----|--1--------2---------3---------4---------5---------6---------7-|

!~       SUBROUTINE ratios2moisg(R,moisg)

!~ !-----|--1--------2---------3---------4---------5---------6---------7-|
!~ !
!~ !       Version implemeting the formulation of Roche, 2011
!~ !
!~ !       Variables d'entree  :
!~ !         R : les ratios isotopiques molaires
!~ !         rmoisg : la quantité d'eau précipitable normale
!~ !       Variables de sortie :
!~ !         rmoisg : la quantité d'eau précipitable isotopique
!~ !
!~ !      Remarque : dues à la symétrie des formulations dérivées pour la
!~ !       quantité d'eau précipitable et l'humidité spécifique vraie,
!~ !       cette routine est valable dans les deux cas.
!~ !
!~ !      Auteur : Didier M. Roche
!~ !      Date   : 15 décembre 2010
!~ !      Derniere modification : 30 janvier 2011
!~ !-----|--1--------2---------3---------4---------5---------6---------7-|

!~        USE iso_param_mod, ONLY : iso_LAT, iso_LON, neauiso, ieau18,     &
!~                           ieau17, ieaud, ieau16, M16, M17, M18, MD, ieau

!~        IMPLICIT NONE

!~        REAL, DIMENSION(iso_LAT,iso_LON,neauiso), INTENT(INOUT) :: moisg
!~        REAL, DIMENSION(iso_LAT,iso_LON,neauiso), INTENT(IN) :: R

!~ !-----|--1--------2---------3---------4---------5---------6---------7-|
!~ !      Variables locales
!~ !-----|--1--------2---------3---------4---------5---------6---------7-|

!~        INTEGER :: i
!~        REAL, DIMENSION(iso_LAT,iso_LON) :: gamma_r, psi_r, rho_r, m_eau

!~        rho_r(:,:) = 2.0d0*R(:,:,ieaud) / (1.0d0 + R(:,:,ieaud))
!~        gamma_r(:,:) = 1.0d0 + R(:,:,ieau17) + R(:,:,ieau18)
!~        psi_r(:,:) = 1.0d0 / gamma_r(:,:) - rho_r(:,:)

!~        CALL M_H2O(R,m_eau)

!~ !-----|--1--------2---------3---------4---------5---------6---------7-|
!~ !       q18
!~ !-----|--1--------2---------3---------4---------5---------6---------7-|

!~       moisg(:,:,ieau18)=moisg(:,:,ieau)*(R(:,:,ieau18)/gamma_r(:,:))&
!~                         *(M18/m_eau(:,:))


!~ !-----|--1--------2---------3---------4---------5---------6---------7-|
!~ !       q17
!~ !-----|--1--------2---------3---------4---------5---------6---------7-|

!~       moisg(:,:,ieau17)=moisg(:,:,ieau)*(R(:,:,ieau17)/gamma_r(:,:))&
!~                         *(M17/m_eau(:,:))

!~ !-----|--1--------2---------3---------4---------5---------6---------7-|
!~ !       qd
!~ !-----|--1--------2---------3---------4---------5---------6---------7-|

!~        moisg(:,:,ieaud) = moisg(:,:,ieau)*rho_r(:,:)*(MD/m_eau(:,:))

!~ !-----|--1--------2---------3---------4---------5---------6---------7-|
!~ !       q16
!~ !-----|--1--------2---------3---------4---------5---------6---------7-|

!~        moisg(:,:,ieau16) = moisg(:,:,ieau)*psi_r(:,:)*(M16/m_eau(:,:))

!~       END SUBROUTINE ratios2moisg

!~ !-----|--1--------2---------3---------4---------5---------6---------7-|
!~ !
!~ !-----|--1--------2---------3---------4---------5---------6---------7-|

!~       SUBROUTINE ratios2moisgMJ(R,moisg)

!~ !-----|--1--------2---------3---------4---------5---------6---------7-|
!~ !
!~ !       Version implemeting the formulation of Merlivat & Jouzel, 1979
!~ !
!~ !       Variables d'entree  :
!~ !         R : les ratios isotopiques molaires
!~ !         rmoisg : la quantité d'eau précipitable normale
!~ !       Variables de sortie :
!~ !         rmoisg : la quantité d'eau précipitable isotopique
!~ !
!~ !      Remarque : dues à la symétrie des formulations dérivées pour la
!~ !       quantité d'eau précipitable et l'humidité spécifique vraie,
!~ !       cette routine est valable dans les deux cas.
!~ !
!~ !      Auteur : Didier M. Roche
!~ !      Date   : 23 mai 2011
!~ !      Derniere modification : 23 mai 2011
!~ !-----|--1--------2---------3---------4---------5---------6---------7-|

!~        USE iso_param_mod, ONLY : iso_LAT, iso_LON, neauiso, ieau18,     &
!~                           ieau17, ieaud, ieau16, ieau

!~        IMPLICIT NONE

!~        REAL, DIMENSION(iso_LAT,iso_LON,neauiso), INTENT(INOUT) :: moisg
!~        REAL, DIMENSION(iso_LAT,iso_LON,neauiso), INTENT(IN) :: R

!~ !-----|--1--------2---------3---------4---------5---------6---------7-|
!~ !      Variables locales
!~ !-----|--1--------2---------3---------4---------5---------6---------7-|

!~        INTEGER :: i
!~        REAL, DIMENSION(iso_LAT,iso_LON) :: gamma_r, psi_r, rho_r, m_eau

!~ !-----|--1--------2---------3---------4---------5---------6---------7-|
!~ !       q18
!~ !-----|--1--------2---------3---------4---------5---------6---------7-|

!~       moisg(:,:,ieau18)=moisg(:,:,ieau)*R(:,:,ieau18)

!~ !-----|--1--------2---------3---------4---------5---------6---------7-|
!~ !       q17
!~ !-----|--1--------2---------3---------4---------5---------6---------7-|

!~       moisg(:,:,ieau17)=moisg(:,:,ieau)*R(:,:,ieau17)

!~ !-----|--1--------2---------3---------4---------5---------6---------7-|
!~ !       qd
!~ !-----|--1--------2---------3---------4---------5---------6---------7-|

!~        moisg(:,:,ieaud) = moisg(:,:,ieau)*R(:,:,ieaud)

!~       END SUBROUTINE ratios2moisgMJ

!~ !-----|--1--------2---------3---------4---------5---------6---------7-|
!~ !
!~ !-----|--1--------2---------3---------4---------5---------6---------7-|

!~       SUBROUTINE M_H2O(ratio,mhdeuxo)

!~ !-----|--1--------2---------3---------4---------5---------6---------7-|
!~ !       Variables d'entree  :
!~ !         ratio : les ratios isotopiques molaires
!~ !       Variables de sortie :
!~ !         mhdeuxo : la masse molaire de l'eau complete ...
!~ !-----|--1--------2---------3---------4---------5---------6---------7-|

!~        USE iso_param_mod, ONLY : iso_LAT, iso_LON, neauiso, ieau17,     &
!~            ieau18, ieaud, M16, M17, M18, MD

!~        IMPLICIT NONE

!~        REAL, DIMENSION(iso_LAT,iso_LON,neauiso), INTENT (IN) :: ratio
!~        REAL, DIMENSION(iso_LAT,iso_LON), INTENT(OUT) :: mhdeuxo

!~ !-----|--1--------2---------3---------4---------5---------6---------7-|
!~ !      Variables locales
!~ !-----|--1--------2---------3---------4---------5---------6---------7-|

!~        REAL, DIMENSION(iso_LAT,iso_LON) :: gamma_r, psi_r, rho_r

!~        gamma_r(:,:) = 1.0d0 + ratio(:,:,ieau17) + ratio(:,:,ieau18)
!~        rho_r(:,:) = 2.0d0*ratio(:,:,ieaud) / (1.0d0 + ratio(:,:,ieaud))
!~        psi_r(:,:) = 1.0d0 / gamma_r(:,:) - rho_r(:,:)

!~        mhdeuxo(:,:) = (ratio(:,:,ieau18)*M18 + ratio(:,:,ieau17)*M17) / &
!~                       gamma_r(:,:) + psi_r(:,:) * M16 + rho_r(:,:) * MD

!~        END SUBROUTINE M_H2O


#endif 
!-----|--1--------2---------3---------4---------5---------6---------7-|
!
!-----|--1--------2---------3---------4---------5---------6---------7-|
      end module iso_funcs
