!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!
!   Copyright 2026 FRATRES-E (https://github.com/FRATRES-E)
!     FRamework for fAst TRansient Earth-system Studies and Education

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

!dmr -- Adding the choice of components through the pre-processing options:

#include "choixcomposantes.h"

      module oceanic_silica_mod

       !! version: v1.0
       !! display: public private protected
       !! proc_internals: true


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: oceanic_silica_mod

!!     @author  Nathan Roquin (nr)
!!     @author  Didier M. Roche  (dmr)
!!     @date Creation date: March, 20th, 2026

!!     @brief This module oceanic_silica_mod is handling the creation of a particuliar type of fried noodles ...

!>
!>     DESCRIPTION : Here add the long_description of the module ...
!>        - Subroutine [NAME] : This subroutine is for blabla
!>        - Formula: $$ \frac{d\lambda}{dt} , \frac{d\phi}{dt},  \frac{dz}{dt} $$
!>     @reference References: papers or other documents to be cited... [Site the website if possible](https://iloveclim.eu)
!>     @date Last modification: SLastChangedDate$
!>     @author Last modified by : tfa, tsa!

!!     REVISION HISTORY:
!!        YYYY-MM-DD - Initial Version
!!        TODO_YYYY-MM-DD - TODO_describe_appropriate_changes to be done or discussed - TODO_name
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       
       use, intrinsic :: iso_fortran_env, only: stdin=>input_unit, stdout=>output_unit, stderr=>error_unit
       use global_constants_mod, only: dblp=>dp, ip
       use declars_mod, only: LT, JT, NOC_CBR
       use loveclim_transfer_mod, only: ZZ

       implicit none
       !private
       public

       ! people -- PUBLIC variables, functions, subroutines.

       !dmr&nr --- Biogenic silica variables are carbon cycle, hence (LT,JT,NOC_CBR)

       REAL(kind=dblp), DIMENSION(LT,JT,NOC_CBR) :: silice   ! traceur silice
       REAL(kind=dblp), DIMENSION(LT,JT,NOC_CBR) :: opal     !traceur opal
       
       REAL(kind=dblp) :: dremopal ! taux de dissolution de silice (à chercher vrai valeur)
       REAL(kind=dblp) :: ropal                      ! ratio Si:P (à chercher valeur à prendre)
       REAL(kind=dblp) :: export                    ! biomasse exporté
       REAL(kind=dblp) :: wopal, wcal                         ! vitesse de chute opal et calcium
       REAL(kind=dblp) :: kremin_carb                         !taux remin car
       REAL(kind=dblp) :: bkopal

       
       !variable à chercher dans iloveclim
       REAL(kind=dblp), DIMENSION(LT,JT,NOC_CBR) :: ptho    ! lambOpal (profil température)
       REAL(kind=dblp) :: ecan                      ! fraction de mortalité réutilisé immédiatement
       REAL(kind=dblp) :: eher                        ! fraction brouté qui est utilisé pour metabo (partie redissoute)
       REAL(kind=dblp) :: phyto_senesc   ! phyto mort naturelle du à l'age (récupérer seulement phyto_senesc)
       REAL(kind=dblp) :: zoo_mortal   !mortalité du zoo
       REAL(kind=dblp) :: zoo_egest    !grazing to POC (ejection direct)
       REAL(kind=dblp), DIMENSION(LT,JT,NOC_CBR) :: caco3_m, caco3_m_b_sh, caco3_dif !concentration de caco3 ajuster
       REAL(kind=dblp), DIMENSION(LT,JT,NOC_CBR) :: SUE_ca_3D      !chute ca dans colonne d'eau (??)
              
      contains


! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  SUBROUTINE PART >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

! **********************************************************************************************************************************
      SUBROUTINE init_oceanic_silica(silice,opal,dremopal,wopal,wcal,bkopal, kremin_carb)
! **********************************************************************************************************************************

!!      AUTHOR : Nathan Roquin (nr)
!!      DESCRIPTION: Subroutine is for initialiser les variable utile pour ce module --> Here add the long_description of the module ...
!!      REFERENCES: papers or other documents to be cited...(including link when possible)
!!      CALL : This subroutine is call in "module_file"

       implicit none

       ! toutes est variables global
       real(kind=dblp), DIMENSION(LT,JT,NOC_CBR), intent(inout)   :: silice     ! variablea [kmol/m3]
       real(kind=dblp), DIMENSION(LT,JT,NOC_CBR), intent(out)   :: opal     ! variableb [kmol/m3]
       real(kind=dblp), intent(out)  :: dremopal          ! variablec [unit] (de 0.01 à 0.0075)
       real(kind=dblp), intent(out)  :: wopal !(m/j * fraction j? )
       real(kind=dblp), intent(out)  :: wcal
       real(kind=dblp), intent(out)  :: bkopal !(cst)
       real(kind=dblp), intent(out)  :: kremin_carb !(cst)


!~        !variables à appeler dans iloveclim
!~        REAL(kind=dblp), DIMENSION(ni,nj,nk),intent(out) :: ptho    ! lambOpal (profil température) (unit?)
!~        REAL(kind=dblp), intent(out) :: ecan                      ! fraction de mortalité réutilisé immédiatement (cst)
!~        REAL(kind=dblp), intent(out) :: eher                        ! fraction brouté qui est utilisé pour metabo (partie redissoute) (cst)
!~        REAL(kind=dblp), intent(out) :: phyto_senesc   ! phyto mort naturelle du à l'age (récupérer seulement phyto_senesc) [unit?]
!~        REAL(kind=dblp), intent(out) :: zoo_mortal   !mortalité du zoo [unit?]
!~        REAL(kind=dblp), intent(out) :: zoo_egest    !grazing to POC (ejection direct) [unit?]
!~        REAL(kind=dblp), DIMENSION(ni,nj,nk), intent(out) :: caco3_m, caco3_m_b_sh, caco3_dif !concentration de caco3 ajuster [kmol/m3]
!~        REAL(kind=dblp), DIMENSION(ni,nj,nk), intent(out) :: SUE_ca_3D      !chute ca dans colonne d'eau (??)


       ! Local variables
       real(kind=dblp) :: sinkspeed_opal, sinkspeed_cal
       real(kind=dblp) :: silica0, silica0surf, opal0
       real(kind=dblp), parameter  :: dtb =1._dblp  ! variablef [unit] [fraction de jours]


       INTEGER :: rc,fu
       CHARACTER(len=18)  :: file_path ="oceanic_silica.nml"




       NAMELIST /silicaInit/silica0, silica0surf, opal0, dremopal, ropal, sinkspeed_opal, sinkspeed_cal, bkopal, kremin_carb
  
   
       INQUIRE (file=file_path, iostat=rc)

       IF (rc /= 0) THEN
                WRITE (stderr, '("Error: input file ", a, " does not exist")') file_path
                STOP
       ENDIF

      ! Open and read Namelist file.
       OPEN (action='read', file=file_path, iostat=rc, newunit=fu)

       IF (rc /= 0) WRITE (stderr, '("Error: invalid Namelist format")')
       READ (nml=silicaInit, iostat=rc, unit=fu)
       
       CLOSE(fu)
      
       silice(:,:,:) = silica0
       silice(:,1,:) = silica0surf
       
       opal(:,:,:)   = opal0
       
       wopal = sinkspeed_opal * dtb   !initialisé car n'évolue pas dans le temps
       wcal  = sinkspeed_cal  * dtb   !initialisé car n'évolue pas dans le temps
      
      END SUBROUTINE init_oceanic_silica
! **********************************************************************************************************************************


! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  FUNCTIONS PART >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!===================================================================================================================================
      SUBROUTINE update_silice_opal_caco3( dremopal, ptho, opal, silice, ropal, export, &
      zoo_mortal, ecan, phyto_senesc, zoo_egest, caco3_m, caco3_m_b_sh, caco3_dif) !result(returnValue)
!===================================================================================================================================

!>      AUTHOR : Nathan Roquin (nr)
!>      DESCRIPTION: Subroutine is for mettre à jour à chaque pas de temps, la quantité de silice, opal et cacO3
!>      --> Here add the long_description of the module
!>      La production de silice et d'opale dépend de la production d'opal (opalrem) mais aussi de la dissolution de l'opal (opalrem)
!>      Opalrem dépend de la quantité d'opal déjà disponible ainsi que de la température ainsi que d'un coeficient de dissolution de silice
!>      delsil dépend de la quantité de silice disponible ou si trop faible, elle suit une loi de Michaelis Mentens 
!>      En deuxième partie nous calculons la production de caco3 en concurence avec la quantité de silice (fac_si)
!>      REFERENCES: papers or other documents to be cited Ilyna et al. 2013(doi:10.1029/2012MS000178)

!>      Input variable :
!>         - inParam1 = blabla
!>         - inParam2 = blabla2
!>      Output variable : GPP_O2.

       implicit none

       !variables globale
       REAL(kind=dblp), intent(in) :: dremopal    ![unit]
       REAL(kind=dblp), DIMENSION(LT,JT,NOC_CBR), intent(in)   :: ptho
       REAL(kind=dblp),DIMENSION(LT,JT,NOC_CBR), intent(inout) :: opal
       REAL(kind=dblp),DIMENSION(LT,JT,NOC_CBR), intent(inout) :: silice
       REAL(kind=dblp), intent(in) :: ropal
       REAL(kind=dblp), intent(inout) :: export
       REAL(kind=dblp), intent(in) :: zoo_mortal
       REAL(kind=dblp), intent(in) :: ecan
       REAL(kind=dblp), intent(in) :: phyto_senesc
       REAL(kind=dblp), intent(in) :: zoo_egest
       REAL(kind=dblp), DIMENSION(LT,JT,NOC_CBR), intent(inout) :: caco3_m
       REAL(kind=dblp), DIMENSION(LT,JT,NOC_CBR), intent(inout) :: caco3_m_b_sh
       REAL(kind=dblp), DIMENSION(LT,JT,NOC_CBR), intent(inout) :: caco3_dif


       !variables locales
       REAL(kind=dblp), DIMENSION(LT,JT,NOC_CBR) :: opalrem     ![unit]
       REAL(kind=dblp) :: fac_si
       REAL(kind=dblp), DIMENSION(LT,JT,NOC_CBR) :: avsil
       REAL(kind=dblp), DIMENSION(LT,JT,NOC_CBR) :: delsil

       INTEGER(kind=ip) :: i,j,k

       export = zoo_mortal * (1._dblp - ecan) + phyto_senesc + zoo_egest

       DO i = 1, LT
         DO k = 1, NOC_CBR
           DO j = 1, JT
             !remi de l'opal
             ! opal a initialisé a 1.e-8_wp (kmol/m3)
             opalrem(i,k,j) = dremopal * 0.1_dblp * (ptho(i,k,j)+3.0_dblp) * opal(i,k,j)

             !silice toujours positif
             avsil(i,k,j) = MAX(0._dblp, silice(i,k,j))

             !Michaelis Mentens
             fac_si = avsil(i,k,j) / (avsil(i,k,j) + bkopal)

             !croissance opale
             delsil(i,k,j) = MIN(ropal * fac_si * export, 0.5_dblp * avsil(i,k,j))

             !mise à jour silice et opal dans la maille
             silice(i,k,j)=silice(i,k,j) - delsil(i,k,j) + opalrem(i,k,j)
             opal(i,k,j)=opal(i,k,j) + delsil(i,k,j) - opalrem(i,k,j)

             !-------------------------------------------------------------------------------------------------------
             !ajout concurence entre caco3 et silice
             !-------------------------------------------------------------------------------------------------------

             !mise à jour de terme caco3 de iloveclim
             caco3_m(i,k,j)      = caco3_m(i,k,j)      * fac_si
             caco3_m_b_sh(i,k,j) = caco3_m_b_sh(i,k,j) * fac_si
             caco3_dif(i,k,j)    = caco3_dif(i,k,j)    * fac_si

             END DO
           END DO
         END DO


      END SUBROUTINE update_silice_opal_caco3

!===================================================================================================================================
      SUBROUTINE sink_opal_caco3(dremopal, ptho, wopal, opal, kremin_carb, SUE_ca_3D, wcal) !result(returnValue)
!===================================================================================================================================

!>      AUTHOR : Nathan Roquin (nr)
!>      DESCRIPTION: Subroutine is for calculer le transport de l'opal et caco3 par chute
!>       --> Here add the long_description of the module ...
!>      Il dépende tous les deux de facteur de reminéralisation et vitesse de chute
!>      REFERENCES: papers or other documents to be cited Ilyna et al. 2013(doi:10.1029/2012MS000178)

      implicit none


      ! variables globales
      REAL(kind=dblp), intent(in) :: dremopal
      REAL(kind=dblp), DIMENSION(LT,JT,NOC_CBR), intent(in) :: ptho
      REAL(kind=dblp), intent(inout) :: wopal
      REAL(kind=dblp), DIMENSION(LT,JT,NOC_CBR), intent(inout) :: opal
      REAL(kind=dblp), intent(inout) :: kremin_carb
      REAL(kind=dblp), DIMENSION(LT,JT,NOC_CBR), intent(inout) :: SUE_ca_3D
      REAL(kind=dblp), intent(out) :: wcal


      ! variable locale
      REAL(kind=dblp), DIMENSION(LT,JT,NOC_CBR) :: k_remin_opal
      INTEGER(kind=ip) :: i,j,k



      DO i = 1, LT
        DO k = 1, NOC_CBR
          DO j = 2, JT
            !conversion en taux de dissolution/remineralisation (lambOpal du papier)
            k_remin_opal(i,k,j) = dremopal * 0.1_dblp * (ptho(i,k,j-1)+3._dblp)

            !transfert de l'opal d'une celule à une autre
            opal(i,k,j) = opal(i,k,j-1) /  (1._dblp + k_remin_opal(i,k,j) * zz(j)/wopal)

            !taux de remin caco3
            kremin_carb = 0.01_dblp

            !transfert caco3 d'une cellule à une autre
            SUE_ca_3D(i,k,j)= SUE_ca_3D(i,k,j-1) / (1._dblp + kremin_carb * zz(j)/wcal)

            !qu'est ce que le i,j,n d'iloveclim ainsi que le JT?

          END DO
        END DO
      END DO


      END SUBROUTINE sink_opal_caco3

!===================================================================================================================================

end module oceanic_silica_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
