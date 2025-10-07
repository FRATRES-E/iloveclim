!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!   Copyright 2023 Didier M. Roche (a.k.a. dmr)

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
#include "choixcomposantes.h"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      MODULE check_content

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       USE global_constants_mod, ONLY: dblp=>dp, silp=>sp, sip, big_dp => alt_olympus_mons, ip
       USE declars_mod, only: LT, JT, NOC_CBR

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      IMPLICIT NONE

      PRIVATE
      
      PUBLIC :: update_CCreservoirs
 
        REAL(KIND=dblp), PARAMETER                :: bigN_dp = -1.0*big_dp, d_zero = 0.0_dblp
        INTEGER(kind=ip), PARAMETER               :: length_keep = 100

        REAL(KIND=dblp), DIMENSION(0:length_keep) :: total_alkK = d_zero, total_cccK = d_zero, total_phoK = d_zero
        INTEGER(kind=ip)                          :: idx_keep = 0

      CONTAINS

! ---

      SUBROUTINE update_CCreservoirs(keep,display)

      USE loveclim_transfer_mod, only: MGT, DVOL
      USE marine_bio_mod,        only: OPO4, OALK, ODIC, ODOC, ODOCS, OPOC, OetaC_POMoxid, OetaC_DOMoxid_1D
      !, Oeta, ODIC, ODOC, ODOCS, OPOC
      USE mbiota_mod,            only: PHYTO_M, ZOO_M, SCALE_M     
       
      INTEGER(kind=ip)              :: i,j,n 
      REAL(kind=dblp)               :: tot_alk, tot_phos, global_dic, vtmp, cav_oc
      LOGICAL, OPTIONAL, INTENT(IN) :: keep, display 
       
      do j=1,JT
       do i=1,LT
        do n=1,NOC_CBR
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

          if (MGT(i,j,n).eq.1) then
!~           global_po4=global_po4+OPO4(i,j,n)*DVOL_prev(i,j,n)
!~           global_no3=global_no3+ONO3(i,j,n)*DVOL_prev(i,j,n)
!~           global_alk=global_alk+OALK(i,j,n)*DVOL_prev(i,j,n)
!~           global_oc13=global_oc13+                                      &
!~                       OC13(i,j,n)*DVOL_prev(i,j,n)
!~           global_oc14=global_oc14+                                      &
!~                       OC14(i,j,n)*DVOL_prev(i,j,n)
!~           global_dic=global_dic+ODIC(i,j,n)*DVOL(i,j,n)


          tot_alk=tot_alk+OALK(i,j,n)*DVOL(i,j,n)*1000_dblp
          
          vtmp=((PHYTO_M(i,j,n)+ZOO_M(i,j,n)+ ODOC(i,j,n)+ODOCS(i,j,n)) &
! [UPDATEG]                      *Oeta(j,4)/Oeta(j,5)+OPOC(i,j,n))                 &
                      *OetaC_POMoxid(i,J,n)/OetaC_DOMoxid_1D(j)+OPOC(i,j,n))       &
                      *DVOL(i,j,n)

          tot_phos=tot_phos+vtmp/OetaC_POMoxid(i,J,n)+OPO4(i,j,n)*DVOL(i,j,n)
          cav_oc=cav_oc+(PHYTO_M(i,J,n)+ZOO_M(i,J,n)+ODOC(i,J,n)+OPOC(i,j,n)+ODOCS(i,J,n)+ODIC(i,J,n)*1.e6)*DVOL(i,J,n)     &
                       *12*SCALE_M*1.028

!~           tot_nit_prev=tot_nit_prev+vtmp/Oeta(j,4)*Oeta(j,1)+ONO3(i,j,n)&
!~                       *DVOL_prev(i,j,n)
           endif

          enddo
        enddo
      enddo

      if (PRESENT(keep).AND.keep) then
         idx_keep = idx_keep + 1
         if (idx_keep.LE.length_keep) then
            total_alkK(idx_keep) = tot_alk
            total_cccK(idx_keep) = cav_oc
            total_phoK(idx_keep) = tot_phos
         endif
      endif
      

      if (PRESENT(display).AND.display) then
        if (idx_keep.LE.length_keep) then
         WRITE(*,'(3(A,G12.5))') "update_CCreservoirs | ALK = ", total_alkK(idx_keep), " | ALKdiff "                            &
                                     , total_alkK(idx_keep)-total_alkK(idx_keep-1), " | % "                                     &
                                     , (total_alkK(idx_keep)-total_alkK(idx_keep-1))/total_alkK(idx_keep)
         WRITE(*,'(3(A,G12.5))') "update_CCreservoirs | CCC = ", total_cccK(idx_keep), " | CCCdiff "                            &
                                     , total_cccK(idx_keep)-total_cccK(idx_keep-1), " | % "                                     &
                                     , (total_cccK(idx_keep)-total_cccK(idx_keep-1))/total_cccK(idx_keep)  
         WRITE(*,'(3(A,G12.5))') "update_CCreservoirs | PHO = ", total_phoK(idx_keep), " | PHOdiff "                            &
                                     , total_phoK(idx_keep)-total_phoK(idx_keep-1), " | % "                                     &
                                     , (total_phoK(idx_keep)-total_phoK(idx_keep-1))/total_phoK(idx_keep) 
        endif
      endif
      

      END SUBROUTINE update_CCreservoirs


      END MODULE check_content

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
