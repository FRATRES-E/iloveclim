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
!
!      Auteur : Didier M. Roche 
!      Date   : 25 mai 2011
!      Derniere modification : 25 mai 2011, 23 juin 2021 (version 21)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       SUBROUTINE ec_isoatminit


       use global_constants_mod, only: ip, dblp=>dp
       
#if ( WISOATM == 1 )       

       use comsurf_mod, only: fractoc
       use iso_param_mod, only: ratio_oceanatm
       use iso_nml_mod, only: ocean_ratio_mode, ocean_d18o_file, r_ocean_fixed, dexc_ocean_frac, o17exc_ocean_frac,       &
                              iso_ratios
       USE lectNC, only: lect_2D_NC
       use comatm, only: nlat, nlon, iwat18

       IMPLICIT NONE

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Variables locales
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       REAL, DIMENSION(nlon, nlat) :: dataread
       REAL, DIMENSION(nlat, nlon) :: dataecb

       character(len=4), parameter :: varia="d18o"
       INTEGER :: i,j

! dmr&clo   wiso:2c:0001 the isotopologues are transported (wiso_vap_mod); until wiso:2c:0002 to 0004 every water source
!           and sink acts on them at the local vapour ratio (alpha = 1, no source composition).
       write(*,*) "ec_isoatminit: WISOATM = 1, isotopologues transported; sources and sinks at the vapour ratio",      &
                  " (alpha = 1, interim, ROADMAP wiso:2c)"

! dmr&clo   wiso:2b:0001 surface ocean ratios from wisoatm.nml: 'fixed' (uniform) or 'file' (d18O map in permil, converted
!           to a fraction on reading; dD and 17O from dexc_ocean and o17exc_ocean, as in iso_ratios)
       if (trim(ocean_ratio_mode) == 'fixed') then
         do j=1,nlon
           do i=1,nlat
             ratio_oceanatm(i,j,:) = r_ocean_fixed(:)
           enddo
         enddo
         write(*,*) "ec_isoatminit: uniform surface ocean ratios (ocean_ratio_mode = fixed)"
         return
       endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Let's read the input file containing surface water d18O
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       CALL lect_2D_NC(TRIM(ocean_d18o_file),dataread,var_nc = TRIM(varia))

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Reading the input file leads to an array indexed differently as ours. We thus reformat it.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       CALL format_nc_to_ec(dataread,dataecb,nlon,nlat)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       The dataset read has negative missing values that are incompatible with our formulation. Change this.
!       (missing values become very large ratios, which fix_input_field replaces where fractoc /= 0)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       WHERE(dataecb.LT.-1E10_dblp)
         dataecb = ABS(dataecb)
       ENDWHERE

! dmr&clo   wiso:2b:0001 the 18O ratio is filled first and dD, 17O are derived from the filled field: in the legacy order
!           (filling after R17 = R18**0.528) missing cells gave R17 ~ 5e3, below the fill threshold, and stayed wrong
       ratio_oceanatm(:,:,iwat18) = 1.0_dblp + dataecb(:,:) * 1.0e-3_dblp
       call fix_input_field(fractoc,ratio_oceanatm,iwat18,iwat18)

       do j=1,nlon
         do i=1,nlat
           ratio_oceanatm(i,j,:) = iso_ratios(ratio_oceanatm(i,j,iwat18) - 1.0_dblp, dexc_ocean_frac, o17exc_ocean_frac)
         enddo
       enddo


#endif

      contains
      
      
       SUBROUTINE fix_input_field(data_in,data_out,nzi,nzj)

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree  : 
!                            data_in, le champs de reference
!                            nx,ny ; la taille de tableau d'entree
!       Variables de sortie : 
!                            data_out, le champs a corriger
!-----|--1--------2---------3---------4---------5---------6---------7-|

!   INSERER ICI LES EVENTUELS "USE MODULE"

       IMPLICIT NONE

!   INSERER ICI LES DECLARATIONS DES VARIABLES D'ENTREE / SORTIE

       REAL, DIMENSION(:,:),   INTENT(IN)    :: data_in  ! nx, ny
       REAL, DIMENSION(:,:,:), INTENT(INOUT) :: data_out ! nx, ny, nz
 
!~        INTEGER, INTENT(IN) :: nx, ny, nz, nzi, nzj
       INTEGER, INTENT(IN) :: nzi, nzj

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Variables locales
!-----|--1--------2---------3---------4---------5---------6---------7-|

        REAL, PARAMETER :: undef = 10000.0D0
        REAL :: val_moy
        INTEGER :: nb_moy
        INTEGER :: i, j, k, l, ii, jj, kk, ll, m, n, mm, nn, z
        INTEGER :: nx, ny, nz
        INTEGER, PARAMETER :: s = 1
!-----|--1--------2---------3---------4---------5---------6---------7-|
!-----|--1--------2---------3---------4---------5---------6---------7-|

       nx = UBOUND(data_out,DIM=1)
       ny = UBOUND(data_out,DIM=2)
       nz = UBOUND(data_out,DIM=3)
       
       DO z = nzi, nzj

       DO i=1,nx
        DO j=1,ny

         IF ((data_in(i,j).NE.0.0d0).AND.(data_out(i,j,z).GE.undef))    &
         THEN

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       On effectue la moyenne des plus proches avec une valeur
!       pertinente (code repris de routageEAU.f90)
!-----|--1--------2---------3---------4---------5---------6---------7-|

          val_moy = 0.0d0
          nb_moy = 0

          ii = i
          jj = j
          m = -s
          n = s

          DO k=m, n, n
           DO l=m, n, n

              mm = ii+k
              nn = jj+l

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Attention, il ne faut pas etre circulaire en latitude !!! 
!-----|--1--------2---------3---------4---------5---------6---------7-|
              if (mm.LT.1) mm = 1  
              if (mm.GT.nx) mm = nx

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       En revanche il FAUT etre circulaire en longitude !!! 
!-----|--1--------2---------3---------4---------5---------6---------7-|

              if (nn.LT.1) nn = ny + nn
              if (nn.GT.ny) nn = nn - ny

              if (.NOT.(data_out(mm,nn,z).GT.undef)) then
                val_moy = val_moy + data_out(mm,nn,z)
                nb_moy = nb_moy + 1
              endif
           ENDDO
          ENDDO

          IF(nb_moy.NE.0) THEN
            data_out(i,j,z) = val_moy / nb_moy
          ELSE
            WRITE(*,*) "ec_isoatminit: I can't fix location: ", i, j
            stop 1
          ENDIF

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Fin de la partie "fix" du code
!-----|--1--------2---------3---------4---------5---------6---------7-|
         ENDIF

        ENDDO
       ENDDO

       ENDDO

       END SUBROUTINE fix_input_field
      

       END SUBROUTINE ec_isoatminit



! dmr --- The End of All Things (op. cit.)
