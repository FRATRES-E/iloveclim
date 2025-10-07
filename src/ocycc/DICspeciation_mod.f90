!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:21 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:21 CET 2009

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      module DICspeciation_mod
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      subroutine DICspeciation_surface
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       


      use loveclim_transfer_mod, ONLY: TM, SM, MGT, mid_level
      use const_mod, only: gpes, rho0
      use carbonate_speciation_mod, only: incche
      use marine_bio_mod, only: ODIC, OALK, osCO2, oxpCO2, oxCO2, oxHCO3, oxCO3
      use declars_mod, only: LT, NOC_CBR

      implicit none

      real z_h, p_bar
      integer i,j,n


!at the surface
      J=1
      do n=1,NOC_CBR 
       do i=1,LT 

         if (MGT(i,j,n).eq.1) then
            ! pressure in bar=rho * g * z = Kg/m3 * N/kg * m = N/m2 = Pa
            ! 1 bar = 1e5 Pa
            ! rho = b * rho0 /g ?
            ! en attendant use rho0
            ! beware mid_level is negative...
            !write(*,*) 'computation of p_bar', rho0, gpes,mid_level(j+1)
            p_bar= rho0 * gpes * (-1) * mid_level(j+1)*1e-5
            !write(*,*) ' p_bar', p_bar

!nb            call incche(TM(i,j,n),SM(i,j,n),ODIC(i,j,n),

            call incche(TM(i,j,n)+273.15,SM(i,j,n),p_bar,ODIC(i,j,n),   &
                       OALK(i,j,n),osCO2(i,n),oxpCO2(i,n),oxCO2(i,n),oxHCO3(i,n),oxCO3(i,n),z_h)

        
        endif

      enddo
     enddo

!     write(*,*) osCO2(i,n),oxpCO2(i,n),oxCO2(i,n),oxHCO3(i,n),oxCO3(i,n), z_h, p_bar
!     write(*,*) 'DICspeciation', ODIC(20,1,60), oxCO2(20,60)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|       
       end subroutine DICspeciation_surface
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|       
       end module DICspeciation_mod
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
