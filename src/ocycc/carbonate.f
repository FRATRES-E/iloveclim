c********************************************************************
      SUBROUTINE CARBONATE
c********************************************************************
c Compare CO3-- a CO3--sat
c
c By N. Bouttes
c 2008
c Jansen, GBC, 2002
c********************************************************************
      USE declars_mod
!      USE params_mod
!      USE ocean_mod
      USE marine_bio_mod
      USE carbonate_mod !CO3sat_ar
      USE loveclim_transfer_mod, only: MGT, TM, SM, mid_level
      use carbonate_speciation_mod, only: incche
      use const_mod, only: gpes, rho0

      implicit none

      INTEGER i,n,j

      REAL sCO2
      REAL xpCO2
      REAL xCO2
      REAL xHCO3
      REAL xCO3
      REAL p_bar
      REAL z_h


      do j=1,JT
        do n=1,NOC_CBR
          do i=1,LT

            p_bar= rho0 * gpes * (-1) * mid_level(j+1)*1e-5

            if (MGT(i,j,n).eq.1) then
!nb               call incche(TM(i,j,n),SM(i,j,n),ODIC(i,j,n),
               call incche(TM(i,j,n)+273.15,SM(i,j,n),p_bar, ODIC(i,j,n),
!nb     >                  OALK(i,j,n),sCO2,xpCO2,xCO2,xHCO3,xCO3)
     >                  OALK(i,j,n),sCO2,xpCO2,xCO2,xHCO3,xCO3, z_h)

c              aragonite
               if ((xCO3*1000000).ge.CO3sat_ar(j)) then !supersaturated
                 sat_ar(i,j,n)=0
c                 print*,sat(i,j,n),xCO3*1000000,CO3sat_ar(j)
               else !undersaturated
                 sat_ar(i,j,n)=1
c                 print*,sat(i,j,n),xCO3*1000000,CO3sat_ar(j)
               endif

c              calcite
               if ((xCO3*1000000).ge.CO3sat_ca(j)) then !supersaturated
                 sat_ca(i,j,n)=0
c                 print*,sat(i,j,n),xCO3*1000000,CO3sat_ar(j)
               else !undersaturated
                 sat_ca(i,j,n)=1
c                 print*,sat(i,j,n),xCO3*1000000,CO3sat_ar(j)
               endif


            else !mgtne1 pas ocean
              sat_ar(i,j,n)=-1
              sat_ca(i,j,n)=-1

            endif

          enddo
        enddo
      enddo

      END SUBROUTINE CARBONATE

