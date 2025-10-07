!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:35 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:35 CET 2009
!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Ce sous-programme est uns substitution LUDUS à ec_runoff inclus
!       initialement dans landmodel.f
!
!      Au passage, il a été simplifié par la suppression des inclusifs
!       pour l'ISM == 1 de P. Huybrechts
!
!      Auteur : Didier M. Roche
!      Date   : 04 décembre 2010
!      Derniere modification : Didier M. Roche, 05 juillet 2012
!-----|--1--------2---------3---------4---------5---------6---------7-|

      SUBROUTINE ec_runoff_2(istep)

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree :
!       Variables de sorties :
!-----|--1--------2---------3---------4---------5---------6---------7-|

!-----|--1--------2---------3---------4---------5---------6---------7-|
!
!-----|--1--------2---------3---------4---------5---------6---------7-|

      use comemic_mod, only: nstpyear,fracto
      use comunit,     only:

      use global_constants_mod, only: ip, dp
      USE routageEAU_mod, ONLY: eni, enj, flux_per_basin, river_basin

c~ #if ( ISOATM >= 2 )
c~       USE iso_param_mod, ONLY : ieau, neauiso, ieau18, rsmow, delta
c~ #endif

      use comland_mod, only: runofo, fractl, dareas, runofl, bmoism, bmoisg
     &                     , epsl, dtland
#if ( ICEBERG == 2 && ISM !=2 )
     &                     , sntoicebl, sntoicebo,iceb_toto
      use input_icemask, only: icemask
      use iceberg_mod, only: iceberg_info_out_id         
#endif

      use comatm, only: nlat, nlon, iwater

      IMPLICIT NONE

!-----|--1--------2---------3---------4---------5---------6---------7-|
!
!-----|--1--------2---------3---------4---------5---------6---------7-|


      INTEGER(kind=ip), INTENT(IN) :: istep

!-----|--1--------2---------3---------4---------5---------6---------7-|
!        Variables locales
!-----|--1--------2---------3---------4---------5---------6---------7-|
      integer(kind=ip) ::   i,j

      integer(kind=ip) :: runoff_per_basin_txt_id
      
c~ #if ( ISOATM >= 2 )
c~       INTEGER(kind=ip) :: k
c~       real(kind=dp), dimension(neauiso) :: variso
c~       character*10 varisonm
c~ #endif

!#if ( CALVFLUX >= 1)
!      real :: mass_grisrunl(nlat,nlon)
!      integer :: irunlabelgris
!      character*6 :: frunl
!      logical :: mexist
cdmr --- factor for converting the precipitation from GRISLI
cdmr ---  (meters3/year) to the current runofl format (m.s-1)
cmab --- will be divided by the area in the code
!      real, parameter :: conv_fact = 1.d0/(360.d0*24.d0*3600.d0)
!      REAL(KIND=8) :: sumprecipnegat, watconstest
!#endif

!-----|--1--------2---------3---------4---------5---------6---------7-|
!        Optional stuff
!-----|--1--------2---------3---------4---------5---------6---------7-|
#define DIAGNOSIS_FLAG 1
#if ( DIAGNOSIS_FLAG == 1 )
      REAL(kind=dp) :: totallnd, totaloc
#endif

!-----|--1--------2---------3---------4---------5---------6---------7-|
!        Inits ...
!-----|--1--------2---------3---------4---------5---------6---------7-|

#if ( DIAGNOSIS_FLAG == 1 )
      totallnd = 0.0d0
      runofo = 0.0d0
#endif

#if ( ICEBERG == 2 && ISM != 2)
         sntoicebo = 0.
         iceb_toto = 0.
#endif

!-----|--1--------2---------3---------4---------5---------6---------7-|
! dmr    Main loop
!        We compute the runoff in each given i,j grid point
!-----|--1--------2---------3---------4---------5---------6---------7-|
      do i=1,nlat
        do j=1,nlon

c~ #if ( ISOATM >= 2 )
          runofl(i,j,:)=0.
c~ #else
c~           runofl(i,j)=0.
c~ #endif

! mab: epsl = 1e-10  
     
          if (fractl(i,j).gt.epsl) then

            IF ((eni(i,j).EQ.0).OR.(enj(i,j).EQ.0)) THEN
              WRITE(*,*) "Probleme en : ", i,j, "!!", eni(i,j), enj(i,j)
              READ(*,*)
            ENDIF

c~ #if ( ISOATM >= 2 )
!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Computation of isotopic runoff
!     == No isotopic fractionation is considered here ==
!-----|--1--------2---------3---------4---------5---------6---------7-|
!mab: bmoism = maximum bottom moisture
c~           if (bmoisg(i,j,ieau).gt.bmoism(i,j)) then

c~ #if ( ISOATM == 2 )
c~            DO k=ieau+1,neauiso
c~              bmoisg(i,j,k) = bmoisg(i,j,k) -
c~      &((bmoisg(i,j,ieau)-bmoism(i,j))/bmoisg(i,j,ieau))*bmoisg(i,j,k)
c~              runofl(i,j,k)=(
c~      &((bmoisg(i,j,ieau)-bmoism(i,j))/bmoisg(i,j,ieau))*bmoisg(i,j,k)
c~      &)/dtland

c~            ENDDO
c~ #elif ( ISOATM == 3)
c~            WRITE(*,*) "OPTION non implementee, ec_runoff2"
c~ #endif
c~               runofl(i,j,ieau)=(bmoisg(i,j,ieau)-bmoism(i,j))/dtland
c~               bmoisg(i,j,ieau)=bmoism(i,j)

c~           endif ! (bmoisg(i,j,ieau).gt.bmoism(i,j))

c~ #else /* ISOATM < 2  thus no isotopes in runoff */

! mab/dmr : dtland = 3600*24/(iatm*ilan) = 14
! dmr: the unit of bmois is m
! dmr: runofl is hence there is m.s-1.step-1 (iland)
! dmr: input is thus consistent in using m.s-1
     
         if (bmoisg(i,j,iwater).gt.bmoism(i,j)) then
              runofl(i,j,iwater)=
     >            (bmoisg(i,j,iwater)-bmoism(i,j))/dtland
              bmoisg(i,j,iwater)=bmoism(i,j)
         endif

c~ #endif
#if ( DIAGNOSIS_FLAG == 1 )
c~ #if ( ISOATM >= 2 )
c~           totallnd = totallnd + runofl(i,j,ieau)*dareas(i)*fractl(i,j)
c~ #else
          totallnd = totallnd + runofl(i,j,iwater)*dareas(i)*fractl(i,j)
c~ #endif
#endif

         endif ! (fractl(i,j).gt.epsl)


        enddo ! on j, nlon
      enddo ! on i, nlat

!      if(istp .eq. 1) print*,'ec_runoff, sumprecip: ',sumprecipnegat
!      if(iyear .gt. 1) print*,'iyear!: ',iyear
!      if(iyear .gt. 1) print*,'ec_runoff, sumprecip: ',sumprecipnegat

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Redistribution sur les cases de sortie ad hoc
!        On re-décrit les cases ECBilt pour trouver leurs sorties ...
!-----|--1--------2---------3---------4---------5---------6---------7-|


!!! [DEPRECATED] 2019-09-27
!!!### #if ( WATGRISCONS == 1 )
!!!###        irunlabelgris = irunlabel + iyear-1
!!!###       if (istp .eq. 1) then
!!!###          write(frunl,'(I6.6)') irunlabelgris
!!!###          inquire(file='startdata/mass_grisrunECB'//frunl//'.dat',
!!!###      &   exist=mexist)
!!!###          if(mexist) then
!!!###           write(*,*)'MEXIST IN RUNOFFLAND: ',mexist,frunl,irunlabelgris
!!!###            open(1050,file='startdata/mass_grisrunECB'//frunl//'.dat',
!!!###      &     form='unformatted')
!!!###            read(1050) mass_grisrunl
!!!###          else
!!!###           write(*,*)'MEXIST IN RUNOFFLAND: ',mexist,frunl,irunlabelgris
!!!###            mass_grisrunl(:,:)=0d0
!!!###          endif
!!!###          close(1050)
!!!###        runoflECB(:,:) = runoflECB(:,:) + mass_grisrunl(:,:)
!!!###        endif
!!!### #endif
!!! [DEPRECATED] 2019-09-27

      do i=1,nlat
        do j=1,nlon
          if ((fractl(i,j).gt.epsl)) then
c~ #if ( ISOATM >= 2 )
c~             DO k=ieau,neauiso
c~               runofo(eni(i,j),enj(i,j),k)=runofo(eni(i,j),enj(i,j),k)
c~      &        + runofl(i,j,k) * fractl(i,j)/(1.0-fractl(eni(i,j),enj(i,j)))
c~      &        * dareas(i)/dareas(eni(i,j))
c~             ENDDO
c~ #else /* ISOATM < 2  thus no isotopes in runoff */

            runofo(eni(i,j),enj(i,j),iwater)=runofo(eni(i,j),enj(i,j),iwater)
     &      + runofl(i,j,iwater) * fractl(i,j)/(1.0-fractl(eni(i,j),enj(i,j)))
     &      *dareas(i)/dareas(eni(i,j))

            flux_per_basin(river_basin(i,j))=flux_per_basin(river_basin(i,j))
     &      + runofl(i,j,iwater) * fractl(i,j) * dareas(i)

#if ( ICEBERG == 2 && ISM != 2 )
            if (mod(istep,nstpyear).eq.0) then ! only do conversion at the last timestep of the year
              sntoicebo(eni(i,j),enj(i,j))= sntoicebo(eni(i,j),enj(i,j))
     &        +sntoicebl(i,j) * fractl(i,j)/(1.0-fractl(eni(i,j),enj(i,j)))
     &        *dareas(i)/dareas(eni(i,j))
            endif
#endif

c~ #endif /* On ISOATM <> 2 */



          endif ! on fractl > epsl
        enddo ! on j, nlon
      enddo ! on i, nlat

#if ( ICEBERG == 2 && ISM != 2 )
      if (mod(istep,nstpyear).eq.0) then ! last timestep of the year
        do i=1,nlat
          do j=1,nlon
            iceb_toto = iceb_toto + sntoicebo(i,j)*dareas(i)*fracto(i,j)
          end do
        end do
        WRITE(iceberg_info_out_id,*) "yearly sum iceberg flux ocean-grid [m3]"
     &                               ,iceb_toto,istep
      endif          
#endif

#if ( DIAGNOSIS_FLAG == 1 )
!-----|--1--------2---------3---------4---------5---------6---------7-|
! dmr   Summarising water flux to the ocean cells
!       Warning: do not take into account the flux from mab in
!         RUNOFFGRIS since conditions is on fractl
! dmr   Practically, no impact since totaloc is not used anywhere
! dmr   Edit: Use this code only if DIAGNOSIS_FLAG set to 1
!-----|--1--------2---------3---------4---------5---------6---------7-|
      totaloc = 0
      do i=1,nlat
        do j=1,nlon
        if ((1.0-fractl(i,j)).gt.epsl) then
c~ #if ( ISOATM >= 2 )
c~         totaloc = totaloc + runofo(i,j,ieau)*dareas(i)*(1.0-fractl(i,j))
c~ #else
        totaloc = totaloc + runofo(i,j,iwater)*dareas(i)*(1.0-fractl(i,j))
c~ #endif
        endif
        enddo
      enddo
#endif /* On DIAGNOSIS_FLAG */

!-----|--1--------2---------3---------4---------5---------6---------7-|
! dmr  Output the ruoff per defined land basins (initially for fjd)
!-----|--1--------2---------3---------4---------5---------6---------7-|

      IF (mod(istep,nstpyear).EQ.0) THEN ! end of year
        OPEN(newunit=runoff_per_basin_txt_id,
     &         file='outputdata/land/runoff_per_basin.txt',
     &         form='formatted', position='append', status='unknown')
!dmr&fjd --- Output dimension is Km3/year
        WRITE(runoff_per_basin_txt_id,'(13E20.10)') flux_per_basin/nstpyear*31104000.0E0/1.0E9
        CLOSE(runoff_per_basin_txt_id)
        flux_per_basin = 0.0
      ENDIF

!-----|--1--------2---------3---------4---------5---------6---------7-|
! dmr   Finished business with runoff!
!-----|--1--------2---------3---------4---------5---------6---------7-|

      return

      END SUBROUTINE ec_runoff_2
