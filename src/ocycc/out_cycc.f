!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:58 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:58 CET 2009
!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Cette routine <out_cycc> sert a ecrire les sorties specifiques
!       et globales au cycle du carbone
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Nathaelle Bouttes, Didier M. Roche 
!      Date   : 08 Mars 2010
!      Derniere modification : 08 Mars 2010
!-----|--1--------2---------3---------4---------5---------6---------7-|

       SUBROUTINE out_cycc(KOD,fracgr,darea)

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree  : 
!       Variables de sortie : 
!-----|--1--------2---------3---------4---------5---------6---------7-|

       USE declars_mod

       use file_libs, only: open_f, close_f, fileDescriptor
       use global_constants_mod, only: str_len
       use C_res_mod, only: c13_res_fich, ca13_oc_ini, ca13_la_ini
     >              , dc13at_ini, cav_oc13,cav_la13, c13atm
       USE veget_iso, ONLY: B4T13, B3T13, B2T13, B1T13, B4G13, B3G13
     >             , B2G13, B1G13

       USE loveclim_transfer_mod, only: MGT, DVOL
       USE mod_sync_time, only: NYR
       USE C_res_mod, ONLY: c_res_fich, cav_oc, cav_la, ca_oc_ini
     >              , ca_la_ini, c_ocean_fich, fPOC_fich, fCAL_fich
#if ( KC14 == 1 )
       USE C_res_mod, ONLY: c14_res_fich, c_flux_fich, dc14_fich
     >              , c14_bilan_fich, FC12OA,FC12LA,FC14OA, FC14LA
     >              , cav_oc14, cav_la14
      USE marine_bio_mod, ONLY: FDIC, OC14, FOC14, FOAC14
      USE mbiota_mod, ONLY: SCALE_M
      USE veget_iso, ONLY: B4T14, B3T14, B2T14, B1T14, B4G14, B3G14
     >             , B2G14, B1G14


      USE marine_bio_mod, ONLY: ODIC
      USE carbone_co2, ONLY: C14ATM, C14RSTD, MPRODC14, LIMC14MASSE
     & , N14LIM
#endif
       USE carbone_co2, ONLY: PA_C, PA0_C
#if ( IMSK == 1 )
       USE input_icemask
#endif

#if ( COMATM == 1 )
       use veget_mod
       use comatm, only: nlat, nlon
#endif

#if ( CORAL == 1 )
       use C_res_mod, ONLY: coral_res_fich, c_nino_fich
#endif

       IMPLICIT NONE

#if ( COMATM == 0 )
#include "veget.h"
#endif

!dmr&nb Declaration des parametres de la SUBROUTINE
       INTEGER KOD
       REAL*8 fracgr(nlat,nlon)
       REAL*8 darea(nlat)

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Variables locales
!-----|--1--------2---------3---------4---------5---------6---------7-|

       INTEGER :: fich_num
       INTEGER :: k, i, j, n
       LOGICAL :: existe
!nb - Total carbon (GtC) 
       real carbone_total, carbone_total_ini
!nb - Vegetation carbon reservoirs
       REAL :: b1_tot,b2_tot,b3_tot,b4_tot!,npp_tot
       real :: b1_13tot,b2_13tot,b3_13tot,b4_13tot
!nb - Ocean carbon reservoirs
       real phyto_tot, zoo_tot
       real odoc_tot, odocs_tot
       real odic_tot, caco3_tot
!       real cav_oc_detail(NOC_CBR)

#if ( KC14 == 1 )
cvm - c14 ocean reservoir
        real oc14_tot
cvm - c14 normalise
        real dc14atm, dc14ocn, dc14lnd
        real n14atm, n14ocn, n14lnd, n14tot
        real n12tot, r14tot, alpha14
        real c14_total, carbone14_total, mc14atm   !vm total c14 carbon
cvm - c14 vegetation reservoir
        real b1_14tot,b2_14tot,b3_14tot,b4_14tot,lc14_tot
#endif

       real dum_toprint_un, dum_toprint_de

!-----|--1--------2---------3---------4---------5---------6---------7-|
!
!-----|--1--------2---------3---------4---------5---------6---------7-|
    
      IF (KOD.EQ.-1) THEN     
      
c~       fich_num=274

c~       print*, "c_res_fich", c_res_fich, fich_num

c~       if (c_res_fich.EQ.fich_num) then
c~         existe=.TRUE.

c~         DO WHILE (existe)
c~          INQUIRE(fich_num,OPENED=existe)
c~          print*, "c_res_fich", c_res_fich, fich_num, existe
c~          fich_num=fich_num+1 
c~         END DO 



c~         c_res_fich=fich_num

        OPEN(newunit=c_res_fich
     &   ,file='outputdata/carbon/C_reservoirs.txt',status='unknown')
     
        WRITE(c_res_fich,'((A6,16A10))') "Annee", "PA_C", "PA0_C", 
     &   "C13ATM",
     &   "OC_C", "OC_C_i", "VE_C", "VE_C_i", "TOT_C", "TOT_C_i",
     &   "b1_tot", "b2_tot", "b3_tot", "b4_tot"

c~       endif


      c13_res_fich%isFormatted = .true.
      call open_f(c13_res_fich,'outputdata/carbon/C13_reservoirs.txt'
     & , o_stat="unknown")

      WRITE(c13_res_fich%id,'((A6,10A10))') "Annee", "OC13_i", "OC13_c",
     &   "LA13_i", "LA13_c", "AT13_i", "AT13_c", "B1T13", "B2T13", "B3T13"
     & , "B4T13"


      c_ocean_fich%isFormatted = .true.
      call open_f(c_ocean_fich,'outputdata/carbon/C_ocean.txt'
     &  , o_stat="unknown")

      fPOC_fich%isFormatted = .true.
      call open_f(fPOC_fich,'outputdata/carbon/fPOC_gen.txt'
     &  , o_stat="unknown")

      fCAL_fich%isFormatted = .true.
      call open_f(fCAL_fich,'outputdata/carbon/fCAL_gen.txt'
     &  , o_stat="unknown")


#if ( KC14 == 1 )
c~       c14_res_fich = fich_num
c~       if (c14_res_fich.EQ.fich_num) then
c~         existe=.TRUE.

c~         DO WHILE (existe)
c~          INQUIRE(fich_num,OPENED=existe)
c~          fich_num=fich_num+1 
c~         END DO 

c~         c14_res_fich = fich_num
cvm - C14 reservoirs 
       open (newunit=c14_res_fich
     > ,file='outputdata/carbon/C14_reservoirs.txt', status='unknown')
c~       endif

c~       c_flux_fich = fich_num
c~       if (c_flux_fich.EQ.fich_num) then
c~         existe=.TRUE.

c~         DO WHILE (existe)
c~          INQUIRE(fich_num,OPENED=existe)
c~          fich_num=fich_num+1 
c~         END DO 

c~         c_flux_fich = fich_num

cvm - Flux de carbone entre les reservoirs
       open (newunit=c_flux_fich
     > ,file='outputdata/carbon/C_flux.txt', status='unknown')
c~       endif

c~       dc14_fich = fich_num
c~       if (dc14_fich.EQ.fich_num) then
c~         existe=.TRUE.

c~         DO WHILE (existe)
c~          INQUIRE(fich_num,OPENED=existe)
c~          fich_num=fich_num+1 
c~         END DO 

c~         dc14_fich = fich_num


cvm - DC14
       open (newunit=dc14_fich
     > ,file='outputdata/carbon/DC14.txt', status='unknown')
c~       endif

c~       c14_bilan_fich = fich_num
c~       if (c14_bilan_fich.EQ.fich_num) then
c~         existe=.TRUE.

c~         DO WHILE (existe)
c~          INQUIRE(fich_num,OPENED=existe)
c~          fich_num=fich_num+1 
c~         END DO 

c~         c14_bilan_fich = fich_num



cvm - C14 bilan d atomes
       open (newunit=c14_bilan_fich
     > ,file='outputdata/carbon/C14_bilan.txt', status='unknown')
c~       endif
#endif


#if ( CORAL == 1 )
!nb corals
c~       coral_res_fich = fich_num
c~       if (coral_res_fich.EQ.fich_num) then
c~         existe=.TRUE.

c~         DO WHILE (existe)
c~          INQUIRE(fich_num,OPENED=existe)
c~          fich_num=fich_num+1 
c~         END DO 

c~         coral_res_fich = fich_num
 
       open (newunit=coral_res_fich
     > ,file='outputdata/carbon/Coral_output.txt', status='unknown')

       write (coral_res_fich,'(A6,3A10)')
     >    'NYR  ', ' area_coral ', ' prod_coral ',
     >    ' mass_coral '

c~       endif

       open (newunit=c_nino_fich
     > ,file='outputdata/carbon/C_nino.txt', status='unknown')

       write (c_nino_fich,'(A6,2A10)')
     >    'NYR  ', ' nino3 ', ' nino3_var '


#endif

      RETURN

      ENDIF

!-----|--1--------2---------3---------4---------5---------6---------7-|
!
!-----|--1--------2---------3---------4---------5---------6---------7-|


cnb - Total carbon = ocean carbon + land carbon + atmospheric carbon
       carbone_total = cav_oc + cav_la + PA_C / 0.47            !PgC
       carbone_total_ini = ca_oc_ini + ca_la_ini + PA0_C / 0.47 !PgC
#if ( KC14 == 1 )
       mc14atm = c14atm / 0.40 * 1e15                           !gC14
       print*,'c14atm in out.f', C14ATM
       carbone14_total = cav_oc14 + cav_la14 + mc14atm          !gC14
#endif
cnb - Land reservoir details
        b1_tot=0
        b2_tot=0
        b3_tot=0
        b4_tot=0

        b1_13tot=0
        b2_13tot=0
        b3_13tot=0
        b4_13tot=0
#if ( KC14 == 1 )
cvm - C14 land reservoirs
        b1_14tot=0
        b2_14tot=0
        b3_14tot=0
        b4_14tot=0
#endif
cnb ###        npp_tot=0

!        do k=1,NS
!          do i=1,IT
!            if (FRLND(i,k).gt.0.001) then
!
!            carea(i,k)=SQRA(i,k)*FRLND(i,k)*1.E-12
        do k=1,nlon
          do i=1,nlat
           if (fracgr(i,k).gt.0.001) then
            carea(i,k)=darea(i)*fracgr(i,k)*1.E-12

c            cav_la=cav_la+(b1t(i,k)*st(i,k)+b1g(i,k)*sg(i,k)
c     >  +b2t(i,k)*st(i,k)+b2g(i,k)*sg(i,k)+b3t(i,k)*st(i,k)
c     >  +b3g(i,k)*sg(i,k)+b4t(i,k)*st(i,k)+b4g(i,k)*sg(i,k))
c     >  *carea(i,k)*(1.-icemask(i,k))

#if ( IMSK == 1 )
        b1_tot=b1_tot+(b1t(i,k)*st(i,k)+b1g(i,k)*sg(i,k))
     >  *carea(i,k)*(1.-icemask(i,k))                         !PgC

        b2_tot=b2_tot+(b2t(i,k)*st(i,k)+b2g(i,k)*sg(i,k))
     >  *carea(i,k)*(1.-icemask(i,k))                         !PgC

        b3_tot=b3_tot+(b3t(i,k)*st(i,k)+b3g(i,k)*sg(i,k))
     >  *carea(i,k)*(1.-icemask(i,k))                         !PgC

        b4_tot=b4_tot+(b4t(i,k)*st(i,k)+b4g(i,k)*sg(i,k))
     >  *carea(i,k)*(1.-icemask(i,k))                         !PgC

        b1_13tot=b1_13tot+(b1t13(i,k)*st(i,k)+b1g13(i,k)*sg(i,k))
     >  *carea(i,k)*(1.-icemask(i,k))                         !PgC14              

        b2_13tot=b2_13tot+(b2t13(i,k)*st(i,k)+b2g13(i,k)*sg(i,k))
     >  *carea(i,k)*(1.-icemask(i,k))                         !PgC14

        b3_13tot=b3_13tot+(b3t13(i,k)*st(i,k)+b3g13(i,k)*sg(i,k))
     >  *carea(i,k)*(1.-icemask(i,k))                         !PgC14

        b4_13tot=b4_13tot+(b4t13(i,k)*st(i,k)+b4g13(i,k)*sg(i,k))
     >  *carea(i,k)*(1.-icemask(i,k)) 

#if ( KC14 == 1 )
cvm calcul des reservoirs de c14 dans la biosphere continentale
        if (st(i,k).lt.(tiny(st(i,k))*1000.d0)) then
                st(i,k) = 0.0_dblp
        endif

        b1_14tot=b1_14tot+(b1t14(i,k)*st(i,k)+b1g14(i,k)*sg(i,k))
     >  *carea(i,k)*(1.-icemask(i,k))                         !PgC14              

        b2_14tot=b2_14tot+(b2t14(i,k)*st(i,k)+b2g14(i,k)*sg(i,k))
     >  *carea(i,k)*(1.-icemask(i,k))                         !PgC14

        b3_14tot=b3_14tot+(b3t14(i,k)*st(i,k)+b3g14(i,k)*sg(i,k))
     >  *carea(i,k)*(1.-icemask(i,k))                         !PgC14

        b4_14tot=b4_14tot+(b4t14(i,k)*st(i,k)+b4g14(i,k)*sg(i,k))
     >  *carea(i,k)*(1.-icemask(i,k))                         !PgC14
#endif
#else
        b1_tot=b1_tot+(b1t(i,k)*st(i,k)+b1g(i,k)*sg(i,k))
     >  *carea(i,k)*(1.)                         !PgC

        b2_tot=b2_tot+(b2t(i,k)*st(i,k)+b2g(i,k)*sg(i,k))
     >  *carea(i,k)*(1.)                         !PgC

        b3_tot=b3_tot+(b3t(i,k)*st(i,k)+b3g(i,k)*sg(i,k))
     >  *carea(i,k)*(1.)                         !PgC

        b4_tot=b4_tot+(b4t(i,k)*st(i,k)+b4g(i,k)*sg(i,k))
     >  *carea(i,k)*(1.)                         !PgC
#if ( KC14 == 1 )
cvm calcul des reservoirs de c14 dans la biosphere continentale

        b1_14tot=b1_14tot+(b1t14(i,k)*st(i,k)+b1g14(i,k)*sg(i,k))
     >  *carea(i,k)*(1.)                         !PgC14              

        b2_14tot=b2_14tot+(b2t14(i,k)*st(i,k)+b2g14(i,k)*sg(i,k))
     >  *carea(i,k)*(1.)                         !PgC14

        b3_14tot=b3_14tot+(b3t14(i,k)*st(i,k)+b3g14(i,k)*sg(i,k))
     >  *carea(i,k)*(1.)                         !PgC14

        b4_14tot=b4_14tot+(b4t14(i,k)*st(i,k)+b4g14(i,k)*sg(i,k))
     >  *carea(i,k)*(1.)                         !PgC14
#endif
#endif
           endif
          enddo
        enddo

#if ( KC14 == 1 )
        lc14_tot = b1_14tot + b2_14tot + b3_14tot + b4_14tot  !PgC14
#endif

! Here the same things are missing for the ocean. Could be copied /
! pasted from the code of CLIMBER-svn, file out.f, lines ~917


#if ( KC14 == 1 )
cvm - Ocean c14 reservoir
      oc14_tot=0
      odic_tot=0

      do n=1,NOC_CBR
       do i=1,LT
        do j=1,JT
         if (MGT(i,j,n).eq.1) then
        oc14_tot=oc14_tot+OC14(i,J,n)*1.e6
     >     *DVOL(i,J,n)*14*SCALE_M*1.028

         odic_tot=odic_tot+ODIC(i,J,n)*1.e6
     >     *DVOL(i,J,n)*12*SCALE_M*1.028


         endif
        enddo
       enddo
      enddo


      WRITE(*,*) "Testing ...", KOD
      WRITE(*,*) "1: ", oc14_tot, odic_tot, c14rstd
      WRITE(*,*) "2: ", lc14_tot, cav_la, c14rstd

cvm - Calcul de c14_total : ocean + land 
      dc14atm = ((c14atm/PA_C)/c14rstd - 1.d0) * 1000.d0          !permil
      dc14ocn = ((oc14_tot/odic_tot)/c14rstd - 1.d0) * 1000.d0  !permil     
      !dmr [NOTA] fix for first call when cav_la is 0.0d0 - bad!
      if (cav_la.lt.epsilon(cav_la)) then
          dc14lnd = 0.0d0
      else
          dc14lnd = ((lc14_tot/cav_la)/c14rstd - 1.d0) * 1000.d0    !permil
      endif
      c14_total = oc14_tot + lc14_tot                     !
  
cvm - Calcul de C14 en moles par reservoir
      n14atm = c14atm  / 0.40d0 * 1d15 / 14.d0
      n14ocn = cav_oc14 / 14.d0
      n14lnd = cav_la14 / 14.d0 
      n14tot = n14atm + n14ocn + n14lnd            
      n12tot = carbone_total * 1.d15 / 12.d0 
      r14tot = n14tot / n12tot * 1d12     !vm : *1e12  = pour manipuler des valeurs faibles
      alpha14 = ( r14tot / c14rstd ) / 1d12
#endif


      if (cav_oc.gt.0.0d0) then
        dum_toprint_un = cav_oc13/cav_oc
      else
        dum_toprint_un = 0.0d0
      endif
      if (cav_la.gt.0.0d0) then
        dum_toprint_de = cav_la13/cav_la
      else
        dum_toprint_de = 0.0d0
      endif

      WRITE(c13_res_fich%id,'(i6,10f10.2)') NYR, ca13_oc_ini/ca_oc_ini
     &, dum_toprint_un, ca13_la_ini/ca_la_ini, dum_toprint_de
     &, dc13at_ini-1000.0, c13atm-1000.0
     &, (b1_13tot/b1_tot)-1000.0, (b2_13tot/b2_tot)-1000.0
     &, (b3_13tot/b3_tot)-1000.0, (b4_13tot/b4_tot)-1000.0


!nb - Write in C_reservoir.txt
!nb      print*,NYR,b1,b2,b3,b4
      write (c_res_fich,'(i6,16f10.2)') NYR,PA_C,PA0_C,
     >   c13atm-1000.0,
     >   cav_oc,ca_oc_ini,cav_la,ca_la_ini,
     >   carbone_total,carbone_total_ini,
     >   b1_tot,b2_tot,b3_tot,b4_tot
!     >   ,cav_oc_detail(1),cav_oc_detail(2),
!     >   cav_oc_detail(3)

#if ( KC14 == 1 )
cvm - Write in C14_reservoirs.txt
      write (c14_res_fich,'(i6,3f14.2,6f16.2)') 
     >   NYR,PA_C,PA0_C,carbone_total, 
     >   mc14atm, cav_oc14, cav_la14,
     >   carbone14_total, mPRODC14, limC14masse

cvm - Write in C14_bilan.txt
      write (c14_bilan_fich,'(i6,10f17.2)')
     >  NYR, n14atm, n14ocn, n14lnd, n14tot, n14lim,
     >  r14tot, alpha14

cvm - Write in DC14.txt
      write (dc14_fich,'(i6,6f14.2)') 
     >  NYR, dc14atm, dc14ocn, dc14lnd

cvm - Write in C_flux.txt
      write (c_flux_fich,'(i6,2f11.2,6e16.5)') NYR,PA_C,PA0_C,
     >   FC14OA, FC14LA, FC12OA, FC12LA, FDIC, FOAC14

#endif

!-----|--1--------2---------3---------4---------5---------6---------7-|
!
!-----|--1--------2---------3---------4---------5---------6---------7-|
    
       IF (KOD.EQ.-2) THEN

       CLOSE(c_res_fich)
       call close_f(c13_res_fich)
       call close_f(c_ocean_fich)

#if ( KC14 == 1 )
       CLOSE(c14_res_fich)
       CLOSE(c_flux_fich)
       CLOSE(dc14_fich)
       CLOSE(c14_bilan_fich)
#endif

#if ( CORAL == 1 )
       CLOSE(coral_res_fich)
       CLOSE(c_nino_fich)
#endif
       ENDIF
       
       END SUBROUTINE out_cycc
