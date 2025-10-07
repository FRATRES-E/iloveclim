!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:58 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:58 CET 2009

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8--|
                   SUBROUTINE ECO2(KOD,fracgr,darea)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8--|
!
!                    CO2 ATMOSPHERIC CONCENTRATION
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8--|
!
!dmr&nb On suppose que le CO2 est DEJA initialise avant cette SUBROUTINE
!
!NB : not anymore, CO2 (and C13atm) is init in restart_cc
! only if restart (KLSR=1), if fresh start (KLSR=0) then 
! rest_cc is not read
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8--|
! dmr
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8--|


      use global_constants_mod, only: dblp=>dp, ip

      use veget_iso,      only: b4g13, b4g14, b3g13, b3g14, b2g13, b2g14
     &                        , b1g13, b1g14, b4t13, b4t14, b3t13, b3t14
     &                        , b2t13, b2t14, b1t13, b1t14
      use carbone_co2,    only: c14atm, c14dec, c14rstd, new_run_c, pa0_c, pa_c
     &                        , pa_c_d
#if ( CEMIS == 1 )
     &                        , cemis, nb_emis
      use comunit,        only: iuo
      use newunit_mod,    only: carbon_emission_dat_id
#endif
#if ( CARAIB > 0 )
      use ec_ca2lbm,      only: stock_carbon_caraib, veget_frac
      use comsurf_mod,    only: fractn, nld
      use ec_co2ca,       only: frac_land, j_antarctique
#endif

      use C_res_mod,      only: c13atm, ca13_at_ini, ca13_la_ini
     &             , ca13_la_rest, ca13_oc_ini, ca13_oc_rest, ca14_la_ini
     &             , ca14_la_rest, ca14_oc_ini, ca14_oc_rest, ca_la_ini
     &             , ca_la_rest, ca_oc_ini, ca_oc_rest, ca_oc_vol, cav_la
     &             , cav_la13, cav_la14, cav_la14_b, cav_la_b, cav_la_p, cav_oc
     &             , cav_oc13, cav_oc14, cav_oc14_b, cav_oc2, cav_oc_b, cav_oc_p
     &             , coc_odoc, coc_odoc13, coc_odocs, coc_odocs13, dc13at_ini
     &             , fc14la, fc14oa
     &             ,emis_cum, emis_c13_cum

      use loveclim_transfer_mod, only: KLSR

#if ( IMSK == 1 )
      use input_icemask,  only: icemask
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8--|
! dmr  Marine carbon cycle variables
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8--|

#if ( OCYCC == 1 )
      use declars_mod,    only: LT, JT, NOC_CBR
      use loveclim_transfer_mod, only: dvol, mgt
      use mod_sync_time, only: kendy, nyr
      use mbiota_mod     ,only: scale_m 
      use mbiota_mod,     only: zoo_m, phyto_m, PHYTO_M13, ZOO_M13
      use marine_bio_mod, only: odic, oc13, oc14, odoc, odocs, odoc13, odocs13
     &                        , odic_diff, odoc, opoc
#endif

#if ( CORAL == 1 )
      use coral_mod, only: C_car_a
      use mbiota_mod, only: SCANU
#endif

!dmr&nb [TODO] !!
!tbd #if ( MEDUSA == 1 )
!      use mbiota_mod,     only: summary_flux_O2S_calc, 
!     &                          summary_flux_O2S_orgm,
!     &                          summary_flux_S2O_dic, 
!tbd     &                          riverine_input_dic
!#endif
      
      use comatm,         only: nlat, nlon

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8--|
! dmr   VECODE carbon cycle variables
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8--|

      use veget_mod,      only: carea, sg, b4g, st, b4t, b3g, b3t
     &                        , b2g,b2t,b1g, b1t

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8--|
! dmr   CARAIB carbon cycle variables
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8--|

#if ( CARAIB > 0 )
      use ec_co2ca,  only: j_antarctique, frac_land, n_pix_caraib
#endif

#if ( INTERACT_CYCC == 2 )
      use atmos_composition_mod, only: get_PGA_CO2
#endif

      implicit none


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8--|
! dmr&nb Declaration des parametres de la SUBROUTINE
!>    @param[in]  KOD
!>    @param[in]  fracgr
!>    @param[in]  darea
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8--|

       integer(kind=ip), intent(in) :: KOD
       real(kind=dblp), dimension(nlat,nlon), intent(in)  :: fracgr
       real(kind=dblp), dimension(nlat),      intent(in)  :: darea


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8--|
! dmr   Variables locales
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8--|

       real(kind=dblp) :: ca_beta, d_oc, d_la, PCO2D, PCO2VAR, sumflux
#if ( CEMIS == 1 )
       integer(kind=ip) :: ii
       integer(kind=ip) :: yy
#endif
       integer(kind=ip):: i,j,k,n

#if ( CORAL == 1 )
       real(kind=dblp) :: ca_car_a
#endif
!dmr&nb [TODO] !!
!tbd #if ( MEDUSA == 1 )
!       real(kind=dblp) :: delta_carb_MEDUSA
!     &                  , delta_riverndic_MEDUSA
!       real(kind=dblp), save :: cumul_delta_carbMEDUSA
!     &                        , cumul_rivern_inpMEDUSA
!#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8--|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8--|


cnb
#if ( CARAIB > 0 )
!#include "/home/climwork/textier/caraib-git/com_18/parameter.common"
      !integer pix_car
      !parameter (pix_car=541)
!#include "/home/acclimate/nbouttes/caraib-git_cpl/com_18/
!     >veccarb.common"
      common /veccarb/ stock_ysoilr(n_pix_caraib)
      real*4  stock_ysoilr

      integer(kind=4):: ij
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8--|
!dmr&nb Facteur de conversion p.p.m <=> GigaTonnes Carbone
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8--|

      ca_beta=0.47


C.... initialization of global carbon pool, beginning of simulation

      if (KOD.eq.0) then
      WRITE(*,*)' '
      WRITE(*,*) 'initialisation in eco2'

!cnb  read restart for total carbon stock in each reservoir
!     and init PA0_C (CO2) and c13atm
!     only if restart (KLSR=1)
      if (KLSR.eq.1) then
         call restart_cc(1)
      endif

      ! PA0_C = PA0_C + 940.0
      PA_C=PA0_C
      PA_C_D=PA0_C
      WRITE(*,*) 'atmosphere PA_C, c13atm, c14atm'
      WRITE(*,*) PA_C, c13atm, c13atm-1000, c14atm

!#if ( KC14 == 1 )
!      C14ATM=C14ATM0
!#endif
      ca_oc_ini=0
      ca_la_ini=0
      ca13_oc_ini=0
      ca13_la_ini=0
#if ( KC14 == 1 )
      ca14_oc_ini=0.
      ca14_la_ini=0.
#endif
      ca_oc_vol=0
!      ca13_at_ini=(c13atm-1000.)*pa/ca_beta
      ca13_at_ini=(c13atm-1000.)*PA_C/ca_beta
      WRITE(*,*) 'c13_at_ini dans eco2 ', ca13_at_ini
      dc13at_ini = c13atm
      coc_odoc=0
      coc_odocs=0
      coc_odoc13=0
      coc_odocs13=0
      emis_cum=0
      emis_c13_cum=0

#if ( CEMIS == 1 )
!read input file for carbon emissions
!        open(1,file='inputdata/carbon_emission.dat',form='unformatted',
!     >       access='direct',recl=2*nb_emix)
!         read(iuo+40,*) yy
!         print *, yy
         do ii=1,nb_emis
           !read(iuo+40,*) yy, cemis(ii)
           read(carbon_emission_dat_id,*) yy, cemis(ii)
           print *, 'test emission ', ii,yy, cemis(ii)
         enddo
        !close(iuo+40)
        close(carbon_emission_dat_id)
#endif

#if ( OCYCC == 1 )
       do n=1,NOC_CBR
        do i=1,LT
          do j=1,JT
          if (MGT(i,j,n).eq.1) then
         ca_oc_ini=ca_oc_ini+(PHYTO_M(i,J,n)+ZOO_M(i,J,n)+
     <         ODOC(i,J,n)+OPOC(i,j,n)+ODOCS(i,J,n)+
     <         ODIC(i,J,n)*1.e6)*DVOL(i,J,n)*12*SCALE_M
     <               *1.028
!dmr&nb [AGARDER] ??
!         ca13_oc_ini=ca13_oc_ini+((PHYTO_M(i,J,n)+ZOO_M(i,J,n)+
!     <         OPOC(i,j,n))*OC13(i,J,n)/ODIC(i,J,n)
!     <         +ODOC13(i,J,n)+ODOCS13(i,J,n)+
!     <         OC13(i,J,n)*1.e6)*DVOL(i,J,n)*12*SCALE_M
!     <               *1.028
!nb 
         ca13_oc_ini=ca13_oc_ini+(PHYTO_M13(i,J,n)+ZOO_M13(i,J,n)
     <         +ODOC13(i,J,n)+ODOCS13(i,J,n)+
     <         OC13(i,J,n)*1.e6)*DVOL(i,J,n)*12*SCALE_M
     <               *1.028

#if ( KC14 == 1 )
         ca14_oc_ini=ca14_oc_ini+OC14(i,J,n)*1.e6*DVOL(i,J,n)*
     >               14*SCALE_M*1.028              ! TgC14
#endif
          ca_oc_vol=ca_oc_vol+DVOL(i,J,n)
          coc_odoc=coc_odoc+ODOC(i,J,n)*DVOL(i,J,n)
          coc_odocs=coc_odocs+ODOCS(i,J,n)*DVOL(i,J,n)
          coc_odoc13=coc_odoc13+ODOC13(i,J,n)*DVOL(i,J,n)
          coc_odocs13=coc_odocs13+ODOCS13(i,J,n)*DVOL(i,J,n)

          endif
          enddo
        enddo
       enddo
#if ( KC14 == 1 )
cvm&dmr --- Change unit of ca14_oc_ini to be consistent in use for cav_oc14_b
         ca14_oc_ini=ca14_oc_ini*1.e15                      ! gC14
#endif

#else /* OCYCC != 1 */
       ca_oc_ini=38000.
       ca13_oc_ini=-0.7*ca_oc_ini
#endif


        WRITE(*,*) 'carbon ocean ca_oc_ini' , ca_oc_ini
        WRITE(*,*) 'C13 ocean ca13_oc_ini' ,ca13_oc_ini
        WRITE(*,*) 'd13C ocean ', ca13_oc_ini/ca_oc_ini

!cnb use restart for total carbon stocks
        if (KLSR.eq.1) then
          WRITE(*,*) 'carbon ocean from restart ca_oc_rest' , ca_oc_rest
          ca_oc_ini=ca_oc_rest ! overwrite initial value using restart
          WRITE(*,*) 'C13 ocean from restart ca13_oc_rest', ca13_oc_rest
          ca13_oc_ini=ca13_oc_rest ! overwrite initial value using restart
#if ( KC14 == 1 )
          WRITE(*,*) 'C14 ocean from restart ca14_oc_rest', ca14_oc_rest
          ca14_oc_ini=ca14_oc_rest ! overwrite initial value using restart
#endif
        endif

        do k=1,nlon
            do i=1,nlat
#if (CARAIB == 0)
           if (fracgr(i,k).gt.0.001) then
            carea(i,k)=darea(i)*fracgr(i,k)*1.E-12
#else

            !carea(i,k)=darea(i)*1.E-12
           if (veget_frac(i,k).gt.0.001) then
            carea(i,k)=darea(i)*veget_frac(i,k)*1.E-12
            !write(*,*) 'carea', i,k, carea(i,k), veget_frac(i,k)
#endif
#if ( CARAIB == 0 )
#if ( IMSK == 1 )
            ca_la_ini=ca_la_ini+(b1t(i,k)*st(i,k)+b1g(i,k)*sg(i,k)
     <  +b2t(i,k)*st(i,k)+b2g(i,k)*sg(i,k)+b3t(i,k)*st(i,k)
     <  +b3g(i,k)*sg(i,k)+b4t(i,k)*st(i,k)+b4g(i,k)*sg(i,k))
     <  *carea(i,k)*(1.-icemask(i,k))
            ca13_la_ini=ca13_la_ini+((b1t13(i,k)+b2t13(i,k)+b3t13(i,k)
     >  +b4t13(i,k))*st(i,k)+(b1g13(i,k)+b2g13(i,k)+
     >      b3g13(i,k)+b4g13(i,k))*sg(i,k))*carea(i,k)
     <   *(1.-icemask(i,k))
#if ( KC14 == 1 )
            ca14_la_ini=ca14_la_ini+((b1t14(i,k)+b2t14(i,k)+b3t14(i,k)
     >  +b4t14(i,k))*st(i,k)+(b1g14(i,k)+b2g14(i,k)+
     >      b3g14(i,k)+b4g14(i,k))*sg(i,k))*carea(i,k)
     >   *(1.-icemask(i,k)) * 1e15                    !gC14
#endif

#else /* IMSK != 1 */
            ca_la_ini=ca_la_ini+(b1t(i,k)*st(i,k)+b1g(i,k)*sg(i,k)
     <  +b2t(i,k)*st(i,k)+b2g(i,k)*sg(i,k)+b3t(i,k)*st(i,k)
     <  +b3g(i,k)*sg(i,k)+b4t(i,k)*st(i,k)+b4g(i,k)*sg(i,k))
     <  *carea(i,k)*(1.)                              !Giga Tonnes Carbone -> vm
            ca13_la_ini=ca13_la_ini+((b1t13(i,k)+b2t13(i,k)+b3t13(i,k)
     >  +b4t13(i,k))*st(i,k)+(b1g13(i,k)+b2g13(i,k)+
     >      b3g13(i,k)+b4g13(i,k))*sg(i,k))*carea(i,k)
     <   *(1.)
#if ( KC14 == 1 )
            ca14_la_ini=ca14_la_ini+((b1t14(i,k)+b2t14(i,k)+b3t14(i,k)
     >  +b4t14(i,k))*st(i,k)+(b1g14(i,k)+b2g14(i,k)+
     >      b3g14(i,k)+b4g14(i,k))*sg(i,k))*carea(i,k)
     >   *(1.) * 1e15                                 !gC14
#endif
#endif
#elif ( CARAIB > 0 )
!#if ( IMSK == 1 )
!cnb  stock_carbon_caraib(lat, lon) en gC (ysoilr dans caraib) *1e-9*1e-6 into
! Gtonne C et *1e12 to remove factor of carea
!            ca_la_ini=ca_la_ini+(stock_carbon_caraib(i,k)*1e-9*1e-6)
!     >  *carea(i,k)*(1.-icemask(i,k))*1e12
!cnb alternative to init veget: set to 1500 GtC
             ca_la_ini=1500.0
!            ca13_la_ini=ca13_la_ini+((b1t13(i,k)+b2t13(i,k)+b3t13(i,k)
!     >  +b4t13(i,k))*st(i,k)+(b1g13(i,k)+b2g13(i,k)+
!     >      b3g13(i,k)+b4g13(i,k))*sg(i,k))*carea(i,k)
!     <   *(1.-icemask(i,k))
!#else /* IMSK != 1 */
!            ca_la_ini=ca_la_ini+(b1t(i,k)*st(i,k)+b1g(i,k)*sg(i,k)
!     <  +b2t(i,k)*st(i,k)+b2g(i,k)*sg(i,k)+b3t(i,k)*st(i,k)
!     <  +b3g(i,k)*sg(i,k)+b4t(i,k)*st(i,k)+b4g(i,k)*sg(i,k))
!     <  *carea(i,k)*(1.)                              !Giga Tonnes
Carbone -> vm
!            ca13_la_ini=ca13_la_ini+((b1t13(i,k)+b2t13(i,k)+b3t13(i,k)
!     >  +b4t13(i,k))*st(i,k)+(b1g13(i,k)+b2g13(i,k)+
!     >      b3g13(i,k)+b4g13(i,k))*sg(i,k))*carea(i,k)
!     <   *(1.)
!#endif
#endif

           endif
          enddo
        enddo

        ca13_la_ini=ca13_la_ini-1000.*ca_la_ini

        WRITE(*,*) 'carbon vegetation ca_la_ini' , ca_la_ini
        WRITE(*,*) 'C13 vegetation ca13_la_ini' , ca13_la_ini
        WRITE(*,*) 'dC13 vegetation ' , ca13_la_ini/ca_la_ini

!cnb use restart for total carbon stocks
        if (KLSR.eq.1) then
          WRITE(*,*) 'carbon vegetation from restart ca_la_rest' ,
     &             ca_la_rest
          ca_la_ini=ca_la_rest ! overwrite initial value using restart
          WRITE(*,*) 'C13 vegetation from restart ca13_la_rest' ,
     &    ca13_la_rest
          ca13_la_ini=ca13_la_rest ! overwrite initial value using restart
#if ( KC14 == 1 )
          WRITE(*,*) 'C14 vegetation from restart ca14_la_rest' ,
     &    ca14_la_rest
          ca14_la_ini=ca14_la_rest ! overwrite initial value using restart
#endif
       endif

        CALL out_cycc(-1,fracgr,darea)


#if ( KC14 == 1 )
       if (new_run_c.eq.1) then

         cav_oc14_b=ca14_oc_ini
         cav_oc_b=ca_oc_ini

cvm         if (INI_STEP.eq.0) cav_oc14_b=48498933.056
cvm         if (INI_STEP.ne.0) cav_oc14_b=cav_oc14
cvm         print *,"cav_oc14_b",cav_oc14_b
cvm         if (INI_STEP.eq.0) cav_oc_b=45312.48
cvm         if (INI_STEP.ne.0) cav_oc_b=cav_oc
cvm         print *,"cav_oc_b",cav_oc_b

cvm - on stocke cav_la14 du pas de temps precedent dans cav_la14_b
cvm avec b comme before
cvm cette variable est necessaire dans le calcul du flux de respiration
cvm de la biosphere continentale vers l atmosphere

         cav_la14_b=ca14_la_ini
         cav_la_b=ca_la_ini

cvm         if (INI_STEP.eq.0) cav_la14_b=2838890.264
cvm         if (INI_STEP.ne.0) cav_la14_b=cav_la14
cvm         print *,"cav_la14_b",cav_la14_b
cvm         if (INI_STEP.eq.0) cav_la_b=2388.76
cvm         if (INI_STEP.ne.0) cav_la_b=cav_la
cvm         print *,"cav_la_b",cav_la_b

       endif


#endif
       CALL out_cycc(1,fracgr,darea)

       WRITE(*,*)' '

!dmr&nb [TODO] !!
!tbd #if (MEDUSA == 1)
!       cumul_delta_carbMEDUSA=0.d0
!       cumul_rivern_inpMEDUSA=0.d0
!#endif       


       else ! KOD <> 0

!      WRITE(*,*), ' '
c======== Calculation of new atmospheric pCO2 ==========

C.... total ocean carbon
c 1.028 - reference seawater density (kg/m3)

       cav_oc=0
       cav_oc2=0
       cav_oc13=0
#if ( KC14 == 1 )
       cav_oc14=0
#endif
#if ( OCYCC == 1 )
       do n=1,NOC_CBR
        do i=1,LT
          do j=1,JT
          if (MGT(i,j,n).eq.1) then
         cav_oc=cav_oc+(PHYTO_M(i,J,n)+ZOO_M(i,J,n)+
     <         ODOC(i,J,n)+OPOC(i,j,n)+ODOCS(i,J,n)+
     <         ODIC(i,J,n)*1.e6)*DVOL(i,J,n)*12*SCALE_M
     <               *1.028
! ### Bricolage
         cav_oc2=cav_oc2+(PHYTO_M(i,J,n)+ZOO_M(i,J,n)+
     <         ODOC(i,J,n)+OPOC(i,j,n)+ODOCS(i,J,n)+
     <         (ODIC(i,J,n)+ODIC_diff(i,j,n))*1.e6)*DVOL(i,J,n)
     <          *12*SCALE_M*1.028
! ### Bricolage
!              cav_oc13=cav_oc13+((PHYTO_M(i,J,n)+ZOO_M(i,J,n)+
!     <         OPOC(i,j,n))*OC13(i,J,n)/ODIC(i,J,n)
!     <         +ODOC13(i,J,n)+ODOCS13(i,J,n)+
!     <         OC13(i,J,n)*1.e6)*DVOL(i,J,n)*12*SCALE_M
!     <               *1.028
!nb 
             cav_oc13=cav_oc13+(PHYTO_M13(i,J,n)+ZOO_M13(i,J,n)
     <         +ODOC13(i,J,n)+ODOCS13(i,J,n)+
     <         OC13(i,J,n)*1.e6)*DVOL(i,J,n)*12*SCALE_M
     <               *1.028
               endif

#if ( KC14 == 1 )
              cav_oc14=cav_oc14+(OC14(i,J,n))*1.e6
     >               *DVOL(i,J,n)*14*SCALE_M
     >                  *1.028 * 1e15                               !gC14
#endif


          enddo
!dmr&vm --- Attention code a verifier par Veronique !!

!         if (KTIME.eq.1) then !real time
!          NYR0 = NYRSR + NYR
!         else !not real time : fixed
!          NYR0 = NYRSR
!         endif
!        NYR01 = -NYR0

!dmr&vm --- Attention code a verifier par Veronique !!

        enddo
       enddo
c~        write(*,*) "cav_oc ===", cav_oc
!              print*, "cav_oc ", cav_oc
!              print*, "cav_oc13 ", cav_oc13
#endif


#if ( KC14 == 1 )
cvm - puis reinitialisation des variables de masse totale de la biosphere
cvm continentale pour chaque isotope du carbone
        cav_la14=0
#endif
        cav_la=0
        cav_la13=0

        do k=1,nlon
            do i=1,nlat
#if (CARAIB == 0)
           if (fracgr(i,k).gt.0.001) then
            carea(i,k)=darea(i)*fracgr(i,k)*1.E-12
#else
            !carea(i,k)=darea(i)*1.E-12
           if (veget_frac(i,k).gt.0.001) then
            carea(i,k)=darea(i)*veget_frac(i,k)*1.E-12
            !write(*,*) 'carea', i,k, carea(i,k), veget_frac(i,k)
#endif

#if ( CARAIB == 0 )
#if ( IMSK == 1 )
            cav_la=cav_la+(b1t(i,k)*st(i,k)+b1g(i,k)*sg(i,k)
     <  +b2t(i,k)*st(i,k)+b2g(i,k)*sg(i,k)+b3t(i,k)*st(i,k)
     <  +b3g(i,k)*sg(i,k)+b4t(i,k)*st(i,k)+b4g(i,k)*sg(i,k))
     <  *carea(i,k)*(1.-icemask(i,k))

            cav_la13=cav_la13+((b1t13(i,k)+b2t13(i,k)+b3t13(i,k)
     >  +b4t13(i,k))*st(i,k)+(b1g13(i,k)+b2g13(i,k)+
     >      b3g13(i,k)+b4g13(i,k))*sg(i,k))*carea(i,k)
     <   *(1.-icemask(i,k))
#if ( KC14 == 1 )
            cav_la14=cav_la14+((b1t14(i,k)+b2t14(i,k)+b3t14(i,k)
     <  +b4t14(i,k))*st(i,k)+(b1g14(i,k)+b2g14(i,k)+
     <      b3g14(i,k)+b4g14(i,k))*sg(i,k))*carea(i,k)
     <   *(1.-icemask(i,k)) * 1e15                            !gC14
#endif

#else /* IMSK != 1 */
            cav_la=cav_la+(b1t(i,k)*st(i,k)+b1g(i,k)*sg(i,k)
     <  +b2t(i,k)*st(i,k)+b2g(i,k)*sg(i,k)+b3t(i,k)*st(i,k)
     <  +b3g(i,k)*sg(i,k)+b4t(i,k)*st(i,k)+b4g(i,k)*sg(i,k))
     <  *carea(i,k)*(1.)
            cav_la13=cav_la13+((b1t13(i,k)+b2t13(i,k)+b3t13(i,k)
     >  +b4t13(i,k))*st(i,k)+(b1g13(i,k)+b2g13(i,k)+
     >      b3g13(i,k)+b4g13(i,k))*sg(i,k))*carea(i,k)
     <   *(1.)
#if ( KC14 == 1 )
            cav_la14=cav_la14+((b1t14(i,k)+b2t14(i,k)+b3t14(i,k)
     <  +b4t14(i,k))*st(i,k)+(b1g14(i,k)+b2g14(i,k)+
     <      b3g14(i,k)+b4g14(i,k))*sg(i,k))*carea(i,k)
     <   *(1.) * 1e15                            !gC14
#endif
#endif
!#elif ( CARAIB > 0 )
!#if ( IMSK == 1 )
            !gC to GtC et on enleve le facteur 1e12 de carea
            !write(*,*) 'dans eco2 en gC: ',i,j,stock_carbon_caraib(i,k)
             !write(*,*) 'carea bis', i,k, carea(i,k), 1.-icemask(i,k)
!            cav_la=cav_la+(stock_carbon_caraib(i,k)*1e-9*1e-6)

!            IF (KENDY.EQ.1)
!     >        write(*,*), 'calcul cav_la', stock_carbon_caraib(i,k),
!     >        stock_carbon_caraib(i,k)*1e-9*1e-6
!            cav_la13=cav_la13+((b1t13(i,k)+b2t13(i,k)+b3t13(i,k)
!     >  +b4t13(i,k))*st(i,k)+(b1g13(i,k)+b2g13(i,k)+
!     >      b3g13(i,k)+b4g13(i,k))*sg(i,k))*carea(i,k)
!     <   *(1.-icemask(i,k))
!#else /* IMSK != 1 */
!            cav_la=cav_la+(b1t(i,k)*st(i,k)+b1g(i,k)*sg(i,k)
!     <  +b2t(i,k)*st(i,k)+b2g(i,k)*sg(i,k)+b3t(i,k)*st(i,k)
!     <  +b3g(i,k)*sg(i,k)+b4t(i,k)*st(i,k)+b4g(i,k)*sg(i,k))
!     <  *carea(i,k)*(1.)
!            cav_la13=cav_la13+((b1t13(i,k)+b2t13(i,k)+b3t13(i,k)
!     >  +b4t13(i,k))*st(i,k)+(b1g13(i,k)+b2g13(i,k)+
!     >      b3g13(i,k)+b4g13(i,k))*sg(i,k))*carea(i,k)
!     <   *(1.)
!#endif
#endif



           endif
          enddo
        enddo

#if ( CARAIB > 0 )
      cav_la=stock_carbon_caraib*1e-15 !in PgC
#endif

        IF (KENDY.EQ.1) WRITE(*,*) 'carbon vegetation cav_la ', cav_la

!cnb a laisser !! car la veget utilise c13atm qui est d13C+1000
        cav_la13=cav_la13-1000.*cav_la
!        WRITE(*,*), 'test cav_la13 ', cav_la13

#if ( CEMIS == 1 )
         if (kendy.eq.1) then !emissions at last day of year
           !cemis in MtCO2
           !emis_cum=emis_cum + cemis(NYR)*1e-3 / 3.67
           !cemis in GtC
           emis_cum=emis_cum + cemis(NYR)
ccc        emis_c13_cum=emis_c13_cum + cemis(NYR)*(-25.)
          if (emis_cum.ne.0) print*,'cemis,emis_cum',NYR, cemis(NYR), emis_cum
         endif
#endif


#if ( CORAL == 0)
        PA_C=PA0_C-(cav_oc-ca_oc_ini+cav_la-ca_la_ini-emis_cum)*ca_beta

!nb test PA_C fixe
!        PA_C=284 !ppm

#elif ( CORAL > 0 )
        !nb with corals remove CO2 from carbonate weathering

        !ca_car_a=C_car_a*SCANU/(TYER/TSTOC)*1e6*12*1.028 ! Pmol/day*1e-6/(360*TDAY/TDAY) 
        ca_car_a=C_car_a*SCANU*1e6*12 ! Pmol/day*1e-6*1e6*g/mol=1e15g/day=Pg/day
        PA_C=PA0_C-(cav_oc-ca_oc_ini+cav_la-ca_la_ini-emis_cum+ca_car_a)*ca_beta !GtC *ca_beta
        !write(*,*), 'ca_car_a et PA_C', ca_car_a, PA_C

!nb [NOTA] Fait ailleurs ? 

!        delta_carb_MEDUSA = summary_flux_O2S_calc 
!     &                    + summary_flux_O2S_orgm 
!     &                    + summary_flux_S2O_dic ! in TmolC.yr-1!

                             ! / 360. yr-1 -> day-1
                             ! / 1000.0 Tmol -> Pmol
                             ! * 12     Pmol -> GtC
                             ! * 0.47   GtC  -> ppm
                             ! end result in ppm.day-1

!        delta_carb_MEDUSA = delta_carb_MEDUSA 
!     &                      / 360.0D0 / 1000.D0 *12.D0 *0.47D0 

!        cumul_delta_carbMEDUSA = cumul_delta_carbMEDUSA 
!     &                                       + delta_carb_MEDUSA
     
!tbd        delta_riverndic_MEDUSA = riverine_input_dic
!tbd     &                      / 360.0D0 / 1000.D0 *12.D0 *0.47D0 
          
!tbd        cumul_rivern_inpMEDUSA = cumul_rivern_inpMEDUSA
!tbd     &                                       + delta_riverndic_MEDUSA

!        PA_C=PA0_C-(cav_oc-ca_oc_ini+cav_la-ca_la_ini-emis_cum)*ca_beta
!     &            -cumul_delta_carbMEDUSA+cumul_rivern_inpMEDUSA
!nb est ce que on ne devrait pas plutot utiliser
! iloveclim_calc_lbtot [molC/yr], 
                            ! / 360.    yr-1 -> day-1
                             ! * 12     molC -> gC
                             ! *1e-15   gC -> GtC
                             ! * ca_beta=0.47   GtC  -> ppm
                             ! end result in ppm.day-1 
!     &            - iloveclim_calc_lbtot / 360.0 *12*1e-15*ca_beta ! units ??
!nb test PA_C fixe
!        PA_C=284 !ppm
     
!        IF (KENDY.EQ.1) THEN
!     
!        WRITE(*,'(A,8D12.5)') "Summary_MEDUSA_fluxes", 
!     &            PA_C,
!     &            summary_flux_O2S_calc 
!     &                    + summary_flux_O2S_orgm 
!     &                    + summary_flux_S2O_dic                   ! in TmolC.yr-1
!     &            , riverine_input_dic                             ! in TmolC.yr-1
!     &            , delta_carb_MEDUSA                              ! in ppm.day-1
!     &            , cumul_delta_carbMEDUSA                         ! in ppm
!     &            , delta_riverndic_MEDUSA                         ! in ppm.day-1
!     &            , cumul_rivern_inpMEDUSA                         ! in ppm
!     &            , -cumul_delta_carbMEDUSA+cumul_rivern_inpMEDUSA ! in ppm
!
!        ENDIF
     

#endif

#if ( INTERACT_CYCC == 2 )
!nb CO2 of carbon cycle (PA_C) = CO2 from radiative code (read in
!GHG.dat)
        PA_C=get_PGA_CO2()
       !write(*,*) 'PGACO2, PA_C', get_PGA_CO2(), PA_C
#endif


!        WRITE(*,*), 'stocks carbon in eco2, PA_C', PA_C
!        WRITE(*,*), 'cav_oc ', cav_oc, 'cav_la', cav_la

      PA_C_D=PA0_C-(cav_oc2-ca_oc_ini+cav_la-ca_la_ini-emis_cum)*ca_beta

        ODIC_diff = 0.0


       c13atm=1000+(ca13_at_ini+emis_c13_cum-
     & (cav_oc13-ca13_oc_ini+cav_la13-ca13_la_ini))/(PA_C/ca_beta)

!nb test valeur fixee dans atm
!       c13atm=993.58

!       WRITE(*,*), 'carbon 13 eco2, c13atm ', c13atm, c13atm-1000
!       WRITE(*,*) 'ca13_at_ini, cav_oc13, cav_la13',
!     &   ca13_at_ini, cav_oc13, cav_la13


cvm - ajout des flux de la biosphere continentale vers l atmosphere
cvm et de l'ocean vers l'atmosphere
cvm - Variation de n14atm liee au flux de la biosphere continentale
cvm vers l'atm s'exprime ainsi:
cvm n14atm(t)=n14atm(t-1)-c14dec*n14lnd(t)-(n14lnd(t)-n14lnd(t-1))
cvm d'ou l'expression ci-dessous
cvm        if (KC14.eq.0)  C14ATM=c14rstd*PA

        IF (KENDY.EQ.1) then
#if ( KC14 == 0 )
           C14ATM=c14rstd*PA_C
#elif ( KC14 == 1 )
           WRITE(*,*) " /*/*/*/ "
           WRITE(*,*) "eco2.f before flux, C14ATM = ", C14ATM

           sumflux = ((cav_la14_b-(1+c14dec)*cav_la14)/1e15
     >        + (cav_oc14_b-(1+c14dec)*cav_oc14)/1e15
cvm     <        + FOAC14* 1.e6 * 1.028 * 14 * SCALE_M
     >           ) * 0.40

cvm           sumflux = ((cav_la14_b-cav_la14)/1e15
cvm     >        + (cav_oc14_b-cav_oc14)/1e15
cvmcvm     <        + FOAC14* 1.e6 * 1.028 * 14 * SCALE_M
cvm     >           ) * 0.40


!dmr&vm --- Following quick fix is NOT conservative !!!
           IF (ABS(sumflux).GE.C14ATM) THEN ! we get a negative atmosphere in 14C !!! Not allowed !!

             WRITE(*,*)
             WRITE(*,*)
             WRITE(*,*)
             WRITE(*,*) 'in eco2.f: Fc14->atm not conservative'
             WRITE(*,*) "====================================="
             WRITE(*,*)
             WRITE(*,*) "Diagnostics ..."
             WRITE(*,*)
             WRITE(*,*) cav_la14_b, cav_la14, cav_oc14_b, cav_oc14
             WRITE(*,*) C14ATM, (C14ATM/PA_C/1.176E-12-1.0)*1000.0
             WRITE(*,*)
             WRITE(*,*)
             WRITE(*,*)

             C14ATM=c14rstd*PA_C
           ELSE
             C14ATM=C14ATM+((cav_la14_b-(1+c14dec)*cav_la14)/1e15
     >        + (cav_oc14_b-(1+c14dec)*cav_oc14)/1e15
cvm     <        + FOAC14* 1.e6 * 1.028 * 14 * SCALE_M
     >           ) * 0.40
cvm             C14ATM=C14ATM+((cav_la14_b-cav_la14)/1e15
cvm     >        + (cav_oc14_b-cav_oc14)/1e15
cvmcvm     <        + FOAC14* 1.e6 * 1.028 * 14 * SCALE_M
cvm     >           ) * 0.40
           ENDIF

cvm           FC14OA = (cav_oc14_b-(1+c14dec)*cav_oc14)/1e3   !flux ocn->atm (kgC14)
cvm           FC14LA = (cav_la14_b-(1+c14dec)*cav_la14)/1e3   !flux lnd->atm (kgC14)
cvm           FC14OA = (cav_oc14_b-cav_oc14)/1e3   !flux ocn->atm (kgC14)
cvm           FC14LA = (cav_la14_b-cav_la14)/1e3   !flux lnd->atm (kgC14)


           FC14OA = (cav_oc14_b-(1+c14dec)*cav_oc14)/1e15*0.40
           FC14LA = (cav_la14_b-(1+c14dec)*cav_la14)/1e15*0.40

cvm           FC14OA = (cav_oc14_b-cav_oc14)/1e15*0.40
cvm           FC14LA = (cav_la14_b-cav_la14)/1e15*0.40

           WRITE(*,*) "eco2.f after flux, C14ATM = ", C14ATM
           WRITE(*,*) "FC14OA = ", FC14OA
           WRITE(*,*) "FC14LA = ", FC14LA
           WRITE(*,*) " /*/*/*/ "
           WRITE(*,*)

#endif /* on KC14 */

        ENDIF  !(vm:KENDY)

        if (cav_oc_p.ne.0.) then
        d_oc=cav_oc-cav_oc_p
        d_la=cav_la-cav_la_p
        endif

!#if ( CEMIS == 1 )
!        open (1,file='OUT/carbon.dat',form='unformatted',
!     >       access='direct',recl=3*4)
!        write (1,rec=NYR) cemis(NYR),d_oc,d_la
!
!        close (1)
!#endif

        cav_oc_p=cav_oc
        cav_la_p=cav_la

      IF (KENDY.EQ.1) CALL out_cycc(NYR,fracgr,darea)

#if ( KC14 == 1 )
         IF (KENDY.EQ.1) THEN

cvm         if (INI_STEP.eq.0) cav_oc14_b=48498933.056
cvm         if (INI_STEP.ne.0) cav_oc14_b=cav_oc14
         cav_oc14_b=cav_oc14
cvm         print *,"cav_oc14_b",cav_oc14_b
cvm         if (INI_STEP.eq.0) cav_oc_b=45312.48
cvm         if (INI_STEP.ne.0) cav_oc_b=cav_oc
         cav_oc_b=cav_oc
cvm         print *,"cav_oc_b",cav_oc_b

cvm - on stocke cav_la14 du pas de temps precedent dans cav_la14_b
cvm avec b comme before
cvm cette variable est necessaire dans le calcul du flux de respiration
cvm de la biosphere continentale vers l atmosphere

cvm         if (INI_STEP.eq.0) cav_la14_b=2838890.264
cvm         if (INI_STEP.ne.0) cav_la14_b=cav_la14
         cav_la14_b=cav_la14
cvm         print *,"cav_la14_b",cav_la14_b
cvm         if (INI_STEP.eq.0) cav_la_b=2388.76
cvm         if (INI_STEP.ne.0) cav_la_b=cav_la
         cav_la_b=cav_la
cvm         print *,"cav_la_b",cav_la_b

         ENDIF
#endif

      endif ! KOD == 0


      return
      end subroutine eco2

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8--|
! dmr   The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8--|
