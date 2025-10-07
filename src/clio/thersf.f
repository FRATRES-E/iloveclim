!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:55 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:55 CET 2009

      SUBROUTINE thersf(ntnow,istep)
!dmr&afq      SUBROUTINE thersf(ntnow,ntrmax,istep,irunlabelclio)
! mab: ither=1,ntrmax (ntrmax = 1?!)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!   This routine prepares the computation of ice thermodynamics. It also
!   computes surface heat and salt flux at the surface of the ocean
!  modif : 24/12/99

      use comemic_mod, only: flag_snow_to_flux, iyear, nstpyear

#if ( ISOOCN >= 1 )
      USE iso_param_mod, ONLY : ieau, neauiso, ieau17, ieaud, ieau18
     & ,delta, rsmow

      use para0_mod,     only : owatert, owisostrt, owisostop, ocnw18
     & ,oc2atwisoindx

#if ( ISOBERG == 2)
       USE isoatm_mod, ONLY: dlbmini
#endif
#endif

!dmr&afq --- Commented out since code will have to be re-used in CONSEAU

#if ( SHELFMELT == 1 )
      !use interpol_clio_grisli_mod, only : interpol_one_field_clio
      use input_subgrid2L, only : icliog,jcliog,bmshelf_clionord
      use input_timerCplGRIS, only: timcplgrisday
#endif
#if ( F_PALAEO >=1 && UNCORFLUX >= 1 )
      use palaeo_timer_mod, only: celest_year, equil_flag
      use zonesOrmenfwf_mod, only: update_fwf
#endif

#if ( CONSEAU == 1 )
!--- dmr&afq Units:
!---           [ calgrisCLIO ]     = m**3.s-1
!---           [ calgrisCLIOms ]   = m.s-1
!---           [ fwatconseau ]     = m.s-1
      use varsCONSEAU_mod, only: calgrisCLIO, fwatconseau, calgrisCLIOms
#endif
#if ( F_PALAEO_FWF == 1 )
      use varsCONSEAU_mod, only: fwatconseau
#endif
#if ( F_PALAEO_FWF == 2 )
      use bloc0_mod, only: fwatconseau
#endif

      use const_mod,   only: zero, one, rho0, cpo, gpes

      use para0_mod,   only: nsmax, imax, jmax, kmax
      use para_mod,    only: nkb0
      use bloc0_mod,   only: phiss, phivs, scs, frapsav, scal, is1, is2, tms
     >             , unsdz, kfs, dts, scpme, dz, ust2s, scalr, rappes, js1, js2
     >             , ks1, ks2, jeq, spvr, bheat, zflux0s, zfluxms, vcor, zurfow
#if( APPLY_UNCORFWF == 1)
     >             , w_uncor
#endif
      use bloc_mod,    only: aire
      use thermo_mod,  only: npb, npac, tbqb, albqb, hgbqb, hnbqb, fsolgb
     >              , fbbqb, hnpbqb, thcmb, qlbqb, qcmbqb, qstbqb, tfub
     >              , tsb, dmgbqb, dmgwib, psbqb, tabqb, qabqb, vabqb, ratbqb
     >              , qlbbqb, cldqb, fscbqb, fltbqb, fstbqb, qfvbqb, dmnbqb
     >              , dvsbqb, dvbbqb, dvlbqb, dvnbqb, fratsb, fcsb, fleb, hgcolb
      use ice_mod,     only: hgbq, albq, qcmbq, qlbq, fbbq, hnbq, qstobq, fcm1
     >           , fcm2, fsbbq, ts, fleg, fcsg, firg, alct, thcm, dmgwi, hgcol
     >           , qfvbq, fstrbq, dvonbq, dvolbq, dvosbq, dvobbq, fscmbq, ffltbq
     >           , dmnbq, vabq, dmgbq, hgbqp, hnplbq, fwat, reslum, tfu, ratbqo
     >           , tbq, fcscn, flecn, fsolcn, tenagx, tenagy, fsolcn, tmics
     >           , fsolg, psbq, tabq, qabq, ratbqg, cloud, iicebers1, iicebers2
     >           , fevabq, ddtb, rhon, rhog, rhoesn, uscomi, hglim, xln, parlat
     >           , sglace, jicebers1, jicebers2, ficebergn, ficebergs, areicebn
     >           , areicebs, toticesm, iicebern1,iicebern2, jicebern1,jicebern2
     >           , xlg
#if ( ICEBERG > 0 )
      use iceberg_mod, only: dVol_icb,fw_icb, heat_icb
#endif      
      use dynami_mod,  only: ug, vg, area
!cfc  use trace_mod

      use global_constants_mod, only: dblp=>dp, ip
      use gather_and_scater_mod, only: gather, scater

!#if ( F_PALAEO_FWF == 1 && APPLY_UNCORFWF == 0 ) 
#if ( F_PALAEO_FWF >= 1 || CONSEAU == 1 ) 
      use update_clio_bathy_tools, only: sum_flux_out
#endif

      use newunit_clio_mod, only: clio3_out_id, mouchard_id      

!---
!Ccpl [Ccp0] => ligne specifique a la version avec [sans] couplage .
!---
      implicit none



      integer(kind=ip), intent(in) :: ntnow, istep

!
!dmr @-@ iceb0
!cfc  include 'trace.com'

      real*8 heaticbism(imax,jmax)
!dmr&afq unused      real*8 heaticbs,heaticbs2
!dmr&afq unused      real*8 phistots,phistotn

      real*8 fluxmean

#if ( UNCORFLUX == 1 )
      real*8 oldphiss(imax,jmax,0:nsmax)
#endif

      real   qref,qclim
      real   dzsdts ! ISO2 see below,DUM

!dmr&mab --- debug !
      real :: tot_Vicb_thersf = 0.0d0
!dmr&mab --- debug !

#if ( CONSEAU == 1 && HEATCALV == 1 )
      integer          :: kheatlimit = 12 ! i.e. approx. 300 meters
      double precision :: unsdz_som
#endif

      common/icbism/ heaticbism

      COMMON/CLIO2A/Qref,Qclim
      COMMON/CLIO2A2/fluxmean

!fdtcn     Transit variable used to compute the vertical heat flux
!     at the ocean surface in ice covered regions (local in thersf)
!
#if ( ISOOCN >= 1 )
      real(kind=dblp), dimension(imax,jmax,owatert) :: zhnpbq, zfwat
#else
      real(kind=dblp), dimension(imax,jmax)         :: zhnpbq, zfwat
#endif
      real(kind=dblp), dimension(imax,jmax) :: fdtcn, zhgbqp, ain, zinda, ifvt
     >               , qdtcn, qlbsbq
#if ( ISOOCN >= 1 )
!      dimension DUM(neauiso), pme(neauiso), zisowatsur(neauiso)
!      dimension variso(neauiso)
      real, dimension(owatert) :: DUM, pme, zisowatsur, variso
      integer :: atmisoindx, iz
      real(kind=dblp) :: totisotop, volocean
      character*10 varisonm
#else
      real :: DUM, pme
#endif
      common/lowat/ zfwat
!     dimension rrro(imax,jmax),dd1o(imax,jmax),dd2o(imax,jmax)
!dmr @-@ iceb0
! JONO
! JA heat and FWF 2layers
! plugging melted iceberg volume per layer per day dVol_icb(i,j,k) {m3/day}
! into sea-layer temperature and salinity scal(i,j,k,1) and scal(i,j,k,2)
! WITHOUT going thru a coupler!! BEWARE when changing timestep-sizes!!

      double precision voltohfl
!dmr @-@ iceb0

!
!nb debut
      real*8 fluxbrines(imax,jmax)
      common / common_brines / fluxbrines
!nb fin
#if( UNCORFLUX == 1 )
!dmr&mb ---  Adding an uncorrected freshwater flux ...
!afq ---      real*8 fluxarea, totflxarea, fluxinmeters(imax,jmax), ! flux in meters
      real*8 fluxarea, totflxarea,
     &       lawrflux,hudsonflux,waisflux,rudiflux
      real*8 flux_pourCLIO(imax,jmax)
!dmr&mb ---
#endif

#if ( UNCORFLUX == 1 || APPLY_UNCORFWF == 1 )
      real*8 fluxinmeters(imax,jmax) ! flux in meters
#endif

#if ( APPLY_UNCORFWF == 1 ) /* Apply uncorrected fwf to the ocean */
      real*8 given_fluxinmetres(imax,jmax) ! flux in meters
#endif

#if ( F_PALAEO_FWF >= 1 && APPLY_UNCORFWF == 0 ) 
      real*8 corr_fluxinmetres(imax,jmax) ! flux in meters
#endif

#if ( SHELFMELT == 1 )
      real*8, dimension (imax,jmax,kmax) :: bmshelf
      real*8 :: gammatg,fmeltg,fluxtg
#endif

      integer(kind=ip) :: i, j, nbpb, k, jp1, ip1, ns, kounter, ia, iadv, ial
     >                  , iflt, klayer, nbpac
      real(kind=dblp)  :: zeps0, zeps1, zeps2, ustmax, ustmin, za, zh, raptime
     >                  , zindb, zindg, zfontn, ustsg, fntlat, zfnsol, pareff
     >                  , alfagr, rhoo, rhofr, rhogr, gred, hg0, fac1hg, fac2hg
     >                  , gamafr, tenagm, dumfac, vfrx, vgx, vgy, vfry, vrel
     >                  , salinf, salflux, zqfv, zfdtm, fons, watsur
     >                  , aflux, dtemp1, dtemp2, fnor, fsud, zdepthi, tfreezloc
     >                  , gammat, scs_mean, zcums, zcums2, tmean
     >                  , smean, dztot, fluxt, phivssup, phivssups, ficesh
     >                  , ficess

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1) Initialization of diagnostic variables.                          |
!-----------------------------------------------------------------------
!       write(106,*)'albq,hgbq',albq(imax-2,5),hgbq(imax-2,5)
!
!age    dtj=dts(ks2)/(86400.*365.)
      do 5 j=js1,js2
         do  3 i=is1(j),is2(j)
            dvosbq(i,j) = 0.0
            dvobbq(i,j) = 0.0
            dvolbq(i,j) = 0.0
            dvonbq(i,j) = 0.0
 3       continue
 5    continue
!

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2) Numerical parameters.                                            |
!-----------------------------------------------------------------------
!
      zeps0 = 1.0e-16
      zeps1 = 1.0e-20
      zeps2 = 1.0e-04
      ustmax= 2.0d-02
!     ustmax= 4.0d-02
      ustmin= 5.0d-03

!dmr @-@ iceb0
! JONO JA heat2layers
! factor: iceberg volume to heat to temperature change of a m3 of sea water
! {(roiceb*Cfusheaticb)/(rho0*Cpo)} ((*ddtb(=86400)/86400)) = (see 604:)
!   voltohfl = 916.7*334000 / 1030*4002 ,, = 74.2778611
! MAB!!! : CHANGED RHOICEB TO 910.0 SO THAT IT IS THE SAME AS IN GRISLI!!!!
! leading to temperature change {K} per day
!   of layer k of volume unsdz(k)/area(i,j):
! voltohfl*dVol_icb(i,j,k)*unsdz(k)/area(i,j)
!! MAB as Ro was changed, also voltohfl has to be changed:
!!    voltohfl=74.2778611
      voltohfl=73.734977
!dmr @-@ iceb0

!
!nb-debut
!nb-fin

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  3) Initialization of some arrays.                                   |
!-----------------------------------------------------------------------
!
!dmr & mab: debug
      if (mod(istep,nstpyear).eq.1) then
        tot_Vicb_thersf=0d0
      endif
!dmr & mab: debug

      if (ntnow.eq.1) then

         do j=js1,js2
            do i=is1(j),is2(j)
               alct(i,j) = 0.0
            enddo
         enddo


      endif
      do 40 j=js1,js2
         do 30 i=is1(j),is2(j)
            fstrbq(i,j) = 0.0
            fscmbq(i,j) = 0.0
            ffltbq(i,j) = 0.0
            qfvbq(i,j)  = 0.0
            hgcol(i,j)  = 0.0
            hnbq(i,j)   = hnbq(i,j)*max(zero,
     &           sign(one,hnbq(i,j)-zeps2))
            dmnbq(i,j)  = 0.0
            dmgbq(i,j)  = 0.0
            dmgwi(i,j)  = 0.0
!age      ageg(i,j)  = (1.0-zindg)*ageg(i,j)+zindg*dtj
 30      continue
 40   continue
!
!
!--division by the number of thermodynamic steps
!--hnplbq and fwat are in kg per day and per m2
      raptime=ddtb/86400.
      do j=js1,js2
         do i=is1(j),is2(j)
#if ( ISOOCN >= 1 )
#if ( 0 )
!dmr --- Background checking of fluxes ...
          k = ieau18
          variso(:) = fwat(i,j,:)
          varisonm = "fwat_thsf0"

!          if ((variso(k)*variso(ieau)).LT.0.0d0) then
!           WRITE(*,*) "Inconsistent sign in variso: ", varisonm
!           WRITE(*,*) variso(k),variso(ieau), i, j
!           WRITE(*,*) zhnpbq(i,j,k), zhnpbq(i,j,ieau)
!           WRITE(*,*) dmnbq(i,j)
!          endif

          if (abs(variso(k)).gt.abs(variso(ieau)*1.0E-1)) then
            WRITE(*,*) "Too much isotopes in: ",varisonm, i,j
            WRITE(*,*) variso(k), variso(ieau)
            WRITE(*,*) delta(variso(:),k)
            WRITE(*,*) dmnbq(i,j)
            READ(*,*)
          endif
#endif
            zhnpbq(i,j,:) = raptime * hnplbq(i,j,:)
            zfwat(i,j,:)  = raptime * fwat(i,j,:)
#else
            zhnpbq(i,j) = raptime * hnplbq(i,j)
            zfwat(i,j)  = raptime * fwat(i,j)
#endif
         enddo
      enddo

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  4) Treatment of particular cases.                                   |
!-----------------------------------------------------------------------
!
      do 210 j=js1,js2
         do 200 i=is1(j),is2(j)
            zindb = tms(i,j,ks2)*(1.0-max(zero,
     &           sign(one,-(hnbq(i,j)+hgbq(i,j)))))
!
!---4.1. Snow is transformed into ice if the original
!        ice cover disappears.
!-----------------------------------------------------
!
            zindg      = tms(i,j,ks2)*max(zero,
     &           sign(one,-hgbq(i,j)))
            hgbq(i,j)  = hgbq(i,j)+zindg*rhon*hnbq(i,j)/rho0
            hnbq(i,j)  = (1.0-zindg)*hnbq(i,j)+zindg
     &           *hgbq(i,j)*(rho0-rhog)/rhon
!     dmgbq(i,j) = zindg*(1.0-albq(i,j))*rhog*hgbq(i,j)
            dmgwi(i,j) = zindg*(1.0-albq(i,j))*rhog*hgbq(i,j)
!
!-- 4.2. The lead fraction, albq, must be little than
!        or equal equal to amax (ice ridging).
!-----------------------------------------------------
!
            za    =   zindb*min(one,(1.0-albq(i,j))*uscomi)
            hnbq(i,j)   = hnbq(i,j)*za
            hgbq(i,j)   = hgbq(i,j)*za
            qstobq(i,j) = qstobq(i,j)*za
            albq(i,j)   = 1.0-zindb*(1.0-albq(i,j))/max(za,zeps1)
!age        agen(i,j)   = agen(i,j)*vnbq0(i,j)/
!age &                    max(dvn+vnbq0(i,j),zeps0)+dtj
!age        ageg(i,j)   = ageg(i,j)*vgbq0(i,j)/
!age &                    max(dvg+vgbq0(i,j),zeps0)+dtj
!age        agen(i,j)   = (1.0-zindn)*agen(i,j)
!age        ageg(i,j)   = (1.0-zindg)*ageg(i,j)
!
!-- 4.3. The in situ ice thickness, hgbq, must be equal to
!        or greater than hglim.
!----------------------------------------------------------
!
            zh          = max(one,zindb*hglim/max(hgbq(i,j),zeps1))
            hnbq(i,j)   = hnbq(i,j)*zh
            hgbq(i,j)   = hgbq(i,j)*zh
            qstobq(i,j) = qstobq(i,j)*zh
            albq(i,j)   = (albq(i,j)+(zh-1.0))/zh
 200     continue
 210  continue
!
!

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  5. Thermodynamics of sea ice.                                       |
!-----------------------------------------------------------------------
!
!--5.1. Partial computation of forcing for the thermo-
!       dynamic sea ice model.
!------------------------------------------------------
!
      do 230 j=js1,js2
         do 220 i=is1(j),is2(j)
!
            zindb      = tms(i,j,ks2)*(one-max(zero
     &           ,sign(one,-(hnbq(i,j)+hgbq(i,j)))))
            zinda(i,j) = 1.0-max(zero,sign(one,-(1.0-albq(i,j))))
            ain(i,j)   = albq(i,j)
!
!--Solar irradiance transmission at the mixed layer
!--bottom.
!
            thcm(i,j) = 1.-reslum(i,j,ks2)
     &           *dz(ks2)*(rho0*cpo)
!         thcm(i,j) = 0.0
!
!--Calculate fdtcn and qdtcn. Limitation of mixed
!--layer temperature.
!
!C        fdtcn(i,j) = .5*zindb*rho0*cpo*dz(ks2)*
!C   &                    (scal(i,j,ks2,1)-tfu(i,j))/dts(ks2)
            ustsg = max(min(sqrt(ust2s(i,j)),ustmax),ustmin)
            fdtcn(i,j) = zindb*rho0*cpo*0.006
     &           *ustsg*(scal(i,j,ks2,1)-tfu(i,j))
            qdtcn(i,j)  = zindb*fdtcn(i,j)*albq(i,j)*ddtb
!
!-- Snow accumulation.
!
#if ( ISOOCN >= 1 )
! dmr --- zhnpbq has dimension owatert
            zhnpbq(i,j,:) = zhnpbq(i,j,:)*rhoesn
#else
            zhnpbq(i,j) = zhnpbq(i,j)*rhoesn
#endif
!
!-- Partial computation of the lead energy budget and
!-- determination of qcmbq.
!
#if ( ISOOCN >= 1 )
! dmr --- zhnpbq has dimension owatert
            zfontn     = zhnpbq(i,j,ieau)*xln/ddtb
#else
            zfontn     = zhnpbq(i,j)*xln/ddtb
#endif
!     if (icoupl .eq. 0) then
!     zfnsol     = zemise*ratbqo(i,j)-zemise*stefan*
!     &                 scal(i,j,ks2,1)*scal(i,j,ks2,1)*
!     &                 scal(i,j,ks2,1)*scal(i,j,ks2,1)
!     &                 -fcscn(i,j)-flecn(i,j)
!     else
!     zfnsol     = ratbqo(i,j)-zemise*stefan*
!     &                 scal(i,j,ks2,1)*scal(i,j,ks2,1)*
!     &                 scal(i,j,ks2,1)*scal(i,j,ks2,1)
!     &                 -fcscn(i,j)-flecn(i,j)
!     endif
            zfnsol     = ratbqo(i,j)-fcscn(i,j)-flecn(i,j)
            qlbq(i,j)  = tms(i,j,ks2)*(fsolcn(i,j)*(1.0-thcm(i,j))
     &           +zfnsol+fdtcn(i,j)
     &           -zfontn+(1.0-zindb)*fsbbq(i,j))*albq(i,j)*ddtb
            fntlat      = 1.0-max(zero,sign(one,-qlbq(i,j)))
            pareff      = 1.0+(parlat-1.0)*zinda(i,j)*fntlat
            qlbsbq(i,j) = qlbq(i,j)*(1.0-pareff)
     &           /max(((one-albq(i,j))*ddtb),zeps0)
            qlbq(i,j)   = pareff*qlbq(i,j)
            qdtcn(i,j)  = pareff*qdtcn(i,j)
            qcmbq(i,j)  = rho0*cpo*dz(ks2)*(tfu(i,j)-scal(i,j,ks2,1))
     &           *(1-zinda(i,j))
!     &                 *.5*(1-zinda(i,j))
!
!--Calculate oceanic heat flux.
!
            fbbq(i,j)  = zindb*(fsbbq(i,j)/max((one-albq(i,j))
     &           ,zeps1)+fdtcn(i,j))
!
!--Store ice volume per unit area for computation of
!--daily thermodynamic ice production.
!--This variable is only needed for output.
!
            zhgbqp(i,j) = hgbq(i,j)*(1.0-albq(i,j))
!
 220     continue
 230  continue

!
!-- 5.2. Select icy points and fulfill arrays for the
!        vectorial grid.
!------------------------------------------------------
!
      nbpb = 0
      do 250 j=js1,js2
         do 240 i=is1(j),is2(j)
            if (albq(i,j).ge.1.0) go to 240
            nbpb      = nbpb+1
            npb(nbpb) = (j-1)*imax+i
 240     continue
 250  continue
!
!--If there is no ice, do nothing.
!
      if (nbpb.eq.0) go to 280
!
      call gather(nbpb,albqb,albq,npb)
      call gather(nbpb,hgbqb,hgbq,npb)
      call gather(nbpb,hnbqb,hnbq,npb)
#if ( ISOOCN >= 1 )
      call gather(nbpb,hnpbqb,zhnpbq(:,:,ieau),npb)
#else
      call gather(nbpb,hnpbqb,zhnpbq,npb)
#endif
      call gather(nbpb,fsolgb,fsolg,npb)
      call gather(nbpb,fbbqb,fbbq,npb)
      call gather(nbpb,thcmb,thcm,npb)
      call gather(nbpb,qlbqb,qlbq,npb)
      call gather(nbpb,qcmbqb,qcmbq,npb)
      call gather(nbpb,qstbqb,qstobq,npb)
      call gather(nbpb,tfub,tfu,npb)
      call gather(nbpb,tsb,ts,npb)
      call gather(nbpb,dmgbqb,dmgbq,npb)
      call gather(nbpb,dmgwib,dmgwi,npb)
      call gather(nbpb,psbqb,psbq,npb)
      call gather(nbpb,tabqb,tabq,npb)
      call gather(nbpb,qabqb,qabq,npb)
      call gather(nbpb,vabqb,vabq,npb)
      call gather(nbpb,ratbqb,ratbqg,npb)
      call gather(nbpb,qlbbqb,qlbsbq,npb)
      call gather(nbpb,cldqb,cloud,npb)
      do 260 k=1,nkb0
         call gather(nbpb,tbqb(1,k),tbq(1,1,k),npb)
 260  continue
!cp2  call gather(nbpb,fderb,fder,npb)
!
!--5.3. Call to ice growth routine.
!-----------------------------------
!
      call fontbc(1,nbpb)
!
!--5.4. Back to the geographic grid.
!-----------------------------------
!
      call scater(nbpb,albq,npb,albqb)
      call scater(nbpb,hnbq,npb,hnbqb)
      call scater(nbpb,hgbq,npb,hgbqb)
      call scater(nbpb,fscmbq,npb,fscbqb)
      call scater(nbpb,ffltbq,npb,fltbqb)
      call scater(nbpb,fstrbq,npb,fstbqb)
      call scater(nbpb,qlbq,npb,qlbqb)
      call scater(nbpb,qfvbq,npb,qfvbqb)
      call scater(nbpb,qstobq,npb,qstbqb)
      call scater(nbpb,ts,npb,tsb)
      call scater(nbpb,dmgbq,npb,dmgbqb)
      call scater(nbpb,dmgwi,npb,dmgwib)
      call scater(nbpb,dmnbq,npb,dmnbqb)
      do 270 k=1,nkb0
         call scater(nbpb,tbq(1,1,k),npb,tbqb(1,k))
270   continue
      call scater(nbpb,dvosbq,npb,dvsbqb)
      call scater(nbpb,dvobbq,npb,dvbbqb)
      call scater(nbpb,dvolbq,npb,dvlbqb)
      call scater(nbpb,dvonbq,npb,dvnbqb)
      call scater(nbpb,firg,npb,fratsb)
      call scater(nbpb,fcsg,npb,fcsb)
      call scater(nbpb,fleg,npb,fleb)
!
280   continue
!
!--5.5. Up-date sea ice thickness.
!---------------------------------
!
!
      do 300 j=js1,js2
         do 290 i=is1(j),is2(j)
            ifvt(i,j) = zinda(i,j)*max(zero,sign(one,-hgbq(i,j)))
            hgbq(i,j) = hgbq(i,j)*(one-max(zero,
     &           sign(one,-(1.0-albq(i,j)))))
 290     continue
 300  continue
!
!--Tricky trick: add 2 to albq in the Southern Hemisphere.
!
      do 307 j=js1,jeq-1
         do 305 i=is1(j),is2(j)
            albq(i,j) = albq(i,j)+2.0
 305     continue
 307  continue
!
!-- 5.6. Select points for lateral accretion and fulfill
!--      arrays for the vectorial grid.
!--------------------------------------------------------
!
! collection thickness of frazil ice (with option Cvhg)
!
! hgcol  : collection thickness of consolidated new ice
! hg0    : frazil ice thickness at the leading edge of
!          the lead
! alfagr : fraction of frazil ice in grease ice (the rest
!          is water)
! rhoo   : seawater density
! rhofr  : frazil ice density
! rhogr  : grease ice density
! vafrmn : minimum surface wind velocity for computation
!          of hgcol (in m/s)
! gamafr : ratio between grease ice and wind speeds
!
      alfagr = 0.35
      rhoo   = 1029.0
      rhofr  = 950.0
      rhogr  = rhoo+alfagr*(rhofr-rhoo)
      gred   = gpes*(rhoo-rhogr)/rhoo
      hg0    = 0.1
      fac1hg = sqrt(2.0*hg0/gred)*2.0/acos(-1.0)
      fac2hg = 0.25/gred
      gamafr = 0.06
      do j=js1,js2
         jp1   = j+1
         do i=is1(j),is2(j)
            if (qcmbq(i,j)-qlbq(i,j).gt.0.0) then
               ip1         = i+1
               tenagm      = sqrt( tenagx(i,j)*tenagx(i,j)
     &              +tenagy(i,j)*tenagy(i,j))
!     approx coupled
               vabq(i,j)   = tenagm/0.041
!
               dumfac      = gamafr*vabq(i,j)/max(1.0d-06,tenagm)
               vfrx        = dumfac*tenagx(i,j)
               vfry        = dumfac*tenagy(i,j)
               vgx         = (ug(i,j)+ug(ip1,j)+ug(i,jp1)+ug(ip1,jp1))*0.25
               vgy         = (vg(i,j)+vg(ip1,j)+vg(i,jp1)+vg(ip1,jp1))*0.25
               vrel        = sqrt((vfrx-vgx)**2+(vfry-vgy)**2)
               hgcol(i,j)  = hg0+(fac1hg+fac2hg*vrel)*vrel
          endif
        enddo
      enddo
!     write(mouchard_id,*) 'collec',int(vabq(50,2)),hgcol(50,2),
!    &             int(vabq(104,59)),hgcol(104,59)
!
      nbpac=0
      do 320 j=js1,js2
         do 310 i=is1(j),is2(j)
            if (qcmbq(i,j)-qlbq(i,j).le.0.0) go to 310
            nbpac       = nbpac+1
            npac(nbpac) = (j-1)*imax+i
 310     continue
 320  continue
!
!--If nbpac = 0, do nothing.
!
      if (nbpac.eq.0) go to 350
!
      call gather(nbpac,albqb,albq,npac)
      call gather(nbpac,hnbqb,hnbq,npac)
      call gather(nbpac,hgbqb,hgbq,npac)
      call gather(nbpac,qlbqb,qlbq,npac)
      call gather(nbpac,qcmbqb,qcmbq,npac)
      call gather(nbpac,qstbqb,qstobq,npac)
      call gather(nbpac,dmgbqb,dmgbq,npac)
      call gather(nbpac,dmgwib,dmgwi,npac)
      call gather(nbpac,tfub,tfu,npac)
      call gather(nbpac,tsb,ts,npac)
      call gather(nbpac,hgcolb,hgcol,npac)
      do 330 k=1,nkb0
         call gather(nbpac,tbqb(1,k),tbq(1,1,k),npac)
 330  continue
      call gather(nbpac,dvlbqb,dvolbq,npac)
!
!-- 5.7. Call lateral accretion routine.
!----------------------------------------
!
      call acrlbq(1,nbpac)
!
!-- 5.8. Back to the geographic grid.
!------------------------------------
!
      call scater(nbpac,albq,npac,albqb)
      call scater(nbpac,hnbq,npac,hnbqb)
      call scater(nbpac,hgbq,npac,hgbqb)
      call scater(nbpac,qstobq,npac,qstbqb)
      call scater(nbpac,ts,npac,tsb)
      call scater(nbpac,dmgbq,npac,dmgbqb)
      call scater(nbpac,dmgwi,npac,dmgwib)
      do 340 k=1,nkb0
         call scater(nbpac,tbq(1,1,k),npac,tbqb(1,k))
 340  continue
      call scater(nbpac,dvolbq,npac,dvlbqb)
!
 350  continue
!
      do 370 j=js1,js2
         do 360 i=is1(j),is2(j)
!
!--Tricky trick: recover albq values between
!--0 and 1 in the Southern Hemisphere.
!
            albq(i,j) = min(albq(i,j),abs(albq(i,j)-2.0))
!
! rate of lead creation/destruction due to thermodynamics.
!
            alct(i,j) = alct(i,j)+(albq(i,j)-ain(i,j))/ddtb
!
!-- 5.9. Daily thermodynamic ice production.
!--------------------------------------------
!
            hgbqp(i,j) = hgbq(i,j)*(1.0-albq(i,j))-zhgbqp(i,j)
     &           +hgbqp(i,j)
 360     continue
 370  continue
!
!-- 5.10. Compute ages.
!-----------------------
!
!age   do j=1,jmax
!age    do i=1,imax
!age      zindn      = max(zero,sign(one,-hnbq(i,j)))
!age      zindg      = max(zero,sign(one,-hgbq(i,j)))
!age      dvn        = max(zero,hnbq(i,j)*(1.0-albq(i,j))-vnbq0(i,j))
!age      dvg        = max(zero,hgbq(i,j)*(1.0-albq(i,j))-vgbq0(i,j))
!age      agen(i,j)  = agen(i,j)*vnbq0(i,j)/
!age &                 max(dvn+vnbq0(i,j),zeps0)+dtj
!age      ageg(i,j)  = ageg(i,j)*vgbq0(i,j)/
!age &                 max(dvg+vgbq0(i,j),zeps0)+dtj
!age      agen(i,j)  = (1.0-zindn)*agen(i,j)
!age      ageg(i,j)  = (1.0-zindg)*ageg(i,j)
!age    enddo
!age  enddo
!
!
!-- 5.10. Raccord cyclique
!--------------------------
!
!dmr ---      if (ltest.ge.1) then
#if ( L_TEST >= 1 )
!--raccord cyclique pour hgbq,albq,hnbq,ts,tbq,firg,fcsg,fleg,
!                     fsbbq,fdtcn,qstobq,scal:
         call raccord(hgbq(1,1),zero,1,8)
         call raccord(albq(1,1),zero,1,8)
         call raccord(alct(1,1),zero,1,8)
         call raccord(hnbq(1,1),zero,1,8)
         call raccord(ts(1,1),zero,1,8)
         call raccord(tbq(1,1,1),zero,nkb0,8)
         call raccord(firg(1,1),zero,1,8)
         call raccord(fcsg(1,1),zero,1,8)
         call raccord(fleg(1,1),zero,1,8)
         call raccord(fsbbq(1,1),zero,1,8)
         call raccord(qstobq(1,1),zero,1,8)
!age  call raccord(agen(1,1),zero,1,8)
!age  call raccord(ageg(1,1),zero,1,8)
         call raccord(scal(1,1,ks2,1),zero,1,8)
#endif
!dmr ---      endif

!
!-- 5.11. Sea ice/ocean interface.
!---------------------------------
!
!
!--CALCULATE fcm1, fcm2, fwat, AND fsbbq.
!
      do ns=1,nsmax
         if (scpme(ns).eq.spvr) then
            do j=js1,js2
               do i=is1(j),is2(j)
                  scs(i,j,ns) = scal(i,j,ks2,ns)
               enddo
            enddo
         else
            do j=js1,js2
               do i=is1(j),is2(j)
                  scs(i,j,ns) = scpme(ns)
               enddo
            enddo
         endif
      enddo

!-AM
! In order to correct for lack of consistency between salinity restoring
! fluxes defined in different routines restoring fluxes are set in this
! routine and transmitted through common array "frapsav"
!
!--borne inf. pour salinite (dans terme rappel sur la sal. pour FW Flx)
      salinf=5.0d0
!
      dzsdts=dz(ks2)/dts(ks2)
!-AM

#if ( F_PALAEO_FWF >= 1 && APPLY_UNCORFWF == 0 )
        corr_fluxinmetres(:,:)=0
#endif

      do 390 j=js1,js2
         do 380 i=is1(j),is2(j)
            salflux    =  34.7
!     salflux    =  scal(i,j,ks2,2)
            ia         =  1.0-max(zero,sign(one,-(1.0-albq(i,j))))
            iflt       =  zinda(i,j)*(1-ia)*(1-ifvt(i,j))
            ial        =  ifvt(i,j)*ia+(1-ifvt(i,j))*
     &           (one-max(zero,sign(one,
     &           albq(i,j)-ain(i,j))))
            fcm1(i,j)  =  ain(i,j)*fsolcn(i,j)+(1.-ain(i,j))*fstrbq(i,j)
            iadv       =  (1-ia)*zinda(i,j)
            zqfv       =  qfvbq(i,j)+qdtcn(i,j)
            fcm2(i,j)  = -fcm1(i,j)*(1-thcm(i,j))+
     &           iflt*((1-iadv)*ffltbq(i,j)+fscmbq(i,j))+
     &           (1-ia*(1-ial))*(ial*qcmbq(i,j)+
     &           (1-ial)*(qlbq(i,j)-iadv*zqfv))/ddtb
!C        zdtm       =  (zqfv+ffltbq(i,j)*ddtb*iflt)/
!C   &                  (rho0*cpo*dz(ks2))
!C        scal(i,j,ks2,1) =  scal(i,j,ks2,1)+iadv*zdtm
            zfdtm      = (zqfv/ddtb+ffltbq(i,j)*iflt)*iadv
            fdtcn(i,j) = fdtcn(i,j)-zfdtm
            fsbbq(i,j) =  (1-(ifvt(i,j)+iflt))*fscmbq(i,j)
!     fd0      prs        =  (zfwat(i,j)-zhnpbq(i,j)*
!fd0 &                  (1.-ain(i,j))*rhon)*salflux
! dmgbq variation of ice mass [?] > [kg/m2 step]
! saltflux = 34.7 [psu]
! sglace Salinity of Ice = 6 [psu]
!     [psu kg / m2 step] = [kg/m2 step]* [psu]
            fons       =  dmgbq(i,j)*(salflux-sglace)
!fd0 &                 -dmgwi(i,j)*sglace
!fd0 &                 +dmnbq(i,j)*salflux
!fd0      pmess      =  (prs-fons)/ddtb-fevabq(i,j)
!fd0 &                  *salflux*albq(i,j)

!-----|--1--------2---------3---------4---------5---------6---------7-|
!dmr !!! fevabq is a remnant of old code and is zero in coupled mode
!mab: pme is in kg/m2/s
!-----|--1--------2---------3---------4---------5---------6---------7-|

#if ( ISOOCN >= 1 )
! dmr --- pme has dimension owatert
          pme(:)        =  ! ### see comment above ... -fevabq(i,j)*albq(i,j)
     &            (zfwat(i,j,:)-zhnpbq(i,j,:)*(1.-ain(i,j))*rhon
     &                  -dmnbq(i,j)*rsmow(:)
     &                  )/ddtb

#else
! fevabq Evaporation flux [?] >> [kg /m2 step] > 0
! albq Leads Fraction [1]
! fwat Freshwater flux [kg /m2 step] < ec_co2oc ? sumofwf
! zhnpbq [?] >> [m / step]
! rhon Density of snow [kg/m3]
! dmnbq Variation of snow mass [?] >> [m / step]
! ddtb 86400 [s/step] day!!
! [kg /m2 s] =  [0] + ( [kg /m2 step] - [m / step]*[1]*[kg/m3] - [m /step] )/[s/step]
            pme        = -fevabq(i,j)*albq(i,j)
     &           +(zfwat(i,j)-zhnpbq(i,j)*(1.-ain(i,j))*rhon
     &           -dmnbq(i,j))
     &           /ddtb

#endif

!dmr !!! frapsav is a remnant of old code and is zero in coupled mode
!
!-AM: frapsav
            frapsav(i,j)=rappes(i,j,0)*
     &           (1.D0-scalr(i,j,ks2,2)/max(scal(i,j,ks2,2),salinf))
!
!dmr !!! frapsav is a remnant of old code and is zero in coupled mode
!mab: DUM is in m
#if ( ISOOCN >= 1 )
!dmr --- all variables with dimension owatert ...
          DUM(:) = (pme(:)/rho0
! ### see comment above     & + frapsav(i,j) * dzsdts
     &) *ddtb
          phiss(i,j,0) = -DUM(ieau)
     &                   +phiss(i,j,0)

! dmr --- To add the fwatconseau flux with conservation, just use:
! dmr ---   DUM(:) = DUM(:) + fwatconseau(i,j)*ddtb*isotopic_content_to_determine

#else
! [m / step] =  [kg /m2 s] / [kg/m3] + [1][m/s] ) * [s/step]
            DUM = (pme/rho0 + frapsav(i,j) * dzsdts) *ddtb

#if ( F_PALAEO_FWF >= 1 && APPLY_UNCORFWF == 0 )
! dmr --- To add the fwatconseau flux with conservation, just use:
      DUM = DUM + fwatconseau(i,j)*ddtb
      corr_fluxinmetres(i,j) = corr_fluxinmetres(i,j)
     >                         +fwatconseau(i,j)*ddtb
      sum_flux_out(:,:)=corr_fluxinmetres(:,:)*area(:,:)/ddtb/1e6 ! in Sv
#endif

! phiss0 is water mass flux m for one timestep (day)
! [m / step]

            phiss(i,j,0) = -DUM
     &           +phiss(i,j,0)
#endif

!-dmr&afq --- update-2019-04-17
!- #if ( CONSEAU == 1 )
!-
!-!dmr&afq [NOTA] : fwatconseau is in [kg/(m2 step)]
!-
!-        ! phiss0 is water mass flux m for one timestep (day)
!-        ! [m / step]
!-            phiss(i,j,0) = phiss(i,j,0)-fwatconseau(i,j)/rho0
!-                     ! [m / step] = [m3 / s] / [ m2 ] * [ s / step ]
!-     &               -calgrisCLIO(i,j)*ddtb/area(i,j)
!-#endif

!
!dbug     if (abs(scal(i,j,ks2,2)).lt.epsil) then
!dbug      write(clio3_out_id,*) 'ARRET : thersf, scal(i,j,ks2,2) too small'
!dbug      write(clio3_out_id,*) 'scal(i=',i,',j=',j,')=',scal(i,j,ks2,2)
!dbug      stop
!dbug     endif

! water surface / thickness of the first layer
#if ( ISOOCN >= 1 )
! dmr --- all variables with dimension owatert
          watsur = -DUM(ieau) * unsdz(ks2) ! water flux at the surface that is used for other tracers
          zisowatsur(:) = -DUM(:) * unsdz(ks2) ! same for isotopes ...
#else
!     [1/step] = [m/step] [1/m]
            watsur = -DUM * unsdz(ks2)
#endif

!-----|--1--------2---------3---------4---------5---------6---------7-|
! dmr&afq --- Update of the Temperature flux ...
!-----|--1--------2---------3---------4---------5---------6---------7-|

!
!C        phiss(i,j,1)  =  -(fcm1(i,j)+fcm2(i,j))
!C   &                      *ddtb*unsdz(ks2)/(rho0*cpo)
! temperature flux [K] for one timestep
! phiss modifies scal in scale.f
! water added has same temp as surface layer (scs1)
! fcm2 : Non-solar flux at ocean surface [W/m2]>[kg/s3] (icecom)

! [K/step] = ( )*[s/step]*[1/m]/([kg/m3]*[J / kg K]) + [1/step]*[K]+[K/step]
!            [kg/s3] * [ s /m step ]* ( [ m s2 K /kg ] ) +
!                      [K/step]     +
            phiss(i,j,1) =  -(fcm2(i,j)-fdtcn(i,j))
     &           *ddtb*unsdz(ks2)/(rho0*cpo)
     &           +watsur*scs(i,j,1)
     &           +phiss(i,j,1)
!fd0      phiss(i,j,2) =  pmess
!fd0 &                      *ddtb*unsdz(ks2)/rho0
!fd0 &                      +phiss(i,j,2)

!-dmr&afq --- update-2019-04-17
!-#if ( CONSEAU == 1 )
!-!dmr&afq same as for other type of water(s) ... (sparkling, still ...)
!-            phiss(i,j,1) = phiss(i,j,1)
!-     >           -fwatconseau(i,j)/rho0*unsdz(ks2)*scs(i,j,1)
!-     >           -calgrisCLIO(i,j)*ddtb/area(i,j)*unsdz(ks2)*scs(i,j,1)
!-#endif

!-----|--1--------2---------3---------4---------5---------6---------7-|
! dmr&afq --- Update of the Salinity flux ...
!-----|--1--------2---------3---------4---------5---------6---------7-|

! salt flux [?] for sea ice, because it is floating ontop, no away pushing
! we increase the salinity, by dumping salt (negative salt on sea ice), when it melts negative salt is freed. for transport
! >> scs is zero, because no salt in freshwater. so no impact from watsur (but possible if scs2 not zero)
! fons :

! [psu/step] =  [psu kg/m2 step]*[1/m]*[m3/kg] + [1/step]*[psu] + [psu/step]
            phiss(i,j,2) = -fons*unsdz(ks2)/rho0
     &           +watsur*scs(i,j,2)
     &           +phiss(i,j,2)

!-dmr&afq --- update-2019-04-17
!-#if ( CONSEAU == 1 )
!-!dmr&afq same as for other type of water(s) ... (sparkling, still ... and non salted!)
!-            phiss(i,j,2) = phiss(i,j,2)
!-     >           -fwatconseau(i,j)/rho0*unsdz(ks2)*scs(i,j,2)
!-     >           -calgrisCLIO(i,j)*ddtb/area(i,j)*unsdz(ks2)*scs(i,j,2)
!-#endif

!nb debut
!nb          fluxbrines(i,j) = -fons*unsdz(ks2)/rho0
!     [psu/step] = [psu kg/m2 step]*[1/m]*[m3/kg]
            fluxbrines(i,j) = fons*unsdz(ks2)/rho0 !flux de sel spoitif
!           write(mouchard_id,*) 'fluxbrines thersf',fluxbrines
!nb fin
!cfc      phiss(i,j,3) = watsur*scs(i,j,3)
!cfc &                   +phiss(i,j,3)
!cfc      phiss(i,j,4) = watsur*scs(i,j,4)
!cfc &                   +phiss(i,j,4)

#if ( ISOOCN >= 1 )
! dmr --- Problem here phiss has dimension 0:nsmax !!!
! dmr      but zisowatsur is owatert
! dmr --- Hence conversion needed.


        do iz = owisostrt, owisostop ! array nsmax ...

           call oc2atwisoindx(iz,atmisoindx)
           phiss(i,j,iz) =  phiss(i,j,iz) + zisowatsur(atmisoindx)

        enddo

! dmr ISO2        DO iz = ieau17, ieaud
! dmr ISO2
! dmr ISO2           phiss(i,j,iz) = phiss(i,j,iz) + zisowatsur(iz)
! dmr ISO2
! dmr ISO2        ENDDO
#endif

#if ( CONSEAU == 1 && HEATCALV == 1 && ICEBERG < 2)

! Heat needed to HEAT (a kg of) SEA WATER by dT degrees {K}:
!  Heat{J} = dT{K} * water{kg} * Cpo{J/KgK} -->
!  dT{K} = Heat{J} / (  water{m3} * rhowater{kg/m3} * Cpo{J/KgK} )
!  dT{K} = Heat{J} / ( (area(i,j)/unsdz(k)){m3} * rhowater{kg/m3} * Cpo{J/KgK} )
! ICEBERG MELTING HEAT per GRID (i,j), per timestep DAY:
!  LH_icb = fonte_icb(i,j){m3/day} * roiceb{kg/m3} * latent heat(=Chicb{J/kg})
! or PER LAYER k:
!  LH_icb = dVol_icb(i,j,k){m3/day} * roiceb{kg/m3} * latent heat(=Chicb{J/kg})
! so the sea water cooling dT, by the iceberg melt is:
!
!  dT{K} = voltohfl * dVol_icb(i,j,k) / (area(i,j)/unsdz(k)) per day
         unsdz_som = 0.0
         DO kounter=kheatlimit,(kmax-1)
            klayer=ks2-kounter+1
            unsdz_som = unsdz_som + unsdz(klayer)*tms(i,j,klayer)
         ENDDO

         IF ( unsdz_som .gt. 0. ) THEN      ! afq
         DO kounter=kheatlimit,(kmax-1)
            klayer=ks2-kounter+1
            scal(i,j,klayer,1)=scal(i,j,klayer,1)-
     &          voltohfl*calgrisCLIO(i,j)*ddtb*
     &            (unsdz(klayer)/area(i,j))*(unsdz(klayer)/unsdz_som)
     &            *tms(i,j,klayer)
         ENDDO
         ENDIF
#endif



 380     continue
 390  continue


#if ( F_PALAEO_FWF >= 1 && APPLY_UNCORFWF == 0 )
      write(mouchard_id,'(A,F12.8,A)') 'CorrFluxInMeters: Z',
     &     (sum(corr_fluxinmetres*area))/ddtb/1E6,' [Sv]'
#endif

!
      do 430 k=11,ks2
         do 420 j=js1,js2
            do 410 i=is1(j),is2(j)
               phivs(i,j,k,1)=-reslum(i,j,k)*fcm1(i,j)
     &              *ddtb
     &              +phivs(i,j,k,1)
 410        continue
 420     continue
 430  continue

! geothermic flow
      do 450 j=js1,js2
         do 440 i=is1(j),is2(j)
            aflux=max(zero,dfloat(sign(1,11-kfs(i,j))) )
            phivs(i,j,kfs(i,j),1)=-aflux*bheat*
     &           dts(kfs(i,j))*unsdz(kfs(i,j))/(rho0*cpo)
!          if (i.eq.2) then
!             write(mouchard_id,*) 'bot',j,kfs(i,j),aflux,
!    &              phivs(i,j,kfs(i,j),1),tms(i,j,kfs(i,j))
!          endif
 440     continue
 450  continue



#if ( ICEBERG > 0 )
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! JONO heat of icebergs was added to phiss1 here.
! phiss1 is phiss(i,j,1) is heat flux to atmosphere at ocean surface.
! old style flux-parameterization put in if (flgicb) else construction.
! Presumably, dimension of phiss(i,j,1) is {K/s}
!	degrees of seawater heating per second in a grid.
! fonte_icb(i,j)={m3/day} (unless NICs trafo in iceb_moy)
!---
! JA heat2layers Clarified this bit and separating ocean layers 'k'
!    (eg. separating unsdz(ks2) from the factor voltohfl)
! to explain daily 'fluxes' to
!     sea water temperature (scalar): scal(i,j,k,1)
!     sea water salinity            : scal(i,j,k,2)
!
! Heat needed to HEAT (a kg of) SEA WATER by dT degrees {K}:
!  Heat{J} = dT{K} * water{kg} * Cpo{J/KgK} -->
!  dT{K} = Heat{J} / (  water{m3} * rhowater{kg/m3} * Cpo{J/KgK} )
!  dT{K} = Heat{J} / ( (area(i,j)/unsdz(k)){m3} * rhowater{kg/m3} * Cpo{J/KgK} )
! ICEBERG MELTING HEAT per GRID (i,j), per timestep DAY:
!  LH_icb = fonte_icb(i,j){m3/day} * roiceb{kg/m3} * latent heat(=Chicb{J/kg})
! or PER LAYER k:
!  LH_icb = dVol_icb(i,j,k){m3/day} * roiceb{kg/m3} * latent heat(=Chicb{J/kg})
! so the sea water cooling dT, by the iceberg melt is:
!
!  dT{K} = voltohfl * dVol_icb(i,j,k) / (area(i,j)/unsdz(k)) per day
!
! with:
!  voltohfl = rho_ice*Chicb / rho_water*Cpo
!           = 74.277861069465267
!    with: rho_ice=916.7{kg/m3}; Chicb=334000{J/kg};
! MAB!!! : CHANGED RHOICEB TO 910.0 SO THAT IT IS THE SAME AS IN GRISLI!!!!
!         rho_water=1030.{kg/m3}; cpo =4002.{J/kgK}
!  unsdz = 1/dz
! ps. times ddtb/86400 to get dT per thermodyn timestep ddtb(=86400)
!---
! checked: is the flux per grid or per m2?per m2 cause areicebn+=area(i,j)*tms
!	and area(i,j) is (i presume) in m2.
! note trafo from surface layer ks2=19(or20) to phiss surface =1
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      do i=1,imax
         do j=1,jmax

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! loop over ocean layers from surface (=ks2=20) downwards
            do kounter=1,(kmax-1)
               klayer=ks2-kounter+1

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! Heat from icebergs to different ocean layers
! putting the dVol_icb into heat_icb for outputting and then into scalar temperature (see explanation above)

               heat_icb(i,j,(klayer)) = -voltohfl*dVol_icb(i,j,(klayer))*
     &                 (unsdz(klayer)/area(i,j))
               scal(i,j,klayer,1)=scal(i,j,klayer,1)+heat_icb(i,j,(klayer))
!
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! Freshwater from icebergs to different ocean layers
! salinity_old=salt/layer_volume --> salt=salinity_old*layer_volume
! salinity_new=salt/(layer_volume+meltwater_volume)
!  --->>
! salinity_new=(salinity_old*layer_volume)/(layer_volume+meltwater_volume)
! with:
!   meltwater_volume=dVol_icb*rho_ice/rho_fresh_water
!   rho_ice=910
! MAB!!! : CHANGED RHOICEB TO 910.0 SO THAT IT IS THE SAME AS IN GRISLI!!!!
!   rho_fresh_water=~1000
! NB. shut of fwf in ec_co2oc when using this!
!
#if (WATGRISCONS == 0)
! Put iceberg freshwater volume [m3[ in fw_icb [m] for outputting, then use to change salinity
               fw_icb(i,j,klayer) = area(i,j)*0.910*dVol_icb(i,j,(klayer)) 
               scal(i,j,klayer,2)=scal(i,j,klayer,2)*(area(i,j)*dz(klayer))
     &                               /( (area(i,j)*dz(klayer))+
     &                               0.910*dVol_icb(i,j,(klayer)) )

#else
!### to test what happens if freshwater coming from icebergs is put
!### into the uncorrected flux...
!### m3/day * factor from ice to water in the upper most layer
               if (dVol_icb(i,j,klayer) .gt. 0d0) then
                  mass_day_uncor(i,j,ks2) = mass_day_uncor(i,j,ks2)
     &               + dVol_icb(i,j,(klayer))* 0.910
               endif
#endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!
               tot_Vicb_thersf=tot_Vicb_thersf+dVol_icb(i,j,klayer)

	       if(dVol_icb(i,j,klayer).lt.(-0.0000000001)) then
	          write(mouchard_id,*)'REDALERT, dVol negative '
     &               ,'dVol_icb',dVol_icb(i,j,(klayer))
     &               ,'i,j',i,j
     &               ,'klayer:',klayer
               endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

            enddo ! end do layers loop
         enddo ! end do j
      enddo ! end do i

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      write(mouchard_id,*)
     & 'thersf: total freshwater volume from icebergs',istep,tot_Vicb_thersf

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!PB       if(flag_snow_to_flux) then
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! this is the old style parameterisation:
!PB          fnor=(ficebergn/areicebn)
!PB     &         *ddtb*unsdz(ks2)/(rho0*cpo*86400.)
!PB          do i=iicebern1,iicebern2
!PB             do j=jicebern1,jicebern2
!PB                phiss(i,j,1) =  phiss(i,j,1) +fnor*tms(i,j,ks2)
!PB             enddo
!PB          enddo
!PB          fsud=(ficebergs/areicebs)
!PB     &         *ddtb*unsdz(ks2)/(rho0*cpo*86400.)
!PB          do i=iicebers1,iicebers2
!PB             do j=jicebers1,jicebers2
!PB                phiss(i,j,1) =  phiss(i,j,1) +fsud*tms(i,j,ks2)
!PB             enddo
!PB          enddo
!PB       endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! else (flgicb.FALSE. --> NO DYNAMIC ICEBERGS)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
#else /*  On ICEBERG > 0 */
!dmr [DEPRECATED]      else

#if (CALVFLUX == 2)
!mab: how latent heat related to freshwaterfluxes is considerd!
#if (HEATFWF == 1)
! how latent heat related to freshwaterfluxes is considerd
!##### calving put into the ocean as direct freshwaterflux!!!!! same way #####!
!##### treated as icebergs, therefore put into the latent heat!!
!mab: heat for melting snow is taken into account, yet, when applying iceberg
!mab: calving as direct freshwater flux it is put into the scal variable
!mab: as is done with the icebergs, therefore the old parameterization is "ignored"
!mab: otherwise heat is taken out twice!
      DO kounter=1,(kmax-1)
         klayer=ks2-kounter+1
         DO i=1,imax
            DO j=1,jmax
               scal(i,j,klayer,1)=scal(i,j,klayer,1)-
!mab: mass_day in kg/day has to be transformed into m3/day
     &               voltohfl*(mass_day(i,j)*(1./910.))*
     &               (unsdz(klayer)/area(i,j))
            ENDDO
         ENDDO
       ENDDO
#elif (HEATFWF == 0)
!mab: if FWFflux is applied instead of icebergs but effect should be similar, then
!mab: old style parameterization for take-up of latent heat (homogenously around
!mab: Greenland according to excess snow) is kept
! Iceberg melting
         fnor=(ficebergn/areicebn)
     &        *ddtb*unsdz(ks2)/(rho0*cpo*86400.)
         do i=iicebern1,iicebern2
            do j=jicebern1,jicebern2
               phiss(i,j,1) =  phiss(i,j,1) +fnor*tms(i,j,ks2)
            enddo
         enddo
         fsud=(ficebergs/areicebs)
!dmr @-@ iceb0
     &        *ddtb*unsdz(ks2)/(rho0*cpo*86400.)
         do i=iicebers1,iicebers2
            do j=jicebers1,jicebers2
               phiss(i,j,1) =  phiss(i,j,1) +fsud*tms(i,j,ks2)
            enddo
         enddo
#else
! just a comment for checking
#endif  /* HEATFWF */
#elif (CALVFLUX == 0)
! this is the original old style stuff:
! ficebergn is total heat of northern dsnow melt per day
! clio/sources/defgrid.f:1437:        areicebn=areicebn+area(i,j)*tms(i,j,ks2)
! Iceberg melting
         fnor=(ficebergn/areicebn)
     &        *ddtb*unsdz(ks2)/(rho0*cpo*86400.)
         do i=iicebern1,iicebern2
            do j=jicebern1,jicebern2
               phiss(i,j,1) =  phiss(i,j,1) +fnor*tms(i,j,ks2)
            enddo
         enddo
         fsud=(ficebergs/areicebs)
!dmr @-@ iceb0
     &        *ddtb*unsdz(ks2)/(rho0*cpo*86400.)
         do i=iicebers1,iicebers2
            do j=jicebers1,jicebers2
               phiss(i,j,1) =  phiss(i,j,1) +fsud*tms(i,j,ks2)
            enddo
         enddo
#else
! just a comment for checking
#endif /* on CALVFLUX */
!dmr @-@ iceb0 End of the if(flgicb) loop
!dmr [DEPRECATED]      endif ! on
!dmr @-@ iceb0 End of the if(flgicb) loop
#endif /*  on ICEBERG > 0  */


#if ( ISM == 1 )
      if (flgism) then
         do i=1,imax
            do j=1,jmax
               phiss(i,j,1) =  phiss(i,j,1) + (heaticbism(i,j)*
     &              tms(i,j,ks2))*ddtb*unsdz(ks2)/(rho0*cpo)
            enddo
         enddo
      endif
#endif

!     write(mouchard_id,*) 'delta t ns', fnor, fsud
!     write(mouchard_id,*) 'flux iceb', ficebergn/(areicebn*86400),
!    &                   ficebergs/(areicebs*86400)
!mab: comment for mouchard
!      write(mouchard_id,*) 'flux iceb2',fnor,fsud
#if ( ISOOCN >= 1 )


! dmr --- ISO2
! dmr ---   zflux0s has dimensions owatert
! dmr ---   phiss   has dimensions 0:nsmax

! mean value of phiss(i,j,iz)
        zflux0s(:)=0.0


        iz = 2 !salt !!!

        do j=js1,js2
         do i=is1(j),is2(j)
           zflux0s(iz)=zflux0s(iz)+aire(i,j)*tms(i,j,ks2)
     >             *(phiss(i,j,iz))
         enddo
        enddo

! dmr --- Loop should be on indexes for phiss: owisostrt -> owisostop
        do iz = owisostrt, owisostop

        call oc2atwisoindx(iz,atmisoindx)

        do j=js1,js2
         do i=is1(j),is2(j)
           zflux0s(atmisoindx)=zflux0s(atmisoindx)+aire(i,j)*tms(i,j,ks2)
     >             *(phiss(i,j,iz))
         enddo
        enddo

        enddo

        zflux0s(:)=zflux0s(:)/zurfow

!PB        write(mouchard_id,*) 'zflux salt',zflux0s(2),phiss(15,15,2)

! mean value of phiss2 set to zero on average

        dum(:)=zfluxms(:)*vcor

        iz = 2 ! salt

        do j=js1,js2
         do i=is1(j),is2(j)
           phiss(i,j,iz)=phiss(i,j,iz)-dum(iz)*tms(i,j,ks2)
         enddo
        enddo

!dmr --- Correction as well for the isotopes ...
! dmr --- Loop should be on indexes for phiss: owisostrt -> owisostop

        do iz = owisostrt, owisostop

        call oc2atwisoindx(iz,atmisoindx)

       do j=js1,js2
         do i=is1(j),is2(j)
           phiss(i,j,iz)=phiss(i,j,iz)-dum(atmisoindx)*tms(i,j,ks2)
         enddo
       enddo

      enddo

!dmr --- Calcul du delta 18 global de l ocean
      totisotop=0.0
      volocean=0.0
      do i=1,imax
        do j=1,jmax
          do k=1,kmax
          volocean=volocean+aire(i,j)*dz(k)*tms(i,j,k)
          totisotop=totisotop+scal(i,j,k,ocnw18)*aire(i,j)*dz(k)
     & *tms(i,j,k)
          enddo
        enddo
      enddo
      totisotop=totisotop/volocean
      IF (MOD(istep,360).EQ.0) then
       WRITE(*,*) "Mean delta18O ocean: ", istep,
     &(totisotop/rsmow(ieau18)-1.0)*1000.0

      ENDIF

!#if ( OOISO == 1 )
!      totisotop2=totisotop
!      d18Omar = ((totisotop/rsmow(ieau18))-1.0)*1000.0
!      print *, 'd18Omar', d18Omar
!#endif

#else

! calculate the mean value and substract ( conserve it globally)
! mean value of phiss(i,j,2)
! phiss(i,j,2) : [psu/step]
!
      zflux0s=0.0
!-dmr&afq --- update-2019-04-17
!-#if ( CONSEAU == 1 )
!-      scs_mean = 0.0
!-#endif
      do j=js1,js2
         do i=is1(j),is2(j)
            zflux0s=zflux0s+aire(i,j)*tms(i,j,ks2)
     >           *(phiss(i,j,2))
!-dmr&afq --- update-2019-04-17
!-#if ( CONSEAU == 1 )
!-            scs_mean = scs_mean + scs(i,j,2)
!-#endif
         enddo
      enddo

! zflux0s same unit as phiss(i,j,2) but spread over the planet [m2/m2]
! zurfow is locally the area of the ocean ...
      zflux0s=zflux0s/zurfow
!mab: comment for mouchard
!      write(mouchard_id,*) 'zflux salt',zflux0s,phiss(15,15,2),zurfow
! mean value of phiss2 set to zero on average
      dum=zfluxms*vcor
!-dmr&afq --- update-2019-04-17
!-#if ( CONSEAU == 1 )
!-      dum = dum - trendconseau/rho0*unsdz(ks2)*scs_mean
!-#endif
      do j=js1,js2
         do i=is1(j),is2(j)
           phiss(i,j,2)=phiss(i,j,2)-dum*tms(i,j,ks2)
         enddo
      enddo
#endif


! Iceshelf melting

#if ( LGMSWITCH == 0 )
!     gammat is in m/s
      gammat=1E-4
#if ( BATHY > 0 )
      gammat=0E-4
#endif

#else
      gammat=0E-4
#endif

      fluxmean=0.0
      zcums=0.0
      zcums2=0.0
      do j=2,8
         do i=is1(j),is2(j)
            tmean=0.0
            smean=0.0
            dztot=0.0
            do k=9,14
               tmean=tmean+scal(i,j,k,1)*dz(k)*tms(i,j,k)
               smean=smean+scal(i,j,k,2)*dz(k)*tms(i,j,k)
               dztot=dztot+dz(k)*tms(i,j,k)
            enddo
            dztot=max(dztot,1D-5)
            tmean=tmean/dztot*tms(i,j,ks2)
            smean=smean/dztot*tms(i,j,ks2)
            zdepthi=200.0
            tfreezloc= -0.0575*(smean)
     &           + 1.710523e-3*(smean)**1.5
     &           - 2.154996e-4*(smean)
     &           *(smean)
     &           - 7.53e-4*zdepthi+273.15

!     fluxt is the total amount of heat available to melt the ice shelf (W)
            fluxt=rho0*cpo*gammat*(tmean-tfreezloc)*tmics(i,j)
!     if (tmics(i,j).gt.0.0) then
!     write(mouchard_id,*) i,j,tmean,(tmean-tfreezloc)
!     endif
!     write(mouchard_id,*) fluxt,i,j,tmics(i,j)
!     write(mouchard_id,*) tmean,tfreezloc,smean
            phivssup=fluxt*ddtb/(rho0*cpo*dztot*area(i,j))
            phivssups=fluxt*ddtb/(xlg*dztot*area(i,j))*salflux
!     phivssups=0.0
            do k=9,14
               phivs(i,j,k,1)=phivs(i,j,k,1)+phivssup*tms(i,j,k)
               phivs(i,j,k,2)=phivssups*tms(i,j,k)
            enddo
            fluxmean=fluxmean+fluxt
            zcums=zcums+fluxt/xlg*salflux
         enddo
      enddo
      toticesm=fluxmean/xlg*360*86400
!     redistribution of the flux
      ficesh=(fluxmean/areicebs)
     &     *ddtb*unsdz(ks2)/(rho0*cpo)
      ficess=(zcums/areicebs)
     &     *ddtb*unsdz(ks2)
!     &           *ddtb*unsdz(ks2)/rho0
!     ficess=0.0

      do i=iicebers1,iicebers2
         do j=jicebers1,jicebers2
            phiss(i,j,1) =  phiss(i,j,1)-ficesh*tms(i,j,ks2)
            phiss(i,j,2) =  phiss(i,j,2)-ficess*tms(i,j,ks2)
         enddo
      enddo
#if ( ISM == 1 )
      if (flgism) then
         qclim=qclim+(fluxmean/360.)
      endif
#endif
!mab: comment for mouchard
!      write(mouchard_id,*) 'ice-shelf',ficesh,ficess
!mab: end of comment

#if ( SHELFMELT == 1 )

      gammatg=1E-4
      fmeltg=15E-3 !7.6E-3  !1.6E-3 Goelzer et al. CPD 2016 AppendixA (from 1.7 to 7.6)
      fluxmean=0.0
      zcums=0.0
      zcums2=0.0
      zdepthi=200.0
      do j=1,jmax
        do i=1,imax
            do k=1,kmax
               if (tms(i,j,k).gt.0.) then
                  tmean=scal(i,j,k,1)*tms(i,j,k)
                  smean=scal(i,j,k,2)*tms(i,j,k)
                  tfreezloc= -0.0575*(smean)
     &                 + 1.710523e-3*(smean)**1.5
     &                 - 2.154996e-4*(smean)
     &                 *(smean)
     &                 - 7.53e-4*zdepthi+273.15
!     fluxt is the total amount of heat available to melt the ice shelf (W)
                  fluxtg=max(rho0*cpo*gammatg*(tmean-tfreezloc),0.d0)
                  bmshelf(i,j,k)=fluxtg*fmeltg/xlg*86400. ! / oceanic dt (day)
               else
                  bmshelf(i,j,k)=100.
               endif
            enddo
         enddo
      enddo

      do j=1,ubound(bmshelf_clionord,dim=2)
         do i=1,ubound(bmshelf_clionord,dim=1)
            ! we sum the daily values
            ! bmshelf_clio is reset in GRISLI
            ! taking into account the decoupling factor
            bmshelf_clionord(i,j,:) =  bmshelf_clionord(i,j,:) +
     &        bmshelf(icliog(i,j,1),jcliog(i,j,1),:)/(dble(timCplGRISday)/360.)
         enddo
      enddo

#endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! VERIF FLUXES
!     zsurfoz=0.0
!     zsurfzz=0.0
!     znsolo=0.0
!     zsolo=0.0
!     znsolo2=0.0
!     zsolo2=0.
!     zlph0=0.0
!     zlph2=0.0
!     zvolg=0.0
!     zvoln=0.0
!     zfdt=0.0
!     do j=js1,js2
!     do i=is1(j),is2(j)
!     zsurfoz=zsurfoz+aire(i,j)*ain(i,j)
!     zsurfzz=zsurfzz+aire(i,j)
!     znsolo=znsolo+ratbqo(i,j)*ain(i,j)
!     &                 *aire(i,j)
!     zsolo=zsolo+fsolcn(i,j)*ain(i,j)*aire(i,j)
!         zsolo2=zsolo2+fcm1(i,j)*aire(i,j)
!         znsolo2=znsolo2+fcm2(i,j)*aire(i,j)
!         zvoln=zvoln+(1.0-albq(i,j))*aire(i,j)*hnbq(i,j)
!    &                              *0.33d0*tms(i,j,ks2)
!         zvolg=zvolg+(1.0-albq(i,j))*aire(i,j)
!    &                     *hgbq(i,j)*0.9*tms(i,j,ks2)
!         zlph0=zlph0+phiss(i,j,0)*aire(i,j)
!    &             *360*tms(i,j,ks2)
!    &             *rho0/1000.0
!         zlph2=zlph2+phiss(i,j,2)*aire(i,j)
!    &             *360/(unsdz(ks2)*salflux)*tms(i,j,ks2)
!    &             *rho0/1000.0
!         zfdt=zfdt+fdtcn(i,j)*aire(i,j)
!        enddo
!       enddo
!       write(93,'(4F11.5)') zsolo/zsurfoz,znsolo/zsurfoz,
!    &     zsolo2/zsurfzz,znsolo2/zsurfzz
!a    &    ,zfdt/zsurfzz


!-dmr&afq --- update-2019-04-17
!- #if ( CALVFLUX == 2 )
!- !##### calving put into the ocean as direct freshwaterflux!!!!! same way #####!
!- !##### treated as icebergs, therefore put into the salinity (yet not into ####!
!- !##### latent heat!!
!-        DO i=1,imax
!-           DO j=1,jmax
!-                   scal(i,j,ks2,2)=scal(i,j,ks2,2)*(
!-      &                 (area(i,j)*dz(ks2))/
!-      &                 ( (area(i,j)*dz(ks2))+
!- !mab: mass_day in kg/day has to be transformed into m3/day
!-      &                 0.910*(mass_day(i,j)*(1./910d0) ) ) )
!-           ENDDO
!-        ENDDO
!- #endif

#if ( APPLY_UNCORFWF == 1 ) /* Apply uncorrected fwf to the ocean */
!dmr --- init of the flux
      given_fluxinmetres(:,:) = 0.d0
#endif

#if ( UNCORFLUX == 1 )
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!MB Icesheet melting as an uncorrected freshwater flux to the ocean surface
!MB Laurentide Icesheet melting
!   after D. Roche and H. Renssen
! Calculate the area and the meltflux



! ----------
!dmr&afq --- Partie pour ajouter des flux d'eau douce en anomalie
!        --- Input : fluxarea, rudiflux
!        --- Output : fluxinmeters
! ----------


#define HEINRICH_TYP 0

#if ( HEINRICH_TYP > 0 )
      fluxarea=0.0d0
      fluxinmeters=0.0

      if ((iyear.GT.100).AND.(iyear.LT.400)) then
         rudiflux = 0.3
      else
         rudiflux = 0.0
      endif
#endif

#if ( HEINRICH_TYP == 1 ) /* Ruddimann Belt forcing cf. Roche et al. 2010 */

      j = 44
      do i=100,108
        fluxarea=fluxarea+(area(i,j)*tms(i,j,ks2))
      enddo

      do j = 42,43
        do i=99,108
          fluxarea=fluxarea+(area(i,j)*tms(i,j,ks2))
        enddo
      enddo

! // flux calculation now //

      j = 44
      do i=100,108
!     [m/step] = [m3/s]/[m2]*[s/step]
            fluxinmeters(i,j) =
     &           ((rudiflux*1.E6/fluxarea)*ddtb)
      enddo

      do j = 42,43
        do i=99,108
!     [m/step] = [m3/s]/[m2]*[s/step]
            fluxinmeters(i,j) =
     &           ((rudiflux*1.E6/fluxarea)*ddtb)
         enddo
      enddo

#elif ( HEINRICH_TYP == 2 ) /* Labrador Sea forcing cf. Roche et al. 2010 */

      do j = 48,51
        do i=99,100
          fluxarea=fluxarea+(area(i,j)*tms(i,j,ks2))
        enddo
      enddo

! // flux calculation now //

      do j = 48,51
        do i=99,100

!     [m/step] = [m3/s]/[m2]*[s/step]
            fluxinmeters(i,j) =
     &           ((rudiflux*1.E6/fluxarea)*ddtb)
         enddo
      enddo


#elif ( HEINRICH_TYP == 3 ) /* Fennoscandian Ice Sheet margin cf. Roche et al. 2010 */

        i = 107
        do j=50,53
          fluxarea=fluxarea+(area(i,j)*tms(i,j,ks2))
        enddo

        i = 108
        do j=47,54
          fluxarea=fluxarea+(area(i,j)*tms(i,j,ks2))
        enddo

        i = 109
        do j=47,53
          fluxarea=fluxarea+(area(i,j)*tms(i,j,ks2))
        enddo

! // flux calculation now //

        i = 107
        do j=50,53
!     [m/step] = [m3/s]/[m2]*[s/step]
            fluxinmeters(i,j) =
     &           ((rudiflux*1.E6/fluxarea)*ddtb)
        enddo

        i = 108
        do j=47,54
!     [m/step] = [m3/s]/[m2]*[s/step]
            fluxinmeters(i,j) =
     &           ((rudiflux*1.E6/fluxarea)*ddtb)
        enddo

        i = 109
        do j=47,53
!     [m/step] = [m3/s]/[m2]*[s/step]
            fluxinmeters(i,j) =
     &           ((rudiflux*1.E6/fluxarea)*ddtb)
        enddo

#endif

#define EIGHT_POINT_FIVE 0
#if ( EIGHT_POINT_FIVE == 1 )
! ---hli set freshwater for 8.5ka state from Zhang et al., 2016
! ---hli, line(1654-1670) for Hudson Strait+River
      lawrflux = 0.05   ! St. Lawrence Outlet
      hudsonflux = 0.03  ! Hudson Strait+River
      waisflux = 0.0    ! West Antarc. Ice Sheet
c     HUDSON Strait+River
      do i=99,100
         do j=48,51
            fluxarea=fluxarea+(area(i,j)*tms(i,j,ks2))
         enddo
      enddo
      do i=99,100
         do j=48,51
c     [m/step] = [m3/s]/[m2]*[s/step]
            fluxinmeters(i,j) = fluxinmeters(i,j)+
     &           ((hudsonflux*1.E6/fluxarea)*ddtb)
!         print*,'i=',i,'j=',j,fluxinmeters(i,j)
         enddo
      enddo

! ---hli,line(1673-1686) for St. Lawrence Outlet
      fluxarea=0.0
      do i=97,97
         do j=44,44
            fluxarea=fluxarea+(area(i,j)*tms(i,j,ks2))
         enddo
      enddo
      do i=97,97
         do j=44,44
c     [m/step] = [m3/s]/[m2]*[s/step]
            fluxinmeters(i,j) = fluxinmeters(i,j)+
     &           ((lawrflux*1.E6/fluxarea)*ddtb)
!             print*,'i=',i,'j=',j,fluxinmeters(i,j)
         enddo
      enddo
#endif

! ----------
!dmr&afq --- F_PALAEO + UNCORFLUX = hosing scenario coherent avec le forcage transitoire (2017/04/04)
! ----------
#if ( F_PALAEO >=1 && UNCORFLUX >= 1 )
!     [m/step] = [m3/s]/[m2]*[s/step]

!dmr&afq --- run is an equilibrium, should add the current year to the actual date to allow equilibrium hosing scenarios
       if ( equil_flag.GT.0) then
           call update_fwf(celest_year-(iyear-1),flux_pourCLIO)
           !write(*,*) "Current date in thersf.f", celest_year-iyear
       else
           call update_fwf(celest_year,flux_pourCLIO)
           write(mouchard_id,*) "Current date in thersf.f", celest_year
       endif

       fluxinmeters(:,:) = flux_pourCLIO(:,:) * ddtb
#endif


!--- dmr&afq [TO BE DELETED]
! ----------
!dmr&afq --- UNCORRUNOFF => ajouter le runoff comme flux non compensé ...
! ----------

!- #if ( UNCORRUNOFF == 1 )
!- !MB   Routing of landmodel runoff as uncorrected freshwater flux!
!- !     reset fluxinmeters (remove this to combine)
!-       fluxinmeters = 0.0
!- !
!- !      write(*,*) sum(fluxinmeters*area)/ddtb/1E6,
!- !     &    sum(fwatgis*area)/ddtb/1E6
!-
!-       WHERE(abs(fwatgis).GT.0.0d0)
!- !     [m/step] = [kg/m2 step] [1] / [kg/m3]
!-          fluxinmeters(:,:) = fluxinmeters(:,:) +
!-      &              fwatgis(:,:)*tms(:,:,ks2)/rho0
!-       ENDWHERE
!-
!- #endif

!     debuging
      oldphiss=phiss

#endif /* on UNCORFLUX == 1, above */

#if ( APPLY_UNCORFWF == 1 ) /* Apply uncorrected fwf to the ocean */
!dmr --- added to provide a mean to combine different uncorrected fluxes
      given_fluxinmetres(:,:) = given_fluxinmetres(:,:)
     &        + fluxinmeters(:,:)

#if ( CONSEAU == 1 ) /*   Water conservation between ISM (GRISLI) and CLIO   */

      given_fluxinmetres(:,:) = given_fluxinmetres(:,:)
     &        + (fwatconseau(:,:) + calgrisCLIOms(:,:))*ddtb
      sum_flux_out(:,:)=given_fluxinmetres(:,:)*area(:,:)/ddtb/1e6 ! in Sv

#endif /* On CONSEAU */

#if ( F_PALAEO_FWF >= 1 ) /*   Water conservation between forced ISM scenario */
      given_fluxinmetres(:,:) = given_fluxinmetres(:,:)
     &                        + fwatconseau(:,:)*ddtb
      sum_flux_out(:,:)=given_fluxinmetres(:,:)*area(:,:)/ddtb/1e6 ! in Sv
#endif

#endif /* on APPLY_UNCORFWF == 1 */


!      write(*,*) sum(given_fluxinmetres*area)/ddtb/1E6

#if ( APPLY_UNCORFWF == 1 ) /* Apply uncorrected fwf to the ocean */
!MB Call freshwater flux function
!dmr --- Modif out of the previous UNCORFLUX version
!dmr --- All uncorrected fluxes are assumed to be in fluxinmeter
      call notcompedmeltflux(given_fluxinmetres,imax,jmax,kmax,ks2,phiss,
     &     unsdz,scs,w_uncor)
!     w_uncor will be applied in conti3d to w, after mean w is calculated!

#if ( F_PALAEO_FWF >= 1 || CONSEAU == 1 )
      write(mouchard_id,'(A,F12.8,A)') 'FluxInMeters: Z',
     &     (sum(given_fluxinmetres*area))/ddtb/1E6,' [Sv]'
#endif

#endif

#if ( UNCORFLUX == 1 )

!MB debugging infos
!           [m/step]/[s/step] = [m/s]
      write(mouchard_id,'(A,F12.8,A)') 'W_UNCOR: Z',sum(w_uncor)/dts(ks2),
     &     ' [m/s]'
      write(mouchard_id,'(A,F12.8,A)') 'W_UNCOR: Z',
     &     (sum(w_uncor*area))/ddtb/1E6,' [Sv]'
      write(mouchard_id,'(A,F12.8,A)') 'FluxInMeters: Z',
     &     (sum(given_fluxinmetres*area))/ddtb/1E6,' [Sv]'
      write(mouchard_id,'(A,2F18.8,A)') 'Phiss1 Z',
     &     sum(phiss(:,:,1)),sum(oldphiss(:,:,1)),'[K]'
      write(mouchard_id,'(A,F12.8,A)') 'Phiss2 Z',
     &     sum(phiss(:,:,2)-oldphiss(:,:,2)),'[psu]'
!      write(*,'(i4,122F15.6)') (j,tms(:,j,ks2),j=1,jmax)
!      write(*,'(i4,122F15.6)') (j,scs(:,j,1),j=1,jmax)

!MB End of freshwater flux
#endif

      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine thersf -

      CONTAINS

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      SUBROUTINE notcompedmeltflux(waterheight,im,jm,km,isurf,phissloc,
     &     unsdzed,scsloc,f_uncor)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!     Authors: D. Roche and M. Blaschek
!     Last Change: 13.2.2012, M. Blaschek
!     Apply a not compensated freshwater flux to phiss and w (conti3d.f)
!     INPUTS:
!       waterheight(imax,jmax)   freshwater [m/step]
!     FIXED IN/OUTPUTS:
!       im,jm,km,isurf,unsdzed  >> imax,jmax,kmax,ks2,unsdz
!       phissloc                >> phiss
!       scsloc                  >> scs
!     OUTPUTS:
!       f_uncor                 >> w_uncor
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
#if ( ISOOCN >= 1 || ISOBERG == 2 )
       USE iso_param_mod, ONLY: neauiso,rsmow, ieau17, ieau18, ieaud
       USE isoatm_mod, ONLY: dlbmini
       use para0_mod, only: ocnw17, ocnw18, ocnw2h
#endif

      INTEGER, INTENT(IN) :: im, jm, km, isurf
!mab: for including isotopes in calving flux and runoff
#if ( ISOOCN >= 1 )
      INTEGER :: iz
#endif

      REAL, DIMENSION(im,jm), INTENT(IN) :: waterheight
      REAL, DIMENSION(im,jm,km), INTENT(IN) :: scsloc
      REAL, DIMENSION(km), INTENT(IN) :: unsdzed
      REAL, DIMENSION(im,jm,0:km), INTENT(INOUT) :: phissloc
      REAL, DIMENSION(im,jm), INTENT(OUT) :: f_uncor

      f_uncor(:,:) = 0.0
      WHERE (ABS(waterheight(:,:)).GT.0.0d0)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!     water mass flux in m for one timestep >> apply to w
         f_uncor(:,:) = -waterheight(:,:)
!         f_uncor(:,:) = 0.0  ! no freshwater added
!         phissloc(:,:,0) = -waterheight(:,:)  ! will be corrected in conti3d by mean w
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!     temperature flux (temperature of freshwater)
!     Should it be colder than scs(:,:,1) = temperature at surface
!     [K/step] = [m/step][1/m]*[K] + [K/step]
         phissloc(:,:,1) = -waterheight(:,:)*unsdzed(isurf)
     &        *scsloc(:,:,1)+phissloc(:,:,1)
!     Scale the impact of freshwater to temperature with 0 deg Celsius sea water (artificially)
!         phissloc(:,:,1) = -waterheight(:,:)*unsdzed(isurf)
!     &        *273.15+phissloc(:,:,1)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!     salt flux
!     [psu/step] = [m/step][1/m]*[psu] + [psu/step]
         phissloc(:,:,2) = -waterheight(:,:)*unsdzed(isurf)
     &        *scsloc(:,:,2)+phissloc(:,:,2)

#if ( ISOOCN >= 1 || ISOBERG == 2 ) /* is water isotopes, then specify a water isotopic content */
         phissloc(:,:,ocnw18) = -waterheight(:,:)*unsdzed(isurf)
     &        *(dlbmini(ieau18) + 1.00 ) *rsmow(ieau18)+phissloc(:,:,ocnw18)
         phissloc(:,:,ocnw17) = -waterheight(:,:)*unsdzed(isurf)
     &        *(dlbmini(ieau17) + 1.00 ) *rsmow(ieau17)+phissloc(:,:,ocnw17)
         phissloc(:,:,ocnw2h) = -waterheight(:,:)*unsdzed(isurf)
     &        *(dlbmini(ieaud) + 1.00 ) *rsmow(ieaud)+phissloc(:,:,ocnw2h)
#endif

! Method 2, ala H. Goosse and P. Mathiot as negativ salt flux
! dz, tms and scal are missing in the call!!
!         phissloc(:,:,2) = phissloc(:,:,2) + scal(:,:,isurf,2)
!     &       *tms(:,:,isurf)*( 1 - dz(isurf) / (dz(isurf)+waterheight(:,:)) )
      ENDWHERE
      END SUBROUTINE notcompedmeltflux

      end


