!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:50 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:50 CET 2009

      SUBROUTINE outnet(ja,xjour,xjour1)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  This routine computes the average of some variables and write it
!  on the ouput files.
!  ATTENTION cette routine devrait etre valable meme si le pas de temps est
!  different de 1 jour. Pour plus de generalite il faudrait introduire la
!  variable xjours indiquant le jour au pas de temps suivant pour savoir si
!  on restera dans le meme mois ou la meme annee. De plus il faudra calculer
!  njucum par un cumul lors de chaque passage.
!  modif : 24/12/99
!
!  This is a major revision of outave that implements netcdf-output of
!  monthly and/or annual means following as closely as possible the
!  "NetCDF Climate and Forecast (CF) Metadata Conventions" (refer to
!  http://www.cgd.ucar.edu/cms/eaton/cf-metadata/index.html for more
!  information). The update retains full compatibility with the former
!  output to cresum.dat etc. Outside outave.f only minor changes have
!  been applied to defgrid.f, where true gridpoint longitudes and
!  latitudes of the combined WW- and AA-grids are computed, and to
!  bloc0.com for some additional arrays that needed to be defined.
!
!  This work has been funded by the Deutsche ForschungsGemeinschaft
!  (DFG) within the research framework "Passagen" (please visit
!  http://www.passagen.uni-kiel.de). Extensions of the code have
!  been enclosed by comment lines

!DFG start ...
!DFG end
!
!  Even though this revision of the code is a bit lengthy due to the
!  administrative nature of netcdf, it is fairly well commented and
!  should be easily understood and flexible enough for modifications
!  according to other users needs.
!
!  Oct-2002, Peter Herrmann (ph@phsck.de)
!
!---
! Ccpl [Ccp0] => ligne specifique a la version avec [sans] couplage .
!---
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!! START_OF_USE_SECTION

      use const_mod

#if ( ISOOCN >= 1 )
      USE iso_param_mod, ONLY : ieau
      use para0_mod, only: ocnw17, ocnw18, ocnw2h, owatert
#endif

#if (CALVFLUX > 0 )
      USE iceberg_mod, ONLY : iceflux_in
#endif

#if ( PATH >= 1 )
      use path_mod, only: scalstart_path, particles_fluxes_field, ncaco3
     >                  , npoc
#endif


      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod
      use densit_mod
      use thermo_mod

      use isoslope_mod
      use ice_mod
      use iceberg_mod
      use dynami_mod

! --- BdB 05-2019: added time variable for writing output
      use comemic_mod, only: time_in_years, init_year_clio

      use mchd99_mod, only: nn99

      use stream_mod, only: uuu, vvv, jmtt, kmtt

#if ( F_PALAEO_FWF == 1 || CONSEAU == 1 )
      use update_clio_bathy_tools, only: sum_flux_out
#endif

      use newunit_clio_mod, only: clio3_out_id, mouchard_id

!! END_OF_USE_SECTION

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!! START_OF_INCLUDE_SECTION
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

#include "type.com"

!! END_OF_INCLUDE_SECTION

!- nn99=2 => ecritures auxiliaires sur fichier "mouchard", unit=99
c~       common / mchd99 / nn99
#if ( ISOOCN >= 1 )
      real*8 zfluxmts(owatert)
#endif
      common / zfluc / zfluxmt,zfluxmts
!PB#if ( CALVFLUX == 0 )
      real*8:: calvCLIO(imax,jmax) = 0.0
!PB#endif
!--Output variables and arrays.
!
      integer njcum(0:12)
      integer iwl(30),jwl(30)
! mab: change noumax (and noumax2) to 50 because two lines were added to the output.param
! mab: defining the iceberg distribution (line nr=10), that means to all the other lines
! mab: >= 10 +1 is added
! mab: in outave.f noumax2 = 30
! mab: added IMfice and lmmfice for the output
! mab: added IMfdust and lmmfdist for the output
      parameter (noumax=50,noumax2=50)
!dmr @-@ iceb0
      real chk_Vicb_outave(imax,jmax)
!dmr @-@ iceb0
! Choice noumax :19,35,44,50, noumax2 depending of the choice
      real*4 cmoyan(imax,jmax,noumax),cmoymo(imax,jmax,noumax),
     &       cmoymap(imax,jmax,noumax2,12)
!DFG start real*8 define for lahey compatibility
!DFG original define is plain real
      real cmulti(noumax),cadd(noumax)
!DFG end

      real*4 sort1(kmax),sort2(kmax),sort3(kmax)
      integer ntypou(noumax),nmm(noumax),nmal(noumax),nma(noumax),
     &        ncor(noumax2),nmap(12)
      character*60 titn(noumax)

!dmr --- Ajout d'une declaration propre de spval ..
      REAL*8, PARAMETER :: spval = -1.0d+32

!dmr --- Ajout d'une declaration propre de total_vol ...
!dmr @-@ iceb0
      REAL :: total_vol = 0.0
!dmr @-@ iceb0

! streamFUNCTION part
! options: nstreamupdate=0(no output), 1(monthly) or 2(annually)
!DFG start
! In order to interfere as little as possible with streamfunc.f,
! we here let nstreamout unchanged to yield monthly means computed
! by that original routine. Those transports are exported from
! streamfunc.f through common streamflu. We will then add them up
! here for annual means... Please leave nstreamout untouched!
!DFG end
      parameter (nstreamout=1)

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      common /totalanu/ umm(imax,jmax,kmax),vmm(imax,jmax,kmax),
     &           tmm(imax,jmax,kmax),smm(imax,jmax,kmax),
#if ( ISOOCN >= 1 )
     &  oo17m(imax,jmax,kmax),
     &  oo18m(imax,jmax,kmax),
     &  oohdm(imax,jmax,kmax),
#endif
#if ( OCYCC == 1 )
     &  odocm(imax,jmax,kmax),
     &  odocsm(imax,jmax,kmax),
     &  odicm(imax,jmax,kmax),
     &  opo4m(imax,jmax,kmax),
#if ( OXNITREUX == 1 )
     &  on2om(imax,jmax,kmax),
#else
     &  ono3m(imax,jmax,kmax),
#endif
     &  oalkm(imax,jmax,kmax),
     &  oo2m(imax,jmax,kmax),
     &  oc13m(imax,jmax,kmax),
     &  odoc13m(imax,jmax,kmax),
     &  odocs13m(imax,jmax,kmax),
     &  oc14m(imax,jmax,kmax),
     &  aqpco2m(imax,jmax),
     &  tpp_mam(imax,jmax),
     &  caco3mm(imax,jmax),
#if ( PATH >= 1 )
     &  papartm(imax,jmax,kmax),
     &  padissm(imax,jmax,kmax),
     &  thpartm(imax,jmax,kmax),
     &  thdissm(imax,jmax,kmax),
     &  pflxcam(imax,jmax,kmax),
     &  pflxpom(imax,jmax,kmax),
#endif
#if ( CORAL == 1)
     &  coraream(imax,jmax,kmax),
     &  corprodm(imax,jmax,kmax),
     &  cormassm(imax,jmax,kmax),
     &  omegam(imax,jmax,kmax),
     &  oco3m(imax,jmax,kmax),
     &  tau_bleachm(imax,jmax,kmax),
     &  DHWm(imax,jmax,kmax),
     &  PHm(imax,jmax,kmax),
#endif
#if ( OOISO == 1 )
     &  oo2_2m(imax,jmax,kmax),
     &  oo2_3m(imax,jmax,kmax),
     &  oo2_4m(imax,jmax,kmax),
#endif
#endif
#if ( F_PALAEO_FWF == 1 || CONSEAU == 1 )
     &  outFWFm(imax,jmax),
#endif
     &  uma(imax,jmax,kmax),vma(imax,jmax,kmax),
     &  tma(imax,jmax,kmax),sma(imax,jmax,kmax)
#if ( ISOOCN >= 1 )
     & ,oo17a(imax,jmax,kmax),
     &  oo18a(imax,jmax,kmax),
     &  oohda(imax,jmax,kmax)
#endif
#if ( OCYCC == 1 )
     &  ,odoca(imax,jmax,kmax),
     &  odocsa(imax,jmax,kmax),
     &  odica(imax,jmax,kmax),
     &  opo4a(imax,jmax,kmax),
#if ( OXNITREUX == 1 )
     &  on2oa(imax,jmax,kmax),
#else
     &  ono3a(imax,jmax,kmax),
#endif
     &  oalka(imax,jmax,kmax),
     &  oo2a(imax,jmax,kmax),
     &  oc13a(imax,jmax,kmax),
     &  odoc13a(imax,jmax,kmax),
     &  odocs13a(imax,jmax,kmax),
     &  oc14a(imax,jmax,kmax),
     &  aqpco2a(imax,jmax),
     &  tpp_maa(imax,jmax),
     &  caco3ma(imax,jmax)
#if ( PATH >= 1 )
     &  ,paparta(imax,jmax,kmax),
     &  padissa(imax,jmax,kmax),
     &  thparta(imax,jmax,kmax),
     &  thdissa(imax,jmax,kmax),
     &  pflxcaa(imax,jmax,kmax),
     &  pflxpoa(imax,jmax,kmax)
#endif
#if ( CORAL == 1 )
     &  ,corareaa(imax,jmax,kmax),
     &  corproda(imax,jmax,kmax),
     &  cormassa(imax,jmax,kmax),
     &  omegaa(imax,jmax,kmax),
     &  oco3a(imax,jmax,kmax),
     &  tau_bleacha(imax,jmax,kmax),
     &  DHWa(imax,jmax,kmax),
     &  PHa(imax,jmax,kmax)
#endif
#if ( OOISO == 1 )
     &  ,oo2_2a(imax,jmax,kmax),
     &  oo2_3a(imax,jmax,kmax),
     &  oo2_4a(imax,jmax,kmax)
#endif
#endif
#if ( F_PALAEO_FWF == 1 || CONSEAU == 1 )
     &  ,outFWFa(imax,jmax)
#endif
!DFG start additonal fields for monthly and annual means of w
      common /totalanu/ wma(imax,jmax,kmax+1), wmm(imax,jmax,kmax+1)
!
      common/moyout/ cmoyan,cmoymo,cmoymap
      common/parout/ titn,npack,klevm(4,2),klevv,njm
      common/parout2/ noumef,nmm,nmal,nma,nmmef,nmaef,nmalef
      common/parout3/ ncor,nmap,ntypou,jmois,nptj,iwl,jwl
      common/parout4/ cmulti,cadd,jmoisi,jmoifl
      common/parout5/ znumas(imax,jmax,4),znumav(imax,jmax,5)

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|


!DFG start
!
! declarations for netcdf output
!
! ncidm, mcida: integer IDs for netcdf files
! lcdf[mon,ann]: logical switch for [monthly,annual] mean netcdf dumps
! l[v,m,a]cdf: logical switches for variables to dump to netcdf file
!
! copied from streamfunc.f: jmtt, kmtt, uuu, vvv
!
!     integer ldim, ncidm, ncida, write_2d, write_3d, status,
      integer write_2d, write_3d, status, write_2dnm,
     &        start(4), count(4)
c~       integer jmtt, kmtt
c~       parameter (jmtt=57, kmtt=kmax+1)

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! dmr --- As of July, 12th, 2013, the cumbersome variable declarations
!          are moved to the following include file (for consistency).
!         Beware that ldim has been moved as well
#include "netcdf_out.com"
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      dimension f2d(imax,jmax),etamo(imax,jmax),ubmo(imax,jmax),
     &          vbmo(imax,jmax),etaan(imax,jmax),uban(imax,jmax),
     &          vban(imax,jmax)
c~       real*4 uuu(jmtt,0:kmtt,0:nbsmax),vvv(jmtt,0:nbsmax+4)
c~       common /streamflu/ uuu,vvv
!DFG start define meridional arrays of minimum bathymetry
!    for masking of streamFUNCTION in netcdf output
!
      integer kfloor(jmtt,0:nbsmax)
      common/streambath/ kfloor
!DFG end
      real*8 moc_m(0:nbsmax,jmtt,0:kmtt),mht_m(0:nbsmax,jmtt),
     &       mst_m(0:nbsmax,jmtt),
     &       moc_a(0:nbsmax,jmtt,0:kmtt),mht_a(0:nbsmax,jmtt),
     &       mst_a(0:nbsmax,jmtt)
      character*200 avgname
!dmr BUGfix      character*32 avgname
!DFG end

!dmr --- for newunits
       integer, save :: cresujldat_id ! was "31"
!dmr [UNUSED]       integer, save :: cresujgdat_id ! was "32"
       integer, save :: cresumdat_id  ! was "33"
       integer, save :: cresuadat_id  ! was "34"
       integer, save :: cresaldat_id  ! was "35"
!dmr [UNUSED]       integer, save :: crestmdat_id  ! was "36"
!dmr [UNUSED]       integer, save :: crestadat_id  ! was "37"
       integer, save :: correcwdat_id ! was 42

       integer :: cpointjdat_id, outputparam_id, netcdfoutparam_id

!
!cp0  data njcum/0,31,59,90,120,151,181,212,243,273,304,334,365/
      data njcum/0,30,60,90,120,150,180,210,240,270,300,330,360/
      ijour=int(xjour)
      ijour1=int(xjour1)
!     jours=xjour+dts(ks2)/86400
      undim6=1.D-6
!     write(clio3_out_id,*) 'debut outave.f',ja,jmois,ijour
!      write(clio3_out_id,*) "Numeers : ", numit, nstart, nlast
      if(numit.eq.nstart) then

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1) INITIALIZATIONS.                                                 |
!-----------------------------------------------------------------------
        write(*,*) 'ldim= ', ldim

        zfluxmt=0.0
        zfluxmts=0.0

#if (CALVFLUX > 0 )
!PB
        calvCLIO = iceflux_in
#endif

!--1.1 OPEN FILES FOR INITIALIZATION AND FOR OUTPUTS.
!----------------------------------------------------
!
!--cresujl.dat : CONTAINS ICE LOCAL DAILY OUTPUTS.
!--cresujg.dat : CONTAINS ICE GLOBAL DAILY OUTPUTS.
!--cresum.dat  : CONTAINS MONTHLY OUTPUTS.
!--cresua.dat  : CONTAINS YEARLY OUTPUTS.
!
        open(newunit=cresujldat_id,file='outputdata/ocean/cresujl.dat'
     >      ,form='unformatted')
        nmtot=365
        write(cresujldat_id) nmtot
!dmr [UNUSED]        open(newunit=cresujgdat_id,file='outputdata/ocean/cresujg.dat'
!dmr [UNUSED]     >      ,form='unformatted')
        open(newunit=cresumdat_id,file='outputdata/ocean/cresum.dat'
     >      ,form='unformatted')
        open(newunit=cresuadat_id,file='outputdata/ocean/cresua.dat'
     >       ,form='unformatted')
        open(newunit=cresaldat_id,file='outputdata/ocean/cresal.dat'
     >      ,form='unformatted')
!dmr [UNUSED]        open(newunit=crestmdat_id,file='outputdata/ocean/crestm.dat'
!dmr [UNUSED]     >      ,form='unformatted')
!dmr [UNUSED]        open(newunit=crestadat_id,file='outputdata/ocean/cresta.dat'
!dmr [UNUSED]     >      ,form='unformatted')

        open(newunit=correcwdat_id,file='correcw.dat'
     >      ,form='formatted')


!
!--1.2. SELECTION OF POINTS FOR LOCAL DAILY OUTPUTS.
!---------------------------------------------------
!
!--cpointj.dat: CONTAINS POINTS FOR LOCAL DAILY OUTPUTS.
!
        nptj=0
        open(newunit=cpointjdat_id,file='cpointj.dat',status='old'
     &      ,err=6)
        read(cpointjdat_id,*) nptj
        do npw=1,nptj
          read(cpointjdat_id,*) iwl(npw),jwl(npw)
        enddo
        close(cpointjdat_id)
6      continue
!
!--1.3 Computation of the month.
!---------------------------------
!
        jmois=1
        do ii=1,ijour
           if (ijour.gt.njcum(jmois)) jmois=jmois+1
        enddo
        jmoisi=jmois
        ijouri=ijour-njcum(jmois-1)-1
        jmoifl=0
        jjoufl=0
        njm=0
!
!--1.4. INITIALIZATION OF CUMULATIVE ARRAYS
!-------------------------------------------
!


! streamFUNCTION part
        if (nstreamout.ne.0) call streamfunc(-1,0.d0,nstreamout)

        do n=1,noumax
          do j=1,jmax
            do i=1,imax
              cmoyan(i,j,n)=0.0
              cmoymo(i,j,n)=0.0
            enddo
          enddo
         enddo
        do n=1,noumax2
          do j=1,jmax
            do i=1,imax
              do jj=1,12
                cmoymap(i,j,n,jj)=0.0
              enddo
            enddo
          enddo
        enddo
        do k=1,kmax
          do j=1,jmax
            do i=1,imax
              umm(i,j,k)= 0.0
              vmm(i,j,k)= 0.0
              tmm(i,j,k)= 0.0
              smm(i,j,k)= 0.0
              uma(i,j,k)= 0.0
              vma(i,j,k)= 0.0
              tma(i,j,k)= 0.0
              sma(i,j,k)= 0.0
#if ( ISOOCN >= 1 )
              oo17m(i,j,k)= 0.0
              oo18m(i,j,k)= 0.0
              oohdm(i,j,k)= 0.0
              oo17a(i,j,k)= 0.0
              oo18a(i,j,k)= 0.0
              oohda(i,j,k)= 0.0
#endif
#if ( OCYCC == 1 )
              odocm(i,j,k)= 0.0
              odocsm(i,j,k)= 0.0
              odicm(i,j,k)= 0.0
              opo4m(i,j,k)= 0.0
#if ( OXNITREUX == 1 )
              on2om(i,j,k)= 0.0
#else
              ono3m(i,j,k)= 0.0
#endif
              oalkm(i,j,k)= 0.0
              oo2m(i,j,k)= 0.0
              oc13m(i,j,k)= 0.0
              odoc13m(i,j,k)= 0.0
              odocs13m(i,j,k)= 0.0
              oc14m(i,j,k)= 0.0
#if ( PATH >= 1 )
              papartm(i,j,k) = 0.0
              padissm(i,j,k) = 0.0
              thpartm(i,j,k) = 0.0
              thdissm(i,j,k) = 0.0
              pflxcam(i,j,k) = 0.0
              pflxpom(i,j,k) = 0.0
#endif
#if ( CORAL == 1 )
              coraream(i,j,k)= 0.0
              corprodm(i,j,k)= 0.0
              cormassm(i,j,k)= 0.0
              omegam(i,j,k)= 0.0
              oco3m(i,j,k)= 0.0
              tau_bleachm(i,j,k)=0.0
              DHWm(i,j,k)=0.0
              PHm(i,j,k)=0.0
#endif
#if ( OOISO == 1 )
              oo2_2m(i,j,k)= 0.0
              oo2_3m(i,j,k)= 0.0
              oo2_4m(i,j,k)= 0.0
#endif
!dmr --- End of monthly arrays

              odoca(i,j,k)= 0.0
              odocsa(i,j,k)= 0.0
              odoca(i,j,k)= 0.0
              odocsa(i,j,k)= 0.0
              odica(i,j,k)= 0.0
              opo4a(i,j,k)= 0.0
#if ( OXNITREUX == 1 )
              on2oa(i,j,k)= 0.0
#else
              ono3a(i,j,k)= 0.0
#endif
              oalka(i,j,k)= 0.0
              oo2a(i,j,k)= 0.0
              oc13a(i,j,k)= 0.0
              odoc13a(i,j,k)= 0.0
              odocs13a(i,j,k)= 0.0
              oc14a(i,j,k)= 0.0
#if ( PATH >= 1 )
              paparta(i,j,k) = 0.0
              padissa(i,j,k) = 0.0
              thparta(i,j,k) = 0.0
              thdissa(i,j,k) = 0.0
              pflxcaa(i,j,k) = 0.0
              pflxpoa(i,j,k) = 0.0
#endif
#if ( CORAL == 1 )
              corareaa(i,j,k)= 0.0
              corproda(i,j,k)= 0.0
              cormassa(i,j,k)= 0.0
              omegaa(i,j,k)= 0.0
              oco3a(i,j,k)= 0.0
              tau_bleacha(i,j,k)= 0.0
              DHWa(i,j,k)= 0.0
              PHa(i,j,k)= 0.0
#endif
#if ( OOISO == 1 )
              oo2_2a(i,j,k)= 0.0
              oo2_3a(i,j,k)= 0.0
              oo2_4a(i,j,k)= 0.0
#endif
#endif

            enddo
          enddo
        enddo
!
!DFG start initialisation of arrays for netcdf-output
!
        do j=1,jmax
          do i=1,imax
            f2d(i,j)=0.0
            etamo(i,j)=0.0
            ubmo(i,j)=0.0
            vbmo(i,j)=0.0
            etaan(i,j)=0.0
            uban(i,j)=0.0
            vban(i,j)=0.0
#if ( OCYCC == 1 )
            aqpco2m(i,j)= 0.0
            tpp_mam(i,j)= 0.0
            caco3mm(i,j)= 0.0
            aqpco2a(i,j)= 0.0
            tpp_maa(i,j)= 0.0
            caco3ma(i,j)= 0.0
#endif
#if ( F_PALAEO_FWF == 1 || CONSEAU == 1 )
            outFWFm(i,j)= 0.0
            outFWFa(i,j)= 0.0
#endif
          enddo
        enddo
        do k=1,kmax+1
          do j=1,jmax
            do i=1,imax
              wmm(i,j,k)=0.0
              wma(i,j,k)=0.0
            enddo
          enddo
        enddo
        do nb=0,nbsmax
          do j=1,jmtt
            do k=0,kmtt
              moc_m(nb,j,k)=0.0
              moc_a(nb,j,k)=0.0
            enddo
            mht_m(nb,j)=0.0
            mst_m(nb,j)=0.0
            mht_a(nb,j)=0.0
            mst_a(nb,j)=0.0
          enddo
        enddo

!DFG end


!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--1.5. Reading output.param
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

        open(newunit=outputparam_id,file='output.param')
!
        read(outputparam_id,*)
        read(outputparam_id,*)
        read(outputparam_id,*)
        read(outputparam_id,*)
        read(outputparam_id,*)
!       read(outputparam_id,*) noumef
        read(outputparam_id,*) npack
!
! Correspondance between npack and he number of fields
        noumef = 0
!dmr @-@ iceb0
        if (npack.eq.1) noumef=20
!dmr orig       if (npack.eq.1) noumef=18
!dmr jochem       if (npack.eq.1) noumef=19
!dmr @-@ iceb0
! mab: adding +1 to the definitions set her to adapt it
! mab: (besides npack.eq.1 because this was already adapted to 19 before)
! mab: orig:       if (npack.eq.2/3/4) noumef=34/43/47
        if (npack.eq.2) noumef=36
        if (npack.eq.3) noumef=45
        if (npack.eq.4) noumef=50
!
        if (noumef.gt.noumax) then
          write(clio3_out_id,*)'STOP : Too much fields for monthly mean'
          write(clio3_out_id,*)'noumef=',noumef,'  noumax=',noumax
          STOP
        endif
!DFG        if (noumef.lt.noumax) then
!DFG          write(clio3_out_id,*) 'WARNING : You can reduce the memory requirements'
!DFG          write(clio3_out_id,*) '  by reducing noumax in outave : noumax=',noumax
!DFG          write(clio3_out_id,*) '  and it can de reduced to noumef=', noumef
!DFG        endif
!dmr @-@ iceb0
! JONO_out adding other case so noumef and noumax are printed:
       if (noumef.eq.noumax)
     &  write(clio3_out_id,*) 'noumef=',noumef,'=(good for you)=noumax='
     &                        ,noumax
!dmr @-@ iceb0
!
! Computation of the levels for averages
        do kz=ks2,ks1,-1
          if ((-60).ge.zw(kz)) then
             klevv=kz
             goto 141
          endif
        enddo
141     continue
        do kk= 1,4
          klevm(kk,1)=0
          klevm(kk,2)=0
        enddo
!mab: commented for mouchard
!        if (nn99.eq.2) then
!          write(mouchard_id,*)
!          write(mouchard_id,'(A,I3)')
!     &     ' Vertical averages for outputs(outave) : klevv=', klevv
!       endif
!mab: end of comment
        read(outputparam_id,*)
        do kk=1,4
          read(outputparam_id,*) zlim1,zlim2
          do kz=ks2,ks1,-1
            if ((-zlim1).ge.zw(kz)) then
              klevm(kk,1)=kz
              goto 142
            endif
          enddo
142       continue
          do kz=ks2,ks1,-1
            if ((-zlim2).ge.zw(kz)) then
              klevm(kk,2)=kz+1
              goto 143
            endif
          enddo
143       if (klevm(kk,2).eq.0) klevm(kk,2)=1
!mab: comment for mouchard
!          if (nn99.eq.2) then
!            write(mouchard_id,*) 'Layer number',kk
!            write(mouchard_id,*) 'lim1 and lim2 :',zlim1,zlim2
!            write(mouchard_id,*) 'lev1 and lev2 :',klevm(kk,1),klevm(kk,2)
!          endif
!mab: end of comment
        enddo
!  Number of value for each average
        do j=2,jmax-1
          do i=2,imax-1
            znumav(i,j,5)=0.0
            do kk=klevv,ku2
              znumav(i,j,5)=znumav(i,j,5)+tmu(i,j,kk)
            enddo
            znumav(i,j,5)=max(znumav(i,j,5),undim6)
          enddo
        enddo

        do j=ju1,ju2
          do i=iu1(j),iu2(j)
            do kz=1,4
              znumav(i,j,kz)=0.0
              do kk=klevm(kz,2),klevm(kz,1)
                znumav(i,j,kz)=znumav(i,j,kz)+tmu(i,j,kk)
              enddo
              znumav(i,j,kz)=max(znumav(i,j,kz),undim6)
            enddo
          enddo
        enddo

        do j=js1,js2
          do i=is1(j),is2(j)
            do kz=1,4
              znumas(i,j,kz)=0.0
              do kk=klevm(kz,2),klevm(kz,1)
                znumas(i,j,kz)=znumas(i,j,kz)+tms(i,j,kk)
              enddo
              znumas(i,j,kz)=max(znumas(i,j,kz),undim6)
            enddo
          enddo
        enddo
! mab: here all the necessary lines of output.param are read in
        read(outputparam_id,*)
        read(outputparam_id,*)
        read(outputparam_id,*)
        do npw=1,noumef
          read(outputparam_id,'(A)') titn(npw)
          read(outputparam_id,*) i,ntypou(npw),nmm(npw),nmal(npw)
     &        ,nma(npw), cmulti(npw),cadd(npw)
        enddo
        close(outputparam_id)
        nmmef=0
        nmalef=0
        nmaef=0
        do npw=1,noumef
          if (nmm(npw).eq.1) then
            nmmef=nmmef+1
          endif
          if (nma(npw).eq.1) then
            nmaef=nmaef+1
          endif
          if (nmal(npw).eq.1) then
            nmalef=nmalef+1
            ncor(nmalef)=npw
          endif
        enddo
!        write(clio3_out_id,*) 'nmmef etc'
!        write(clio3_out_id,*) nmmef,nmaef,nmalef,ncor
        if (nmalef.gt.noumax2) then
          write(clio3_out_id,*)
     &                     'STOP : Too much fields for all period mean'
          STOP
        endif
        if (nmalef.lt.noumax2) then
          write(clio3_out_id,*)
     &               'WARNING : You can reduce the memory requirements'
          write(clio3_out_id,*)
     &                     '  by reducing noumax2 in outave : noumax2=',
     &                     noumax2
         write(clio3_out_id,*) '  and it can de reduced to nmalef='
     &, nmalef
        endif
!
!DFG start
!
!-- 1.6. PREPARATION OF NETCDF OUTPUT (IF DESIRED)
!
!-- 1.6.1 read netcdfout.param to decide what needs to be done.
!
!         note:
!
!         when netcdfout.param is not found, or netcdf output is
!         disabled therein, outave will perform binary dumps (to
!         cresum.dat etc.) as usual.
!
! --  preset switches, IDs and record counters
!
        lcdfmon=.false.
        lcdfann=.false.
        do i=1,ldim
          lvcdf(i)=.false.
          lmcdf(i)=.false.
          lacdf(i)=.false.
        enddo
        do i=1,ldim+2
          icdfNvars(i)=0
          icdfMvars(i)=0
          icdfAvars(i)=0
        enddo
!
! -- preset spval with CLIO3 standard missing value of -100.
!
!dmr     spval=spvr
!dmr     spval=-1.0d+32

!
        open(newunit=netcdfoutparam_id,file='netcdfout.param'
     &      ,status='old',err=128)
        write(clio3_out_id,*) 'MESSAGE : file netcdfout.param found,',
     &    '  checking for monthly/annual mean dumps to netcdf file:'
!
! -- switches for monthly/annual mean dumps
!
        read(netcdfoutparam_id,*)
        read(netcdfoutparam_id,*) lcdfmon, maxmrecs
        read(netcdfoutparam_id,*)
        read(netcdfoutparam_id,*) lmcdf
        read(netcdfoutparam_id,*)
        read(netcdfoutparam_id,*) lcdfann, maxarecs
        read(netcdfoutparam_id,*)
        read(netcdfoutparam_id,*) lacdf

        !write(*,*) 'lcdfmon= ', lcdfmon, 'lcdfann= ', lcdfann
        !write(*,*) maxmrecs, maxarecs
        !write(*,*) lacdf

        if (lcdfmon) then
          write(clio3_out_id,*)
     &               '   monthly mean netcdf dumps will be written'//
     &               ' using the native model grid.'
!mdr ---          spval=-1.0d+32
        endif

        if (lcdfann) then
          write(clio3_out_id,*)
     &               '   annual mean netcdf dumps will be written'//
     &               ' using the native model grid.'
!mdr ---          spval=-1.0d+32
        endif

        goto 129
  128   write(clio3_out_id,'(a,/)')
     &    'WARNING : file netcdfout.param not found',
     &    '  outave will perform dumps to cresum etc. as usual.'
  129   close(netcdfoutparam_id)
!DFG end
!
! end of preparation if numit=nstart
!
      endif
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--2. Computation of averages                                          |
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!
!--2.1. WRITE LOCAL DAILY OUTPUTS.
!---------------------------------
!
      if (mod(ja,nwjl).eq.nwtest) then
!
            nom2=4
            nom3=3
            write(cresujldat_id) nom2,nom3
            write(cresujldat_id) ja,ijour
            do np=1,nptj
              write(cresujldat_id) 'Epaisseur de neige   (71)'
              write(cresujldat_id) hnbq(iwl(np),jwl(np))
              write(cresujldat_id) 'Epaisseur de glace   (72)'
              write(cresujldat_id) hgbq(iwl(np),jwl(np))
              write(cresujldat_id) 'Epaisseur cree       (73)'
              write(cresujldat_id) hgbqp(iwl(np),jwl(np))
              write(cresujldat_id) 'Aire des leads       (74)'
              write(cresujldat_id) albq(iwl(np),jwl(np))
!             write(cresujldat_id) 'Temp de surface      (75)'
!             write(cresujldat_id) ts(iwl(np),jwl(np))
!             write(cresujldat_id) 'Temp de neige        (76)'
!             write(cresujldat_id) tbq(iwl(np),jwl(np),1)
!             write(cresujldat_id) 'Temp de glace 1      (77)'
!             write(cresujldat_id) tbq(iwl(np),jwl(np),2)
!             write(cresujldat_id) 'Temp de glace 2      (78)'
!             write(cresujldat_id) tbq(iwl(np),jwl(np),3)
!             write(cresujldat_id) 'var volume surf      (79)'
!             write(cresujldat_id) dvosbq(iwl(np),jwl(np))
!             write(cresujldat_id) 'var volume bott      (80)'
!             write(cresujldat_id) dvobbq(iwl(np),jwl(np))
!             write(cresujldat_id) 'var volume lead      (81)'
!             write(cresujldat_id) dvolbq(iwl(np),jwl(np))
!             write(cresujldat_id) 'var volume neige     (82)'
!             write(cresujldat_id) dvonbq(iwl(np),jwl(np))
!             write(cresujldat_id) 'flux base glace      (83)'
!             write(cresujldat_id) fbbq(iwl(np),jwl(np))
!             write(cresujldat_id) 'temp air             (84)'
!             write(cresujldat_id) tabq(iwl(np),jwl(np))
!             write(cresujldat_id) 'humid air            (85)'
!             write(cresujldat_id) qabq(iwl(np),jwl(np))
!             write(cresujldat_id) 'vit vent             (86)'
!             write(cresujldat_id) vabq(iwl(np),jwl(np))
!             write(cresujldat_id) 'flux sol             (87)'
!             write(cresujldat_id) fsolg(iwl(np),jwl(np))
!             write(cresujldat_id) 'flux IR              (88)'
!             write(cresujldat_id) firg(iwl(np),jwl(np))
!             write(cresujldat_id) 'flux sensible        (89)'
!             write(cresujldat_id) fcsg(iwl(np),jwl(np))
!             write(cresujldat_id) 'flux latent          (90)'
!             write(cresujldat_id) fleg(iwl(np),jwl(np))

              do kk=1,kmax
                 zms=tms(iwl(np),jwl(np),kk)
                 sort1(kk)=zms*(scal(iwl(np),jwl(np),kk,1)
     &                    -273.15) + (1.0-zms)*(-100.0)
                 sort2(kk)=zms*scal(iwl(np),jwl(np),kk,2)
     &                     + (1.0-zms)*(-100.0)
                 sort3(kk)=zms*q2turb(iwl(np),jwl(np),kk)
     &                     + (1.0-zms)*(-100.0)
!tk0             sort3(kk)=zms*bvf(iwl(np),jwl(np),kk)
!tk0 &                     + (1.0-zms)*(-100.0)
              enddo
              write(cresujldat_id) 'temperature ocean    (91)'
!             write(cresujldat_id) (scal(iwl(np),jwl(np),kv,1),kv=1,kmax)
              write(cresujldat_id) (sort1(kv),kv=1,kmax)
              write(cresujldat_id) 'Salinite ocean       (92)'
!             write(cresujldat_id) (scal(iwl(np),jwl(np),kv,2),kv=1,kmax)
              write(cresujldat_id) (sort2(kv),kv=1,kmax)
              write(cresujldat_id) 'Q2turb               (93)'
!             write(cresujldat_id) 'n2                   (93)'
!             write(cresujldat_id) (q2turb(iwl(np),jwl(np),kv),kv=2,kmax+1)
              write(cresujldat_id) (sort3(kv),kv=1,kmax)

            enddo
!
      endif
!
! computation of annula mean surface velocity
      zfluxmt=zfluxmt+zflux0
      zfluxmts=zfluxmts+zflux0s
!
      if (nwtal.eq.1.or.mod(ja,nwm).eq.nwtest) then
        njm=njm+1
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--2.2. CUMULATION OF GLOBAL MONTHLY OUTPUTS.
!--------------------------------------------
! streamFUNCTION part
! cumulation for monthly and/or yearly outputs !
        if (nstreamout.ne.0) call streamfunc(0,0.d0,nstreamout)
!
!DFG start accumulation of ssh, ubar, vbar
!
        do j=1,jmax
          do i=1,imax
            etamo(i,j) = etamo(i,j) + eta(i,j)
             ubmo(i,j) =  ubmo(i,j) +  ub(i,j)
             vbmo(i,j) =  vbmo(i,j) +  vb(i,j)
          enddo
        enddo
!DFG end
!
!DFG start
!       do k=1,kmax
!         do j=2,jmax-1
!           do i=2,imax-1
!DFG end
        do k=1,kmax
          do j=1,jmax
            do i=1,imax
              umm(i,j,k) = umm(i,j,k)+u(i,j,k)
              vmm(i,j,k) = vmm(i,j,k)+v(i,j,k)
              tmm(i,j,k) = tmm(i,j,k)+scal(i,j,k,1)
              smm(i,j,k) = smm(i,j,k)+scal(i,j,k,2)
              uma(i,j,k) = uma(i,j,k)+u(i,j,k)
              vma(i,j,k) = vma(i,j,k)+v(i,j,k)
              tma(i,j,k) = tma(i,j,k)+scal(i,j,k,1)
              sma(i,j,k) = sma(i,j,k)+scal(i,j,k,2)
#if ( ISOOCN >= 1 )
              oo17m(i,j,k)= oo17m(i,j,k)+scal(i,j,k,ocnw17)
              oo18m(i,j,k)= oo18m(i,j,k)+scal(i,j,k,ocnw18)
              oohdm(i,j,k)= oohdm(i,j,k)+scal(i,j,k,ocnw2h)
              oo17a(i,j,k)= oo17a(i,j,k)+scal(i,j,k,ocnw17)
              oo18a(i,j,k)= oo18a(i,j,k)+scal(i,j,k,ocnw18)
              oohda(i,j,k)= oohda(i,j,k)+scal(i,j,k,ocnw2h)
#endif
#if ( OCYCC == 1 )
              odocm(i,j,k)= odocm(i,j,k)+scal(i,j,k,3)
              odocsm(i,j,k)= odocsm(i,j,k)+scal(i,j,k,4)
              odicm(i,j,k)= odicm(i,j,k)+scal(i,j,k,5)
              opo4m(i,j,k)= opo4m(i,j,k)+scal(i,j,k,6)
              oalkm(i,j,k)= oalkm(i,j,k)+scal(i,j,k,8)
              oo2m(i,j,k)= oo2m(i,j,k)+scal(i,j,k,9)
              oc13m(i,j,k)= oc13m(i,j,k)+scal(i,j,k,10)
              odoc13m(i,j,k)= odoc13m(i,j,k)+scal(i,j,k,11)
              odocs13m(i,j,k)= odocs13m(i,j,k)+scal(i,j,k,12)
              oc14m(i,j,k)= oc14m(i,j,k)+scal(i,j,k,13)
#if ( PATH >= 1 )
            padissm(i,j,k) = padissm(i,j,k)+scal(i,j,k,scalstart_path+0)
            papartm(i,j,k) = papartm(i,j,k)+scal(i,j,k,scalstart_path+1)
            thdissm(i,j,k) = thdissm(i,j,k)+scal(i,j,k,scalstart_path+2)
            thpartm(i,j,k) = thpartm(i,j,k)+scal(i,j,k,scalstart_path+3)
            pflxcam(i,j,k) = pflxcam(i,j,k)
     &                        +particles_fluxes_field(i,j,k,ncaco3)
            pflxpom(i,j,k) = pflxpom(i,j,k)
     &                        +particles_fluxes_field(i,j,k,npoc)
#endif
#if ( CORAL == 1 )
           coraream(i,j,k)= coraream(i,j,k)+coral_area_clio(i,j,k)
           corprodm(i,j,k)= corprodm(i,j,k)+coral_prod_clio(i,j,k)
           cormassm(i,j,k)= cormassm(i,j,k)+coral_mass_clio(i,j,k)
           omegam(i,j,k)= omegam(i,j,k)+omega_arag_clio(i,j,k)
           oco3m(i,j,k)= oco3m(i,j,k)+oco3_clio(i,j,k)
           tau_bleachm(i,j,k)= tau_bleachm(i,j,k)+tau_bleach_clio(i,j,k)
           DHWm(i,j,k)= DHWm(i,j,k)+DHW_clio(i,j,k)
           PHm(i,j,k)= PHm(i,j,k)+PH_clio(i,j,k)
#endif
#if ( OOISO == 1 )
              oo2_2m(i,j,k)= oo2_2m(i,j,k)+scal(i,j,k,oo2iso16)
              oo2_3m(i,j,k)= oo2_3m(i,j,k)+scal(i,j,k,oo2iso17)
              oo2_4m(i,j,k)= oo2_4m(i,j,k)+scal(i,j,k,oo2iso18)
#endif
              odoca(i,j,k)= odoca(i,j,k)+scal(i,j,k,3)
              odocsa(i,j,k)= odocsa(i,j,k)+scal(i,j,k,4)
              odica(i,j,k)= odica(i,j,k)+scal(i,j,k,5)
              opo4a(i,j,k)= opo4a(i,j,k)+scal(i,j,k,6)
              oalka(i,j,k)= oalka(i,j,k)+scal(i,j,k,8)
              oo2a(i,j,k)= oo2a(i,j,k)+scal(i,j,k,9)
              oc13a(i,j,k)= oc13a(i,j,k)+scal(i,j,k,10)
              odoc13a(i,j,k)= odoc13a(i,j,k)+scal(i,j,k,11)
              odocs13a(i,j,k)= odocs13a(i,j,k)+scal(i,j,k,12)
              oc14a(i,j,k)= oc14a(i,j,k)+scal(i,j,k,13)
#if ( PATH >= 1 )
            padissa(i,j,k) = padissa(i,j,k)+scal(i,j,k,scalstart_path+0)
            paparta(i,j,k) = paparta(i,j,k)+scal(i,j,k,scalstart_path+1)
            thdissa(i,j,k) = thdissa(i,j,k)+scal(i,j,k,scalstart_path+2)
            thparta(i,j,k) = thparta(i,j,k)+scal(i,j,k,scalstart_path+3)
            pflxcaa(i,j,k) = pflxcaa(i,j,k)
     &                        +particles_fluxes_field(i,j,k,ncaco3)
            pflxpoa(i,j,k) = pflxpoa(i,j,k)
     &                        +particles_fluxes_field(i,j,k,npoc)
#endif
#if ( CORAL == 1 )
              corareaa(i,j,k)= corareaa(i,j,k)+coral_area_clio(i,j,k)
              corproda(i,j,k)= corproda(i,j,k)+coral_prod_clio(i,j,k)
              cormassa(i,j,k)= cormassa(i,j,k)+coral_mass_clio(i,j,k)
              !if (coral_prod_clio(i,j,k).ne.0) then
              !if (corproda(i,j,k).ne.0) then
              !  write(*,*) 'in outnet.f',coral_prod_clio(i,j,k)
              !  write(*,*) 'in outnet.f',corproda(i,j,k)
              !endif
              omegaa(i,j,k)= omegaa(i,j,k)+omega_arag_clio(i,j,k)
              oco3a(i,j,k)= oco3a(i,j,k)+oco3_clio(i,j,k)
              tau_bleacha(i,j,k)=tau_bleacha(i,j,k)
     &                          +tau_bleach_clio(i,j,k)
              DHWa(i,j,k)=DHWa(i,j,k)+DHW_clio(i,j,k)
              PHa(i,j,k)=PHa(i,j,k)+PH_clio(i,j,k)
#endif
#if ( OOISO == 1 )
              oo2_2a(i,j,k)= oo2_2a(i,j,k)+scal(i,j,k,oo2iso16)
              oo2_3a(i,j,k)= oo2_3a(i,j,k)+scal(i,j,k,oo2iso17)
              oo2_4a(i,j,k)= oo2_4a(i,j,k)+scal(i,j,k,oo2iso18)
#endif
#if ( OXNITREUX == 1 )
              on2om(i,j,k)= on2om(i,j,k)+scal(i,j,k,7)
              on2oa(i,j,k)= on2oa(i,j,k)+scal(i,j,k,7)
#else
              ono3m(i,j,k)= ono3m(i,j,k)+scal(i,j,k,7)
              ono3a(i,j,k)= ono3a(i,j,k)+scal(i,j,k,7)
#endif
#endif
            enddo
          enddo
        enddo

#if ( OCYCC == 1 )
!dmr --- La variable aqpco2 est un ajout a 2D !!
          do j=1,jmax
            do i=1,imax
             aqpco2m(i,j)= aqpco2m(i,j)+oxpco2_clio(i,j)
             aqpco2a(i,j)= aqpco2a(i,j)+oxpco2_clio(i,j)

             tpp_maa(i,j)= tpp_maa(i,j)+TPPma_clio(i,j)
             tpp_mam(i,j)= tpp_mam(i,j)+TPPma_clio(i,j)

             caco3ma(i,j)= caco3ma(i,j)+CACO3ma_clio(i,j)
             caco3mm(i,j)= caco3mm(i,j)+CACO3ma_clio(i,j)
          enddo
        enddo
#endif

#if ( F_PALAEO_FWF == 1 || CONSEAU == 1 )
          !print*, "GET FWF output"
          do j=1,jmax
            do i=1,imax
             outFWFm(i,j)= outFWFm(i,j)+sum_flux_out(i,j)
             outFWFa(i,j)= outFWFa(i,j)+sum_flux_out(i,j)
          enddo
        enddo
#endif

!DFG start
        do k=1,kmax+1
          do j=1,jmax
            do i=1,imax
              wmm(i,j,k) = wmm(i,j,k)+w(i,j,k)
              wma(i,j,k) = wma(i,j,k)+w(i,j,k)
            enddo
          enddo
        enddo
!DFG end

!DFG start
!        if (npack.ge.1) then
        if ((npack.ge.1).or.lcdfmon.or.lcdfann) then
!DFG end
          do j=2,jmax-1
            do i=2,imax-1
!dmr @-@ iceb0
               chk_Vicb_outave(i,j)=0.
#if (CALVFLUX == 2)
                  chk_Vicb_outave(i,j)=chk_Vicb_outave(i,j)+
     &              dVol_calv(i,j,ks2)
#else
               do klyr= 1, ks2
                  chk_Vicb_outave(i,j)=chk_Vicb_outave(i,j)+
     &              dVol_icb(i,j,klyr)
               enddo
#endif
               total_vol=total_vol+chk_Vicb_outave(i,j)
!dmr @-@ iceb0
                zindh       = max(zero,sign(one,hgbq(i,j)*
     &                        (1.0-albq(i,j))-0.10))
                zinda       = max(zero,sign(one,
     &                         (1.0-albq(i,j))-0.10))
                zindb       = zindh*zinda
                cmoymo(i,j,1)= cmoymo(i,j,1)+hnbq(i,j)
                cmoymo(i,j,2)= cmoymo(i,j,2)+hgbq(i,j)
                cmoymo(i,j,3)= cmoymo(i,j,3)+hgbqp(i,j)
                cmoymo(i,j,4)= cmoymo(i,j,4)+albq(i,j)
                cmoymo(i,j,5)= cmoymo(i,j,5)+ts(i,j)
                cmoymo(i,j,6)= cmoymo(i,j,6)+fbbq(i,j)
                cmoymo(i,j,7)= cmoymo(i,j,7)+
     &                ug(i,j)*tmu(i,j,ks2)*zindb
                cmoymo(i,j,8)= cmoymo(i,j,8)+
     &                vg(i,j)*tmu(i,j,ks2)*zindb
! mab: here +1 is added to the cmoymo(i,j,>= 9(IMF) and now >=10) (DISTICB)
! mab: orig:
!                cmoymo(i,j,9) = cmoymo(i,j,9)+scal(i,j,ks2,1)
!                cmoymo(i,j,10)= cmoymo(i,j,10)+scal(i,j,ks2,2)
!                cmoymo(i,j,15)= cmoymo(i,j,15) - phiss(i,j,1)*
!     &                         (rho0*cpo)/(dts(ks2)*unsdz(ks2))
!     &                         +fcm1(i,j)
!                cmoymo(i,j,16)= cmoymo(i,j,16)+86400./dts(ks2)*rho0*
!     &                  ( phiss(i,j,2)/(unsdz(ks2)*34.7)
!     &                    -phiss(i,j,0)
!C    &                   +rappes(i,j,0)*dz(ks2)
!C    &             *(scal(i,j,ks2,2)-scalr(i,j,ks2,2))/scal(i,j,ks2,2)
!     &                  )
!                cmoymo(i,j,17) = cmoymo(i,j,17) + ai(ks2) * sqrt(
!     &          c4x(i,j,ks2)*c4x(i,j,ks2)+c4y(i,j,ks2)*c4y(i,j,ks2) )
!                cmoymo(i,j,18)=cmoymo(i,j,18)+ eta(i,j)
! mab: new:

!age            cmoymo(i,j,18)=cmoymo(i,j,18)+ageg(i,j)
!dmr @-@ iceb0
! JONO WIEA cmoymo thinks he is getting a 'volume'(m3) PER SEC per day,
! so we have to multiply with daysec=86400
! JA heat2layers            cmoymo(i,j,18)=cmoymo(i,j,18)+fonte_icb(i,j)*86400.
!dmr --- From JONO, WIEA, moved 18 ---> 19, 18 already used in LOVECLIM
!dmr --- Back change the maximum number!

!dmr&mab --- 18th may, 2011
!dmr&mab --- The flux chk_Vicb_outave is in m3.day-1 so we need to divide
!dmr&mab ---  the flux by 86400. to get m3.s-1, not multiply it!
!dmr&mab orig                cmoymo(i,j,19)=cmoymo(i,j,19)+chk_Vicb_outave(i,j)*86400.
! mab: added new line, 9 for icebergs!!and adapted  all the other lines
! mab:          cmoymo(i,j,19)=cmoymo(i,j,19)+chk_Vicb_outave(i,j)/86400.

              cmoymo(i,j,9)=cmoymo(i,j,9)+chk_Vicb_outave(i,j)/86400.
!mab: add iceberg distribution
!PB              cmoymo(i,j,10)=cmoymo(i,j,10)+disticeb(i,j)
! JONO added fonte aug03, flag this output in output.param
!dmr @-@ iceb0


                cmoymo(i,j,11) = cmoymo(i,j,11)+scal(i,j,ks2,1)
                cmoymo(i,j,12)= cmoymo(i,j,12)+scal(i,j,ks2,2)
                cmoymo(i,j,17)= cmoymo(i,j,17) - phiss(i,j,1)*
     &                         (rho0*cpo)/(dts(ks2)*unsdz(ks2))
     &                         +fcm1(i,j)
                cmoymo(i,j,18)= cmoymo(i,j,18)+86400./dts(ks2)*rho0*
     &                  ( phiss(i,j,2)/(unsdz(ks2)*34.7)
     &                    -phiss(i,j,0)
!    &                   +rappes(i,j,0)*dz(ks2)
!    &             *(scal(i,j,ks2,2)-scalr(i,j,ks2,2))/scal(i,j,ks2,2)
     &                  )
              cmoymo(i,j,19) = cmoymo(i,j,19) + ai(ks2) * sqrt(
     &          c4x(i,j,ks2)*c4x(i,j,ks2)+c4y(i,j,ks2)*c4y(i,j,ks2) )
                cmoymo(i,j,20)=cmoymo(i,j,20)+ eta(i,j)

            enddo
          enddo
!dmr @-@ iceb0
!mab: too much output in mouchard
!                write(mouchard_id,*)'outave: total volume from bergs',total_vol
!dmr @-@ iceb0
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
          do j=2,jmax-1
            do i=2,imax-1
              zztmp1=0.0
              zztmp2=0.0
              do kk=klevv,ku2
                zztmp1=zztmp1+tmu(i,j,kk)*u(i,j,kk)
                zztmp2=zztmp2+tmu(i,j,kk)*v(i,j,kk)
              enddo
! mab: added +1 (IMF) +1 (DISTICB)
              cmoymo(i,j,13)= cmoymo(i,j,13) +
     &                      zztmp1/znumav(i,j,5)
              cmoymo(i,j,14)= cmoymo(i,j,14) +
     &                      zztmp2/znumav(i,j,5)
            enddo
          enddo
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
          do j=js1,js2
            do  i=is1(j),is2(j)
              kmin3 = ks2
              ztest3 = 0.0
              kmin4 = ks2
              ztest4 = 0.0

              tloc=scal(i,j,ks2,1)-273.15
              sloc=scal(i,j,ks2,2)
              ccb1 = cstrho(4)*sloc
     &         + (cstrho(2)-cstrho(3)*tloc)*tloc
              ccb2 = cstrho(5)
     &              + (cstrho(6) - cstrho(7)*tloc)*tloc
     &              - (cstrho(8) + cstrho(9)*tloc)*sloc
              sigsurf=1000.0*(1.0/(cstrho(0)+
     &                 ccb2/(ccb1+cfb1z4(ks2+1)) ) - 1.0 )
!             write(155,*) sigsurf*tms(i,j,ks2)

              do k=ks2,kfs(i,j),-1
                ztest3= ztest3+max(zero,sign(one,
     &                    bvf(i,j,k) ) )
                ztest3=min(ztest3,one)
                kmin3 = int(1.0-ztest3)*k+int(ztest3)*kmin3

                tloc=scal(i,j,k,1)-273.15
                sloc=scal(i,j,k,2)
                ccb1 = cstrho(4)*sloc
     &         +  (cstrho(2)-cstrho(3)*tloc)*tloc
                ccb2 = cstrho(5)
     &               + (cstrho(6) - cstrho(7)*tloc)*tloc
     &               - (cstrho(8) + cstrho(9)*tloc)*sloc
                sigz=1000.0*(1.0/(cstrho(0)+ ccb2
     &              /(ccb1+cfb1z4(ks2+1)) ) -1.0 )
                ztest4= ztest4+max(zero,sign(one,
     &                     -0.02+abs(sigz-sigsurf) ) )
                ztest4=min(ztest4,one)
                kmin4 = int(1.0-ztest4)*k+int(ztest4)*kmin4
               enddo
!mb: added +1(IMF) +1 (DISTICB)
              cmoymo(i,j,16)= cmoymo(i,j,16)-zw(kmin3-1)
              cmoymo(i,j,15)= cmoymo(i,j,15)-zw(kmin4)
              zmix_CC(i,j)  = -zw(kmin4)
            enddo
          enddo
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
        endif
!-------
!DFG start
!        if (npack.ge.2) then
        if ((npack.ge.2).or.lcdfmon.or.lcdfann) then
!DFG end
          do kz=1,4
! mab: added +1  (orig:21/22 (IMF)/now 23(DIST ICB))
            kji=23+(kz-1)*4
            do kk=klevm(kz,2),klevm(kz,1)
              do j=ju1,ju2
                do i=iu1(j),iu2(j)
                  cmoymo(i,j,kji)=cmoymo(i,j,kji)+
     &                tmu(i,j,kk)*u(i,j,kk)/znumav(i,j,kz)
                  cmoymo(i,j,kji+1)=cmoymo(i,j,kji+1)+
     &                tmu(i,j,kk)*v(i,j,kk)/znumav(i,j,kz)
                enddo
              enddo
            enddo
          enddo
          do kz=1,4
! mab: added +1 (orig: 19/20(IMF)/21(DISTICB))
            kji=21+(kz-1)*4
            do kk=klevm(kz,2),klevm(kz,1)
              do j=js1,js2
                do i=is1(j),is2(j)
                  cmoymo(i,j,kji)=cmoymo(i,j,kji)+
     &               tms(i,j,kk)*scal(i,j,kk,1)/znumas(i,j,kz)
                  cmoymo(i,j,kji+1)=cmoymo(i,j,kji+1)+
     &               tms(i,j,kk)*scal(i,j,kk,2)/znumas(i,j,kz)
                enddo
              enddo
            enddo
          enddo
        endif
!-------
!DFG start
!        if (npack.ge.3) then
        if ((npack.ge.3).or.lcdfmon.or.lcdfann) then
!DFG end
          do j=js1,js2
            do  i=is1(j),is2(j)
! mab: added +1 (IMF) +1 (DISTICB)
              cmoymo(i,j,37)= cmoymo(i,j,37)+fcm1(i,j)
              cmoymo(i,j,38)= cmoymo(i,j,38)+tenagx(i,j)
              cmoymo(i,j,39)= cmoymo(i,j,39)+tenagy(i,j)
              cmoymo(i,j,40)= cmoymo(i,j,40)+fsolg(i,j)*
     &                   (1.0-albq(i,j))+albq(i,j)*fsolcn(i,j)
              cmoymo(i,j,41)= cmoymo(i,j,41)+(1.0-albq(i,j))
     &               *firg(i,j)+albq(i,j)*ratbqo(i,j)
              cmoymo(i,j,42)= cmoymo(i,j,42)+(1.0-albq(i,j))
     &               *fcsg(i,j)+albq(i,j)*fcscn(i,j)
              cmoymo(i,j,43)= cmoymo(i,j,43)+(1.0-albq(i,j))
     &               *fleg(i,j)+albq(i,j)*flecn(i,j)
              cmoymo(i,j,44)=cmoymo(i,j,44)+phiss(i,j,2)*rho0*
     &                  86400./(unsdz(ks2)*34.7*dts(ks2))
              cmoymo(i,j,45)= cmoymo(i,j,45)+albege(i,j)*
     &               (1-albq(i,j))+ albecn(i,j)*albq(i,j)
            enddo
          enddo
        endif
!-------
!DFG start
!        if (npack.ge.4) then
        if ((npack.ge.4).or.lcdfmon.or.lcdfann) then
!DFG end
          do j=js1,js2
            do  i=is1(j),is2(j)
! mab: added +1(IMF) +1 (DISTICB)
              cmoymo(i,j,46)= cmoymo(i,j,46)+tabq(i,j)
              cmoymo(i,j,47)= cmoymo(i,j,47)+qabq(i,j)
              cmoymo(i,j,48)= cmoymo(i,j,48)+vabq(i,j)
#if ( ISOOCN >= 1 )
              cmoymo(i,j,49)= cmoymo(i,j,49)+hnplbq(i,j,ieau)
#else
              cmoymo(i,j,49)= cmoymo(i,j,49)+hnplbq(i,j)
#endif
!mab: added calving flux
              cmoymo(i,j,50) = (calvCLIO(i,j)*0.91)/(86400.*360.)
            enddo
          enddo
        endif
!-------
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!-- COMPUTES MONTHLY MEAN
!
!         if (mod(ja,nwm).eq.nwtest.and.ijour.eq.njcum(jmois)) then
!         if (ijour.eq.njcum(jmois)) then

      if (ijour.eq.njcum(jmois).and.ijour.ne.ijour1) then
!
!           njm = (njcum(jmois)-njcum(jmois-1))-ijouri*(1-jjoufl)
        jjoufl = 1
        usnjm = 1.0/dble(njm)
        njm=0
!
        do n=1,noumef
          do j=1,jmax
            do i=1,imax
!             cmoymo(i,j,n)=cmoymo(i,j,n)*usnjm
              cmoymo(i,j,n)=(cmoymo(i,j,n)*usnjm)*cmulti(n)+cadd(n)
            enddo
          enddo
        enddo
!
        do jj=1,jmax
          do ij=1,imax
            zindp=max(0.0,sign(1.0,cmoymo(ij,jj,2)-0.2))
            cmoymo(ij,jj,6)=cmoymo(ij,jj,6)*zindp
          enddo
        enddo
!
        do n=1,noumef
          do j=1,jmax
            do i=1,imax
              zindo       = tms(i,j,ks2)
              cmoymo(i,j,n)=zindo*cmoymo(i,j,n)+(1.0-zindo)*spvr
            enddo
          enddo
        enddo
        if (npack.ge.2) then
          do kk=1,4
            do j=1,jmax
              do i=1,imax
                zindo1      = tms(i,j,klevm(kk,1))
                zindo2      = tmu(i,j,klevm(kk,1))
! mab: added +1 (orig: 19)(IMF) +1 (DISTICB)
                n=21+(kk-1)*4
                cmoymo(i,j,n)=zindo1*cmoymo(i,j,n)+(1.0-zindo1)*spvr
                cmoymo(i,j,n+1)=zindo1*cmoymo(i,j,n+1)+(1.0-zindo1)*spvr
                cmoymo(i,j,n+2)=zindo2*cmoymo(i,j,n+2)+(1.0-zindo2)*spvr
                cmoymo(i,j,n+3)=zindo2*cmoymo(i,j,n+3)+(1.0-zindo2)*spvr
              enddo
            enddo
          enddo
        endif
!
!DFG start computing of rest of 2d monthly means not yet contained
!    in cmoymo
!
        do j=1,jmax
          do i=1,imax
            etamo(i,j) = etamo(i,j)*usnjm
             ubmo(i,j) =  ubmo(i,j)*usnjm*0.01
             vbmo(i,j) =  vbmo(i,j)*usnjm*0.01
          enddo
        enddo
!DFG end
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! WRITES MONTHLY MEAN
!DFG start
!
        if (lcdfmon) then
!
! Do we need to create a new NetCDF-file?
!
          if (iMtimrec.eq.0) then
!
! generate a unique filename
!
! --- BdB 05-2019: adapted name equivalent to atmosphere output
            write(avgname,'(a,a,a,a,i6.6,a)')
     &        'outputdata/ocean/',
     &        'CLIO3_',TRIM(refexp),'_mon_avrg_',
     &         int(time_in_years(ja-init_year_clio)),'.nc'
!
! -- pass monthly mean var switches to prep_out
!
            do i=1,ldim
              lvcdf(i)=lmcdf(i)
            enddo
            do i=1,ldim+2
              icdfNvars(i)=0
              icdfMvars(i)=0
            enddo
            call prep_out(ncidm,'m',TRIM(avgname))
!
! -- save IDs of desired variables
!
            do i=1,ldim+2
              icdfMvars(i)=icdfNvars(i)
            enddo
          endif
!
! Update time record. Data are dumped at the end of the month, the
! timestamp represents the entire averaging period. See also the
! timerecord definition with the 360 day calendar and the 30 day
! averaging period attributes.
!
          !rtim=dfloat((ja-1)*360+ijour)/30.0
! ---     BdB 05-2019: replaced with global time value
          rtim=time_in_years(ja-init_year_clio) + dfloat(ijour)/360.
          iMtimrec=iMtimrec+1
          status=nf_put_var1_double(ncidm,IMtime,iMtimrec,rtim)
        endif
!DFG end
        if (mod(ja,nwm).eq.nwtest) then
!mab: comment for mouchard
!          if (nn99.eq.2) write(mouchard_id,'(2A,3I6)') 'Write monthly mean',
!     &      ' on file cresum.out ; year,month,day=', ja, jmois, ijour
!mab: end of comment
          write(cresujldat_id) ja,jmois
          write(cresujldat_id) nmmef
          do n=1,noumef
            if (nmm(n).eq.1) then
              write(cresujldat_id) n,ntypou(n),titn(n)
              do j=1,jmax
                write(cresujldat_id) (cmoymo(i,j,n),i=1,imax)
              enddo
            endif
          enddo

! streamFUNCTION part
          if (nstreamout.eq.1) call streamfunc(1,usnjm,nstreamout)
!DFG start streamFUNCTION part
          do nb=0,nbsmax
            do j=1,jmtt
              do k=0,kmtt
                moc_m(nb,j,k) = uuu(j,k,nb)
                moc_a(nb,j,k) = moc_a(nb,j,k)+  moc_m(nb,j,k)
                uuu(j,k,nb) = 0.0
              enddo
              do k=0,kfloor(j,nb)
                moc_m(nb,j,k) = spval
                moc_a(nb,j,k) = spval
              enddo
              mht_m(nb,j) = vvv(j,nb)
              mst_m(nb,j) = vvv(j,nb+4)
              mht_a(nb,j) = mht_a(nb,j)+ mht_m(nb,j)
              mst_a(nb,j) = mst_a(nb,j)+ mst_m(nb,j)
              vvv(j,nb) = 0.0
              vvv(j,nb+4) = 0.0
            enddo
          enddo
!
! -- apply mask
          do nb=1,nbsmax
            do j=1,16
              do k=0,kmtt
                moc_m(nb,j,k) = spval
              enddo
              mht_m(nb,j) = spval
              mst_m(nb,j) = spval
            enddo
          enddo
!
! dump of moc fields
!
          if (lcdfmon) then
            start(1)=1
            start(2)=1
            start(3)=1
            start(4)=iMtimrec
            count(1)=nbsmax+1
            count(2)=jmtt
            count(3)=kmtt
            count(4)=1
            if (lmmmoc ) status=nf_put_vara_double(ncidm,IMmoc ,
     &                     start,count, moc_m)
            start(3)=iMtimrec
            start(4)=0
            count(1)=nbsmax+1
            count(2)=jmtt
            count(3)=1
            count(4)=0
            if (lmmmht) status=nf_put_vara_double(ncidm,IMmht,
     &                     start,count,mht_m)
            if (lmmmst) status=nf_put_vara_double(ncidm,IMmst,
     &                     start,count,mst_m)
          endif
!
! -- clean up
          do j=1,jmtt
            do nb=0,nbsmax
              do k=0,kmtt
                moc_m(nb,j,k) = 0.0
              enddo
              mht_m(nb,j) = 0.0
              mst_m(nb,j) = 0.0
            enddo
          enddo
!DFG end
          if ((nwtom.eq.1).or.lcdfmon) then
!DFG start
!            do j=2,jmax-1
!             do i=2,imax-1
!              do k=1,kmax
            do k=1,kmax
              do j=1,jmax
                do i=1,imax
!                  umm(i,j,k)= tmu(i,j,k)*umm(i,j,k)*usnjm
!     &                     +(1-tmu(i,j,k))*spvr
!                  vmm(i,j,k)= tmu(i,j,k)*vmm(i,j,k)*usnjm
!     &                     +(1-tmu(i,j,k))*spvr
!                  tmm(i,j,k)= tms(i,j,k)*(tmm(i,j,k)*usnjm-273.15)
!     &                     +(1-tms(i,j,k))*spvr
!                  smm(i,j,k)= tms(i,j,k)*smm(i,j,k)*usnjm
!     &                     +(1-tms(i,j,k))*spvr
                  umm(i,j,k)= tmu(i,j,k)*umm(i,j,k)*usnjm
     &                     +(1-tmu(i,j,k))*spval
                  vmm(i,j,k)= tmu(i,j,k)*vmm(i,j,k)*usnjm
     &                     +(1-tmu(i,j,k))*spval
                  tmm(i,j,k)= tms(i,j,k)*(tmm(i,j,k)*usnjm-273.15)
     &                     +(1-tms(i,j,k))*spval
                  smm(i,j,k)= tms(i,j,k)*smm(i,j,k)*usnjm
     &                     +(1-tms(i,j,k))*spval
                  wmm(i,j,k)= tms(i,j,k)*wmm(i,j,k)*usnjm*86400.
     &                   +(1-tms(i,j,k))*spval
#if ( ISOOCN >= 1)
                  oo17m(i,j,k)= tms(i,j,k)*oo17m(i,j,k)*usnjm
     &                     +(1-tms(i,j,k))*spval
                  oo18m(i,j,k)= tms(i,j,k)*oo18m(i,j,k)*usnjm
     &                     +(1-tms(i,j,k))*spval
                  oohdm(i,j,k)= tms(i,j,k)*oohdm(i,j,k)*usnjm
     &                     +(1-tms(i,j,k))*spval
#endif
#if ( OCYCC == 1 )
                  odocm(i,j,k)= tms(i,j,k)*odocm(i,j,k)*usnjm
     &                     +(1-tms(i,j,k))*spval
                  odocsm(i,j,k)= tms(i,j,k)*odocsm(i,j,k)*usnjm
     &                     +(1-tms(i,j,k))*spval
                  odicm(i,j,k)= tms(i,j,k)*odicm(i,j,k)*usnjm
     &                     +(1-tms(i,j,k))*spval
                  opo4m(i,j,k)= tms(i,j,k)*opo4m(i,j,k)*usnjm
     &                     +(1-tms(i,j,k))*spval
#if ( OXNITREUX == 1 )
                  on2om(i,j,k)= tms(i,j,k)*on2om(i,j,k)*usnjm
     &                     +(1-tms(i,j,k))*spval
#else
                  ono3m(i,j,k)= tms(i,j,k)*ono3m(i,j,k)*usnjm
     &                     +(1-tms(i,j,k))*spval
#endif
                  oalkm(i,j,k)= tms(i,j,k)*oalkm(i,j,k)*usnjm
     &                     +(1-tms(i,j,k))*spval
                  oo2m(i,j,k)= tms(i,j,k)*oo2m(i,j,k)*usnjm
     &                     +(1-tms(i,j,k))*spval
                  oc13m(i,j,k)= tms(i,j,k)*oc13m(i,j,k)*usnjm
     &                     +(1-tms(i,j,k))*spval
                  odoc13m(i,j,k)= tms(i,j,k)*odoc13m(i,j,k)*usnjm
     &                     +(1-tms(i,j,k))*spval
                  odocs13m(i,j,k)= tms(i,j,k)*odocs13m(i,j,k)*usnjm
     &                     +(1-tms(i,j,k))*spval
                  oc14m(i,j,k)= tms(i,j,k)*oc14m(i,j,k)*usnjm
     &                     +(1-tms(i,j,k))*spval

#if ( PATH >= 1 )
                  papartm(i,j,k)= tms(i,j,k)*papartm(i,j,k)*usnjm
     &                     +(1-tms(i,j,k))*spval
                  padissm(i,j,k)= tms(i,j,k)*padissm(i,j,k)*usnjm
     &                     +(1-tms(i,j,k))*spval
                  thpartm(i,j,k)= tms(i,j,k)*thpartm(i,j,k)*usnjm
     &                     +(1-tms(i,j,k))*spval
                  thdissm(i,j,k)= tms(i,j,k)*thdissm(i,j,k)*usnjm
     &                     +(1-tms(i,j,k))*spval
                  pflxcam(i,j,k)= tms(i,j,k)*pflxcam(i,j,k)*usnjm
     &                     +(1-tms(i,j,k))*spval
                  pflxpom(i,j,k)= tms(i,j,k)*pflxpom(i,j,k)*usnjm
     &                     +(1-tms(i,j,k))*spval
#endif
#if ( CORAL == 1 )
                  coraream(i,j,k)= tms(i,j,k)*coraream(i,j,k)*usnjm
     &                     +(1-tms(i,j,k))*spval
                  corprodm(i,j,k)= tms(i,j,k)*corprodm(i,j,k)*usnjm
     &                     +(1-tms(i,j,k))*spval
                  cormassm(i,j,k)= tms(i,j,k)*cormassm(i,j,k)*usnjm
     &                     +(1-tms(i,j,k))*spval
                  omegam(i,j,k)= tms(i,j,k)*omegam(i,j,k)*usnjm
     &                     +(1-tms(i,j,k))*spval
                  oco3m(i,j,k)= tms(i,j,k)*oco3m(i,j,k)*usnjm
     &                     +(1-tms(i,j,k))*spval
                  tau_bleachm(i,j,k)=tms(i,j,k)*tau_bleachm(i,j,k)*usnjm
     &                     +(1-tms(i,j,k))*spval
                  DHWm(i,j,k)=tms(i,j,k)*DHWm(i,j,k)*usnjm
     &                     +(1-tms(i,j,k))*spval
                  PHm(i,j,k)=tms(i,j,k)*PHm(i,j,k)*usnjm
     &                     +(1-tms(i,j,k))*spval
#endif
#if ( OOISO == 1 )
                  oo2_2m(i,j,k)= tms(i,j,k)*oo2_2m(i,j,k)*usnjm
     &                     +(1-tms(i,j,k))*spval
                  oo2_3m(i,j,k)= tms(i,j,k)*oo2_3m(i,j,k)*usnjm
     &                     +(1-tms(i,j,k))*spval
                  oo2_4m(i,j,k)= tms(i,j,k)*oo2_4m(i,j,k)*usnjm
     &                     +(1-tms(i,j,k))*spval
#endif


#endif

                enddo
              enddo
            enddo
            k=kmax+1
            do j=1,jmax
              do i=1,imax
                wmm(i,j,k)= tms(i,j,k-1)*wmm(i,j,k)*usnjm*86400.
     &                 +(1-tms(i,j,k-1))*spval
#if ( OCYCC == 1 )
                 aqpco2m(i,j)= tms(i,j,k-1)*aqpco2m(i,j)*usnjm
     &                 +(1-tms(i,j,k-1))*spval
                 tpp_mam(i,j)= tms(i,j,k-1)*tpp_mam(i,j)*usnjm
     &                 +(1-tms(i,j,k-1))*spval
                 caco3mm(i,j)= tms(i,j,k-1)*caco3mm(i,j)*usnjm
     &                 +(1-tms(i,j,k-1))*spval
#endif
#if ( F_PALAEO_FWF == 1 || CONSEAU == 1 )
                 outFWFm(i,j)= tms(i,j,k-1)*outFWFm(i,j)*usnjm
     &                 +(1-tms(i,j,k-1))*spval
#endif
              enddo
            enddo
!
! dump of 3d monthly means
!
            if (lcdfmon) then
              if (lmmtemp) then
                status=write_3d(ncidm,IMtemp,tmm,kmax  ,iMtimrec)
              endif
              if (lmmsalt)
     &          status=write_3d(ncidm,IMsalt,smm,kmax  ,iMtimrec)
              if (lmmu   )
     &          status=write_3d(ncidm,IMu   ,umm,kmax  ,iMtimrec)
              if (lmmv   )
     &          status=write_3d(ncidm,IMv   ,vmm,kmax  ,iMtimrec)
              if (lmmw   )
     &          status=write_3d(ncidm,IMw   ,wmm,kmax+1,iMtimrec)
#if ( ISOOCN >= 1)
              if (lmmoo17)
     &          status=write_3d(ncidm,IMoo17,oo17m,kmax  ,iMtimrec)
              if (lmmoo18)
     &          status=write_3d(ncidm,IMoo18,oo18m,kmax  ,iMtimrec)
              if (lmmoohd)
     &          status=write_3d(ncidm,IMoohd,oohdm,kmax  ,iMtimrec)
#endif
#if ( OCYCC == 1 )
              if (lmmodoc)
     &          status=write_3d(ncidm,IModoc,odocm,kmax  ,iMtimrec)
              if (lmmodocs)
     &          status=write_3d(ncidm,IModocs,odocsm,kmax  ,iMtimrec)
              if (lmmodic)
     &          status=write_3d(ncidm,IModic,odicm,kmax  ,iMtimrec)
              if (lmmopo4)
     &          status=write_3d(ncidm,IMopo4,opo4m,kmax  ,iMtimrec)
#if ( OXNITREUX == 1 )
              if (lmmon2o)
     &          status=write_3d(ncidm,IMon2o,on2om,kmax  ,iMtimrec)
#else
              if (lmmono3)
     &          status=write_3d(ncidm,IMono3,ono3m,kmax  ,iMtimrec)
#endif
              if (lmmoalk)
     &          status=write_3d(ncidm,IMoalk,oalkm,kmax  ,iMtimrec)
              if (lmmoo2)
     &          status=write_3d(ncidm,IMoo2,oo2m,kmax  ,iMtimrec)
              if (lmmoc13) then
                status=write_3d(ncidm,IMoc13,oc13m,kmax  ,iMtimrec)
              endif
              if (lmmodoc13)
     &          status=write_3d(ncidm,IModoc13,odoc13m,kmax  ,iMtimrec)
              if (lmmodocs13)
     &          status=write_3d(ncidm,IModocs13,odocs13m,kmax,iMtimrec)
              if (lmmoc14)
     &          status=write_3d(ncidm,IMoc14,oc14m,kmax  ,iMtimrec)

              if (lmmaqpco2) then
                status=write_2dnm(ncidm,IMaqpco2,aqpco2m,iMtimrec)
              endif
              if (lmmtpp_ma)
     &          status=write_2dnm(ncidm,IMtpp_ma,tpp_mam,iMtimrec)
              if (lmmcaco3m)
     &          status=write_2dnm(ncidm,IMcaco3m,caco3mm,iMtimrec)

#if ( PATH >= 1 )
              if (lmmpapart) then
                status=write_3d(ncidm,IMpapart,papartm,kmax  ,iMtimrec)
              endif
              if (lmmpadiss) then
                status=write_3d(ncidm,IMpadiss,padissm,kmax  ,iMtimrec)
              endif
              if (lmmthpart) then
                status=write_3d(ncidm,IMthpart,thpartm,kmax  ,iMtimrec)
              endif
              if (lmmthdiss) then
                status=write_3d(ncidm,IMthdiss,thdissm,kmax  ,iMtimrec)
              endif
              if (lmmpflxca) then
                status=write_3d(ncidm,IMpflxca,pflxcam,kmax  ,iMtimrec)
              endif
              if (lmmpflxpo) then
                status=write_3d(ncidm,IMpflxpo,pflxpom,kmax  ,iMtimrec)
              endif
#endif
#if ( CORAL == 1 )
              if (lmmcorarea)
     &         status=write_3d(ncidm,IMcorarea,coraream,kmax  ,iMtimrec)
              if (lmmcorprod)
     &         status=write_3d(ncidm,IMcorprod,corprodm,kmax  ,iMtimrec)
              if (lmmcormass)
     &         status=write_3d(ncidm,IMcormass,cormassm,kmax  ,iMtimrec)
              if (lmmomega)
     &          status=write_3d(ncidm,IMomega,omegam,kmax     ,iMtimrec)
              if (lmmoco3)
     &          status=write_3d(ncidm,IMoco3,oco3m,kmax       ,iMtimrec)
              if (lmmtaub)
     &           status=write_3d(ncidm,IMtaub,tau_bleachm,kmax,iMtimrec)
              if (lmmdhw)
     &           status=write_3d(ncidm,IMdhw,DHWm,kmax        ,iMtimrec)
              if (lmmdhw)
     &           status=write_3d(ncidm,IMdhw,PHm,kmax         ,iMtimrec)
#endif
#if ( OOISO == 1 )
              if (lmmoo2_2)
     &          status=write_3d(ncidm,IMoo2_2,oo2_2m,kmax  ,iMtimrec)
              if (lmmoo2_3)
     &          status=write_3d(ncidm,IMoo2_3,oo2_3m,kmax  ,iMtimrec)
              if (lmmoo2_4)
     &          status=write_3d(ncidm,IMoo2_4,oo2_4m,kmax  ,iMtimrec)
#endif
#endif
#if ( F_PALAEO_FWF == 1 || CONSEAU == 1 )
              if (lmmoutFWF)
     &          status=write_2dnm(ncidm,IMoutFWF,outFWFm,iMtimrec)
#endif
            endif
!
!            ntotoc=999
!            write(cresaldat_id) ja,jmois
!            write(cresaldat_id) ntotoc
!            write(cresaldat_id) umm
!            write(cresaldat_id) vmm
!            write(cresaldat_id) tmm
!            write(cresaldat_id) smm
!
!            do j=2,jmax-1
!             do i=2,imax-1
!              do k=1,kmax
            do k=1,kmax
              do j=1,jmax
                do i=1,imax
!DFG end
                  umm(i,j,k)= 0.0
                  vmm(i,j,k)= 0.0
                  tmm(i,j,k)= 0.0
                  smm(i,j,k)= 0.0
#if ( ISOOCN >= 1)
                  oo17m(i,j,k)= 0.0
                  oo18m(i,j,k)= 0.0
                  oohdm(i,j,k)= 0.0
#endif
#if ( OCYCC == 1 )
                  odocm(i,j,k)= 0.0
                  odocsm(i,j,k)= 0.0
                  odicm(i,j,k)= 0.0
                  opo4m(i,j,k)= 0.0
#if ( OXNITREUX == 1 )
                  on2om(i,j,k)= 0.0
#else
                  ono3m(i,j,k)= 0.0
#endif
                  oalkm(i,j,k)= 0.0
                  oo2m(i,j,k)= 0.0
                  oc13m(i,j,k)= 0.0
                  odoc13m(i,j,k)= 0.0
                  odocs13m(i,j,k)= 0.0
                  oc14m(i,j,k)= 0.0
#if ( PATH >= 1 )
                  papartm(i,j,k) = 0.0
                  padissm(i,j,k) = 0.0
                  thpartm(i,j,k) = 0.0
                  thdissm(i,j,k) = 0.0
                  pflxcam(i,j,k) = 0.0
                  pflxpom(i,j,k) = 0.0
#endif
#if ( CORAL == 1 )
                  coraream(i,j,k)= 0.0
                  corprodm(i,j,k)= 0.0
                  cormassm(i,j,k)= 0.0
                  omegam(i,j,k)= 0.0
                  oco3m(i,j,k)= 0.0
                  tau_bleachm(i,j,k)= 0.0
                  DHWm(i,j,k)= 0.0
                  PHm(i,j,k)= 0.0
#endif
#if ( OOISO == 1 )
                  oo2_2m(i,j,k)= 0.0
                  oo2_3m(i,j,k)= 0.0
                  oo2_4m(i,j,k)= 0.0
#endif
#endif
                enddo
              enddo
            enddo

#if ( OCYCC == 1 )
            aqpco2m(:,:)= 0.0
            tpp_mam(:,:) = 0.0
            caco3mm(:,:) = 0.0
#endif
#if ( F_PALAEO_FWF == 1 || CONSEAU == 1 )
            outFWFm(:,:) = 0.0
#endif

            do k=1,kmax+1
              do j=1,jmax
                do i=1,imax
                  wmm(i,j,k)= 0.0
                enddo
              enddo
            enddo
          endif
!
!DFG start output of 2d monthly means
!
! We need to care about the size of the floats defined for
! netcdf-output. Whereas everything in the netcdf-file is
! defined as double precision, cmoymo is single (even that
! may sometimes yield problems of accuracy...). Nevertheless,
! we someimes need to convert units anyway, so we simply
! copy all the data to the double precision work arry f2d
! first...
!
          if (lcdfmon) then
!
! ice model data:
!
! - snow cover above ice
!
            if (lmmhsn) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=cmoymo(i,j, 1)
                enddo
              enddo
              status=write_2d(ncidm,IMhsn  ,
     &                        f2d,spval,.true. ,.false.,iMtimrec)
            endif
!
! - snow precipitation
!mab: +1 (DISTICB)
            if (lmmsnow) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=cmoymo(i,j,49)
                enddo
              enddo
              status=write_2d(ncidm,IMsnow ,
     &                        f2d,spval,.true. ,.false.,iMtimrec)
            endif
!ccccc
!
! - ice thickness
!
            if (lmmhic) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=cmoymo(i,j, 2)
                enddo
              enddo
              status=write_2d(ncidm,IMhic  ,
     &                        f2d,spval,.true. ,.false.,iMtimrec)
            endif
!
! - ice production
!
            if (lmmhicp) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=cmoymo(i,j, 3)
                enddo
              enddo
              status=write_2d(ncidm,IMhicp ,
     &                        f2d,spval,.true. ,.false.,iMtimrec)
            endif
!
! - lead fraction
!
            if (lmmalbq) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=cmoymo(i,j, 4)
                enddo
              enddo
              status=write_2d(ncidm,IMalbq ,
     &                        f2d,spval,.true. ,.false.,iMtimrec)
            endif
!
! - ice temperature
!
            if (lmmtice) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=cmoymo(i,j, 5)
                enddo
              enddo
              status=write_2d(ncidm,IMtice ,
     &                        f2d,spval,.true. ,.false.,iMtimrec)
            endif
!
! - oceanic heat flux at ice base
!
            if (lmmfb) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=cmoymo(i,j, 6)
                enddo
              enddo
              status=write_2d(ncidm,IMfb ,
     &                        f2d,spval,.true. ,.false.,iMtimrec)
            endif
!
! - zonal ice velocity
!
            if (lmmuice) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=cmoymo(i,j, 7)
                enddo
              enddo
              status=write_2d(ncidm,IMuice ,
     &                        f2d,spval,.false.,.true. ,iMtimrec)
            endif
!
! -meridional ice velocity
!
            if (lmmvice) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=cmoymo(i,j, 8)
                enddo
              enddo
              status=write_2d(ncidm,IMvice ,
     &                        f2d,spval,.false.,.true. ,iMtimrec)
            endif
!
! - meltflux of icebergs
!
            if (lmmfice) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=cmoymo(i,j, 9)
                enddo
              enddo
              status=write_2d(ncidm,IMfice ,
     &                        f2d,spval,.false.,.true. ,iMtimrec)
            endif
!c - distribution of icebergs
!
            if (lmmfdist) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=cmoymo(i,j, 10)
                enddo
              enddo
              status=write_2d(ncidm,IMfdist ,
     &                        f2d,spval,.false.,.true. ,iMtimrec)
            endif
!
! 2-dimensional ocean (or ice/ocean) data:
!
! - sea surface height
!
            if (lmmssh ) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=etamo(i,j)
                enddo
              enddo
              status=write_2d(ncidm,IMssh  ,
     &                        f2d,spval,.true. ,.false.,iMtimrec)
            endif
!
! - ubar
!
            if (lmmubar) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=ubmo(i,j)
                enddo
              enddo
              status=write_2d(ncidm,IMubar ,
     &                        f2d,spval,.false.,.true. ,iMtimrec)
            endif
!
! - vbar
!
            if (lmmvbar) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=vbmo(i,j)
                enddo
              enddo
              status=write_2d(ncidm,IMvbar ,
     &                        f2d,spval,.false.,.true. ,iMtimrec)
            endif
!
! - sea surface temperature
!
! mab: added +1 to all the calculations made below  (IMF) +1 (DISTICB)
            if (lmmsst ) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=cmoymo(i,j,11)
                enddo
              enddo
              status=write_2d(ncidm,IMsst  ,
     &                        f2d,spval,.true. ,.false.,iMtimrec)
            endif
!
! - sea surface salinity
!
            if (lmmsss ) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=cmoymo(i,j,12)
                enddo
              enddo
              status=write_2d(ncidm,IMsss  ,
     &                        f2d,spval,.true. ,.false.,iMtimrec)
            endif
!
! - depth of mixed layer
!
            if (lmmzmix) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=cmoymo(i,j,15)
                enddo
              enddo
              status=write_2d(ncidm,IMzmix ,
     &                        f2d,spval,.true. ,.false.,iMtimrec)
            endif
!
! - depth of convection
!
            if (lmmzcnv) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=cmoymo(i,j,16)
                enddo
              enddo
              status=write_2d(ncidm,IMzcnv ,
     &                        f2d,spval,.true. ,.false.,iMtimrec)
            endif
!
! - surface heat flux
!
            if (lmmshflx) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=cmoymo(i,j,17)
                enddo
              enddo
              status=write_2d(ncidm,IMshflx,
     &                        f2d,spval,.true. ,.false.,iMtimrec)
            endif
!
! - surface freshwater flux
!
            if (lmmsfflx) then
!
! cmoymo(i,j,18) is in cm/day = 3.6m/year
!
              fac=1.0/3.6
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=cmoymo(i,j,18)*fac
                enddo
              enddo
              status=write_2d(ncidm,IMsfflx,
     &                        f2d,spval,.true. ,.false.,iMtimrec)
            endif
!
! - G-M slope parameter
!
            if (lmmmsl) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=cmoymo(i,j,19)
                enddo
              enddo
              status=write_2d(ncidm,IMmsl,
     &                        f2d,spval,.true. ,.false.,iMtimrec)
            endif
!
! - zonal windstress component
!   CLIO3 uses dynes/cm**2 for the stresses, but we want N/m**2 !
!
            if (lmmtx) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=cmoymo(i,j,38)*0.1
                enddo
              enddo
              status=write_2d(ncidm,IMtx,
     &                        f2d,spval,.false.,.true. ,iMtimrec)
            endif
!
! - meridional windstress component
!
            if (lmmty) then
              do j=1,jmax
                do i=1,imax
                  f2d(i,j)=cmoymo(i,j,39)*0.1
                enddo
              enddo
              status=write_2d(ncidm,IMty,
     &                        f2d,spval,.false.,.true. ,iMtimrec)
            endif
!
! sync file (such that user may access it with analysis tools while
! model is still running)
!
            status=nf_sync(ncidm)
!
! close netcdf-file at end of run or when record limit is reached
!
            if (numit.eq.nlast) status=nf_close(ncidm)
            if (iMtimrec.eq.maxmrecs) then
              status=nf_close(ncidm)
              iMtimrec=0
            endif
          endif
!
        endif
!
!                       UPDATE jmois.
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--2.4. CUMULATION FOR GLOBAL YEARLY OUTPUTS AND ALL PERIOD MEAN
!---------------------------------------------------------------
!
        do n=1,nmalef
          do j=1,jmax
            do i=1,imax
              cmoymap(i,j,n,jmois)=cmoymap(i,j,n,jmois)+
     &                                 cmoymo(i,j,ncor(n))
            enddo
          enddo
        enddo
        nmap(jmois)=nmap(jmois)+1

        if (mod(ja,nwa).eq.nwtest) then
!
          do n=1,noumef
            do j=1,jmax
              do i=1,imax
                cmoyan(i,j,n)=cmoyan(i,j,n)+cmoymo(i,j,n)
              enddo
            enddo
          enddo
!
!DFG start accumulation of 2d momentum and SSH
!
          do j=1,jmax
            do i=1,imax
              etaan(i,j) = etaan(i,j) + etamo(i,j)
               uban(i,j) =  uban(i,j) +  ubmo(i,j)
               vban(i,j) =  vban(i,j) +  vbmo(i,j)
            enddo
          enddo
!DFG end
        endif
!
!                       UPDATE jmois.
        jmois = jmois+1
        if (jmois.gt.12) jmois = 1
!DFG start
!
! now that we have updated all annual sums, we can clear all
! monthly means
!
!DFG end
        do n=1,noumef
          do j=1,jmax
            do i=1,imax
              cmoymo(i,j,n) =0.0
            enddo
          enddo
        enddo
!-----
!DFG start cleanup
!
        do j=1,jmax
          do i=1,imax
            etamo(i,j) = 0.0
             ubmo(i,j) = 0.0
             vbmo(i,j) = 0.0
          enddo
        enddo
!DFG end

!-----
      endif

      if (mod(int(xjour),int(yeaday)).eq.0) then

      zfluxm=zfluxmt/yeaday
      zfluxms=zfluxmts/yeaday
#if ( ISOOCN >= 1 )
      write(correcwdat_id,'(6E25.18)') zfluxm,zfluxms
#else
      write(correcwdat_id,*) zfluxm,zfluxms
#endif
      zfluxmt=0.0
      zfluxmts=0.0

      endif


!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!-- WRITE GLOBAL YEARLY OUTPUTS.
!
      if (mod(ja,nwa).eq.nwtest.and.ijour.eq.int(yeaday)
     &               .and.ijour.ne.ijour1) then
!
        usanms = 1.0/real((12-(jmoisi-1)*(1-jmoifl)))
        jmoifl = 1
        do n=1,noumef
          do j=1,jmax
            do i=1,imax
              zindo        = tms(i,j,ks2)
              cmoyan(i,j,n)=zindo*cmoyan(i,j,n)*usanms+(1.0-zindo)*spvr
            enddo
          enddo
        enddo
!
!mab: commented, otherwise too much mouchard output
!        if (nn99.eq.2) write(mouchard_id,'(2A,3I6)') 'Write annual mean ',
!     &    ' on file cresua.out ; year,month,day=', ja, jmois, ijour
!mab: end of commented part
        write(cresuadat_id) ja,jmois
        write(cresuadat_id) nmaef
        do n=1,noumef
          if (nma(n).eq.1) then
            write(cresuadat_id) n,ntypou(n),titn(n)
            do j=1,jmax
              write(cresuadat_id) (cmoyan(i,j,n),i=1,imax)
            enddo
          endif
        enddo
!
!DFG start netcdf output of 2d annual fields
!
        do j=1,jmax
          do i=1,imax
            etaan(i,j) = etaan(i,j)*usanms
             uban(i,j) =  uban(i,j)*usanms
             vban(i,j) =  vban(i,j)*usanms
          enddo
        enddo
        if (lcdfann) then
!
! Do we need to create a new NetCDF-file?
!
          if (iAtimrec.eq.0) then
!
! generate a unique filename
!
! --- BdB 05-2019: adapted name equivalent to atmosphere output
            write(avgname,'(a,a,a,a,i6.6,a)')
     &        'outputdata/ocean/','CLIO3_',
     &        TRIM(refexp),'_ann_avrg_',
     &        int(time_in_years(ja-init_year_clio)),'.nc'
!
! -- pass annual mean var switches to prep_out
!
            do i=1,ldim
              lvcdf(i)=lacdf(i)
            enddo
            do i=1,ldim+2
              icdfNvars(i)=0
              icdfAvars(i)=0
            enddo
            call prep_out(ncida,'a',TRIM(avgname))
!
! -- save IDs of desired variables
!
            do i=1,ldim+2
              icdfAvars(i)=icdfNvars(i)
            enddo
          endif
!
! Update time record. Data are dumped at the end of the year, the
! timestamp represents the entire averaging period. See also the
! timerecord definition with the 360 day calendar and the
! 360 day averaging period attributes.
!
          !rtim=dfloat(ja)
! ---     BdB 05-2019: replaced with global time value
          rtim=time_in_years(ja-init_year_clio+1)

          iAtimrec=iAtimrec+1
          status=nf_put_var1_double(ncida,IAtime,iAtimrec,rtim)
!
! ice model data:
!
! - snow cover above ice
!
          if (lamhsn) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=cmoyan(i,j, 1)
              enddo
            enddo
            status=write_2d(ncida,IAhsn  ,
     &                      f2d,spval,.true. ,.false.,iAtimrec)
          endif
!
! - snow precipitation
!
! mab: add +1 to all the calculations made below (if cmoyan(i,j,>=9)(IMF) +1(DISTICB)
          if (lamsnow) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=cmoyan(i,j,49)
              enddo
            enddo
            status=write_2d(ncida,IAsnow ,
     &                      f2d,spval,.true. ,.false.,iAtimrec)
          endif
!
! - ice thickness
!
          if (lamhic) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=cmoyan(i,j, 2)
              enddo
            enddo
            status=write_2d(ncida,IAhic  ,
     &                      f2d,spval,.true. ,.false.,iAtimrec)
          endif
!
! - ice production
!
          if (lamhicp) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=cmoyan(i,j, 3)
              enddo
            enddo
            status=write_2d(ncida,IAhicp ,
     &                      f2d,spval,.true. ,.false.,iAtimrec)
          endif
!
! - lead fraction
!
          if (lamalbq) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=cmoyan(i,j, 4)
              enddo
            enddo
            status=write_2d(ncida,IAalbq ,
     &                      f2d,spval,.true. ,.false.,iAtimrec)
          endif
!
! - ice temperature
!
          if (lamtice) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=cmoyan(i,j, 5)
              enddo
            enddo
            status=write_2d(ncida,IAtice ,
     &                      f2d,spval,.true. ,.false.,iAtimrec)
          endif
!
! - oceanic heat flux at ice base
!
          if (lamfb) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=cmoyan(i,j, 6)
              enddo
            enddo
            status=write_2d(ncida,IAfb ,
     &                      f2d,spval,.true. ,.false.,iAtimrec)
          endif
!
! - zonal ice velocity
!
          if (lamuice) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=cmoyan(i,j, 7)
              enddo
            enddo
            status=write_2d(ncida,IAuice ,
     &                      f2d,spval,.false.,.true. ,iAtimrec)
          endif
!
! -meridional ice velocity
!
          if (lamvice) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=cmoyan(i,j, 8)
              enddo
            enddo
            status=write_2d(ncida,IAvice ,
     &                      f2d,spval,.false.,.true. ,iAtimrec)
          endif
!
! - meltflux of icebergs
!
          if (lamfice) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=cmoyan(i,j, 9)
              enddo
            enddo
            status=write_2d(ncida,IAfice ,
     &                      f2d,spval,.false.,.true. ,iAtimrec)
          endif
! -spatial distribution of icebergs
!
          if (lamfdist) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=cmoyan(i,j,10)
              enddo
            enddo
            status=write_2d(ncida,IAfdist ,
     &                      f2d,spval,.false.,.true. ,iAtimrec)
          endif

! 2-dimensional ocean (or ice/ocean) data:
!
! - sea surface height
!
          if (lamssh ) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=etaan(i,j)
              enddo
            enddo
            status=write_2d(ncida,IAssh  ,
     &                      f2d,spval,.true. ,.false.,iAtimrec)
          endif
!
! - ubar
!
          if (lamubar) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=uban(i,j)
              enddo
            enddo
            status=write_2d(ncida,IAubar ,
     &                      f2d,spval,.false.,.true. ,iAtimrec)
          endif
!
! - vbar
!
          if (lamvbar) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=vban(i,j)
              enddo
            enddo
            status=write_2d(ncida,IAvbar ,
     &                      f2d,spval,.false.,.true. ,iAtimrec)
          endif

!
!mab: added +1 in all the following to take into account icbdist
! - sea surface temperature
!
          if (lamsst ) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=cmoyan(i,j,11)
              enddo
            enddo
            status=write_2d(ncida,IAsst  ,
     &                      f2d,spval,.true. ,.false.,iAtimrec)
          endif
!
! - sea surface salinity
!
          if (lamsss ) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=cmoyan(i,j,12)
              enddo
            enddo
            status=write_2d(ncida,IAsss  ,
     &                      f2d,spval,.true. ,.false.,iAtimrec)
          endif
!
! - depth of mixed layer
!
          if (lamzmix) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=cmoyan(i,j,15)
              enddo
            enddo
            status=write_2d(ncida,IAzmix ,
     &                      f2d,spval,.true. ,.false.,iAtimrec)
          endif
!
! - depth of convection
!
          if (lamzcnv) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=cmoyan(i,j,16)
              enddo
            enddo
            status=write_2d(ncida,IAzcnv ,
     &                      f2d,spval,.true. ,.false.,iAtimrec)
          endif
!
! - surface heat flux
!
          if (lamshflx) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=cmoyan(i,j,17)
              enddo
            enddo
            status=write_2d(ncida,IAshflx,
     &                      f2d,spval,.true. ,.false.,iAtimrec)
          endif
!
! - surface freshwater flux
!
          if (lamsfflx) then
            fac=1.0/3.6
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=cmoyan(i,j,18)*fac
              enddo
            enddo
            status=write_2d(ncida,IAsfflx,
     &                      f2d,spval,.true. ,.false.,iAtimrec)
          endif
!
! - G-M slope parameter
!
          if (lammsl) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=cmoyan(i,j,19)
              enddo
            enddo
            status=write_2d(ncida,IAmsl,
     &                      f2d,spval,.true. ,.false.,iAtimrec)
          endif
!
! - zonal windstress component
!
          if (lamtx) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=cmoyan(i,j,38)
              enddo
            enddo
            status=write_2d(ncida,IAtx,
     &                      f2d,spval,.false.,.true. ,iAtimrec)
          endif
!
! - meridional windstress component
!
          if (lamty) then
            do j=1,jmax
              do i=1,imax
                f2d(i,j)=cmoyan(i,j,39)
              enddo
            enddo
            status=write_2d(ncida,IAty,
     &                      f2d,spval,.false.,.true. ,iAtimrec)
          endif
        endif
!
! cleanup of annual files that have been written to disk meanwhile
!
        do j=1,jmax
          do i=1,imax
            etaan(i,j) = 0.0
             uban(i,j) = 0.0
             vban(i,j) = 0.0
          enddo
        enddo
!DFG end netcdf output
!
        do n=1,noumef
          do j=1,jmax
            do i=1,imax
              cmoyan(i,j,n)=0.0
            enddo
          enddo
        enddo

! streamFUNCTION part
        if (nstreamout.eq.2) call streamfunc(1,usanms*usnjm,nstreamout)
!DFG start streamFUNCTION part
!
! -- time average
        usnja=1./12.
        do j=1,jmtt
          do nb=0,nbsmax
            do k=0,kmtt
              moc_a(nb,j,k) =  moc_a(nb,j,k)*usnja
            enddo
            do k=0,kfloor(j,nb)
              moc_a(nb,j,k) = spval
            enddo
            mht_a(nb,j) = mht_a(nb,j)*usnja
            mst_a(nb,j) = mst_a(nb,j)*usnja
          enddo
        enddo
!
! -- apply mask
        do nb=1,nbsmax
          do j=1,16
            do k=0,kmtt
              moc_a(nb,j,k) = spval
            enddo
            mht_a(nb,j) = spval
            mst_a(nb,j) = spval
          enddo
        enddo
!
! --write to netcdf file
        if (lcdfann) then
          start(1)=1
          start(2)=1
          start(3)=1
          start(4)=iAtimrec
          count(1)=nbsmax+1
          count(2)=jmtt
          count(3)=kmtt
          count(4)=1
          if (lammoc ) status=nf_put_vara_double(ncida,IAmoc ,
     &                   start,count, moc_a)
          start(3)=iAtimrec
          start(4)=0
          count(3)=1
          count(4)=0
          if (lammht) status=nf_put_vara_double(ncida,IAmht,
     &                   start,count,mht_a)
          if (lammst) status=nf_put_vara_double(ncida,IAmst,
     &                   start,count,mst_a)
        endif
!
! -- wipe arrays
        do j=1,jmtt
          do nb=0,nbsmax
            do k=0,kmtt
               moc_a(nb,j,k) = 0.0
            enddo
            mht_a(nb,j) = 0.0
            mst_a(nb,j) = 0.0
          enddo
        enddo
!DFG end

!DFG start computing annual mean 3d fields
!        if (nwtoa.eq.1) then
        if ((nwtoa.eq.1).or.lcdfann) then
          usnja=1/yeaday
!          do j=2,jmax-1
!           do i=2,imax-1
!            do k=1,kmax
          do k=1,kmax
            do j=1,jmax
              do i=1,imax
!              uma(i,j,k)= tmu(i,j,k)*uma(i,j,k)*usnja
!     &                   +(1-tmu(i,j,k))*spvr
!              vma(i,j,k)= tmu(i,j,k)*vma(i,j,k)*usnja
!     &                   +(1-tmu(i,j,k))*spvr
!              tma(i,j,k)= tms(i,j,k)*tma(i,j,k)*usnja
!     &                   +(1-tms(i,j,k))*spvr
!              sma(i,j,k)= tms(i,j,k)*sma(i,j,k)*usnja
!     &                   +(1-tms(i,j,k))*spvr
                uma(i,j,k)= tmu(i,j,k)*uma(i,j,k)*usnja
     &                     +(1-tmu(i,j,k))*spval
                vma(i,j,k)= tmu(i,j,k)*vma(i,j,k)*usnja
     &                     +(1-tmu(i,j,k))*spval
                tma(i,j,k)= tms(i,j,k)*(tma(i,j,k)*usnja-273.15)
     &                     +(1-tms(i,j,k))*spval
                sma(i,j,k)= tms(i,j,k)*sma(i,j,k)*usnja
     &                     +(1-tms(i,j,k))*spval
#if ( ISOOCN >= 1)
                  oo17a(i,j,k)= tms(i,j,k)*oo17a(i,j,k)*usnja
     &                     +(1-tms(i,j,k))*spval
                  oo18a(i,j,k)= tms(i,j,k)*oo18a(i,j,k)*usnja
     &                     +(1-tms(i,j,k))*spval
                  oohda(i,j,k)= tms(i,j,k)*oohda(i,j,k)*usnja
     &                     +(1-tms(i,j,k))*spval
#endif
#if ( OCYCC == 1 )
                  odoca(i,j,k)= tms(i,j,k)*odoca(i,j,k)*usnja
     &                     +(1-tms(i,j,k))*spval
                  odocsa(i,j,k)= tms(i,j,k)*odocsa(i,j,k)*usnja
     &                     +(1-tms(i,j,k))*spval
                  odica(i,j,k)= tms(i,j,k)*odica(i,j,k)*usnja
     &                     +(1-tms(i,j,k))*spval
                  opo4a(i,j,k)= tms(i,j,k)*opo4a(i,j,k)*usnja
     &                     +(1-tms(i,j,k))*spval
#if ( OXNITREUX == 1 )
                  on2oa(i,j,k)= tms(i,j,k)*on2oa(i,j,k)*usnja
     &                     +(1-tms(i,j,k))*spval
#else
                  ono3a(i,j,k)= tms(i,j,k)*ono3a(i,j,k)*usnja
     &                     +(1-tms(i,j,k))*spval
#endif
                  oalka(i,j,k)= tms(i,j,k)*oalka(i,j,k)*usnja
     &                     +(1-tms(i,j,k))*spval
                  oo2a(i,j,k)= tms(i,j,k)*oo2a(i,j,k)*usnja
     &                     +(1-tms(i,j,k))*spval
                  oc13a(i,j,k)= tms(i,j,k)*oc13a(i,j,k)*usnja
     &                     +(1-tms(i,j,k))*spval
                  odoc13a(i,j,k)= tms(i,j,k)*odoc13a(i,j,k)*usnja
     &                     +(1-tms(i,j,k))*spval
                  odocs13a(i,j,k)= tms(i,j,k)*odocs13a(i,j,k)*usnja
     &                     +(1-tms(i,j,k))*spval
                  oc14a(i,j,k)= tms(i,j,k)*oc14a(i,j,k)*usnja
     &                     +(1-tms(i,j,k))*spval
#if ( PATH >= 1 )
                  paparta(i,j,k)= tms(i,j,k)*paparta(i,j,k)*usnja
     &                     +(1-tms(i,j,k))*spval
                  padissa(i,j,k)= tms(i,j,k)*padissa(i,j,k)*usnja
     &                     +(1-tms(i,j,k))*spval
                  thparta(i,j,k)= tms(i,j,k)*thparta(i,j,k)*usnja
     &                     +(1-tms(i,j,k))*spval
                  thdissa(i,j,k)= tms(i,j,k)*thdissa(i,j,k)*usnja
     &                     +(1-tms(i,j,k))*spval
                  pflxcaa(i,j,k)= tms(i,j,k)*pflxcaa(i,j,k)*usnja
     &                     +(1-tms(i,j,k))*spval
                  pflxpoa(i,j,k)= tms(i,j,k)*pflxpoa(i,j,k)*usnja
     &                     +(1-tms(i,j,k))*spval
#endif
#if ( CORAL == 1 )
                  !for area no need to divide by number of day as
                  !already annual
                  corareaa(i,j,k)= tms(i,j,k)*corareaa(i,j,k)*usnja
                  !corareaa(i,j,k)= tms(i,j,k)*corareaa(i,j,k)
     &                     +(1-tms(i,j,k))*spval
                  corproda(i,j,k)= tms(i,j,k)*corproda(i,j,k)*usnja
                  !corproda(i,j,k)= tms(i,j,k)*corproda(i,j,k)
     &                     +(1-tms(i,j,k))*spval
!              write(*,*) 'ICI outnet'
!              if ((corproda(i,j,k).ne.0)
!     >           .and. (corproda(i,j,k).ne.spval)) then
!                write(*,*) 'corproda in outnet.f',corproda(i,j,k)
!              endif
                  cormassa(i,j,k)= tms(i,j,k)*cormassa(i,j,k)*usnja
     &                     +(1-tms(i,j,k))*spval
                  omegaa(i,j,k)= tms(i,j,k)*omegaa(i,j,k)*usnja
     &                     +(1-tms(i,j,k))*spval
                  oco3a(i,j,k)= tms(i,j,k)*oco3a(i,j,k)*usnja
     &                     +(1-tms(i,j,k))*spval
                  tau_bleacha(i,j,k)= tms(i,j,k)*tau_bleacha(i,j,k)
     &                               *usnja
     &                     +(1-tms(i,j,k))*spval
                  DHWa(i,j,k)= tms(i,j,k)*DHWa(i,j,k)*usnja
     &                     +(1-tms(i,j,k))*spval
                  PHa(i,j,k)= tms(i,j,k)*PHa(i,j,k)*usnja
     &                     +(1-tms(i,j,k))*spval


#endif
#if ( OOISO == 1 )
                  oo2_2a(i,j,k)= tms(i,j,k)*oo2_2a(i,j,k)*usnja
     &                     +(1-tms(i,j,k))*spval
                  oo2_3a(i,j,k)= tms(i,j,k)*oo2_3a(i,j,k)*usnja
     &                     +(1-tms(i,j,k))*spval
                  oo2_4a(i,j,k)= tms(i,j,k)*oo2_4a(i,j,k)*usnja
     &                     +(1-tms(i,j,k))*spval
#endif
#endif
                wma(i,j,k)= tms(i,j,k)*wma(i,j,k)*usnja*86400.
     &                       +(1-tms(i,j,k))*spval
              enddo
            enddo
          enddo
          k=kmax+1
          do j=1,jmax
            do i=1,imax
              wma(i,j,k)= tms(i,j,k-1)*wma(i,j,k)*usnja*86400.
     &                   +(1-tms(i,j,k-1))*spval
#if ( OCYCC == 1 )

              aqpco2a(i,j)= tms(i,j,k-1)*aqpco2a(i,j)*usnja
     &                     +(1-tms(i,j,k-1))*spval
              tpp_maa(i,j)= tms(i,j,k-1)*tpp_maa(i,j)*usnja
     &                     +(1-tms(i,j,k-1))*spval
              caco3ma(i,j)= tms(i,j,k-1)*caco3ma(i,j)*usnja
     &                     +(1-tms(i,j,k-1))*spval
#endif
#if ( F_PALAEO_FWF == 1 || CONSEAU == 1 )
              outFWFa(i,j)= tms(i,j,k-1)*outFWFa(i,j)*usnja
     &                     +(1-tms(i,j,k-1))*spval
#endif
            enddo
          enddo
!
! netcdf output of 3d annual means
!
#define debuggy 0
          if (lcdfann) then
            if (lamtemp) then
#if ( debuggy == 1 )
         WRITE(*,*) "Writing variable lmatemp: IAtemp = ", IAtemp
#endif
              status=write_3d(ncida,IAtemp,tma,kmax  ,iAtimrec)
            endif
            if (lamsalt) then
#if ( debuggy == 1 )
         WRITE(*,*) "Writing variable lmasalt: IAsalt = ", IAsalt
#endif
              status=write_3d(ncida,IAsalt,sma,kmax  ,iAtimrec)
            endif
            if (lamu   )
     &        status=write_3d(ncida,IAu   ,uma,kmax  ,iAtimrec)
            if (lamv   )
     &        status=write_3d(ncida,IAv   ,vma,kmax  ,iAtimrec)
            if (lamw   )
     &        status=write_3d(ncida,IAw   ,wma,kmax+1,iAtimrec)
#if ( ISOOCN >= 1)
              if (lamoo17)
     &          status=write_3d(ncida,IAoo17,oo17a,kmax  ,iAtimrec)
              if (lamoo18)
     &          status=write_3d(ncida,IAoo18,oo18a,kmax  ,iAtimrec)
              if (lamoohd)
     &          status=write_3d(ncida,IAoohd,oohda,kmax  ,iAtimrec)
#endif
#if ( OCYCC == 1 )
              if (lamodoc) then
#if ( debuggy == 1 )
         WRITE(*,*) "Writing variable lamodoc: IAodoc = ", IAodoc
#endif
                status=write_3d(ncida,IAodoc,odoca,kmax  ,iAtimrec)
              endif
              if (lamodocs) then
#if ( debuggy == 1 )
         WRITE(*,*) "Writing variable lamodocs: IAodocs = ", IAodocs
#endif
                status=write_3d(ncida,IAodocs,odocsa,kmax  ,iAtimrec)
              endif
              if (lamodic) then
#if ( debuggy == 1 )
         WRITE(*,*) "Writing variable lamodic: IAodic = ", IAodic
#endif
                status=write_3d(ncida,IAodic,odica,kmax  ,iAtimrec)
              endif
              if (lamopo4) then
#if ( debuggy == 1 )
         WRITE(*,*) "Writing variable lamopo4: IAopo4 = ", IAopo4
#endif
                status=write_3d(ncida,IAopo4,opo4a,kmax  ,iAtimrec)
              endif
#if ( OXNITREUX == 1 )
              if (lamon2o) then
#if ( debuggy == 1 )
         WRITE(*,*) "Writing variable lamon2o: IAono2 = ", IAono2
#endif
                status=write_3d(ncida,IAon2o,on2oa,kmax  ,iAtimrec)
              endif
#else
              if (lamono3) then
#if ( debuggy == 1 )
         WRITE(*,*) "Writing variable lamono3: IAono3 = ", IAono3
#endif
                status=write_3d(ncida,IAono3,ono3a,kmax  ,iAtimrec)
              endif
#endif
              if (lamoalk) then
#if ( debuggy == 1 )
         WRITE(*,*) "Writing variable lamoalk: IAoalk = ", IAoalk
#endif
                status=write_3d(ncida,IAoalk,oalka,kmax  ,iAtimrec)
              endif
              if (lamoo2) then
#if ( debuggy == 1 )
         WRITE(*,*) "Writing variable lamoo2: IAoo2 = ", IAoo2
#endif
                status=write_3d(ncida,IAoo2,oo2a,kmax  ,iAtimrec)
              endif
              if (lamoc13) then
#if ( debuggy == 1 )
         WRITE(*,*) "Writing variable lamoc13: IAoc13 = ", IAoc13
#endif
                status=write_3d(ncida,IAoc13,oc13a,kmax  ,iAtimrec)
              endif
              if (lamodoc13) then
#if ( debuggy == 1 )
         WRITE(*,*)"Writing variable lamodoc13: IAodoc13 = ",IAodoc13
#endif
                status=write_3d(ncida,IAodoc13,odoc13a,kmax  ,iAtimrec)
              endif
              if (lamodocs13) then
#if ( debuggy == 1 )
         WRITE(*,*)"Writing variable lamodocs13: IAodocs13 = ",IAodocs13
#endif
               status=write_3d(ncida,IAodocs13,odocs13a,kmax  ,iAtimrec)
              endif
              if (lamoc14) then
#if ( debuggy == 1 )
         WRITE(*,*) "Writing variable lamoc14: IAoc14 = ", IAoc14
#endif
                status=write_3d(ncida,IAoc14,oc14a,kmax  ,iAtimrec)
              endif
              if (lamaqpco2) then
#if ( debuggy == 1 )
         WRITE(*,*) "Writing variable lamaqpco2: IAaqpco2 = ",IAaqpco2
#endif
                status=write_2dnm(ncida,IAaqpco2,aqpco2a,iAtimrec)
              endif
              if (lamtpp_ma) then
#if ( debuggy == 1 )
         WRITE(*,*) "Writing variable lamtpp_ma: IAtpp_ma = ",IAtpp_ma
#endif
                status=write_2dnm(ncida,IAtpp_ma,tpp_maa,iAtimrec)
              endif
              if (lamcaco3m) then
#if ( debuggy == 1 )
         WRITE(*,*) "Writing variable lamcaco3m: IAcaco3m = ",IAcaco3m
#endif
                status=write_2dnm(ncida,IAcaco3m,caco3ma,iAtimrec)
              endif

#if ( PATH >= 1 )
              if (lampapart) then
#if ( debuggy == 1 )
         WRITE(*,*) "Writing variable lampapart: IApapart = ", IApapart
#endif
                status=write_3d(ncida,IApapart,paparta,kmax  ,iAtimrec)
              endif
              if (lampadiss) then
#if ( debuggy == 1 )
         WRITE(*,*) "Writing variable lampadiss: IApadiss = ", IApadiss
#endif
                status=write_3d(ncida,IApadiss,padissa,kmax  ,iAtimrec)
              endif
              if (lamthpart) then
#if ( debuggy == 1 )
         WRITE(*,*) "Writing variable lamthpart: IAthpart = ", IAthpart
#endif
                status=write_3d(ncida,IAthpart,thparta,kmax  ,iAtimrec)
              endif
              if (lamthdiss) then
#if ( debuggy == 1 )
         WRITE(*,*) "Writing variable lamthdiss: IAthdiss = ", IAthdiss
#endif
                status=write_3d(ncida,IAthdiss,thdissa,kmax  ,iAtimrec)
              endif
              if (lampflxca) then
#if ( debuggy == 1 )
         WRITE(*,*) "Writing variable lampflxca: IApflxca = ", IApflxca
#endif
                status=write_3d(ncida,IApflxca,pflxcaa,kmax  ,iAtimrec)
              endif
              if (lampflxpo) then
#if ( debuggy == 1 )
         WRITE(*,*) "Writing variable lampflxpo: IApflxpo = ", IApflxpo
#endif
                status=write_3d(ncida,IApflxpo,pflxpoa,kmax  ,iAtimrec)
              endif
#endif

#if ( CORAL == 1 )
              if (lamcorarea) then
#if ( debuggy == 1 )
        WRITE(*,*) "Writing variable lamcorarea: IAcorarea = ",IAcorarea
#endif
               status=write_3d(ncida,IAcorarea,corareaa,kmax  ,iAtimrec)
              endif
              if (lamcorprod) then
#if ( debuggy == 1 )
        WRITE(*,*) "Writing variable lamcorprod: IAcorprod = ",IAcorprod
#endif
               status=write_3d(ncida,IAcorprod,corproda,kmax  ,iAtimrec)
               endif

              if (lamcormass) then
#if ( debuggy == 1 )
        WRITE(*,*) "Writing variable lamcormass: IAcormass = ",IAcormass
#endif
               status=write_3d(ncida,IAcormass,cormassa,kmax  ,iAtimrec)
               endif

              if (lamomega) then
#if ( debuggy == 1 )
        WRITE(*,*) "Writing variable lamomega: IAomega = ",IAomega
#endif
               status=write_3d(ncida,IAomega,omegaa,kmax  ,iAtimrec)
               endif

              if (lamoco3) then
#if ( debuggy == 1 )
        WRITE(*,*) "Writing variable lamoco3: IAoco3 = ",IAoco3
#endif
               print*,'Write OCO3 a'
               status=write_3d(ncida,IAoco3,oco3a,kmax  ,iAtimrec)
              endif
              if (lamtaub) then
#if ( debuggy == 1 )
        WRITE(*,*) "Writing variable lamtau: IAtau_bleach =",IAtaub
#endif
               print*, 'Write Tau bleaching a'
               status=write_3d(ncida,IAtaub,tau_bleacha,kmax  ,iAtimrec)
              endif
              if (lamdhw) then
#if ( debuggy == 1 )
        WRITE(*,*) "Writing variable lamdhw: IAdhw =",IAdhw
#endif
               print*, 'Write DHW a'
               status=write_3d(ncida,IAdhw,DHWa,kmax  ,iAtimrec)
              endif
              if (lamph) then
#if ( debuggy == 1 )
        WRITE(*,*) "Writing variable lamph: IAph =",IAph
#endif
               print*, 'Write PH a'
               status=write_3d(ncida,IAph,PHa,kmax  ,iAtimrec)
              endif

#endif


#if ( OOISO == 1 )
              if (lamoo2_2) then
#if ( debuggy == 1 )
         WRITE(*,*) "Writing variable lamoo2_2: IAoo2_2 = ", IAoo2_2
#endif
                status=write_3d(ncida,IAoo2_2,oo2_2a,kmax  ,iAtimrec)
              endif

              if (lamoo2_3) then
#if ( debuggy == 1 )
         WRITE(*,*) "Writing variable lamoo2_3: IAoo2_3 = ", IAoo2_3
#endif
                status=write_3d(ncida,IAoo2_3,oo2_3a,kmax  ,iAtimrec)
              endif

              if (lamoo2_4) then
#if ( debuggy == 1 )
         WRITE(*,*) "Writing variable lamoo2_4: IAoo2_4 = ", IAoo2_4
#endif

                status=write_3d(ncida,IAoo2_4,oo2_4a,kmax  ,iAtimrec)
              endif
#endif

#endif

#if ( F_PALAEO_FWF == 1 || CONSEAU == 1 )
              if (lamoutFWF) then
         WRITE(*,*) "Writing variable lamout_FWF: IAoutFWF =",IAoutFWF
#if ( debuggy == 1 )
         WRITE(*,*) "Writing variable lamout_FWF: IAoutFWF =",IAoutFWF
#endif
                status=write_2dnm(ncida,IAoutFWF,outFWFa,iAtimrec)
              endif

#endif
          endif
!
!          ntotoc=999
!          write(cresaldat_id) ja,jmois
!          write(cresaldat_id) ntotoc
!          write(cresaldat_id) uma
!          write(cresaldat_id) vma
!          write(cresaldat_id) tma
!          write(cresaldat_id) sma
!
!            do j=2,jmax-1
!             do i=2,imax-1
!              do k=1,kmax
          do k=1,kmax
            do j=1,jmax
              do i=1,imax
!DFG end
                uma(i,j,k)= 0.0
                vma(i,j,k)= 0.0
                tma(i,j,k)= 0.0
                sma(i,j,k)= 0.0
#if ( ISOOCN >= 1 )
              oo17a(i,j,k)= 0.0
              oo18a(i,j,k)= 0.0
              oohda(i,j,k)= 0.0
#endif
#if ( OCYCC == 1 )
              odoca(i,j,k)= 0.0
              odocsa(i,j,k)= 0.0

              odoca(i,j,k)= 0.0
              odocsa(i,j,k)= 0.0
              odica(i,j,k)= 0.0
              opo4a(i,j,k)= 0.0
#if ( OXNITREUX == 1 )
              on2oa(i,j,k)= 0.0
#else
              ono3a(i,j,k)= 0.0
#endif
              oalka(i,j,k)= 0.0
              oo2a(i,j,k)= 0.0
              oc13a(i,j,k)= 0.0
              odoc13a(i,j,k)= 0.0
              odocs13a(i,j,k)= 0.0
              oc14a(i,j,k)= 0.0
#if ( PATH >= 1 )
              paparta(i,j,k) = 0.0
              padissa(i,j,k) = 0.0
              thparta(i,j,k) = 0.0
              thdissa(i,j,k) = 0.0
              pflxcaa(i,j,k) = 0.0
              pflxpoa(i,j,k) = 0.0
#endif
#if ( CORAL == 1 )
              corareaa(i,j,k)= 0.0
              corproda(i,j,k)= 0.0
              cormassa(i,j,k)= 0.0
              omegaa(i,j,k)= 0.0
              oco3a(i,j,k)= 0.0
              tau_bleacha(i,j,k)= 0.0
              DHWa(i,j,k)= 0.0
              PHa(i,j,k)= 0.0
#endif
#if ( OOISO )
              oo2_2a(i,j,k)= 0.0
              oo2_3a(i,j,k)= 0.0
              oo2_4a(i,j,k)= 0.0
#endif
#endif
              enddo
            enddo
          enddo
!DFG start
          do k=1,kmax+1
            do j=1,jmax
              do i=1,imax
                wma(i,j,k)= 0.0
              enddo
            enddo
          enddo
!
          status=nf_sync(ncida)
!
! close netcdf-file at end of run or when record limit is reached
!
          if (numit.eq.nlast) status=nf_close(ncida)
          if (iAtimrec.eq.maxarecs) then
            status=nf_close(ncida)
            iAtimrec=0
          endif
!DFG end
        endif
!-----
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!-- WRITE GLOBAL ALL PERIOD OUTPUTS.
      if (numit.eq.nlast.and.nwtal.eq.1) then

        do k=1,12
          do n=1,nmalef
            do j=1,jmax
              do i=1,imax
                zindo        = tms(i,j,ks2)
                cmoymap(i,j,n,k)=zindo*cmoymap(i,j,n,k)
     &                          /max(dfloat(nmap(k)),undim6)
     &                          +(1.0-zindo)*spvr
              enddo
            enddo
          enddo
!mab: commented because of mouchard
!         if (nn99.eq.2) write(mouchard_id,'(2A,3I6)') 'Write monthly mean',
!     &       ' on all period on cresal.out ; yy,mm,dd=',ja,jmois,ijour
!mab: end of commented part
          write(cresaldat_id) ja,k
          write(cresaldat_id) nmalef
          do n=1,nmalef
            write(cresaldat_id) ncor(n),ntypou(ncor(n)),titn(ncor(n))
            do j=1,jmax
              write(cresaldat_id) (cmoymap(i,j,n,k),i=1,imax)
            enddo
          enddo
        enddo
      endif
!
      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine outave -
      end

!DFG start
      SUBROUTINE prep_out(ncid,cflag,avgname)
!
!=======================================================================
!
! This routine creates a new averages NetCDF file, defines all the
! necessary variables and their respective grids. All grid-related
! items and time-recordless data will then be dumped to that file.
!
! input:
!         cflag = 'm' for monthly, 'a' for annual means
!       avgname = name of respective netcdf dataset
!
!       ...and all the logical switches vor requested variables
!       that are passed to this routine via common 'logics'
!
! output:
!          ncid = integer ID of generated netcdf file
!
!        ...and all the variable IDs of the requested fields that are
!        returned to the calling routine via common 'icdfout'
!
!=======================================================================
!

!! START_OF_USE_SECTION

      use const_mod

      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod
      use densit_mod
      use thermo_mod

      use isoslope_mod
      use ice_mod
      use dynami_mod
      use reper_mod

      use newunit_clio_mod, only: clio3_out_id

!! END_OF_USE_SECTION

!! START_OF_INCLUDE_SECTION

#include "type.com"
! [SCRPTCOM] #include "para.com"
! [SCRPTCOM] #include "bloc.com"
! [SCRPTCOM] #include "ice.com"
! [SCRPTCOM] #include "dynami.com"
! [SCRPTCOM] #include "thermo.com"
! [SCRPTCOM] #include "densit.com"
! [SCRPTCOM] #include "isoslope.com"
! [SCRPTCOM] #include "reper.com"
#include "netcdf.inc"

!! END_OF_INCLUDE_SECTION

!
      character*1   cflag
      character*(*) avgname
!
!     integer ldim, ncid, write_2d, write_3d, status
      integer ncid, write_2d, write_3d, status

!-----|--1--------2---------3---------4---------5---------6---------7-|
! dmr --- As of July, 12th, 2013, the cumbersome variable declarations
!          are moved to the following include file (for consistency).
!         Beware that ldim has been moved as well
#include "netcdf_out.com"
!-----|--1--------2---------3---------4---------5---------6---------7-|
!
!-----------------------------------------------------------------------
!
      integer recdim, idim, jdim, kdim,
     &                        idimp1, jdimp1, kdimp1
      integer DimIDs(17), dims(4), start(3), count(3)
!
      dimension rmisc(max(imax,jmax,kmax))
      character*14 cbasins(4)
!
      data cbasins /'World Ocean   ','Indian Ocean  ','Pacific Ocean ',
     &              'Atlantic Ocean'/
!
      spval=-1.0d+32
!
!-----------------------------------------------------------------------
! Create new averages NetCDF files.
!-----------------------------------------------------------------------
!
      status=nf_create(avgname,nf_clobber,ncid)
      if (status.ne.nf_noerr) then
        write(clio3_out_id,*) 'ERROR: unable to create netcdf file '
     &                      ,avgname
        stop
      endif
!
!-----------------------------------------------------------------------
! Define the dimensions of staggered fields. Please note that the data
! are stored in the netcdf dataset on their native gridis, unlike those
! generated by the postprocessing routine cresuint, which performs
! horizontal interpolation onto a regular grid with a spacing of 2.5
! degrees longitude. CLIO is discretized on the Arakawa B-grid, and thus
! we deal with both tracer and momentum point longitudes and latitudes
! in the horizontal, and tracer (and horizontal momentum) depths at the
! centers of the cells and vertical momentum at their bottom in the
! vertical. This adds up to 6 dimensions, plus one for the number of
! tracers, and, finally, one for time. CLIO uses a grid combined of a
! regular longitude latitude grid all over the globe except the North
! Atlantic and Arctic Ocean, where a rotated subgrid is used. The two
! meshes are being merged in the equatorial Atlantic. We thus define
! PSEUDO LONgitudes and LATitudes here...
!-----------------------------------------------------------------------
!
      idim=imax-2
      idimp1=idim+1
      jdim=jmax
      jdimp1=jdim+1
      kdim=kmax
      kdimp1=kdim+1
      kdimp2=kdim+2
      status=nf_def_dim(ncid,'ptlon' ,idim        ,DimIDs( 1))
      status=nf_def_dim(ncid,'pulon' ,idim        ,DimIDs( 2))
      status=nf_def_dim(ncid,'ptlat' ,jdim        ,DimIDs( 3))
      status=nf_def_dim(ncid,'pulat' ,jdim        ,DimIDs( 4))
      status=nf_def_dim(ncid,'tdepth',kdim        ,DimIDs( 5))
      status=nf_def_dim(ncid,'wdepth',kdimp1      ,DimIDs( 6))
      status=nf_def_dim(ncid,'wedges',kdimp2      ,DimIDs( 7))
      status=nf_def_dim(ncid,'corners',4          ,DimIDs( 8))
      status=nf_def_dim(ncid,'sflat'  ,57         ,DimIDs( 9))
      status=nf_def_dim(ncid,'sfdepth',kdimp2     ,DimIDs(10))
      status=nf_def_dim(ncid,'sfedges',kdimp2+1   ,DimIDs(11))
      status=nf_def_dim(ncid,'basidx' ,4          ,DimIDs(12))
      status=nf_def_dim(ncid,'ptlonp',idimp1      ,DimIDs(13))
      status=nf_def_dim(ncid,'pulonp',idimp1      ,DimIDs(14))
      status=nf_def_dim(ncid,'ptlatp',jdimp1      ,DimIDs(15))
      status=nf_def_dim(ncid,'pulatp',jdimp1      ,DimIDs(16))
      status=nf_def_dim(ncid,'time'  ,nf_unlimited,DimIDs(17))
      recdim=DimIDs(17)
!
!-----------------------------------------------------------------------
! Define time-recordless information variables.
!-----------------------------------------------------------------------
!
      if ((cflag.eq.'m').or.(cflag.eq.'M')) status=
     &  nf_put_att_text(ncid,NF_GLOBAL,'title',19,'CLIO3 monthly means')
      if ((cflag.eq.'a').or.(cflag.eq.'A')) status=
     &  nf_put_att_text(ncid,NF_GLOBAL,'title',18,'CLIO3 annual means')
      status=nf_put_att_text(ncid,NF_GLOBAL,'institution',45,
     &                 'Vrije Universiteit Amsterdam, The Netherlands')
      status=nf_put_att_text(ncid,NF_GLOBAL,'source',16,
     &                       'iLOVECLIM, v1.0')
!
! pseudo tracer grid longitudes:
!
      status=nf_def_var(ncid,'ptlon',NF_DOUBLE,1,DimIDs(1),IDptlon)
      status=nf_put_att_text(ncid,IDptlon,'units',12,'degrees_east')
      status=nf_put_att_text(ncid,IDptlon,'long_name',28,
     &                       'pseudo tracer grid longitude')
      status=nf_put_att_text(ncid,IDptlon,'modulo',1,'m')
      status=nf_put_att_text(ncid,IDptlon,'topology',8,'circular')
!
! pseudo tracer grid latitudes:
!
      status=nf_def_var(ncid,'ptlat',NF_DOUBLE,1,DimIDs(3),IDptlat)
      status=nf_put_att_text(ncid,IDptlat,'units',13,'degrees_north')
      status=nf_put_att_text(ncid,IDptlat,'long_name',27,
     &                       'pseudo tracer grid latitude')
!
! pseudo momentum grid longitudes:
!
      status=nf_def_var(ncid,'pulon',NF_DOUBLE,1,DimIDs(2),IDpulon)
      status=nf_put_att_text(ncid,IDpulon,'units',12,'degrees_east')
      status=nf_put_att_text(ncid,IDpulon,'long_name',30,
     &                       'pseudo momentum grid longitude')
      status=nf_put_att_text(ncid,IDpulon,'modulo',1,'m')
      status=nf_put_att_text(ncid,IDpulon,'topology',8,'circular')
!
! pseudo momentum grid latitudes:
!
      status=nf_def_var(ncid,'pulat',NF_DOUBLE,1,DimIDs(4),IDpulat)
      status=nf_put_att_text(ncid,IDpulat,'units',13,'degrees_north')
      status=nf_put_att_text(ncid,IDpulat,'long_name',29,
     &                       'pseudo momentum grid latitude')
!
! tracer (and horizontal momentum) grid depths:
!
      status=nf_def_var(ncid,'tdepth',NF_DOUBLE,1,DimIDs(5),IDtdepth)
      status=nf_put_att_text(ncid,IDtdepth,'units',6,'meters')
      status=nf_put_att_text(ncid,IDtdepth,'long_name',22,
     &                       'depth of tracer points')
      status=nf_put_att_text(ncid,IDtdepth,'point_spacing',6,'uneven')
      status=nf_put_att_text(ncid,IDtdepth,'positive',2,'up')
      status=nf_put_att_text(ncid,IDtdepth,'edges',6,'wdepth')
!
! vertical momentum grid depths (also edges of tracer grid depth:
!
      status=nf_def_var(ncid,'wdepth',NF_DOUBLE,1,DimIDs(6),IDwdepth)
      status=nf_put_att_text(ncid,IDwdepth,'units',6,'meters')
      status=nf_put_att_text(ncid,IDwdepth,'long_name',33,
     &                       'depth of vertical momentum points')
      status=nf_put_att_text(ncid,IDwdepth,'point_spacing',6,'uneven')
      status=nf_put_att_text(ncid,IDwdepth,'positive',2,'up')
      status=nf_put_att_text(ncid,IDwdepth,'edges',6,'wedges')
!
! edges of vertical momentum grid depths:
!
      status=nf_def_var(ncid,'wedges',NF_DOUBLE,1,DimIDs(7),IDwedges)
      status=nf_put_att_text(ncid,IDwedges,'units',6,'meters')
      status=nf_put_att_text(ncid,IDwedges,'long_name',32,
     &                       'edges of vertical momentum boxes')
      status=nf_put_att_text(ncid,IDwedges,'point_spacing',6,'uneven')
      status=nf_put_att_text(ncid,IDwedges,'positive',2,'up')
!
! latitudes at which overturning streamFUNCTIONs are evaluated:
!
      status=nf_def_var(ncid,'sflat',NF_DOUBLE,1,DimIDs(9),IDsflat)
      status=nf_put_att_text(ncid,IDsflat,'units',13,'degrees_north')
      status=nf_put_att_text(ncid,IDsflat,'long_name',8,'latitude')
!
! depths for streamFUNCTION computation:
!
      status=nf_def_var(ncid,'sfdepth',NF_DOUBLE,1,DimIDs(10),IDsfdepth)
      status=nf_put_att_text(ncid,IDsfdepth,'units',6,'meters')
      status=nf_put_att_text(ncid,IDsfdepth,'long_name',20,
     &                       'streamFUNCTION depth')
      status=nf_put_att_text(ncid,IDsfdepth,'point_spacing',6,'uneven')
      status=nf_put_att_text(ncid,IDsfdepth,'positive',2,'up')
      status=nf_put_att_text(ncid,IDsfdepth,'edges',7,'sfedges')
!
! corresponding edges
!
      status=nf_def_var(ncid,'sfedges',NF_DOUBLE,1,DimIDs(11),IDsfedges)
      status=nf_put_att_text(ncid,IDsfedges,'units',6,'meters')
      status=nf_put_att_text(ncid,IDsfedges,'long_name',30,
     &                       'edges of streamFUNCTION depths')
      status=nf_put_att_text(ncid,IDsfedges,'point_spacing',6,'uneven')
      status=nf_put_att_text(ncid,IDsfedges,'positive',2,'up')
!
! basin index axis
!
      status=nf_def_var(ncid,'basidx',NF_INT,1,DimIDs(12),IDbasidx)
      status=nf_put_att_text(ncid,IDbasidx,'long_name',9,'longitude')
!
! basin names:
!
!      status=nf_def_var(ncid,'basins',NF_CHAR,2,14,DimIDs(12),IDbasnam)
!      status=nf_put_att_text(ncid,IDbasnam,'long_name',11,'basin names')
!
!
! true tracer grid longitudes:
!
      dims(1)=DimIDs(1)
      dims(2)=DimIDs(3)
      dims(3)=0
      status=nf_def_var(ncid,'tlon',NF_DOUBLE,2,dims,IDtlon)
      status=nf_put_att_text(ncid,IDtlon,'units',12,'degrees_east')
      status=nf_put_att_text(ncid,IDtlon,'long_name',9,'longitude')
      status=nf_put_att_text(ncid,IDtlon,'point_spacing',6,'uneven')
      status=nf_put_att_text(ncid,IDtlon,'bounds',11,'tlon_bounds')
      status=nf_put_att_double(ncid,IDtlon,
     &                         'missing_value',NF_DOUBLE,1,spval)
!
! corresponding edges
!
      dims(3)=DimIDs(8)
      status=nf_def_var(ncid,'tlon_bounds',NF_DOUBLE,3,dims,IDtlon_e)
      status=nf_put_att_text(ncid,IDtlon_e,'units',12,'degrees_east')
      status=nf_put_att_text(ncid,IDtlon_e,'long_name',15,
     &                       'longitude edges')
      status=nf_put_att_text(ncid,IDtlon_e,'point_spacing',6,'uneven')
      status=nf_put_att_double(ncid,IDtlon_e,
     &                         'missing_value',NF_DOUBLE,1,spval)
!
! true tracer grid longitude bounds:
!
      dims(1)=DimIDs(13)
      dims(2)=DimIDs(15)
      dims(3)=0
      status=nf_def_var(ncid,'tlonp',NF_DOUBLE,2,dims,IDtlonp)
      status=nf_put_att_text(ncid,IDtlonp,'units',12,'degrees_east')
      status=nf_put_att_text(ncid,IDtlonp,'long_name',16,
     &                       'longitude bounds')
      status=nf_put_att_text(ncid,IDtlonp,'point_spacing',6,'uneven')
      status=nf_put_att_double(ncid,IDtlonp,
     &                         'missing_value',NF_DOUBLE,1,spval)
!
! true tracer grid latitudes:
!
      dims(1)=DimIDs(1)
      dims(2)=DimIDs(3)
      dims(3)=0
      status=nf_def_var(ncid,'tlat',NF_DOUBLE,2,dims,IDtlat)
      status=nf_put_att_text(ncid,IDtlat,'units',13,'degrees_north')
      status=nf_put_att_text(ncid,IDtlat,'long_name',8,'latitude')
      status=nf_put_att_text(ncid,IDtlat,'point_spacing',6,'uneven')
      status=nf_put_att_text(ncid,IDtlat,'bounds',11,'tlat_bounds')
      status=nf_put_att_double(ncid,IDtlat,
     &                         'missing_value',NF_DOUBLE,1,spval)
!
! corresponding edges
!
      dims(3)=DimIDs(8)
      status=nf_def_var(ncid,'tlat_bounds',NF_DOUBLE,3,dims,IDtlat_e)
      status=nf_put_att_text(ncid,IDtlat_e,'units',13,'degrees_north')
      status=nf_put_att_text(ncid,IDtlat_e,'long_name',14,
     &                       'latitude edges')
      status=nf_put_att_text(ncid,IDtlat_e,'point_spacing',6,'uneven')
      status=nf_put_att_double(ncid,IDtlat_e,
     &                         'missing_value',NF_DOUBLE,1,spval)
!
! true tracer grid latitude bounds:
!
      dims(1)=DimIDs(13)
      dims(2)=DimIDs(15)
      dims(3)=0
      status=nf_def_var(ncid,'tlatp',NF_DOUBLE,2,dims,IDtlatp)
      status=nf_put_att_text(ncid,IDtlatp,'units',13,'degrees_north')
      status=nf_put_att_text(ncid,IDtlatp,'long_name',15,
     &                       'latitude bounds')
      status=nf_put_att_text(ncid,IDtlatp,'point_spacing',6,'uneven')
      status=nf_put_att_double(ncid,IDtlatp,
     &                         'missing_value',NF_DOUBLE,1,spval)
!
! true momentum grid longitudes:
!
      dims(1)=DimIDs(2)
      dims(2)=DimIDs(4)
      dims(3)=0
      status=nf_def_var(ncid,'ulon',NF_DOUBLE,2,dims,IDulon)
      status=nf_put_att_text(ncid,IDulon,'units',12,'degrees_east')
      status=nf_put_att_text(ncid,IDulon,'long_name',9,'longitude')
      status=nf_put_att_text(ncid,IDulon,'point_spacing',6,'uneven')
      status=nf_put_att_text(ncid,IDulon,'bounds',11,'ulon_bounds')
      status=nf_put_att_double(ncid,IDulon,
     &                         'missing_value',NF_DOUBLE,1,spval)
!
! corresponding edges
!
      dims(3)=DimIDs(8)
      status=nf_def_var(ncid,'ulon_bounds',NF_DOUBLE,3,dims,IDulon_e)
      status=nf_put_att_text(ncid,IDulon_e,'units',12,'degrees_east')
      status=nf_put_att_text(ncid,IDulon_e,'long_name',15,
     &                       'longitude edges')
      status=nf_put_att_text(ncid,IDulon_e,'point_spacing',6,'uneven')
      status=nf_put_att_double(ncid,IDulon_e,
     &                         'missing_value',NF_DOUBLE,1,spval)
!
! true momentum grid longitude bounds:
!
      dims(1)=DimIDs(14)
      dims(2)=DimIDs(16)
      dims(3)=0
      status=nf_def_var(ncid,'ulonp',NF_DOUBLE,2,dims,IDulonp)
      status=nf_put_att_text(ncid,IDulonp,'units',12,'degrees_east')
      status=nf_put_att_text(ncid,IDulonp,'long_name',16,
     &                       'longitude bounds')
      status=nf_put_att_text(ncid,IDulonp,'point_spacing',6,'uneven')
      status=nf_put_att_double(ncid,IDulonp,
     &                         'missing_value',NF_DOUBLE,1,spval)
!
! true momentum grid latitudes:
!
      dims(1)=DimIDs(2)
      dims(2)=DimIDs(4)
      dims(3)=0
      status=nf_def_var(ncid,'ulat',NF_DOUBLE,2,dims,IDulat)
      status=nf_put_att_text(ncid,IDulat,'units',13,'degrees_north')
      status=nf_put_att_text(ncid,IDulat,'long_name',8,'latitude')
      status=nf_put_att_text(ncid,IDulat,'point_spacing',6,'uneven')
      status=nf_put_att_text(ncid,IDulat,'bounds',11,'ulat_bounds')
      status=nf_put_att_double(ncid,IDulat,
     &                         'missing_value',NF_DOUBLE,1,spval)
!
! corresponding edges:
!
      dims(3)=DimIDs(8)
      status=nf_def_var(ncid,'ulat_bounds',NF_DOUBLE,3,dims,IDulat_e)
      status=nf_put_att_text(ncid,IDulat_e,'units',13,'degrees_north')
      status=nf_put_att_text(ncid,IDulat_e,'long_name',14,
     &                       'latitude edges')
      status=nf_put_att_text(ncid,IDulat_e,'point_spacing',6,'uneven')
      status=nf_put_att_double(ncid,IDulat_e,
     &                         'missing_value',NF_DOUBLE,1,spval)
      dims(3)=0
!
! true momentum grid latitude bounds:
!
      dims(1)=DimIDs(14)
      dims(2)=DimIDs(16)
      dims(3)=0
      status=nf_def_var(ncid,'ulatp',NF_DOUBLE,2,dims,IDulatp)
      status=nf_put_att_text(ncid,IDulatp,'units',13,'degrees_north')
      status=nf_put_att_text(ncid,IDulatp,'long_name',15,
     &                       'latitude bounds')
      status=nf_put_att_text(ncid,IDulatp,'point_spacing',6,'uneven')
      status=nf_put_att_double(ncid,IDulatp,
     &                         'missing_value',NF_DOUBLE,1,spval)
!
! angle of rotation of true grid agains pseudo lon-lat grid
!
      dims(1)=DimIDs(2)
      dims(2)=DimIDs(4)
      status=nf_def_var(ncid,'angle',NF_DOUBLE,2,dims,IDangle)
      status=nf_put_att_text(ncid,IDangle,'units',7,'radians')
      status=nf_put_att_text(ncid,IDangle,'long_name',28,
     &                       'clockwise rotation of x-axis')
      status=nf_put_att_text(ncid,IDangle,'point_spacing',6,'uneven')
      status=nf_put_att_double(ncid,IDangle,
     &                         'missing_value',NF_DOUBLE,1,spval)
!      status=nf_put_att_text(ncid,IDangle,'coordinates',9,'ulon ulat')
!
! zonal (x) length of the grid cell sides:
!
      dims(1)=DimIDs(1)
      dims(2)=DimIDs(3)
      status=nf_def_var(ncid,'dxs1',NF_DOUBLE,2,dims,IDdxs1)
      status=nf_put_att_text(ncid,IDdxs1,'units',6,'meters')
      status=nf_put_att_text(ncid,IDdxs1,'long_name',30,
     &                       'zonal length of grid box sides')
      status=nf_put_att_text(ncid,IDdxs1,'point_spacing',6,'uneven')
      status=nf_put_att_double(ncid,IDdxs1,
     &                         'missing_value',NF_DOUBLE,1,spval)
!      status=nf_put_att_text(ncid,IDdxs1,'coordinates',9,'tlon tlat')
!
! meridional (y) length of the grid cell sides:
!
      status=nf_def_var(ncid,'dxs2',NF_DOUBLE,2,dims,IDdxs2)
      status=nf_put_att_text(ncid,IDdxs2,'units',6,'meters')
      status=nf_put_att_text(ncid,IDdxs2,'long_name',35,
     &                       'meridional length of grid box sides')
      status=nf_put_att_text(ncid,IDdxs2,'point_spacing',6,'uneven')
      status=nf_put_att_double(ncid,IDdxs2,
     &                         'missing_value',NF_DOUBLE,1,spval)
!      status=nf_put_att_text(ncid,IDdxs2,'coordinates',9,'tlon tlat')
!
! zonal (x) width at the grid cell centres:
!
      status=nf_def_var(ncid,'dxc1',NF_DOUBLE,2,dims,IDdxc1)
      status=nf_put_att_text(ncid,IDdxc1,'units',6,'meters')
      status=nf_put_att_text(ncid,IDdxc1,'long_name',31,
     &                       'zonal width at grid box centres')
      status=nf_put_att_text(ncid,IDdxc1,'point_spacing',6,'uneven')
      status=nf_put_att_double(ncid,IDdxc1,
     &                         'missing_value',NF_DOUBLE,1,spval)
!      status=nf_put_att_text(ncid,IDdxc1,'coordinates',9,'tlon tlat')
!
! meridional (y) width at the grid cell centres:
!
      status=nf_def_var(ncid,'dxc2',NF_DOUBLE,2,dims,IDdxc2)
      status=nf_put_att_text(ncid,IDdxc2,'units',6,'meters')
      status=nf_put_att_text(ncid,IDdxc2,'long_name',36,
     &                       'meridional width at grid box centres')
      status=nf_put_att_text(ncid,IDdxc2,'point_spacing',6,'uneven')
      status=nf_put_att_double(ncid,IDdxc2,
     &                         'missing_value',NF_DOUBLE,1,spval)
!      status=nf_put_att_text(ncid,IDdxc2,'coordinates',9,'tlon tlat')
!
! horizontal area of grid cells:
!
      status=nf_def_var(ncid,'area',NF_DOUBLE,2,dims,IDarea)
      status=nf_put_att_text(ncid,IDarea,'units',9,'meters**2')
      status=nf_put_att_text(ncid,IDarea,'long_name',31,
     &                       'surface area of tracer grid box')
      status=nf_put_att_text(ncid,IDarea,'point_spacing',6,'uneven')
      status=nf_put_att_double(ncid,IDarea,
     &                         'missing_value',NF_DOUBLE,1,spval)
!      status=nf_put_att_text(ncid,IDarea,'coordinates',9,'tlon tlat')
!
! tracer grid mask:
!
      status=nf_def_var(ncid,'tmask',NF_INT,2,dims,IDtmask)
      status=nf_put_att_text(ncid,IDtmask,'units',4,'none')
      status=nf_put_att_text(ncid,IDtmask,'long_name',27,
     &                       'horizontal tracer grid mask')
      status=nf_put_att_int(ncid,IDtmask,
     &                         'missing_value',NF_INT,1,99)
!      status=nf_put_att_text(ncid,IDtmask,'coordinates',9,'tlon tlat')
!
! momentum grid mask:
!
      dims(1)=DimIDs(2)
      dims(2)=DimIDs(4)
      status=nf_def_var(ncid,'umask',NF_INT,2,dims,IDumask)
      status=nf_put_att_text(ncid,IDumask,'units',4,'none')
      status=nf_put_att_text(ncid,IDumask,'long_name',29,
     &                       'horizontal momentum grid mask')
      status=nf_put_att_int(ncid,IDumask,
     &                         'missing_value',NF_INT,1,99)
!      status=nf_put_att_text(ncid,IDumask,'coordinates',9,'ulon ulat')
!
! bathymetry:
!
      dims(1)=DimIDs(1)
      dims(2)=DimIDs(3)
      status=nf_def_var(ncid,'h',NF_DOUBLE,2,dims,IDbathy)
      status=nf_put_att_text(ncid,IDbathy,'units',6,'meters')
      status=nf_put_att_text(ncid,IDbathy,'long_name',31,
     &                       'bathymetry on tracer gridpoints')
      status=nf_put_att_double(ncid,IDbathy,
     &                         'missing_value',NF_DOUBLE,1,spval)
!      status=nf_put_att_text(ncid,IDbathy,'coordinates',9,'tlon tlat')
!
! Coriolis parameter:
!
      dims(1)=DimIDs(2)
      dims(2)=DimIDs(4)
      status=nf_def_var(ncid,'fcor',NF_DOUBLE,2,dims,IDfcor)
      status=nf_put_att_text(ncid,IDfcor,'units',9,'1/seconds')
      status=nf_put_att_text(ncid,IDfcor,'long_name',35,
     &                       'Coriolis parameter on momentum grid')
      status=nf_put_att_double(ncid,IDfcor,
     &                         'missing_value',NF_DOUBLE,1,spval)
!      status=nf_put_att_text(ncid,IDfcor,'coordinates',9,'ulon ulat')
!
!-----------------------------------------------------------------------
! Define variables and their attributes.
!-----------------------------------------------------------------------
!
! time:
!
      status=nf_def_var(ncid,'time',NF_DOUBLE,1,recdim,IDtime)
      status=nf_put_att_text(ncid,IDtime,'calendar',7,'360_day')
      if ((cflag.eq.'m').or.(cflag.eq.'M')) then
        status=nf_put_att_text(ncid,IDtime,'units',32,
     &                         'months since 0000-00-00 00:00:00')
        status=nf_put_att_text(ncid,IDtime,'delta_t',19,
     &                            '0000-01-00 00:00:00')
        status=nf_put_att_text(ncid,IDtime,'avg_period',19,
     &                            '0000-01-00 00:00:00')
      endif
      if ((cflag.eq.'a').or.(cflag.eq.'A')) then
        status=nf_put_att_text(ncid,IDtime,'units',31,
     &                         'years since 0000-00-00 00:00:00')
        status=nf_put_att_text(ncid,IDtime,'delta_t',19,
     &                            '0001-00-00 00:00:00')
        status=nf_put_att_text(ncid,IDtime,'avg_period',19,
     &                            '0001-00-00 00:00:00')
      endif
      itimrec=0
!
! fields on 3 d tracer grid
!
      dims(1)=DimIDs(1)
      dims(2)=DimIDs(3)
      dims(3)=DimIDs(5)
      dims(4)=DimIDs(17)
!
! potential temperature:
!
      if (lnctemp) then
        status=nf_def_var(ncid,'temp',NF_DOUBLE,4,dims,IDtemp)
        status=nf_put_att_text(ncid,IDtemp,'units',7,'Celsius')
        status=nf_put_att_text(ncid,IDtemp,'long_name',30,
     &                              'averaged potential temperature')
        status=nf_put_att_double(ncid,IDtemp,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDtemp,'coordinates',21,
!     &                         'tlon tlat tdepth time')
      endif
!
! salinity:
!
      if (lncsalt) then
        status=nf_def_var(ncid,'salt',NF_DOUBLE,4,dims,IDsalt)
        status=nf_put_att_text(ncid,IDsalt,'units',3,'PSU')
        status=nf_put_att_text(ncid,IDsalt,'long_name',17,
     &                              'averaged salinity')
        status=nf_put_att_double(ncid,IDsalt,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDsalt,'coordinates',21,
!     &                         'tlon tlat tdepth time')
      endif
#if ( ISOOCN >= 1 )
!
! ocean OO17:
!
      if (lncoo17) then
        status=nf_def_var(ncid,'oo17',NF_DOUBLE,4,dims,IDoo17)
        status=nf_put_att_text(ncid,IDoo17,'units',4,'UNIT')
        status=nf_put_att_text(ncid,IDoo17,'long_name',13,
     &                              'averaged OO17')
        status=nf_put_att_double(ncid,IDoo17,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDoo17,'coordinates',21,
!     &                         'tlon tlat tdepth time')
      endif
!
! ocean OO18:
!
      if (lncoo18) then
        status=nf_def_var(ncid,'oo18',NF_DOUBLE,4,dims,IDoo18)
        status=nf_put_att_text(ncid,IDoo18,'units',4,'UNIT')
        status=nf_put_att_text(ncid,IDoo18,'long_name',13,
     &                              'averaged OO18')
        status=nf_put_att_double(ncid,IDoo18,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDoo18,'coordinates',21,
!     &                         'tlon tlat tdepth time')
      endif
!
! ocean OOHD:
!
      if (lncoohd) then
        status=nf_def_var(ncid,'oohd',NF_DOUBLE,4,dims,IDoohd)
        status=nf_put_att_text(ncid,IDoohd,'units',4,'UNIT')
        status=nf_put_att_text(ncid,IDoohd,'long_name',13,
     &                              'averaged OOHD')
        status=nf_put_att_double(ncid,IDoohd,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDoohd,'coordinates',21,
!     &                         'tlon tlat tdepth time')
      endif
#endif

#if ( OCYCC == 1 )
!
! ocean ODOC:
!
      if (lncodoc) then
        status=nf_def_var(ncid,'odoc',NF_DOUBLE,4,dims,IDodoc)
        status=nf_put_att_text(ncid,IDodoc,'units',4,'UNIT')
        status=nf_put_att_text(ncid,IDodoc,'long_name',13,
     &                              'averaged ODOC')
        status=nf_put_att_double(ncid,IDodoc,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDodoc,'coordinates',21,
!     &                         'tlon tlat tdepth time')
      endif
!
! ocean ODOCS:
!
      if (lncodocs) then
        status=nf_def_var(ncid,'odocs',NF_DOUBLE,4,dims,IDodocs)
        status=nf_put_att_text(ncid,IDodocs,'units',4,'UNIT')
        status=nf_put_att_text(ncid,IDodocs,'long_name',13,
     &                              'averaged ODOCS')
        status=nf_put_att_double(ncid,IDodocs,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDodocs,'coordinates',21,
!     &                         'tlon tlat tdepth time')
      endif
!
! ocean ODIC:
!
      if (lncodic) then
        status=nf_def_var(ncid,'odic',NF_DOUBLE,4,dims,IDodic)
        status=nf_put_att_text(ncid,IDodic,'units',4,'UNIT')
        status=nf_put_att_text(ncid,IDodic,'long_name',13,
     &                              'averaged ODIC')
        status=nf_put_att_double(ncid,IDodic,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDodic,'coordinates',21,
!     &                         'tlon tlat tdepth time')
      endif
!
! ocean OPO4:
!
      if (lncopo4) then
        status=nf_def_var(ncid,'opo4',NF_DOUBLE,4,dims,IDopo4)
        status=nf_put_att_text(ncid,IDopo4,'units',4,'UNIT')
        status=nf_put_att_text(ncid,IDopo4,'long_name',13,
     &                              'averaged OPO4')
        status=nf_put_att_double(ncid,IDopo4,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDopo4,'coordinates',21,
!     &                         'tlon tlat tdepth time')
      endif
#if ( OXNITREUX == 1 )
!
! ocean ON2O:
!
      if (lncon2o) then
        status=nf_def_var(ncid,'on2o',NF_DOUBLE,4,dims,IDon2o)
        status=nf_put_att_text(ncid,IDon2o,'units',4,'UNIT')
        status=nf_put_att_text(ncid,IDon2o,'long_name',13,
     &                              'averaged ON2O')
        status=nf_put_att_double(ncid,IDon2o,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDon2o,'coordinates',21,
!     &                         'tlon tlat tdepth time')
      endif
#else
!
! ocean ONO3:
!
      if (lncono3) then
        status=nf_def_var(ncid,'ono3',NF_DOUBLE,4,dims,IDono3)
        status=nf_put_att_text(ncid,IDono3,'units',4,'UNIT')
        status=nf_put_att_text(ncid,IDono3,'long_name',13,
     &                              'averaged ONO3')
        status=nf_put_att_double(ncid,IDono3,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDono3,'coordinates',21,
!     &                         'tlon tlat tdepth time')
      endif

#endif
!
! ocean OALK:
!
      if (lncoalk) then
        status=nf_def_var(ncid,'oalk',NF_DOUBLE,4,dims,IDoalk)
        status=nf_put_att_text(ncid,IDoalk,'units',4,'UNIT')
        status=nf_put_att_text(ncid,IDoalk,'long_name',13,
     &                              'averaged OALK')
        status=nf_put_att_double(ncid,IDoalk,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDoalk,'coordinates',21,
!     &                         'tlon tlat tdepth time')
      endif
!
! ocean OO2:
!
      if (lncoo2) then
        status=nf_def_var(ncid,'oo2',NF_DOUBLE,4,dims,IDoo2)
        status=nf_put_att_text(ncid,IDoo2,'units',4,'UNIT')
        status=nf_put_att_text(ncid,IDoo2,'long_name',13,
     &                              'averaged OO2')
        status=nf_put_att_double(ncid,IDoo2,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDoo2,'coordinates',21,
!     &                         'tlon tlat tdepth time')
      endif

!
! ocean OC13:
!
      if (lncoc13) then
        status=nf_def_var(ncid,'oc13',NF_DOUBLE,4,dims,IDoc13)
        status=nf_put_att_text(ncid,IDoc13,'units',4,'UNIT')
        status=nf_put_att_text(ncid,IDoc13,'long_name',13,
     &                              'averaged OC13')
        status=nf_put_att_double(ncid,IDoc13,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDoc13,'coordinates',21,
!     &                         'tlon tlat tdepth time')
      endif
!
! ocean ODOC13:
!
      if (lncodoc13) then
        status=nf_def_var(ncid,'odoc13',NF_DOUBLE,4,dims,IDodoc13)
        status=nf_put_att_text(ncid,IDodoc13,'units',4,'UNIT')
        status=nf_put_att_text(ncid,IDodoc13,'long_name',13,
     &                              'averaged ODOC13')
        status=nf_put_att_double(ncid,IDodoc13,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDodoc13,'coordinates',21,
!     &                         'tlon tlat tdepth time')
      endif
!
! ocean ODOCS13:
!
      if (lncodocs13) then
        status=nf_def_var(ncid,'odocs13',NF_DOUBLE,4,dims,IDodocs13)
        status=nf_put_att_text(ncid,IDodocs13,'units',4,'UNIT')
        status=nf_put_att_text(ncid,IDodocs13,'long_name',13,
     &                              'averaged ODOCS13')
        status=nf_put_att_double(ncid,IDodocs13,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDodocs13,'coordinates',21,
!     &                         'tlon tlat tdepth time')
      endif
!
! ocean OC14:
!
      if (lncoc14) then
        status=nf_def_var(ncid,'oc14',NF_DOUBLE,4,dims,IDoc14)
        status=nf_put_att_text(ncid,IDoc14,'units',4,'UNIT')
        status=nf_put_att_text(ncid,IDoc14,'long_name',13,
     &                              'averaged OC14')
        status=nf_put_att_double(ncid,IDoc14,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDoc14,'coordinates',21,
!     &                         'tlon tlat tdepth time')
      endif

#if ( PATH >= 1 )

!
! ocean PaPart:
!
      if (lncpapart) then
        status=nf_def_var(ncid,'papart',NF_DOUBLE,4,dims,IDpapart)
        status=nf_put_att_text(ncid,IDpapart,'units',4,'UNIT')
        status=nf_put_att_text(ncid,IDpapart,'long_name',13,
     &                              'averaged PaPart')
        status=nf_put_att_double(ncid,IDpapart,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDpapart,'coordinates',21,
!     &                         'tlon tlat tdepth time')
      endif
!
! ocean PaDiss:
!
      if (lncpadiss) then
        status=nf_def_var(ncid,'padiss',NF_DOUBLE,4,dims,IDpadiss)
        status=nf_put_att_text(ncid,IDpadiss,'units',4,'UNIT')
        status=nf_put_att_text(ncid,IDpadiss,'long_name',13,
     &                              'averaged PaDiss')
        status=nf_put_att_double(ncid,IDpadiss,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDpadiss,'coordinates',21,
!     &                         'tlon tlat tdepth time')
      endif
!
! ocean ThPart:
!
      if (lncthpart) then
        status=nf_def_var(ncid,'thpart',NF_DOUBLE,4,dims,IDthpart)
        status=nf_put_att_text(ncid,IDthpart,'units',4,'UNIT')
        status=nf_put_att_text(ncid,IDthpart,'long_name',13,
     &                              'averaged ThPart')
        status=nf_put_att_double(ncid,IDthpart,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDthpart,'coordinates',21,
!     &                         'tlon tlat tdepth time')
      endif
!
! ocean ThDiss:
!
      if (lncthdiss) then
        status=nf_def_var(ncid,'thdiss',NF_DOUBLE,4,dims,IDthdiss)
        status=nf_put_att_text(ncid,IDthdiss,'units',4,'UNIT')
        status=nf_put_att_text(ncid,IDthdiss,'long_name',13,
     &                              'averaged ThDiss')
        status=nf_put_att_double(ncid,IDthdiss,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDthdiss,'coordinates',21,
!     &                         'tlon tlat tdepth time')
      endif

!
! ocean PFluxCaCO3:
!
      if (lncpflxca) then
        status=nf_def_var(ncid,'pflxca',NF_DOUBLE,4,dims,IDpflxca)
        status=nf_put_att_text(ncid,IDpflxca,'units',4,'UNIT')
        status=nf_put_att_text(ncid,IDpflxca,'long_name',13,
     &                              'averaged PFlxCa')
        status=nf_put_att_double(ncid,IDpflxca,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDthdiss,'coordinates',21,
!     &                         'tlon tlat tdepth time')
      endif

!
! ocean PFluxPOC:
!
      if (lncpflxpo) then
        status=nf_def_var(ncid,'pflxpo',NF_DOUBLE,4,dims,IDpflxpo)
        status=nf_put_att_text(ncid,IDpflxpo,'units',4,'UNIT')
        status=nf_put_att_text(ncid,IDpflxpo,'long_name',13,
     &                              'averaged PFlxPOC')
        status=nf_put_att_double(ncid,IDpflxpo,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDthdiss,'coordinates',21,
!     &                         'tlon tlat tdepth time')
      endif

#endif

#if ( CORAL == 1 )
!!
!nb coral area
!
      if (lnccorarea) then
        status=nf_def_var(ncid,'corarea',NF_DOUBLE,4,dims,IDcorarea)
        status=nf_put_att_text(ncid,IDcorarea,'units',4,'UNIT')
        status=nf_put_att_text(ncid,IDcorarea,'long_name',13,
     &                              'coral area')
        status=nf_put_att_double(ncid,IDcorarea,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDodic,'coordinates',21,
!     &                         'tlon tlat tdepth time')
      endif

!nb coral production
!
      if (lnccorprod) then
        status=nf_def_var(ncid,'corprod',NF_DOUBLE,4,dims,IDcorprod)
        status=nf_put_att_text(ncid,IDcorprod,'units',4,'UNIT')
        status=nf_put_att_text(ncid,IDcorprod,'long_name',13,
     &                              'coral production')
        status=nf_put_att_double(ncid,IDcorprod,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDodic,'coordinates',21,
!     &                         'tlon tlat tdepth time')
      endif

!nb coral mass
!
      if (lnccormass) then
        status=nf_def_var(ncid,'cormass',NF_DOUBLE,4,dims,IDcormass)
        status=nf_put_att_text(ncid,IDcormass,'units',4,'UNIT')
        status=nf_put_att_text(ncid,IDcormass,'long_name',10,
     &                              'coral mass')
        status=nf_put_att_double(ncid,IDcormass,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDodic,'coordinates',21,
!     &                         'tlon tlat tdepth time')
      endif

!nb aragonite omega
!
      if (lncomega) then
        status=nf_def_var(ncid,'omega',NF_DOUBLE,4,dims,IDomega)
        status=nf_put_att_text(ncid,IDomega,'units',4,'UNIT')
        status=nf_put_att_text(ncid,IDomega,'long_name',13,
     &                              'omega arag')
        status=nf_put_att_double(ncid,IDomega,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDodic,'coordinates',21,
!     &                         'tlon tlat tdepth time')
      endif

!nb CO32- concentration
!
      if (lncoco3) then
        status=nf_def_var(ncid,'oco3',NF_DOUBLE,4,dims,IDoco3)
        status=nf_put_att_text(ncid,IDoco3,'units',4,'UNIT')
        status=nf_put_att_text(ncid,IDoco3,'long_name',13,
     &                              'OCO3')
        status=nf_put_att_double(ncid,IDoco3,
     &                           'missing_value',NF_DOUBLE,1,spval)
      endif
!nb tau bleaching
!
      if (lnctaub) then
        status=nf_def_var(ncid,'tau_bleach',NF_DOUBLE,4,dims,IDtaub)
        status=nf_put_att_text(ncid,IDtaub,'units',4,'UNIT')
        status=nf_put_att_text(ncid,IDtaub,'long_name',13,
     &                              'Tau bleaching')
        status=nf_put_att_double(ncid,IDtaub,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDodic,'coordinates',21,
!     &                         'tlon tlat tdepth time')
      endif
!nb DHW
!
      if (lncdhw) then
        status=nf_def_var(ncid,'DHW',NF_DOUBLE,4,dims,IDdhw)
        status=nf_put_att_text(ncid,IDdhw,'units',4,'UNIT')
        status=nf_put_att_text(ncid,IDdhw,'long_name',13,
     &                              'DHW')
        status=nf_put_att_double(ncid,IDdhw,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDodic,'coordinates',21,
!     &                         'tlon tlat tdepth time')
      endif
!nb PH
!
      if (lncph) then
        status=nf_def_var(ncid,'PH',NF_DOUBLE,4,dims,IDph)
        status=nf_put_att_text(ncid,IDph,'units',4,'UNIT')
        status=nf_put_att_text(ncid,IDph,'long_name',13,
     &                              'PH')
        status=nf_put_att_double(ncid,IDph,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDodic,'coordinates',21,
!     &                         'tlon tlat tdepth time')
      endif

#endif

!ocean O2 isotopes
#if ( OOISO == 1 )
      if (lncoo2_2) then
        status=nf_def_var(ncid,'oo2_16',NF_DOUBLE,4,dims,IDoo2_2)
        status=nf_put_att_text(ncid,IDoo2_2,'units',4,'UNIT')
        status=nf_put_att_text(ncid,IDoo2_2,'long_name',13,
     &                              'averaged OO2')
        status=nf_put_att_double(ncid,IDoo2_2,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDoo2_2,'coordinates',21,
!     &                         'tlon tlat tdepth time')
      endif
      if (lncoo2_3) then
        status=nf_def_var(ncid,'oo2_17',NF_DOUBLE,4,dims,IDoo2_3)
        status=nf_put_att_text(ncid,IDoo2_3,'units',4,'UNIT')
        status=nf_put_att_text(ncid,IDoo2_3,'long_name',13,
     &                              'averaged OO2')
        status=nf_put_att_double(ncid,IDoo2_3,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDoo2,'coordinates',21,
!     &                         'tlon tlat tdepth time')
      endif
      if (lncoo2_4) then
        status=nf_def_var(ncid,'oo2_18',NF_DOUBLE,4,dims,IDoo2_4)
        status=nf_put_att_text(ncid,IDoo2_4,'units',4,'UNIT')
        status=nf_put_att_text(ncid,IDoo2_4,'long_name',13,
     &                              'averaged OO2')
        status=nf_put_att_double(ncid,IDoo2_4,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDoo2_4,'coordinates',21,
!     &                         'tlon tlat tdepth time')
      endif
#endif

#endif
!
! fields on 3 d momentum grid
!
      dims(1)=DimIDs(2)
      dims(2)=DimIDs(4)
      dims(3)=DimIDs(5)
      dims(4)=DimIDs(17)
!
! u:
!
      if (lncu) then
        status=nf_def_var(ncid,'u',NF_DOUBLE,4,dims,IDu)
        status=nf_put_att_text(ncid,IDu,'units',13,'meters/second')
        status=nf_put_att_text(ncid,IDu,'long_name',33,
     &                          'averaged zonal velocity component')
        status=nf_put_att_double(ncid,IDu,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDu,'coordinates',21,
!     &                         'ulon ulat tdepth time')
      endif
!
! v:
!
      if (lncv) then
        status=nf_def_var(ncid,'v',NF_DOUBLE,4,dims,IDv)
        status=nf_put_att_text(ncid,IDv,'units',13,'meters/second')
        status=nf_put_att_text(ncid,IDv,'long_name',38,
     &                      'averaged meridional velocity component')
        status=nf_put_att_double(ncid,IDv,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDv,'coordinates',21,
!     &                          'ulon ulat tdepth time')
      endif
!
! w:
!
      dims(1)=DimIDs(1)
      dims(2)=DimIDs(3)
      dims(3)=DimIDs(6)
      if (lncw) then
        status=nf_def_var(ncid,'w',NF_DOUBLE,4,dims,IDw)
        status=nf_put_att_text(ncid,IDw,'units',10,'meters/day')
        status=nf_put_att_text(ncid,IDw,'long_name',36,
     &                      'averaged vertical velocity component')
        status=nf_put_att_double(ncid,IDw,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDw,'coordinates',21,
!     &                         'tlon tlat wdepth time')
      endif
!
! fields on horizontal momentum grid
!
        dims(1)=DimIDs(2)
        dims(2)=DimIDs(4)
        dims(3)=DimIDs(17)
!
! ubar:
!
      if (lncubar) then
        status=nf_def_var(ncid,'ubar',NF_DOUBLE,3,dims,IDubar)
        status=nf_put_att_text(ncid,IDubar,'units',13,'meters/second')
        status=nf_put_att_text(ncid,IDubar,'long_name',34,
     &                         'averaged zonal barotropic momentum')
        status=nf_put_att_double(ncid,IDubar,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDubar,'coordinates',14,
!     &                         'ulon ulat time')
      endif
!
! vbar:
!
      if (lncvbar) then
        status=nf_def_var(ncid,'vbar',NF_DOUBLE,3,dims,IDvbar)
        status=nf_put_att_text(ncid,IDvbar,'units',13,'meters/second')
        status=nf_put_att_text(ncid,IDvbar,'long_name',39,
     &                        'averaged meridional barotropic momentum')
        status=nf_put_att_double(ncid,IDvbar,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDvbar,'coordinates',14,
!     &                         'ulon ulat time')
      endif
!
! fields on horizontal tracer grid
!
        dims(1)=DimIDs(1)
        dims(2)=DimIDs(3)
        dims(3)=DimIDs(17)
!
! sea surface height:
!
      if (lncssh) then
        status=nf_def_var(ncid,'ssh',NF_DOUBLE,3,dims,IDssh)
        status=nf_put_att_text(ncid,IDssh,'units',6,'meters')
        status=nf_put_att_text(ncid,IDssh,'long_name',27,
     &                         'averaged sea surface height')
        status=nf_put_att_double(ncid,IDssh,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDssh,'coordinates',14,
!     &                         'tlon tlat time')
      endif
!
! sea surface temperature
!
      if (lncsst) then
        status=nf_def_var(ncid,'sst',NF_DOUBLE,3,dims,IDsst)
        status=nf_put_att_text(ncid,IDsst,'units',7,'Celsius')
        status=nf_put_att_text(ncid,IDsst,'long_name',12,
     &                      'averaged SST')
        status=nf_put_att_double(ncid,IDsst,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDsst,'coordinates',14,
!     &                         'tlon tlat time')
      endif
!dmr --- added the following forgotten in revision 184
!
! fice: flux of water from icebergs ...
!
      if (lncfice) then
        status=nf_def_var(ncid,'fice',NF_DOUBLE,3,dims,IDfice)
        status=nf_put_att_text(ncid,IDfice,'units',4,'m3/s')
        status=nf_put_att_text(ncid,IDfice,'long_name',19,
     &                              'averaged fice-bergs')
        status=nf_put_att_double(ncid,IDfice,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDoo17,'coordinates',21,
!     &                         'tlon tlat tdepth time')
      endif
! fdist: distribution of icebergs per gridcell...
!
      if (lncfdist) then
        status=nf_def_var(ncid,'fdist',NF_DOUBLE,3,dims,IDfdist)
        status=nf_put_att_text(ncid,IDfdist,'units',4,'nr')
        status=nf_put_att_text(ncid,IDfdist,'long_name',19,
     &                              'averaged numer bergs')
        status=nf_put_att_double(ncid,IDfdist,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDoo17,'coordinates',21,
!     &                         'tlon tlat tdepth time')
      endif

#if ( OCYCC == 1 )
!
! surface partial pressure CO2 in the ocean
!
      if (lncaqpco2) then
        status=nf_def_var(ncid,'aqpco2',NF_DOUBLE,3,dims,IDaqpco2)
        status=nf_put_att_text(ncid,IDaqpco2,'units',3,'ppm')
        status=nf_put_att_text(ncid,IDaqpco2,'long_name',34,
     &                      'CO2 partial pressure ocean surface')
        status=nf_put_att_double(ncid,IDaqpco2,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDaqpco2,'coordinates',14,
!     &                         'tlon tlat time')
      endif
!
! surface partial pressure CO2 in the ocean
!
      if (lnctpp_ma) then
        status=nf_def_var(ncid,'tpp_ma',NF_DOUBLE,3,dims,IDtpp_ma)
        status=nf_put_att_text(ncid,IDtpp_ma,'units',3,'unit')
        status=nf_put_att_text(ncid,IDtpp_ma,'long_name',34,
     &                      'total particulate prod. export')
        status=nf_put_att_double(ncid,IDtpp_ma,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDtpp_ma,'coordinates',14,
!     &                         'tlon tlat time')
      endif
!
! surface partial pressure CO2 in the ocean
!
      if (lnccaco3m) then
        status=nf_def_var(ncid,'caco3m',NF_DOUBLE,3,dims,IDcaco3m)
        status=nf_put_att_text(ncid,IDcaco3m,'units',3,'unit')
        status=nf_put_att_text(ncid,IDcaco3m,'long_name',34,
     &                      'CACO3 prod. export')
        status=nf_put_att_double(ncid,IDcaco3m,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDcaco3m,'coordinates',14,
!     &                         'tlon tlat time')
      endif
#endif

#if ( F_PALAEO_FWF == 1 || CONSEAU == 1 )
!FWF to ocean
      if (lncoutFWF) then
        write(*,*) 'write status FWF'
        status=nf_def_var(ncid,'outFWF',NF_DOUBLE,3,dims,IDoutFWF)
        status=nf_put_att_text(ncid,IDoutFWF,'units',3,'Sv')
        status=nf_put_att_text(ncid,IDoutFWF,'long_name',34,
     &                      'FWF in ocean')
        status=nf_put_att_double(ncid,IDoutFWF,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDoutFWF,'coordinates',14,
!     &                         'tlon tlat time')
      endif

#endif
!
! sea surface salinity
!
      if (lncsss) then
        status=nf_def_var(ncid,'sss',NF_DOUBLE,3,dims,IDsss)
        status=nf_put_att_text(ncid,IDsss,'units',3,'PSU')
        status=nf_put_att_text(ncid,IDsss,'long_name',29,
     &                      'averaged sea surface salinity')
        status=nf_put_att_double(ncid,IDsss,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDsss,'coordinates',14,
!     &                         'tlon tlat time')
      endif
!
! surface heat flux:
!
      if (lncshflx) then
        status=nf_def_var(ncid,'shflx',NF_DOUBLE,3,dims,IDshflx)
        status=nf_put_att_text(ncid,IDshflx,'units',6,'W/m**2')
        status=nf_put_att_text(ncid,IDshflx,'long_name',26,
     &                         'averaged surface heat flux')
        status=nf_put_att_double(ncid,IDshflx,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDshflx,'coordinates',14,
!     &                         'tlon tlat time')
      endif
!
! surface freshwater flux:
!
      if (lncsfflx) then
        status=nf_def_var(ncid,'sfflx',NF_DOUBLE,3,dims,IDsfflx)
        status=nf_put_att_text(ncid,IDsfflx,'units',11,'meters/year')
        status=nf_put_att_text(ncid,IDsfflx,'long_name',32,
     &                         'averaged surface freshwater flux')
        status=nf_put_att_double(ncid,IDsfflx,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDsfflx,'coordinates',14,
!     &                         'tlon tlat time')
      endif
!
! mixed layer depth:
!
      if (lnczmix) then
        status=nf_def_var(ncid,'zmix',NF_DOUBLE,3,dims,IDzmix)
        status=nf_put_att_text(ncid,IDzmix,'units',6,'meters')
        status=nf_put_att_text(ncid,IDzmix,'long_name',43,
     &                 'averaged depth of ocean surface mixed layer')
        status=nf_put_att_double(ncid,IDzmix,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDzmix,'coordinates',14,
!     &                         'tlon tlat time')
      endif
!
! maximum convection depth:
!
      if (lnczcnv) then
        status=nf_def_var(ncid,'zcnv',NF_DOUBLE,3,dims,IDzcnv)
        status=nf_put_att_text(ncid,IDzcnv,'units',6,'meters')
        status=nf_put_att_text(ncid,IDzcnv,'long_name',28,
     &                 'averaged depth of convection')
        status=nf_put_att_double(ncid,IDzcnv,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDzcnv,'coordinates',14,
!     &                         'tlon tlat time')
      endif
!
! G-M slope parameter
!
      if (lncmsl) then
        status=nf_def_var(ncid,'msl',NF_DOUBLE,3,dims,IDmsl)
        status=nf_put_att_text(ncid,IDmsl,'units',4,'none')
        status=nf_put_att_text(ncid,IDmsl,'long_name',18,
     &                      'averaged G-M slope')
        status=nf_put_att_double(ncid,IDmsl,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDmsl,'coordinates',14,
!     &                         'tlon tlat time')
      endif
!
! ice model data:
!
      if (lnchic) then
        status=nf_def_var(ncid,'hice',NF_DOUBLE,3,dims,IDhic)
        status=nf_put_att_text(ncid,IDhic,'units',6,'meters')
        status=nf_put_att_text(ncid,IDhic,'long_name',22,
     &                      'averaged ice thickness')
        status=nf_put_att_double(ncid,IDhic,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDhic,'coordinates',14,
!     &                         'tlon tlat time')
      endif
!
      if (lnchicp) then
        status=nf_def_var(ncid,'hicp',NF_DOUBLE,3,dims,IDhicp)
        status=nf_put_att_text(ncid,IDhicp,'units',6,'meters')
        status=nf_put_att_text(ncid,IDhicp,'long_name',23,
     &                      'averaged ice production')
        status=nf_put_att_double(ncid,IDhicp,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDhicp,'coordinates',14,
!     &                         'tlon tlat time')
      endif
!
      if (lncalbq) then
        status=nf_def_var(ncid,'albq',NF_DOUBLE,3,dims,IDalbq)
        status=nf_put_att_text(ncid,IDalbq,'units',4,'none')
        status=nf_put_att_text(ncid,IDalbq,'long_name',22,
     &                      'averaged lead fraction')
        status=nf_put_att_double(ncid,IDalbq,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDalbq,'coordinates',14,
!     &                         'tlon tlat time')
      endif
!
      if (lnchsn) then
        status=nf_def_var(ncid,'hsn',NF_DOUBLE,3,dims,IDhsn)
        status=nf_put_att_text(ncid,IDhsn,'units',6,'meters')
        status=nf_put_att_text(ncid,IDhsn,'long_name',23,
     &                      'averaged snow thickness')
        status=nf_put_att_double(ncid,IDhsn,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDhsn,'coordinates',14,
!     &                         'tlon tlat time')
      endif
!
      if (lncsnow) then
        status=nf_def_var(ncid,'snow',NF_DOUBLE,3,dims,IDsnow)
        status=nf_put_att_text(ncid,IDsnow,'units',6,'meters')
        status=nf_put_att_text(ncid,IDsnow,'long_name',27,
     &                      'averaged snow precipitation')
        status=nf_put_att_double(ncid,IDsnow,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDsnow,'coordinates',14,
!     &                         'tlon tlat time')
      endif
!
      if (lnctice) then
        status=nf_def_var(ncid,'tice',NF_DOUBLE,3,dims,IDtice)
        status=nf_put_att_text(ncid,IDtice,'units',7,'Celsius')
        status=nf_put_att_text(ncid,IDtice,'long_name',24,
     &                      'averaged ice temperature')
        status=nf_put_att_double(ncid,IDtice,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDtice,'coordinates',14,
!     &                         'tlon tlat time')
      endif
!
      if (lncfb) then
        status=nf_def_var(ncid,'fb',NF_DOUBLE,3,dims,IDfb)
        status=nf_put_att_text(ncid,IDfb,'units',6,'W/m**2')
        status=nf_put_att_text(ncid,IDfb,'long_name',30,
     &                      'averaged heat flux at ice base')
        status=nf_put_att_double(ncid,IDfb,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDfb,'coordinates',14,
!     &                         'tlon tlat time')
      endif
!
      dims(1)=DimIDs(2)
      dims(2)=DimIDs(4)
      if (lncuice) then
        status=nf_def_var(ncid,'uice',NF_DOUBLE,3,dims,IDuice)
        status=nf_put_att_text(ncid,IDuice,'units',13,'meters/second')
        status=nf_put_att_text(ncid,IDuice,'long_name',27,
     &                      'averaged zonal ice velocity')
        status=nf_put_att_double(ncid,IDuice,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDuice,'coordinates',14,
!     &                         'ulon ulat time')
      endif
!
      if (lncvice) then
        status=nf_def_var(ncid,'vice',NF_DOUBLE,3,dims,IDvice)
        status=nf_put_att_text(ncid,IDvice,'units',13,'meters/second')
        status=nf_put_att_text(ncid,IDvice,'long_name',32,
     &                      'averaged meridional ice velocity')
        status=nf_put_att_double(ncid,IDvice,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDvice,'coordinates',14,
!     &                         'ulon ulat time')
      endif
!
      if (lnctx) then
        status=nf_def_var(ncid,'wsx',NF_DOUBLE,3,dims,IDtx)
        status=nf_put_att_text(ncid,IDtx,'units',6,'N/m**2')
        status=nf_put_att_text(ncid,IDtx,'long_name',26,
     &                      'averaged zonal wind stress')
        status=nf_put_att_double(ncid,IDtx,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDtx,'coordinates',14,
!     &                         'ulon ulat time')
      endif
!
      if (lncty) then
        status=nf_def_var(ncid,'wsy',NF_DOUBLE,3,dims,IDty)
        status=nf_put_att_text(ncid,IDty,'units',6,'N/m**2')
        status=nf_put_att_text(ncid,IDty,'long_name',31,
     &                      'averaged meridional wind stress')
        status=nf_put_att_double(ncid,IDty,
     &                           'missing_value',NF_DOUBLE,1,spval)
!        status=nf_put_att_text(ncid,IDty,'coordinates',14,
!     &                         'ulon ulat time')
      endif
!
! meridional streamFUNCTIONs
!
      dims(1)=DimIDs(12)
      dims(2)=DimIDs(9)
      dims(3)=DimIDs(10)
      dims(4)=DimIDs(17)
!
      if (lncmoc) then
        status=nf_def_var(ncid,'moc',NF_DOUBLE,4,dims,IDmoc)
        status=nf_put_att_text(ncid,IDmoc,'units',2,'Sv')
        status=nf_put_att_text(ncid,IDmoc,'long_name',37,
     &                      'meridional overturning streamFUNCTION')
        status=nf_put_att_double(ncid,IDmoc,
     &                           'missing_value',NF_DOUBLE,1,spval)
        status=nf_put_att_text(ncid,IDmoc,'comment',60,
     & 'i-index: 1=global 2=Indian 3=Pacific 4=Atlantic/Arctic ocean')
      endif
!
      dims(1)=DimIDs(12)
      dims(2)=DimIDs(9)
      dims(3)=DimIDs(17)
      if (lncmht) then
        status=nf_def_var(ncid,'mht',NF_DOUBLE,3,dims,IDmht)
        status=nf_put_att_text(ncid,IDmht,'units',2,'PW')
        status=nf_put_att_text(ncid,IDmht,'long_name',25,
     &                      'meridional heat transport')
        status=nf_put_att_double(ncid,IDmht,
     &                           'missing_value',NF_DOUBLE,1,spval)
        status=nf_put_att_text(ncid,IDmht,'comment',60,
     & 'i-index: 1=global 2=Indian 3=Pacific 4=Atlantic/Arctic ocean')
      endif
!
      if (lncmst) then
        status=nf_def_var(ncid,'mst',NF_DOUBLE,3,dims,IDmst)
        status=nf_put_att_text(ncid,IDmst,'units',6,'Sv*PSU')
        status=nf_put_att_text(ncid,IDmst,'long_name',25,
     &                'meridional salt transport')
        status=nf_put_att_double(ncid,IDmst,
     &                           'missing_value',NF_DOUBLE,1,spval)
        status=nf_put_att_text(ncid,IDmst,'comment',60,
     & 'i-index: 1=global 2=Indian 3=Pacific 4=Atlantic/Arctic ocean')
      endif
!
!-----------------------------------------------------------------------
! Leave definition mode.
!-----------------------------------------------------------------------
!
        status=nf_enddef(ncid)
!
!-----------------------------------------------------------------------
! Write out grid, time-recordless and information variables.
!-----------------------------------------------------------------------
!
! Standard netcdf applications like ferret require a strictly
! monotonic axis. We need to specify pseudo latitudes here
! because of the AA grid nested within the entire domain partly
! using grid points that are mapped to latidues higher than
! 90 degrees north on the WW grid. Pseudo longitudes are
! defined here to assure that the axis is strictly monotonic,
! which is not neccessarily guarantedd when the Greenwich meridian
! intersects the model domain and all latitudes are modulo 360.
!
      rmisc(1)=xslon(2,1)
      rdelta=xslon(3,1)-xslon(2,1)
      do i=2,idim
        rmisc(i)=rmisc(i-1)+rdelta
      enddo
      status=nf_put_var_double(ncid,IDptlon,rmisc)
!
      rmisc(1)=yslat(1,1)
      rdelta=yslat(1,2)-yslat(1,1)
      do j=2,jdim
        rmisc(j)=rmisc(j-1)+rdelta
      enddo
      status=nf_put_var_double(ncid,IDptlat,rmisc)
!
      rmisc(1)=xulon(2,1)
      rdelta=xulon(3,1)-xulon(2,1)
      do i=2,idim
        rmisc(i)=rmisc(i-1)+rdelta
      enddo
      status=nf_put_var_double(ncid,IDpulon,rmisc)
!
      rmisc(1)=yulat(1,1)
      rdelta=yulat(1,2)-yulat(1,1)
      do j=2,jdim
        rmisc(j)=rmisc(j-1)+rdelta
      enddo
      status=nf_put_var_double(ncid,IDpulat,rmisc)
!
! pseudo momentum latitudes up to the North Pole
! are used for meridional streamFUNCTION calculations
!
      status=nf_put_var_double(ncid,IDsflat,rmisc)
!
      status=nf_put_var_double(ncid,IDtdepth,z)
      status=nf_put_var_double(ncid,IDwdepth,zw)
      rmisc(1)=zw(1)
      do k=2,kdimp1
        rmisc(k)=z(k-1)
      enddo
      rmisc(kdimp2)=zw(kdimp1)
      status=nf_put_var_double(ncid,IDwedges,rmisc)
!
      rmisc(1)=zw(1)
      do k=1,kmax
        rmisc(k+1)=z(k)
      enddo
      rmisc(kmax+2)=0.0
      status=nf_put_var_double(ncid,IDsfdepth,rmisc)
!
      rmisc(1)=zw(1)
      rmisc(2)=zw(1)+1.0d-6
      do k=2,kmax
        rmisc(k+1)=zw(k)
      enddo
      rmisc(kmax+2)=-1.0d-6
      rmisc(kmax+3)=0.0
      status=nf_put_var_double(ncid,IDsfedges,rmisc)
!
      do k=1,4
        status=nf_put_var1_int(ncid,IDbasidx,k,k)
      enddo
!
!      start(1)=1
!      start(2)=1
!      count(1)=14
!      count(2)=4
!      status=nf_put_vara_text(ncid,IDbasnam,start,count,cbasins)
!
      start(1)=1
      count(1)=imax-2
      count(2)=1
      start(3)=0
      count(3)=0
      do j=1,jmax
        start(2)=j
        status=nf_put_vara_double(ncid,IDtlon,start,count,xslon(2,j))
        status=nf_put_vara_double(ncid,IDtlat,start,count,yslat(2,j))
        status=nf_put_vara_double(ncid,IDulon,start,count,xulon(2,j))
        status=nf_put_vara_double(ncid,IDulat,start,count,yulat(2,j))
        status=nf_put_vara_double(ncid,IDangle,start,count,angle(2,j))
      enddo
!
      start(1)=1
      count(1)=imax-1
      count(2)=1
      start(3)=0
      count(3)=0
      do j=1,jmax+1
        start(2)=j
        status=nf_put_vara_double(ncid,IDtlonp,start,count,xslonp(2,j))
        status=nf_put_vara_double(ncid,IDtlatp,start,count,yslatp(2,j))
        status=nf_put_vara_double(ncid,IDulonp,start,count,xulonp(2,j))
        status=nf_put_vara_double(ncid,IDulatp,start,count,yulatp(2,j))
      enddo
!
      do n=1,4
        start(1)=1
        count(1)=imax-2
        start(2)=1
        count(2)=1
        start(3)=n
        count(3)=1
#if ( HRCLIO == 1 )
      open (61,file='xsedg.dat')
      do j=1,jmax
!L15     write(61,'(122( I1))' ) (msks(i,j),i=1,imax)
         write(61,'(242( F10.5))' ) (ysedg(i,j,4),i=1,imax)
      enddo
      write(61,*)
      close (61)
#endif

        do j=1,jmax
          start(2)=j
          status=nf_put_vara_double(ncid,IDtlon_e,start,count,
     &                              xsedg(2,j,n))
          status=nf_put_vara_double(ncid,IDtlat_e,start,count,
     &                              ysedg(2,j,n))
          status=nf_put_vara_double(ncid,IDulon_e,start,count,
     &                              xuedg(2,j,n))
          status=nf_put_vara_double(ncid,IDulat_e,start,count,
     &                              yuedg(2,j,n))
        enddo
      enddo
      start(3)=0
      count(3)=0
!
      status=write_2d(ncid,IDbathy,hs    ,spval,.true. ,.false.,0)
!
      do j=1,jmax
        do i=1,imax
          fs2cor(i,j)=2.0*fs2cor(i,j)
        enddo
      enddo
      status=write_2d(ncid,IDfcor ,fs2cor,spval,.true. ,.false.,0)
      do j=1,jmax
        do i=1,imax
          fs2cor(i,j)=0.5*fs2cor(i,j)
        enddo
      enddo
!
      status=write_2d(ncid,IDdxs1 ,dxs1  ,spval,.true. ,.false.,0)
      status=write_2d(ncid,IDdxs2 ,dxs2  ,spval,.true. ,.false.,0)
      status=write_2d(ncid,IDdxc1 ,dxc1  ,spval,.true. ,.false.,0)
      status=write_2d(ncid,IDdxc2 ,dxc2  ,spval,.true. ,.false.,0)
      status=write_2d(ncid,IDarea ,area  ,spval,.true. ,.false.,0)
      start(1)=1
      count(1)=imax-2
      count(2)=1
      do j=1,jmax
        start(2)=j
        status=nf_put_vara_int(ncid,IDtmask,start,count,msks(2,j))
        status=nf_put_vara_int(ncid,IDumask,start,count,msku(2,j))
      enddo
#if ( HRCLIO == 1 )
      open (61,file='mask2.dat')
      do j=1,jmax
!L15     write(61,'(122( I1))' ) (msks(i,j),i=1,imax)
         write(61,'(242( I1))' ) (msks(i,j),i=1,imax)
      enddo
      write(61,*)
      close(61)
#endif
!
      status=nf_sync(ncid)
!
      return
      end
!
      integer FUNCTION write_2d(ncid,ID,field,spval,msk_s,msk_u,itim)

!! START_OF_USE_SECTION

      use const_mod

      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod
      use reper_mod

      use newunit_clio_mod, only: clio3_out_id
!! END_OF_USE_SECTION

!! START_OF_INCLUDE_SECTION

#include "type.com"
! [SCRPTCOM] #include "para.com"
! [SCRPTCOM] #include "bloc.com"
! [SCRPTCOM] #include "reper.com"

#include "netcdf.inc"

!! END_OF_INCLUDE_SECTION

      logical msk_s,msk_u
      integer ncid,ID,itim,start(3),count(3),status
      dimension f2d(imax,jmax),field(imax,jmax)
!
! fill work array
!
      do j=1,jmax
        do i=1,imax
          f2d(i,j)=field(i,j)
        enddo
      enddo
!
! Perform combined WW/AA-grid masking. Indices and lon/lat ranges for
! the two distinct subgrid have been taken from netdatL.f. We undefine
! gridpoints (by setting the spval flag) that separate the two partial
! grids, but we do apply the tracer or momentum grid mask here only
! upon request.
!
      do j=1,jmax
        xxAA = xaj1 + dxaj * DFLOAT(j-1)
        yyWW = ywj1 + dywj * DFLOAT(j-1)
        do i=1,imax
          yyAA = yai1 + dyai * DFLOAT(i-1)
          inAA = 1
          if ( yyAA.le.-47. .or. yyAA.ge.68. ) inAA = 0
          if ( xxAA.lt.(2.*yyAA-12.) ) inAA = 0
          if ( xxAA.gt.(284.-2*yyAA) ) inAA = 0
          if ( xxAA.gt.190. .and. xxAA.gt.(220.-yyAA) ) inAA = 0
          if ( xxAA.gt.(2.*yyAA+237.) ) inAA = 0
          xxWW = xwi1 + dxwi * DFLOAT(i-1)
          xxWW = xxWW - 20. + untour
          xxWW = mod(xxWW,untour) + 20.
          inWW = 1
          if ( yyWW.ge.70. ) inWW = 0
          if ( (yyWW.ge.39.) .and. (xxWW.lt.42.) ) inWW = 0
          if ( yyWW.ge.30. .and. yyWW.lt.39. .and. xxWW.lt.36. ) inWW=0
          if ( (yyWW.ge.20.) .and. (xxWW.gt.263.) ) inWW = 0
          if ( (yyWW.ge.dyai) .and. (yyWW.lt.20.) .and.
     &            xxWW.gt.(305.-2.*yyWW) ) inWW = 0
!
          if ( yyWW.ge.30. .and. yyWW.le.70. .and. xxWW.lt.264.
     &              .and. xxWW.gt.(219.+(70.-yyWW)*45./40.) ) inWW=0
!
          if ((inWW.eq.0).and.(inAA.eq.0)) f2d(i,j) = spval
!
! upon request: tracer-/momentum-grid masking
!
          if (msk_s.and.(msks(i,j).eq.0)) f2d(i,j) = spval
          if (msk_u.and.(msku(i,j).eq.0)) f2d(i,j) = spval
        enddo
      enddo
!
! dump to netcdf file
!
      start(1)=1
      count(1)=imax-2
      start(2)=1
      count(2)=1
      start(3)=itim
      count(3)=0
      if (itim.gt.0) count(3)=1
      do j=1,jmax
        start(2)=j

        write_2d=nf_put_vara_double(ncid,ID,start,count,f2d(2,j))
        if (write_2d.ne.nf_noerr) then
          write(clio3_out_id,*) '  ==== '
          write(clio3_out_id,*) '//2D ',write_2d,' writing netcdf file!'
          write(clio3_out_id,*) ' start:',start
          write(clio3_out_id,*) ' count:',count
          write(clio3_out_id,*) '  ncid:', ncid
          write(clio3_out_id,*) '    ID:', ID
          write(clio3_out_id,*) '  ==== '
          status=nf_close(ncid)
          stop
        endif
      enddo
!
      return
      end

      integer FUNCTION write_3d(ncid,ID,f3d,kdim,itim)
!
! no masking needs to be applied here, this is all done when
! averaging data!
!

!! START_OF_USE_SECTION

      use const_mod

      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod
      use reper_mod

      use newunit_clio_mod, only: clio3_out_id
!! END_OF_USE_SECTION

!! START_OF_INCLUDE_SECTION

#include "type.com"
! [SCRPTCOM] #include "para.com"
! [SCRPTCOM] #include "bloc.com"
! [SCRPTCOM] #include "reper.com"

#include "netcdf.inc"

!! END_OF_INCLUDE_SECTION

      logical msk_s,msk_u
      integer ncid,ID,itim,start(4),count(4),status
      dimension f3d(imax,jmax,kdim)
!
      do k=1,kdim
        start(1)=1
        count(1)=imax-2
        count(2)=1
        start(3)=k
        count(3)=1
        start(4)=itim
        count(4)=1
        do j=1,jmax
          start(2)=j
          write_3d=nf_put_vara_double(ncid,ID,start,count,f3d(2,j,k))
          if (write_3d.ne.nf_noerr) then
            status=nf_close(ncid)
            write(clio3_out_id,*) 'ERROR3D ',write_3d,
     &                            ' writing netcdf file!'
            write(clio3_out_id,*) ' start:',start
            write(clio3_out_id,*) ' count:',count
            write(clio3_out_id,*) ' En plus : ', ID, ncid, itim
            stop
          endif
        enddo
      enddo
!
      return
      end
!DFG end
!dmr --- Start a new write_2d function
      integer FUNCTION write_2dnm(ncid,ID,f2d,itim)
!
! no masking needs to be applied here, this is all done when
! averaging data!
!

!! START_OF_USE_SECTION

      use const_mod

      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod
      use reper_mod

      use newunit_clio_mod, only: clio3_out_id
!! END_OF_USE_SECTION

!! START_OF_INCLUDE_SECTION

#include "type.com"
! [SCRPTCOM] #include "para.com"
! [SCRPTCOM] #include "bloc.com"
! [SCRPTCOM] #include "reper.com"

#include "netcdf.inc"

!! END_OF_INCLUDE_SECTION

      logical msk_s,msk_u
      integer ncid,ID,itim,start(4),count(4),status
      dimension f2d(imax,jmax)
!
!      do k=1,kdim
        start(1)=1
        count(1)=imax-2
        count(2)=1
!        start(3)=k
!        count(3)=1
        start(3)=itim
        count(3)=1
        do j=1,jmax
          start(2)=j
          write_2dnm=nf_put_vara_double(ncid,ID,start,count,f2d(2,j))
          if (write_2dnm.ne.nf_noerr) then
            status=nf_close(ncid)
            write(clio3_out_id,*) 'ERROR2Dnm ',write_2dnm,
     &                            ' writing netcdf file!'
            write(clio3_out_id,*) ' start:',start
            write(clio3_out_id,*) ' count:',count
            write(clio3_out_id,*) ' En plus : ', ID, ncid, itim
            stop
          endif
        enddo
!      enddo
!
      return
      end function write_2dnm
!DFG end



