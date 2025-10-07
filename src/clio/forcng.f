!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:46 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:46 CET 2009

      SUBROUTINE forcng(ja,xjour)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--This routine reads the and atmospheric forcing (NCAR)
!  needed by the model. It computes surface fluxes and inter-
!  polates every quantity in time (daily) by means of a cubic spline.
!  Derniere modification 30/05/95
!  difference entre version glace et version CLIO:
!  include 'oceacn.com',tcn-scal
!-----
!  modif : 04/10/99

!! START_OF_USE_SECTION


      use global_constants_mod, only: dblp=>dp, ip      

!dmr [DEPRECATED]      use comemic_mod, only: flgicb

#if ( ISOOCN >= 1 )
      USE iso_param_mod, ONLY : ieau
#endif

      use const_mod

      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod
      use ice_mod
      use dynami_mod
!cfc  use trace_mod

!      use mchd99_mod, only: nn99
      use newunit_clio_mod, only: clio3_out_id      

!! END_OF_USE_SECTION



!cfc  include 'trace.com'

!- nn99=2 => ecritures auxiliaires sur fichier "mouchard", unit=99
c~       common / mchd99 / nn99

      integer(kind=ip) :: nsamt
      parameter(nsamt=12)

! iwind.eq.1, use already computed geostrophic winds.
! iwind.ne.1, compute wind stress from derived geostrophic winds.
! uwind     Wind velocity in the x direction
! vwind     Wind velocity in the y direction
!
        integer(kind=ip) :: iwind
        parameter(iwind=1)

        real(kind=4), save :: rtt(imax,jmax,2,2),rtr(imax,jmax,2,2),
     &         rpr(imax,jmax,2,2),rdw(imax,jmax,2,2,2),
     &         rrn(imax,jmax,2,2),rcd(imax,jmax,2,2),
     &         rhr(imax,jmax,2,2,2),rvv(imax,jmax,2,2),
     &         rru(imax,jmax,2,2),rdv(imax,jmax,2,2)



        real(kind=dblp), dimension(nsamt,2) :: rhcg2, theta
        real(kind=dblp), dimension(imax,jmax) :: psl, ugw, vgw         

        real(kind=dblp), dimension(19) :: budyko
!        real(kind=dblp), dimension(19) :: bunker, budyko

        real(kind=dblp), dimension(imax,jmax),save :: trans, zmodfo
        real(kind=dblp), dimension(imax,jmax) :: tenhrx, tenhry, vabqec
        

!-uwind et vwind en common pour les icebergs
!       dimension uwind(imax,jmax),vwind(imax,jmax)
!


        real(kind=dblp), dimension(imax,jmax), save :: amtt,amtr,ampr,
     &                 amwx,amwy,amrn,
     &                 amcd,amhrx,amhry
     &                ,amvv,amru,amdv
     
        integer(kind=ip), save :: jnddt, mit,mft
!        real(kind=dblp),  save :: cmarsh
!
        data rhcg2  /0.00085,0.00085,0.000850,0.000933,0.001016,0.00110,
     &               0.00110,0.00110,0.001016,0.000933,0.000850,0.00085,
     &               0.00110,0.00110,0.001016,0.000933,0.000850,0.00085,
     &               0.00085,0.00085,0.000850,0.000933,0.001016,0.00110/
        data theta  /33.0,33.0,33.0,29.6,26.3,23.0,
     &               23.0,23.0,26.3,29.6,33.0,33.0,
     &               23.0,23.0,26.3,29.6,33.0,33.0,
     &               33.0,33.0,33.0,29.6,26.3,23.0/
!
!  MARSHUNOVA s cp coefficient in the central Arctic:
!               source: Shine and Crane, 1984.
!
!       data cmartc /0.310,0.285,0.260,0.235,0.210,0.185,
!    &                     0.160,0.185,0.210,0.235,0.260,0.285,
!    &                     0.160,0.185,0.210,0.235,0.260,0.285,
!    &                     0.310,0.285,0.260,0.235,0.210,0.185/
!
!               source: Doronin, 1969.
!
!       data cmartc /0.300,0.300,0.300,0.280,0.270,0.240,
!    &                     0.220,0.230,0.270,0.290,0.300,0.300,
!    &                     0.220,0.230,0.270,0.290,0.300,0.300,
!    &                     0.300,0.300,0.300,0.280,0.270,0.240/
!
!  BUNKER s coefficient (cloudiness effect on LW radiation):
!
!       data bunker /1.00,1.00,0.95,0.90,0.86,0.81,0.75,0.70,0.62,0.60,
!    &                    0.62,0.70,0.75,0.81,0.86,0.90,0.95,1.00,1.00/
!
!  BUDYKO s coefficient (cloudiness effect on LW radiation):
!
!       data budyko /0.82,0.82,0.82,0.82,0.80,0.78,0.76,0.74,0.72,0.70,
!    &               0.68,0.65,0.63,0.61,0.59,0.57,0.55,0.52,0.50/
!       data budyko /1.00,0.98,0.95,0.92,0.89,0.86,0.83,0.80,0.78,0.75,
!    &               0.72,0.69,0.67,0.64,0.61,0.58,0.56,0.53,0.50/
        data budyko /1.00,0.95,0.89,0.83,0.78,0.72,0.67,0.61,0.56,0.50,
     &               0.57,0.62,0.68,0.75,0.90,1.00,1.00,1.00,1.00/




!--- loads of locales


       integer(kind=ip) :: ja,  jndt, jm, jmm1, jmp1, jjt, j, i, l
     >                   , nzmdif, n, i1, i2, j1, j2, k, ihm,  indx, jm1
     >                   , im1, njh, jp1, ip1

       real(kind=dblp) :: xjour, dpi, dtt, dtts6, datet, ttbt, ttat
     >                  , zmopr, xlati, hem, alat, clat, es, evg, evo
     >                  , rhlim, rhoa, rhoafn, usp, gam, zmod, angt
     >                  , costt, sintt, costb, sintb

       integer(kind=ip) :: wind_int_dat_id, wind_int_d_dat_id
     >                  , pres_int_dat_id, pres_int_d_dat_id
     >                  , prcp_int_dat_id, prcp_int_d_dat_id
     >                  , cldsa_int_dat_id, cldsa_int_d_dat_id
     >                  , ten5_int_dat_id, ten5_int_d_dat_id
     >                  , vvent_int_dat_id, vvent_int_d_dat_id
     >                  , temp_int_dat_id, temp_int_d_dat_id
     >                  , humid_int_dat_id, humid_int_d_dat_id
     >                  , runo_int_dat_id, runo_int_d_dat_id
     >                  , ust_int_dat_id, ust_int_d_dat_id
     >                  , modif_pre_id
!
!  Begin.
!
      dpi   = 2.0*pi
      dtt=yeaday/float(nsamt)
      dtts6=dtt/6.0
!     jour  = int(xjour)
!     datet = (dble((ja-1)*yeaday+jour)-0.5)/dtt
      datet = ((ja-1)*yeaday+xjour-1.0)/dtt
      ttbt  = datet-int(datet)
      ttat  = 1.0-ttbt
      jndt  = int(datet)
      jm    = mod(jndt,nsamt)+1
      jmm1  = jm-1+nsamt*(1/jm)
      jmp1  = jm+1-nsamt*(jm/nsamt)
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!       READ NCAR DATA.                                                |
!-----------------------------------------------------------------------
!
      if(numit.eq.nstart) then
#if ( ICEBERG_MODEL > 0 )      
!dmr [DEPRECATED]        if (flgicb) then
          open(newunit=wind_int_dat_id,file='wind.int.dat',
     &        form='unformatted')
          open(newunit=wind_int_d_dat_id,file='wind.int.d.dat',
     &        form='unformatted')
!dmr [DEPRECATED]        endif
#endif
!
!       open(71,file='/u14/grpastr/hgs/data/ncar/ncar.int.dat',
!    &        form='unformatted')
!       open(72,file='/u14/grpastr/hgs/data/ncar/ncar.int.d.dat',
!    &        form='unformatted')
        open(newunit=pres_int_dat_id,file='pres.int.dat',
     &  form='unformatted')
        open(newunit=pres_int_d_dat_id,file='pres.int.d.dat',
     &  form='unformatted')
        open(newunit=prcp_int_dat_id,file='prcp.int.dat',
     &  form='unformatted')
        open(newunit=prcp_int_d_dat_id,file='prcp.int.d.dat',
     &  form='unformatted')
!       open(73,file='/u14/grpastr/hgs/data/precip/pleg.int.dat',
!    &  form='unformatted')
!       open(74,file='/u14/grpastr/hgs/data/precip/pleg.int.d.dat',
!    &  form='unformatted')
        open(newunit=cldsa_int_dat_id,file='cldsa.int.dat',
     &  form='unformatted')
        open(newunit=cldsa_int_d_dat_id,file='cldsa.int.d.dat',
     &  form='unformatted')
!       open(77,file='tenhr.int.dat',
!    &  form='unformatted')
!       open(78,file='tenhr.int.d.dat',
!    &  form='unformatted')
        open(newunit=ten5_int_dat_id,file='ten5.int.dat',
     &  form='unformatted')
        open(newunit=ten5_int_d_dat_id,file='ten5.int.d.dat',
     &  form='unformatted')
!       open(79,file='runoff.anu2.dat',
!    &  form='unformatted')
        open(newunit=vvent_int_dat_id,file='vvent.int.dat',
     &  form='unformatted')
        open(newunit=vvent_int_d_dat_id,file='vvent.int.d.dat',
     &  form='unformatted')
        open(newunit=temp_int_dat_id,file='temp.int.dat',
     &  form='unformatted')
        open(newunit=temp_int_d_dat_id,file='temp.int.d.dat',
     &  form='unformatted')
!       open(84,file='/u14/grpastr/hgs/data/ncar/tdew.int.dat',
!    &  form='unformatted')
!       open(85,file='/u14/grpastr/hgs/data/ncar/tdew.int.d.dat',
!    &  form='unformatted')
        open(newunit=humid_int_dat_id,file='humid.int.dat',
     &  form='unformatted')
        open(newunit=humid_int_d_dat_id,file='humid.int.d.dat',
     &  form='unformatted')
        open(newunit=runo_int_dat_id,file='runo.int.dat',
     &  form='unformatted')
        open(newunit=runo_int_d_dat_id,file='runo.int.d.dat',
     &  form='unformatted')
        open(newunit=ust_int_dat_id,file='ust.int.dat',
     &  form='unformatted')
        open(newunit=ust_int_d_dat_id,file='ust.int.d.dat',
     &  form='unformatted')

        rewind pres_int_dat_id
        rewind pres_int_d_dat_id
        rewind prcp_int_dat_id
        rewind prcp_int_d_dat_id
        rewind cldsa_int_dat_id
        rewind cldsa_int_d_dat_id
        rewind ten5_int_dat_id
        rewind ten5_int_d_dat_id
        rewind vvent_int_dat_id
        rewind vvent_int_d_dat_id
        rewind temp_int_dat_id
        rewind temp_int_d_dat_id
        rewind humid_int_dat_id
        rewind humid_int_d_dat_id
        rewind runo_int_dat_id
        rewind runo_int_d_dat_id
        rewind ust_int_dat_id
        rewind ust_int_d_dat_id
        rewind wind_int_dat_id
        rewind wind_int_d_dat_id
!
        jjt=mod(jndt,nsamt)
        do j=1,jmax
        read(pres_int_dat_id) (rpr(i,j,1,1),i=1,imax)
!       read(pres_int_dat_id) (rtt(i,j,1,1),i=1,imax)
!       read(pres_int_dat_id) (rtr(i,j,1,1),i=1,imax)
!       read(pres_int_dat_id) (rdw(i,j,1,1,1),i=1,imax)
!       read(pres_int_dat_id) (rdw(i,j,2,1,1),i=1,imax)
        read(pres_int_d_dat_id) (rpr(i,j,2,1),i=1,imax)
!       read(pres_int_d_dat_id) (rtt(i,j,2,1),i=1,imax)
!       read(pres_int_d_dat_id) (rtr(i,j,2,1),i=1,imax)
!       read(pres_int_d_dat_id) (rdw(i,j,1,2,1),i=1,imax)
!       read(pres_int_d_dat_id) (rdw(i,j,2,2,1),i=1,imax)
#if( ICEBERG_MODEL > 0 )
!dmr [DEPRECATED]        if (flgicb) then
          read(wind_int_dat_id) (rdw(i,j,1,1,1),i=1,imax)
          read(wind_int_dat_id) (rdw(i,j,2,1,1),i=1,imax)
          read(wind_int_d_dat_id) (rdw(i,j,1,2,1),i=1,imax)
          read(wind_int_d_dat_id) (rdw(i,j,2,2,1),i=1,imax)
!dmr [DEPRECATED]        endif
#endif
        read(prcp_int_dat_id) (rrn(i,j,1,1),i=1,imax)
        read(prcp_int_d_dat_id) (rrn(i,j,2,1),i=1,imax)
        read(cldsa_int_dat_id) (rcd(i,j,1,1),i=1,imax)
        read(cldsa_int_d_dat_id) (rcd(i,j,2,1),i=1,imax)
        read(ten5_int_dat_id) (rhr(i,j,1,1,1),i=1,imax)
        read(ten5_int_dat_id) (rhr(i,j,2,1,1),i=1,imax)
        read(ten5_int_d_dat_id) (rhr(i,j,1,2,1),i=1,imax)
        read(ten5_int_d_dat_id) (rhr(i,j,2,2,1),i=1,imax)
        read(vvent_int_dat_id) (rvv(i,j,1,1),i=1,imax)
        read(vvent_int_d_dat_id) (rvv(i,j,1,2),i=1,imax)
        read(temp_int_dat_id) (rtt(i,j,1,1),i=1,imax)
        read(temp_int_d_dat_id) (rtt(i,j,2,1),i=1,imax)
        read(humid_int_dat_id) (rtr(i,j,1,1),i=1,imax)
        read(humid_int_d_dat_id) (rtr(i,j,2,1),i=1,imax)
        read(runo_int_dat_id) (rru(i,j,1,1),i=1,imax)
        read(runo_int_d_dat_id) (rru(i,j,2,1),i=1,imax)
        read(ust_int_dat_id) (rdv(i,j,1,1),i=1,imax)
        read(ust_int_d_dat_id) (rdv(i,j,2,1),i=1,imax)
        do i=1,imax
        amtt(i,j)=rtt(i,j,1,1)
        amtr(i,j)=rtr(i,j,1,1)
        ampr(i,j)=rpr(i,j,1,1)
!       amwx(i,j)=rdw(i,j,1,1,1)
!       amwy(i,j)=rdw(i,j,2,1,1)
#if ( ICEBERG_MODEL > 0 )
!dmr [DEPRECATED]        if (flgicb) then
          amwx(i,j)=rdw(i,j,1,1,1)
          amwy(i,j)=rdw(i,j,2,1,1)
!dmr [DEPRECATED]        endif
#endif
        amrn(i,j)=rrn(i,j,1,1)
        amcd(i,j)=rcd(i,j,1,1)
        amhrx(i,j)=rhr(i,j,1,1,1)
        amhry(i,j)=rhr(i,j,2,1,1)
        amvv(i,j)=rvv(i,j,1,1)
        amru(i,j)=rru(i,j,1,1)
        amdv(i,j)=rdv(i,j,1,1)
        enddo
        enddo
        do l=1,jjt
        do j=1,jmax
        read(pres_int_dat_id) (rpr(i,j,1,1),i=1,imax)
!       read(pres_int_dat_id) (rtt(i,j,1,1),i=1,imax)
!       read(pres_int_dat_id) (rtr(i,j,1,1),i=1,imax)
!       read(pres_int_dat_id) (rdw(i,j,1,1,1),i=1,imax)
!       read(pres_int_dat_id) (rdw(i,j,2,1,1),i=1,imax)
#if ( ICEBERG_MODEL > 0 )
!dmr [DEPRECATED]        if (flgicb) then
          read(wind_int_dat_id) (rdw(i,j,1,1,1),i=1,imax)
          read(wind_int_dat_id) (rdw(i,j,2,1,1),i=1,imax)
          read(wind_int_d_dat_id) (rdw(i,j,1,2,1),i=1,imax)
          read(wind_int_d_dat_id) (rdw(i,j,2,2,1),i=1,imax)
!dmr [DEPRECATED]        endif
#endif
        read(pres_int_d_dat_id) (rpr(i,j,2,1),i=1,imax)
!       read(pres_int_d_dat_id) (rtt(i,j,2,1),i=1,imax)
!       read(pres_int_d_dat_id) (rtr(i,j,2,1),i=1,imax)
!       read(pres_int_d_dat_id) (rdw(i,j,1,2,1),i=1,imax)
!       read(pres_int_d_dat_id) (rdw(i,j,2,2,1),i=1,imax)
        read(prcp_int_dat_id) (rrn(i,j,1,1),i=1,imax)
        read(prcp_int_d_dat_id) (rrn(i,j,2,1),i=1,imax)
        read(cldsa_int_dat_id) (rcd(i,j,1,1),i=1,imax)
        read(cldsa_int_d_dat_id) (rcd(i,j,2,1),i=1,imax)
        read(ten5_int_dat_id) (rhr(i,j,1,1,1),i=1,imax)
        read(ten5_int_dat_id) (rhr(i,j,2,1,1),i=1,imax)
        read(ten5_int_d_dat_id) (rhr(i,j,1,2,1),i=1,imax)
        read(ten5_int_d_dat_id) (rhr(i,j,2,2,1),i=1,imax)
        read(vvent_int_dat_id) (rvv(i,j,1,1),i=1,imax)
        read(vvent_int_d_dat_id) (rvv(i,j,2,1),i=1,imax)
        read(temp_int_dat_id) (rtt(i,j,1,1),i=1,imax)
        read(temp_int_d_dat_id) (rtt(i,j,2,1),i=1,imax)
        read(humid_int_dat_id) (rtr(i,j,1,1),i=1,imax)
        read(humid_int_d_dat_id) (rtr(i,j,2,1),i=1,imax)
        read(runo_int_dat_id) (rru(i,j,1,1),i=1,imax)
        read(runo_int_d_dat_id) (rru(i,j,2,1),i=1,imax)
        read(ust_int_dat_id) (rdv(i,j,1,1),i=1,imax)
        read(ust_int_d_dat_id) (rdv(i,j,2,1),i=1,imax)
        enddo
        enddo
!       write(95,*) ampr(25,32)
!       write(95,*) amtt(20,40),amtr(20,40),ampr(20,40),amwx(20,40),
!    &              amwy(20,40),amrn(20,40),amcd(20,40),amhrx(20,40),
!    &              amhry(20,40),runoff(20,40)
! modif precip
        nzmdif = 0
        open(newunit=modif_pre_id,file='modif.pre',status='old',err=5)
        read(modif_pre_id,*)
        read(modif_pre_id,*)
        read(modif_pre_id,*)
        read(modif_pre_id,*) nzmdif
        do n=1,nzmdif
          read(modif_pre_id,*)
          read(modif_pre_id,*) i1,i2
          read(modif_pre_id,*) j1,j2
          read(modif_pre_id,*) zmopr
          zmopr=zmopr/(100.*365.0*86400.0)
          do j=j1-1,j2+1
           do i=i1-1,i2+1
            zmodfo(i,j)=zmopr/2.0
           enddo
          enddo
          do j=j1,j2
           do i=i1,i2
            zmodfo(i,j)=zmopr
           enddo
          enddo
        enddo
        close(modif_pre_id)
 5      continue                    ! Landing point when there is an error
                                    ! during the opening of file 'modif.pre'
!
!  Data is given for 66 years, and for the SH (1) and NH (2).
!  Note that year 1989 data is my own extrapolation to ensure last
!    few months of 1988 are forced properly.  STOP PRESS: Now have CFC
!    data going until 1995 based on Elkins et al., 1993 Nature Vol 364
!    pp 780-783.

!=======================================================================
!    Read in atmospheric CFC histories
!=======================================================================
!cfc   if (nn99.eq.2) write(99,977)
!cfc   open(91,form='formatted',file='CFC_atmos_data')
!cfc   do nnyear=1930,1995
!cfc       read(91,979) mmyr, cfc11(nnyear-1929,2),
!cfc &             cfc12(nnyear-1929,2), xcfcn,
!cfc &             cfc11(nnyear-1929,1), cfc12(nnyear-1929,1), xcfcs
!cfc       if (nn99.eq.2) write(99,978) mmyr, cfc11(nnyear-1929,2),
!cfc & cfc12(nnyear-1929,2),cfc11(nnyear-1929,1),cfc12(nnyear-1929,1)
!cfc   enddo
!cfc   close(unit=91)
!cfc   if (nn99.eq.2) write(99,*)
!cfc979    format(i4, 2f8.2, f7.3, 2f8.2, f7.3)
!cfc978    format(i5, 4f11.2)
!cfc977    format('Year:   CFC-11(N): CFC-12(N): CFC-11(S): CFC-12(S):')

!
!   Coefficients for geostrophic derivation of wind stress.
!   (Overland & Colony, 1994; Hibler, 1972).
!
        do k=1,nsamt
          theta(k,1) = radian*theta(k,1)
          theta(k,2) = radian*theta(k,2)
        enddo
!
!  Coefficient for cloudiness effect on LW radiation.
!
!       do j=1,jmax
!         do i=1,imax
!           alat    = asin(covrai(i,j))/radian
!           clat    = (90.0-alat)/10.0
!           indx    = 1+int(clat)
!           ihm     = 1+max(0,sign(1,indx-10))
!           indxp1  = indx+1
!           ihmp1   = 1+max(0,sign(1,indxp1-10))
!           dl      = clat-int(clat)
!           dr      = 1.0-dl
!           do k=1,nsamt
!             cmarsh(k,i,j) = dr*bunker(indx)*cmartc(k,ihm)
!    &                       +dl*bunker(indxp1)*cmartc(k,ihmp1)
!             cmarsh(k,i,j) = dr*cmartc(k,ihm)+dl*cmartc(k,ihmp1)
!           enddo
!         enddo
!       enddo

!      coef for the transition between world ocean and the poles
       do j=1,jmax
         do i=1,imax
!        do 11 i=is1(j),is2(j)
           xlati=asin(covrai(i,j))/radian
           if (xlati.lt.(-65.0)) then
             trans(i,j)= 1.0
           else
             if (xlati.lt.(-50.0)) then
               trans(i,j)=1.0 -(xlati+65.0)/15.0
             else
               if (xlati.lt.55.0) then
                 trans(i,j)=0.0
               else
                 if (xlati.lt.70.0) then
                   trans(i,j)=0.0+(xlati-55)/15.0
                 else
                   if (xlati.lt.90.0) then
                     trans(i,j)=1.0
                   else
                     write(clio3_out_id,*) 'arret dans latitude',i,j,xlati
                   endif
                 endif
               endif
             endif
           endif
        enddo
       enddo
!
        mit=2
        mft=1
        jnddt=jndt
!
      endif
!
      do l=jnddt,jndt
      mit=3-mit
      mft=3-mft
      if (mod(jndt+1,nsamt).eq.0) then
        rewind pres_int_dat_id
        rewind pres_int_d_dat_id
        rewind prcp_int_dat_id
        rewind prcp_int_d_dat_id
        rewind cldsa_int_dat_id
        rewind cldsa_int_d_dat_id
        rewind ten5_int_dat_id
        rewind ten5_int_d_dat_id
        rewind vvent_int_dat_id
        rewind vvent_int_d_dat_id
        rewind temp_int_dat_id
        rewind temp_int_d_dat_id
        rewind humid_int_dat_id
        rewind humid_int_d_dat_id
        rewind runo_int_dat_id
        rewind runo_int_d_dat_id
        rewind ust_int_dat_id
        rewind ust_int_d_dat_id
        rewind wind_int_dat_id
        rewind wind_int_d_dat_id
      endif
      do j=1,jmax
        read(pres_int_dat_id) (rpr(i,j,1,mft),i=1,imax)
!       read(pres_int_dat_id) (rtt(i,j,1,mft),i=1,imax)
!       read(pres_int_dat_id) (rtr(i,j,1,mft),i=1,imax)
!       read(pres_int_dat_id) (rdw(i,j,1,1,mft),i=1,imax)
!       read(pres_int_dat_id) (rdw(i,j,2,1,mft),i=1,imax)
#if ( ICEBERG_MODEL > 0 )
!dmr [DEPRECATED]          if (flgicb) then
          read(wind_int_dat_id) (rdw(i,j,1,1,mft),i=1,imax)
          read(wind_int_dat_id) (rdw(i,j,2,1,mft),i=1,imax)
          read(wind_int_d_dat_id) (rdw(i,j,1,2,mft),i=1,imax)
          read(wind_int_d_dat_id) (rdw(i,j,2,2,mft),i=1,imax)
!dmr [DEPRECATED]          endif
#endif
        read(pres_int_d_dat_id) (rpr(i,j,2,mft),i=1,imax)
!       read(pres_int_d_dat_id) (rtt(i,j,2,mft),i=1,imax)
!       read(pres_int_d_dat_id) (rtr(i,j,2,mft),i=1,imax)
!       read(pres_int_d_dat_id) (rdw(i,j,1,2,mft),i=1,imax)
!       read(pres_int_d_dat_id) (rdw(i,j,2,2,mft),i=1,imax)
        read(prcp_int_dat_id) (rrn(i,j,1,mft),i=1,imax)
        read(prcp_int_d_dat_id) (rrn(i,j,2,mft),i=1,imax)
        read(cldsa_int_dat_id) (rcd(i,j,1,mft),i=1,imax)
        read(cldsa_int_d_dat_id) (rcd(i,j,2,mft),i=1,imax)
        read(ten5_int_dat_id) (rhr(i,j,1,1,mft),i=1,imax)
        read(ten5_int_dat_id) (rhr(i,j,2,1,mft),i=1,imax)
        read(ten5_int_d_dat_id) (rhr(i,j,1,2,mft),i=1,imax)
        read(ten5_int_d_dat_id) (rhr(i,j,2,2,mft),i=1,imax)
        read(vvent_int_dat_id) (rvv(i,j,1,mft),i=1,imax)
        read(vvent_int_d_dat_id) (rvv(i,j,2,mft),i=1,imax)
        read(temp_int_dat_id) (rtt(i,j,1,mft),i=1,imax)
        read(temp_int_d_dat_id) (rtt(i,j,2,mft),i=1,imax)
        read(humid_int_dat_id) (rtr(i,j,1,mft),i=1,imax)
        read(humid_int_d_dat_id) (rtr(i,j,2,mft),i=1,imax)
        read(runo_int_dat_id) (rru(i,j,1,mft),i=1,imax)
        read(runo_int_d_dat_id) (rru(i,j,2,mft),i=1,imax)
        read(ust_int_dat_id) (rdv(i,j,1,mft),i=1,imax)
        read(ust_int_d_dat_id) (rdv(i,j,2,mft),i=1,imax)
        enddo
      jnddt=jndt+1
      enddo
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!       INTERPOLATE DATA.                                              |
!-----------------------------------------------------------------------
!
      do j=1,jmax
        ihm = max(0,sign(1,j-jeq))
        hem  = real(ihm)
      do i=1,imax
!       alat = asin(covrai(i,j))/radian
!       indx = 1+int((90.0-abs(alat))/5.0)
        alat = asin(covrai(i,j))/radian
        clat = (95.0-alat)/10.0
        indx = 1+int(clat)

!     fwat(i,j)=max(zero,((rrn(i,j,1,mft)-rrn(i,j,1,mit))/dtt-
!    &         ((3.0*ttat*ttat-1.0)*rrn(i,j,2,mit)-
!    &          (3.0*ttbt*ttbt-1.0)*rrn(i,j,2,mft))*dtts6+
!    &               amrn(i,j))/dtt)
#if ( ISOOCN >= 1 )
      fwat(i,j,ieau)=max(zero,(rrn(i,j,1,mft)-rrn(i,j,1,mit))/dtt-
#else
      fwat(i,j)=max(zero,(rrn(i,j,1,mft)-rrn(i,j,1,mit))/dtt-
#endif
     &         ((3.0*ttat*ttat-1.0)*rrn(i,j,2,mit)-
     &          (3.0*ttbt*ttbt-1.0)*rrn(i,j,2,mft))*dtts6+
     &         amrn(i,j))
      cloud(i,j)=max(zero,(rcd(i,j,1,mft)-rcd(i,j,1,mit))/dtt-
     &         ((3.0*ttat*ttat-1.0)*rcd(i,j,2,mit)-
     &          (3.0*ttbt*ttbt-1.0)*rcd(i,j,2,mft))*dtts6+
     &         amcd(i,j))
      tabq(i,j)=(rtt(i,j,1,mft)-rtt(i,j,1,mit))/dtt-
     &          ((3.0*ttat*ttat-1.0)*rtt(i,j,2,mit)-
     &          (3.0*ttbt*ttbt-1.0)*rtt(i,j,2,mft))*dtts6+
     &          amtt(i,j)+273.15
!     tdew(i,j)=(rtr(i,j,1,mft)-rtr(i,j,1,mit))/dtt-
!    &          ((3.0*ttat*ttat-1.0)*rtr(i,j,2,mit)-
!    &          (3.0*ttbt*ttbt-1.0)*rtr(i,j,2,mft))*dtts6+
!    &                amtr(i,j)+273.15
      tdew(i,j)=(rtr(i,j,1,mft)-rtr(i,j,1,mit))/dtt-
     &          ((3.0*ttat*ttat-1.0)*rtr(i,j,2,mit)-
     &          (3.0*ttbt*ttbt-1.0)*rtr(i,j,2,mft))*dtts6+
     &          amtr(i,j)
      psl(i,j)=((rpr(i,j,1,mft)-rpr(i,j,1,mit))/dtt-
     &          ((3.0*ttat*ttat-1.0)*rpr(i,j,2,mit)-
     &          (3.0*ttbt*ttbt-1.0)*rpr(i,j,2,mft))*dtts6+
     &          ampr(i,j))*100.0
!     psbq(i,j)=hem*101400.0+(1.0-hem)*98800.0
      psbq(i,j)=psl(i,j)
      vabqec(i,j)=(rvv(i,j,1,mft)-rvv(i,j,1,mit))/dtt-
     &          ((3.0*ttat*ttat-1.0)*rvv(i,j,2,mit)-
     &          (3.0*ttbt*ttbt-1.0)*rvv(i,j,2,mft))*dtts6+
     &          amvv(i,j)
      runoff(i,j)=(rru(i,j,1,mft)-rru(i,j,1,mit))/dtt-
     &          ((3.0*ttat*ttat-1.0)*rru(i,j,2,mit)-
     &          (3.0*ttbt*ttbt-1.0)*rru(i,j,2,mft))*dtts6+
     &          amru(i,j)
      sdvt(i,j)=((rdv(i,j,1,mft)-rdv(i,j,1,mit))/dtt-
     &          ((3.0*ttat*ttat-1.0)*rdv(i,j,2,mit)-
     &          (3.0*ttbt*ttbt-1.0)*rdv(i,j,2,mft))*dtts6+
     &          amdv(i,j))
! modification du forcing P-E
!     if (j.eq.6) then
!     write(112,*) i,runoff(i,6),zmodfo(i,6)
!     endif
      runoff(i,j)=runoff(i,j)+zmodfo(i,j)
!
!
!  Longwave radiation:
!----------------------
!
!  vapor pressures (in millibars).
!
!     evg=611.0*10.0**(9.5*(tdew(i,j)-273.16)/(tdew(i,j)-7.660))*0.01
!     evo=611.0*10.0**(7.5*(tdew(i,j)-273.16)/(tdew(i,j)-35.86))*0.01
!     es1=611.0*10.0**(7.5*(tabq(i,j)-273.16)/(tabq(i,j)-7.660))*0.01
!     es2=611.0*10.0**(7.5*(tabq(i,j)-273.16)/(tabq(i,j)-35.86))*0.01
!     evg=tdew(i,j)*es1
!     evo=tdew(i,j)*es2
      es =611.0*exp(min(sign(17.269*one,tabq(i,j)-too),sign(21.875*one,
     &    tabq(i,j)-too))*abs(tabq(i,j)-too)/(tabq(i,j)-35.86+
     &    max(zero,sign(28.2*one,too-tabq(i,j)))))
      evg=tdew(i,j)*es*0.01
      evo=tdew(i,j)*es*0.01
!
!  1. IDSO and JACKSON (all latitudes). (from Parkinson and Washington).
!
!     ratbqg(i,j)=(1.0+cloud(i,j)*0.275)*
!    &            stefan*(tabq(i,j)*tabq(i,j)*tabq(i,j)*tabq(i,j))*
!    &            (1.0-0.261*exp(-7.77e-04*(273.0-tabq(i,j))
!    &             *(273.0-tabq(i,j))))
!     ratbqo(i,j)=(1.0+cloud(i,j)*0.275)*
!    &            stefan*(tabq(i,j)*tabq(i,j)*tabq(i,j)*tabq(i,j))*
!    &            (1.0-0.261*exp(-7.77e-04*(273.0-tabq(i,j))
!    &             *(273.0-tabq(i,j))))
!
!  1.b IDSO, 1981 (all latitudes).
!
!     ratbqg(i,j)=(1.0+cloud(i,j)*0.275)*
!    &            stefan*(tabq(i,j)*tabq(i,j)*tabq(i,j)*tabq(i,j))*
!    &            (0.70+5.95e-05*evg*exp(1500.0/tabq(i,j)))
!     ratbqo(i,j)=(1.0+cloud(i,j)*0.275)*
!    &            stefan*(tabq(i,j)*tabq(i,j)*tabq(i,j)*tabq(i,j))*
!    &            (0.70+5.95e-05*evo*exp(1500.0/tabq(i,j)))
!
!  2. MAYKUT AND CHURCH (Arctic).
!
!     ratbqg(i,j)=0.7855*(1.0+0.2232*cloud(i,j)**2.75)*
!    &            stefan*(tabq(i,j)*tabq(i,j)*tabq(i,j)*tabq(i,j))
!     ratbqo(i,j)=0.7855*(1.0+0.2232*cloud(i,j)**2.75)*
!    &            stefan*(tabq(i,j)*tabq(i,j)*tabq(i,j)*tabq(i,j))
!
!  3. MARSHUNOVA (Arctic).
!
!     cmar        = 0.5*((cmarsh(jmm1,i,j)+cmarsh(jm,i,j))*ttat+
!    &                   (cmarsh(jm,i,j)+cmarsh(jmp1,i,j))*ttbt)
!     ratbqg(i,j) = stefan*(tabq(i,j)*tabq(i,j)*tabq(i,j)*tabq(i,j))*
!    &              (0.67+0.05*sqrt(evg))*(1.0+cloud(i,j)*cmar)
!     ratbqo(i,j) = stefan*(tabq(i,j)*tabq(i,j)*tabq(i,j)*tabq(i,j))*
!    &              (0.67+0.05*sqrt(evo))*(1.0+cloud(i,j)*cmar)
!
!  3.b MARSHUNOVA with cloud correction depending on latitude.
!
!     ratbqg(i,j)=stefan*(tabq(i,j)*tabq(i,j)*tabq(i,j)*tabq(i,j))*
!    &            (1.0-(0.33-0.05*sqrt(evg))*
!    &                 (1.0-budyko(indx)*cloud(i,j)*cloud(i,j)))
!     ratbqo(i,j)=stefan*(tabq(i,j)*tabq(i,j)*tabq(i,j)*tabq(i,j))*
!    &            (1.0-(0.33-0.05*sqrt(evo))*
!    &                 (1.0-budyko(indx)*cloud(i,j)*cloud(i,j)))
!
!  4. BERLIAND (all latitudes).
!
!     ratbqg(i,j)=stefan*tabq(i,j)*tabq(i,j)*tabq(i,j)*tabq(i,j)*
!    &            (1.0-(0.39-0.05*sqrt(evg))*(1.0-budyko(indx)*
!    &              cloud(i,j)*cloud(i,j)))
!     ratbqo(i,j)=stefan*tabq(i,j)*tabq(i,j)*tabq(i,j)*tabq(i,j)*
!    &            (1.0-(0.39-0.05*sqrt(evo))*(1.0-budyko(indx)
!    &             *cloud(i,j)*cloud(i,j)))
      ratbqg(i,j)=-stefan*tabq(i,j)*tabq(i,j)*tabq(i,j)*tabq(i,j)*
     &            (0.39-0.05*sqrt(evg))*(1.0-budyko(indx)*
     &             cloud(i,j)*cloud(i,j))-4.*stefan*tabq(i,j)*
     &             tabq(i,j)*tabq(i,j)*(ts(i,j)-tabq(i,j))+
     &             stefan*ts(i,j)*ts(i,j)*ts(i,j)*ts(i,j)
!     ratbqo(i,j)=-stefan*tabq(i,j)*tabq(i,j)*tabq(i,j)*tabq(i,j)*
!    &            (0.39-0.05*sqrt(evo))*(1.0-budyko(indx)*
!    &             cloud(i,j)*cloud(i,j))-4.*stefan*tabq(i,j)*
!    &             tabq(i,j)*tabq(i,j)*(tcn(i,j,1)-tabq(i,j))+
!    &             stefan*tcn(i,j,1)*tcn(i,j,1)*tcn(i,j,1)*tcn(i,j,1)
!     ratbqo(i,j)=-stefan*tabq(i,j)*tabq(i,j)*tabq(i,j)*tabq(i,j)*
!    &            (0.39-0.05*sqrt(evo))*(1.0-budyko(indx)*
!    &             cloud(i,j)*cloud(i,j))-4.*stefan*tabq(i,j)*
!    &             tabq(i,j)*tabq(i,j)*(scal(i,j,ks2,1)-tabq(i,j))+
!    &             stefan*scal(i,j,ks2,1)*scal(i,j,ks2,1)*
!    &             scal(i,j,ks2,1)*scal(i,j,ks2,1)
      ratbqo(i,j)=zemise*(-stefan*tabq(i,j)*tabq(i,j)*tabq(i,j)
     &             *tabq(i,j)*(0.39-0.05*sqrt(evo))*(1.0-budyko(indx)*
     &             cloud(i,j)*cloud(i,j))-4.*stefan*tabq(i,j)*
     &             tabq(i,j)*tabq(i,j)*(scal(i,j,ks2,1)-tabq(i,j)) )
!  Comparison with the param of Andreas
!     if (i.eq.50.and.j.lt.12) then
!      toto=(0.601+5.95e-5*evo*exp(1500./tabq(i,j)))*
!    &        stefan*tabq(i,j)*tabq(i,j)*tabq(i,j)*tabq(i,j)
!      titi=stefan*scal(i,j,ks2,1)*scal(i,j,ks2,1)*
!    &             scal(i,j,ks2,1)*scal(i,j,ks2,1)
!      fnet=(toto-titi)*(1-.11*8.0*cloud(i,j))
!      write(clio3_out_id,*) 'IR',j,ratbqo(i,j),fnet
!     endif
      uwind(i,j)=(rdw(i,j,1,1,mft)-rdw(i,j,1,1,mit))/dtt-
     &            ((3.0*ttat*ttat-1.0)*rdw(i,j,1,2,mit)-
     &            (3.0*ttbt*ttbt-1.0)*rdw(i,j,1,2,mft))*dtts6+
     &            amwx(i,j)
      write(*,*)'forcng',uwind(i,j)
      vwind(i,j)=(rdw(i,j,2,1,mft)-rdw(i,j,2,1,mit))/dtt-
     &            ((3.0*ttat*ttat-1.0)*rdw(i,j,2,2,mit)-
     &            (3.0*ttbt*ttbt-1.0)*rdw(i,j,2,2,mft))*dtts6+
     &            amwy(i,j)

      tenhrx(i,j)=(rhr(i,j,1,1,mft)-rhr(i,j,1,1,mit))/dtt-
     &             ((3.0*ttat*ttat-1.0)*rhr(i,j,1,2,mit)-
     &             (3.0*ttbt*ttbt-1.0)*rhr(i,j,1,2,mft))*dtts6+
     &             amhrx(i,j)
      tenhry(i,j)=(rhr(i,j,2,1,mft)-rhr(i,j,2,1,mit))/dtt-
     &             ((3.0*ttat*ttat-1.0)*rhr(i,j,2,2,mit)-
     &             (3.0*ttbt*ttbt-1.0)*rhr(i,j,2,2,mft))*dtts6+
     &             amhry(i,j)

      enddo
      enddo
!
! Set up CFC-11 and CFC-12 atmospheric concentrations as a FUNCTION
!   of hemisphere and time (linearly interpolate yearly data).
! 1st some year variables and time weights for atmospheric cfc conc'ns:

!cfc    cfcyear = 1929.0 + float(ja)
!cfc    intyear = ja
!cfc
!cfc    if(cfcyear.eq.1930) then
!cfc       atmcfc11(1) = 0.0
!cfc       atmcfc11(2) = 0.0
!cfc       atmcfc12(1) = 0.0
!cfc       atmcfc12(2) = 0.0
!cfc       goto 41
!cfc    endif
!cfc
!cfc    if(cfcyear.ge.1931 .and. cfcyear.le.1978) then
!cfc! During and before 1978 data is for 31st December:
!cfc       cfcweighta = (xjour)/yeaday
!cfc       cfcweightb = 1.0 - cfcweighta
!cfc       goto 51
!cfc    endif
!cfc
!cfc    if(cfcyear.eq.1979 .and. xjour.le.180.0) then
!cfc!  In first part of 1979, interpolation is with only 6 months gap!
!cfc       cfcweighta =  (xjour)/182.5
!cfc       cfcweightb = 1.0 - cfcweighta
!cfc       goto 51
!cfc    endif
!cfc
!cfc    if(cfcyear.ge.1979) then
!cfc!  Later than July 1st 1979, data is 12 monthly centred on July 1:
!cfc        if(xjour.le.183.0) then
!cfc            cfcweighta = (xjour + 181.0)/yeaday
!cfc            cfcweightb = 1.0 - cfcweighta
!cfc        endif
!cfc        if(xjour.gt.183.0) then
!cfc!  For months after July, interpolation is between year *ahead* and tyear:
!cfc            intyear = intyear + 1
!cfc            cfcweightb = (548.0 - xjour)/yeaday
!cfc            cfcweighta = 1.0 - cfcweightb
!cfc        endif
!cfc    endif
!cfc
!cfc51      continue
!cfc
!cfc! Atmospheric CFC computed as a FUNCTION of time and hemisphere:
!cfc
!cfc       ihemcfc = 1
!cfc       atmcfc11(ihemcfc) =
!cfc $          cfc11(intyear,ihemcfc)*cfcweighta +
!cfc $          cfc11(intyear-1,ihemcfc)*cfcweightb
!cfc       atmcfc12(ihemcfc) =
!cfc $          cfc12(intyear,ihemcfc)*cfcweighta +
!cfc $          cfc12(intyear-1,ihemcfc)*cfcweightb
!cfc
!cfc       ihemcfc = 2
!cfc       atmcfc11(ihemcfc) =
!cfc $          cfc11(intyear,ihemcfc)*cfcweighta +
!cfc $          cfc11(intyear-1,ihemcfc)*cfcweightb
!cfc       atmcfc12(ihemcfc) =
!cfc $          cfc12(intyear,ihemcfc)*cfcweighta +
!cfc $          cfc12(intyear-1,ihemcfc)*cfcweightb

! Test output some CFC weightings to check time interpolation:
!       if(j.eq.20 .and. prntsi .and. eots) then
!          write(stdout,888) intyear, ihemcfc,
!     $ cfc12(intyear,ihemcfc), cfc12(intyear-1,ihemcfc),
!     $ cfcweighta, cfcweightb,tmonth, atmcfc12
!888       format('Year+1=',i3,' hem=',i2,'cfc12(t1)=',f6.2,
!     $ 'cfc12(t0)=',f6.2,'weight to t1=',f4.2,'weight to t0=',
!     $ f4.2,'since tmonth=',f4.1,'==> CFC-12=',f6.2)
!       endif

!41    continue

!  Set some CFC variables for printing:
!       if(j.eq.15 .and. prntsi .and. eots) then
!          cfc11atm = atmcfc11
!          cfc12atm = atmcfc12
!          write(stdout,889) atmcfc11, atmcfc12
!889        format('Calculated values for Southern Hemisphere CFCs:',
!    $ ' CFC-11 =',f8.3,'  CFC-12=',f8.3)
!       endif
!      if ((xjour/50.0).eq.int(xjour/50.)) then
!      write(117,*) ja,xjour,atmcfc11
!      write(118,*) ja,xjour,atmcfc12
!      endif

      if (iwind.eq.1) go to 400
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  3) Wind velocity and wind stress (OVERLAND and COLONY, 1994).       |
!-----------------------------------------------------------------------
!
!--3.1. Calculate geostrophic wind at the corners of the grid squares.
!---------------------------------------------------------------------
!
      rhlim=2.0*omega*sin(15.0*radian)
      do j=2,jmax
          jm1 = j-1
        do i=1,imax
          im1      =  (i-1)+(imax-2)*(1/i)
          rhoa     =  psbq(i,j)/(287.04*tabq(i,j))
!         rhoafn   =  rhoa*zfn(i,j)+1d-18
          rhoafn   =  rhoa*sign(one,zfn(i,j))*max(abs(zfn(i,j)),rhlim)
!         write(89,*) i,j,rhoa,zfn(i,j)
          ugw(i,j) = -(-(alambd(i,j,1,1,2,1)+
     &                       alambd(i,j,1,2,2,1))*psl(i,jm1)-
     &                 (alambd(i,j,1,1,1,1)
     &                      +alambd(i,j,1,2,1,1))*psl(im1,jm1)+
     &                 (alambd(i,j,1,1,2,2)
     &                      -alambd(i,j,1,2,2,2))*psl(i,j)+
     &                 (alambd(i,j,1,1,1,2)
     &                      -alambd(i,j,1,2,1,2))*psl(im1,j))/
     &                rhoafn

          vgw(i,j) =  ((alambd(i,j,2,2,2,1)
     &                     -alambd(i,j,2,1,2,1))*psl(i,jm1)+
     &                 (alambd(i,j,2,2,2,2)
     &                     -alambd(i,j,2,1,2,2))*psl(i,j)-
     &                 (alambd(i,j,2,2,1,1)
     &                     +alambd(i,j,2,1,1,1))*psl(im1,jm1)-
     &                 (alambd(i,j,2,2,1,2)
     &                     +alambd(i,j,2,1,1,2))*psl(im1,j))/
     &                rhoafn
        enddo
        enddo
!  South pole: never used.
!
      do i=1,imax
         ugw(i,1) = 0.0
         vgw(i,1) = 0.0
      enddo
!
!  Equator.
!
      njh = jeq-1
      do i=1,imax
        ugw(i,njh+1) =0.5*(ugw(i,njh)+ugw(i,njh+2))
        vgw(i,njh+1) =0.5*(vgw(i,njh)+vgw(i,njh+2))
      enddo
!
!  Does not use coastal points for determining wind
!  at oceanic grid.
!
      do j=js1,js2
        jp1 = j+1
        do i=is1(j),is2(j)
          ip1         = (i+1)
          usp         = 1.0/
     &    max(1.0e-12*one,tmu(i,j,ks2)+tmu(ip1,j,ks2)
     &                   +tmu(i,jp1,ks2)+tmu(ip1,jp1,ks2))
          uwind(i,j)  = (ugw(i,j)*tmu(i,j,ks2)+
     &                   ugw(ip1,j)*tmu(ip1,j,ks2)+
     &                   ugw(i,jp1)*tmu(i,jp1,ks2)+
     &                   ugw(ip1,jp1)*tmu(ip1,jp1,ks2))*usp
          vwind(i,j)  = (vgw(i,j)*tmu(i,j,ks2)+
     &                   vgw(ip1,j)*tmu(ip1,j,ks2)+
     &                   vgw(i,jp1)*tmu(i,jp1,ks2)+
     &                   vgw(ip1,jp1)*tmu(ip1,jp1,ks2))*usp
        enddo
      enddo
!
400   continue
!
!--3.2. Calculate wind and wind stress at the center of the grid squares.
!------------------------------------------------------------------------
!
!  gam:  enhancing factor of OVERLAND and COLONY, 1994.
!
      gam = 1.3
      do j=js1,js2
!       ihm   = max(0,sign(1,j-jeq))
!       hem   = 2.0*real(ihm)-1.0
!       thet  = 0.5*((theta(jmm1,2-ihm)+theta(jm,2-ihm))*ttat+
!    &               (theta(jm,2-ihm)+theta(jmp1,2-ihm))*ttbt)
!       rhc2  = 0.5*((rhcg2(jmm1,2-ihm)+rhcg2(jm,2-ihm))*ttat+
!    &               (rhcg2(jm,2-ihm)+rhcg2(jmp1,2-ihm))*ttbt)
!       sint  = sin(thet)
!       cost  = cos(thet)
!       costa = cos(15.0*radian)
!       sinta = sin(15.0*radian)
        do i=is1(j),is2(j)
!         vabq(i,j)   = max(0.1*one,sqrt(uwind(i,j)*uwind(i,j)
!    &                       +vwind(i,j)*vwind(i,j)))
!         uwind(i,j)  = gam*uwind(i,j)
!         vwind(i,j)  = gam*vwind(i,j)
!         vabq(i,j)   = gam*vabq(i,j)
!         tenagx(i,j) = rhc2*vabq(i,j)*
!    &                  (uwind(i,j)*cost-vwind(i,j)*sint*hem)
!         tenagy(i,j) = rhc2*vabq(i,j)*
!    &                  (vwind(i,j)*cost+uwind(i,j)*sint*hem)
!
!         rhoa     =  psbq(i,j)/(287.04*tabq(i,j))
!         tairox(i,j) = cdb(vabq(i,j),tabq(i,j)
!    &                  -scal(i,j,ks2,1))*rhoa*vabq(i,j)*
!    &                  (uwind(i,j)*costa-vwind(i,j)*sinta*hem)
!         tairoy(i,j) = cdb(vabq(i,j),tabq(i,j)
!    &                  -scal(i,j,ks2,1))*rhoa*vabq(i,j)*
!    &                  (vwind(i,j)*costa+uwind(i,j)*sinta*hem)
!         tairox(i,j) = 1.3e-03*rhoa*vabq(i,j)*
!    &                  (uwind(i,j)*costa-vwind(i,j)*sinta*hem)
!         tairoy(i,j) = 1.3e-03*rhoa*vabq(i,j)*
!    &                  (vwind(i,j)*costa+uwind(i,j)*sinta*hem)
!
!         tairox(i,j) = trans(i,j)*tairox(i,j)
!    &                 +(1-trans(i,j))*tenhrx(i,j)
!         tairoy(i,j) = trans(i,j)*tairoy(i,j)
!    &                 +(1-trans(i,j))*tenhry(i,j)
! vent de Esbensen et Kushnir

          vabq(i,j)=vabqec(i,j)
!
! tension Trenberth
          zmod=1.0
          angt=0.*radian
          costt=cos(angt)
          sintt=sin(angt)
          tairox(i,j) = zmod*(tenhrx(i,j)*costt-tenhry(i,j)*sintt*hem)
          tairoy(i,j) = zmod*(tenhrx(i,j)*sintt*hem+tenhry(i,j)*costt)
          angt=0.*radian
          costb=cos(angt)
          sintb=sin(angt)
          tenagx(i,j) = zmod*(tenhrx(i,j)*costb-tenhry(i,j)*sintb*hem)
          tenagy(i,j) = zmod*(tenhrx(i,j)*sintb*hem+tenhry(i,j)*costb)

        enddo
        enddo
!
      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine forcng -
      end
