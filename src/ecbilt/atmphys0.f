!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:35 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:35 CET 2009
!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_iatmphys
!-----------------------------------------------------------------------
! *** initializes variables used in SUBROUTINEs of atmphys.f
!-----------------------------------------------------------------------

      use comatm, only: iwater,phi,p0,initfield
      use comdyn, only: rmount,geopg
      use comphys, only: nlon,nlat,irn,tzero,tempm,temp0g,temp2g,temp4g,
     &                   relhum,tcc,rmoisg,ccisccp,torain,tosnow,ipl,
     &                   tncep,pncep,qancep,z500ncep,ghgipcc,pisccp,
     &                   qmount,lwrt,lwrts,lwrqa,lwrqts,lwrref,lwrghg,
     &                   ampwir,ampeqir,ampanir,ampanir2,hprofw,hprofan,
     &                   hprofan2,hproftrop,hprofeq,sulopt,nyscenmax,
     &                   y1scenghg,iscenghg,isghgstrt,ghgscen,iscensul,
     &                   isceno3,solarvol,solartsi,solarcl,iscentsi,
     &                   iscenvol,solarm,isvolstrt,iso3strt,ghg,costref,
     &                   rlogtl,rlogts,tncep12,facttsi,istsistrt,swrref,
     &                   swrcost,swrsalb,salbref,dragan,dragane,dragla
      use comemic_mod, only: fini, imonth, iyear
      use comrunlabel_mod, only: irunlabelf
      use comsurf_mod, only: nse,noc,nld,epss,ntyps,tsurf,tsurfn,tempsg,
     &                       tempsgn,q10,q10n,qsurf,qsurfn,rmountn,
     &                       fractn,pground,pgroundn
      use comunit, only: iuo
      use global_constants_mod, only: silp=>sp,  dblp=>dp, ip
      use newunit_mod, only: lwrref_dat_id,lwrcoef_dat_id,swrref_dat_id,
     &     swrcoef_dat_id,GHG_dat_id,TSI_RM_dat_id,VOLC_dat_id,
     &     SUL_dat_id,OZONE_dat_id,scenario2Xco2_dat_id,error_id,info_id
      
#if ( IMSK == 1 )
      USE input_icemask, ONLY: icemask
#endif
#if ( CYCC == 2 )
      USE carbone_co2
#endif

#if (WISOATM == 1 )
       USE comatm, ONLY : iwater, nwisos, iwat18, iwat17, iwat2h
       USE iso_param_mod, only: datmini, delta_inv, isoatm_restart
#endif

#if ( DOWNSTS == 1 )
      use vertDownsc_mod,    only: tsurfn_d, tsurf_d
#endif

#if ( ISM == 2 )
      use comemic_mod,       only: fracto
#endif

#define MODIS_CLDS 0
#if ( MODIS_CLDS == 1 )
      use global_constants_mod, only: str_len,dblp=>dp,ip,months_year_i
      use ncio, only: nc_read
#endif

      implicit none

      real*8 fvolcan
      integer ios
      integer(kind=ip):: inatphy_dat_id, wisoatm_restart_dat_id
     &                   , orography_ctl_id, orography_dat_id
     &                   , landfrac_ctl_id, landfrac_dat_id 


      integer i,j,k,l,ireg,im,nn,is,j1,i1,ii,jj,ism
      real*8  beta,draganr,draglar,dum(2),asum,spv
      integer jyear,ilat,jmonth,m,mfin,indxvol,indxtsi
      real*8 tsi,ksw

comphys.f90:      real(silp), dimension(0:iqmtab,0:jqmtab,0:kqmtab):: qmtabel
      real(silp), dimension(27,12)        :: costref_silp, salbref_silp
      real(silp), dimension(8,27,12,0:1)  :: swrref_silp, swrcost_silp
      real(silp), dimension(8,27,12,0:3)  :: swrsalb_silp
      real(silp), dimension(19,27,12)     :: tncep_silp
      real(silp), dimension(27,12)        :: qancep_silp
      real(silp), dimension(19)           :: ghgipcc_silp
      real(silp), dimension(27,12)        :: z500ncep_silp
      real(silp), dimension(27)           :: pisccp_silp
      real(silp), dimension(17)           :: pncep_silp
      real(silp), dimension(32,64,12)     :: ccisccp_silp
      real(silp), dimension(7,27,4,0:1)   :: lwrref_silp
      real(silp), dimension(7,18,27,4,0:1):: lwrt_silp
      real(silp), dimension(7,4,27,4,0:1) :: lwrts_silp
      real(silp), dimension(7,27,4,0:1)   :: lwrqa_silp
      real(silp), dimension(7,4,27,4,0:1) :: lwrqts_silp
      real(silp), dimension(7,19,27,4,0:1):: lwrghg_silp


#if (WISOATM == 1 )
c~       REAL, DIMENSION(nlat,nlon,neauiso) :: rations
c~       REAL, DIMENSION(nlat,nlon,neauiso) :: ration_qatm
      real(kind=dblp), dimension(nlat,nlon), parameter :: Ones_Ecbilt = 1.0

#endif

#if ( ISM == 2 )
!dmr Added for orography output !
      real*4  outdata(nlon,nlat)
#endif

#if ( MODIS_CLDS == 1 )
      character(str_len), parameter ::
     >    file_clouds="inputdata/clt_MODIS_000001-000012_ac-T21.nc"
      real(dblp), dimension(nlon,nlat,months_year_i) :: modis_clds
      integer(ip) :: t
#endif
!
! *** initial atmospheric temperatures and moisture
!
      if (initfield.eq.0) then
        tempm(0)=tzero-35d0
        tempm(1)=tzero-35d0
        tempm(2)=tzero-8d0
        do j=1,nlon
          do i=1,nlat
            temp0g(i,j)=tempm(0)
            temp2g(i,j)=tempm(1)
            temp4g(i,j)=tempm(2)
            tempsg(i,j)=290d0
c~ #if ( ISOATM >= 1 )
c~             rmoisg(i,j,:)=0d0
c~ #else
            rmoisg(i,j,:)=0d0
c~ #endif
            relhum(i,j)=0d0
            q10(i,j)=0d0
            qsurf(i,j)=0d0
            do nn=1,ntyps
              q10n(i,j,nn)=0.0d0
              qsurfn(i,j,nn)=0.0d0
              tempsgn(i,j,nn)=290d0
              pgroundn(i,j,nn)=p0
            enddo
            pground(i,j)=p0
            geopg(i,j,2)=0d0
            tcc(i,j)=ccisccp(i,j,1)
          enddo
        enddo


        do j=1,nlon
          do i=1,nlat
#if ( IMSK == 1 )
            if ( icemask(i,j).gt.0.9) then
              ! do not set temperatures above freezing on ice
              tempsg(i,j) = min(tzero,tempsg(i,j))
            endif
#endif
#if ( DOWNSTS == 1 )
            tsurfn_d(i,j,noc,:)= tempsg(i,j)
            tsurfn_d(i,j,nse,:)= tzero
            tsurfn_d(i,j,nld,:)= tempsg(i,j)
            tsurf_d (i,j,:)    = tempsg(i,j)
#endif
!dmr ajouter icemask et tzero limitation here ...
            tsurfn (i,j,noc)= tempsg(i,j)
            tsurfn (i,j,nse)= tzero
            tsurfn (i,j,nld)= tempsg(i,j)
            tsurf (i,j)= tempsg(i,j)
          enddo
        enddo
        call ec_atmphyszero
      else
        ios=0
        open(newunit=inatphy_dat_id,
     &     file='startdata/inatphy'//fini//'.dat',form='unformatted')
        read(inatphy_dat_id) tsurfn,tempm,temp0g
c~ #if ( ISOATM >= 1 )
c~ ! Petite astuce pour avoir un point de redémarrage compatible avec toutes les
c~ !  versions du modèle : on ne mets pas le redémarrage des isotopes dans le
c~ !  même fichier. En attendant, on peut initialiser à une valeur fixée au rapport
c~ !  que l'on veut, en proportion par rapport à l'eau
c~         read(inatphy_dat_id) rmoisg(:,:,ieau),torain(:,:,ieau),
c~      &      tosnow(:,:,ieau)


c~ ! Même si l'on passe l'intégralité du tableau à la routine (i.e. y compris l'eau totale)
c~ !  la première dimension du tableau n'est pas modifiée par cet appel. D'où l'initialisation
c~ !  de "rations" pour toutes dimentions en ratios (i.e. hors les deux premières)

c~ ! En outre, delta_i = R_i / R_ismow - 1.0 donc :
c~ !   R_i = ( delta_i + 1.0 ) * R_ismow
c~ #if ( ISOATM == 3 )
c~ ! Ici a mettre une initialisation basee sur le calcul des R "cursifs"
c~ !        rations(:,:,ieau18) = ( -10.0D-3 + 1.0D0 ) * r18smow
c~ !        rations(:,:,ieau17) = ( -100.0D-3 + 1.0D0 ) * r17smow
c~ !        rations(:,:,ieaud)  = ( -100.0D-3 + 1.0D0 ) * rdsmow
c~ !        rations(:,:,ieaud)  = 0.0 ! dummy type
c~ !
c~         WRITE(*,*) "Option non implementee !!! , atmphys0.f"
c~ #elif ( ISOATM == 1 || ISOATM == 2 )
c~ ! dmr
c~ ! Initialization to fixed initial values for the version with "R" in MJ79
c~ ! 17O and deuterium are initialized with d-excess and 17Oexcess to zero
c~ ! w.r.t. 18O
c~ ! dmr
c~         FORALL (k=ieau+1:neauiso)
c~           rations(:,:,k) = (datmini(k)+1.0d0) * rsmow(k)
c~         ENDFORALL

c~ #endif
c~ #if ( ISOATM == 3 )
c~         WRITE(*,*) "Option non implementee !!! , atmphys0.f"

        read(inatphy_dat_id) rmoisg(:,:,iwater)
     >                      ,torain(:,:,iwater),tosnow(:,:,iwater)
     

#if ( WISOATM == 1 )
        if (isoatm_restart.EQ.0) then

        do i=iwat17,iwat2h
          rmoisg(:,:,i) = 
     >       delta_inv(rmoisg(:,:,iwater),Ones_Ecbilt(:,:)*datmini(i),i)
          torain(:,:,i) = 
     >       delta_inv(torain(:,:,iwater),Ones_Ecbilt(:,:)*datmini(i),i)
          tosnow(:,:,i) = 
     >       delta_inv(tosnow(:,:,iwater),Ones_Ecbilt(:,:)*datmini(i),i)     
        enddo

        else ! let's read the restart file

      OPEN(newunit=wisoatm_restart_dat_id,
     & FILE='startdata/wisoatm_restart.dat',STATUS='old'
     &,FORM='unformatted')

           READ(wisoatm_restart_dat_id)  rmoisg(:,:,iwater+1:nwisos)
     &              ,torain(:,:,iwater+1:nwisos)
     &              ,tosnow(:,:,iwater+1:nwisos)

         CLOSE(wisoatm_restart_dat_id)
        endif

#endif
! Fin de l'initialisation des isotopes sans/avec relecture de fichier


        close(iuo+95)
      endif

#if ( ISM == 2 )
! output of orography and land-fraction
      open(newunit=orography_ctl_id,
     &  file='outputdata/atmos/orography.ctl')
      write(orography_ctl_id,fmt="('dset   ^orography.dat')")
      write(orography_ctl_id,fmt="('options big_endian')")
      write(orography_ctl_id,fmt="('undef ',1p,e12.4)") -1.0e20
      write(orography_ctl_id,fmt="('title ECBILT orography')")
      write(orography_ctl_id,
     &  fmt="('xdef ',i3,' linear ',2f7.2)") 64,0.00,5.625
      write(orography_ctl_id,fmt="('ydef ',i3,' levels')") 32
      write(orography_ctl_id,
     &  fmt="(' -85.7606 -80.2688 -74.7445 -69.2130 -63.6786')")
      write(orography_ctl_id,
     &  fmt="(' -58.1430 -52.6065 -47.0696 -41.5325 -35.9951')")
      write(orography_ctl_id,
     &  fmt="(' -30.4576 -24.9199 -19.3822 -13.8445 -8.30670')")
      write(orography_ctl_id,
     &  fmt="(' -2.76890 2.76890 8.30670 13.8445 19.3822')")
      write(orography_ctl_id,
     &  fmt="('  24.9199 30.4576 35.9951 41.5325 47.0696')")
      write(orography_ctl_id,
     &  fmt="('  52.6065 58.1430 63.6786 69.2130 74.7445')")
      write(orography_ctl_id,fmt="('  80.2688 85.7606')")
      write(orography_ctl_id,
     &  fmt="('zdef ',i3,' linear ',2f7.2)") 1,0.00,1.00
      write(orography_ctl_id,
     &  fmt="('tdef ',i4,' linear 1jan0001  1YR')") 1
      write(orography_ctl_id,fmt="('vars 1')")
      write(orography_ctl_id,
     &  fmt="('orog       1  99 orography ECBILT')")
      write(orography_ctl_id,fmt="('endvars')")
      close(orography_ctl_id)
      do i=1,nlon
        do j=1,nlat
          outdata(i,j)=rmount(j,i)
        enddo
      enddo
      open(newunit=orography_dat_id,
     &  CONVERT='BIG_ENDIAN',file='outputdata/atmos/orography.dat
     &         ',form='unformatted',
     &         access='direct',recl=Size(outdata)*Kind(outdata(1,1)))
      write(orography_dat_id,REC=1) outdata
      close(orography_dat_id)
      open(newunit=landfrac_ctl_id,file='outputdata/atmos/landfrac.ctl')
      write(landfrac_ctl_id,fmt="('dset   ^landfrac.dat')")
      write(landfrac_ctl_id,fmt="('options big_endian')")
      write(landfrac_ctl_id,fmt="('undef ',1p,e12.4)") -1.0e20
      write(landfrac_ctl_id,fmt="('title ECBILT land fraction')")
      write(landfrac_ctl_id,
     &  fmt="('xdef ',i3,' linear ',2f7.2)") 64,0.00,5.625
      write(landfrac_ctl_id,fmt="('ydef ',i3,' levels')") 32
      write(landfrac_ctl_id,
     &  fmt="(' -85.7606 -80.2688 -74.7445 -69.2130 -63.6786')")
      write(landfrac_ctl_id,
     &  fmt="(' -58.1430 -52.6065 -47.0696 -41.5325 -35.9951')")
      write(landfrac_ctl_id,
     &  fmt="(' -30.4576 -24.9199 -19.3822 -13.8445 -8.30670')")
      write(landfrac_ctl_id,
     &  fmt="(' -2.76890 2.76890 8.30670 13.8445 19.3822')")
      write(landfrac_ctl_id,
     &  fmt="('  24.9199 30.4576 35.9951 41.5325 47.0696')")
      write(landfrac_ctl_id,
     &  fmt="('  52.6065 58.1430 63.6786 69.2130 74.7445')")
      write(landfrac_ctl_id,fmt="('  80.2688 85.7606')")
      write(landfrac_ctl_id,
     &  fmt="('zdef ',i3,' linear ',2f7.2)") 1,0.00,1.00
      write(landfrac_ctl_id,
     &  fmt="('tdef ',i4,' linear 1jan0001  1YR')") 1
      write(landfrac_ctl_id,fmt="('vars 1')")
      write(landfrac_ctl_id,
     &  fmt="('landfrac       1  99 land fraction ECBILT')")
      write(landfrac_ctl_id,fmt="('endvars')")
      close(landfrac_ctl_id)
      do i=1,nlon
        do j=1,nlat
          outdata(i,j)=1.0-fracto(j,i)
        enddo
      enddo
      open(newunit=landfrac_dat_id,
     &  CONVERT='BIG_ENDIAN',file='outputdata/atmos/landfrac.dat'
     & ,form='unformatted',
     &         access='direct',recl=Size(outdata)*Kind(outdata(1,1)))
      write(landfrac_dat_id,REC=1) outdata
      close(landfrac_dat_id)
#endif
      do j=1,nlon
        do i=1,nlat
          rmountn(i,j,nld)=rmount(i,j)
          rmountn(i,j,nse)=0d0
          rmountn(i,j,noc)=0d0
          if (rmountn(i,j,nld).lt.0d0) rmountn(i,j,nld)=0d0
          if (fractn(i,j,nld).lt.epss) then
            rmountn(i,j,nld)=0d0
            qmount(i,j)=0d0
          else
            qmount(i,j)=rmountn(i,j,nld)
          endif
        enddo
      enddo

      do j=1,nlon
        do i=1,nlat
          asum=0d0
          do i1=-1,1
            do j1=-1,1
              ii=i+i1
              jj=j+j1
              if (ii.lt.1) then
                ii=1
                jj=jj+nlon/2
              endif
              if (ii.gt.nlat) then
                ii=nlat
                jj=jj+nlon/2
              endif
              if (jj.lt.1) jj=jj+nlon
              if (jj.gt.nlon) jj=jj-nlon
              asum=asum+rmountn(ii,jj,nld)
            enddo
          enddo
          qmount(i,j)=asum/9d0
        enddo
      enddo


!***  longwave radiation parameterisation, based on a linearisation
!***  of KRCM with respect to reference T and q profiles and other
!***  variables (greenhousegases)


      read(lwrref_dat_id) irn,ipl,pisccp_silp,pncep_silp,z500ncep_silp
      read(lwrref_dat_id) tncep_silp,qancep_silp,ghgipcc_silp,
     &                      ccisccp_silp
      read(lwrref_dat_id) lwrref_silp
      
      pisccp   = dble(pisccp_silp)
      pncep    = dble(pncep_silp)
      z500ncep = dble(z500ncep_silp)
      tncep    = dble(tncep_silp)
      qancep   = dble(qancep_silp)
      ghgipcc  = dble(ghgipcc_silp)
      ccisccp  = dble(ccisccp_silp)      
      lwrref   = dble(lwrref_silp)


      read(lwrcoef_dat_id) lwrt_silp,lwrts_silp,lwrqts_silp,lwrqa_silp
     &            ,lwrghg_silp
     
     
      lwrt   = dble(lwrt_silp)
      lwrts  = dble(lwrts_silp)
      lwrqts = dble(lwrqts_silp)
      lwrqa  = dble(lwrqa_silp)
      lwrghg = dble(lwrghg_silp)

#if ( MODIS_CLDS == 1 )

      call nc_read(file_clouds,"clt_CdmsRegrid",modis_clds)

      do j=1,nlon
        do i=1,nlat
          do t=1,months_year_i
            ccisccp(i,j,t)  = modis_clds(j,i,t)/100.0
          enddo
        enddo
      enddo

#endif


!**   amplifcation of freshwater feedback
      lwrqa(:,:,:,:)=lwrqa(:,:,:,:)*AMPWIR
!**   amplifcation of feedback at the equator
      lwrqa(:,9,:,:)=lwrqa(:,9,:,:)*(1.+(AMPEQIR-1.0)/2.0)
      lwrqa(:,10,:,:)=lwrqa(:,10,:,:)*(1.+(AMPEQIR-1.0)/2.0)
      lwrqa(:,11,:,:)=lwrqa(:,11,:,:)*AMPEQIR
      lwrqa(:,12,:,:)=lwrqa(:,12,:,:)*AMPEQIR
      lwrqa(:,13,:,:)=lwrqa(:,13,:,:)*(1.+(AMPEQIR-1.0)/2.0)
      lwrqa(:,14,:,:)=lwrqa(:,14,:,:)*(1.+(AMPEQIR-1.0)/2.0)
      lwrqa(:,23,:,:)=lwrqa(:,23,:,:)*AMPANIR
      lwrqa(:,1,:,:)=lwrqa(:,1,:,:)*AMPANIR2
      lwrqa(:,2,:,:)=lwrqa(:,2,:,:)*AMPANIR2
      lwrqa(:,3,:,:)=lwrqa(:,3,:,:)*AMPANIR2

!     write(info_id,*) "Modif IR scheme"
!     write(info_id,*) AMPWIR,AMPEQIR,expIR
!     write(info_id,*) HPROFW,HPROFTROP,1.0+(HPROFTROP-1.0)/2.0,HPROFEQ

!***  update moisture profile used in the linearization of the
!***  radiative scheme in the tropics:easy surrogate to a change
!***  mean IR flux in the model without affecting sensitivity
      do ism=1,12
        qancep(1,ism)=qancep(1,ism)*HPROFAN2
        qancep(2,ism)=qancep(2,ism)*HPROFAN2
        qancep(3,ism)=qancep(3,ism)*HPROFAN2
        qancep(4,ism)=qancep(4,ism)*HPROFW
        qancep(5,ism)=qancep(5,ism)*(1.0+(HPROFTROP-1.0)/2.0)
        qancep(6,ism)=qancep(6,ism)*(1.0+(HPROFTROP-1.0)/2.0)
        qancep(7,ism)=qancep(7,ism)*HPROFTROP
        qancep(8,ism)=qancep(8,ism)*HPROFTROP
        qancep(9,ism)=qancep(9,ism)*HPROFTROP
        qancep(10,ism)=qancep(10,ism)*HPROFTROP

        qancep(11,ism)=qancep(11,ism)*HPROFEQ
        qancep(12,ism)=qancep(12,ism)*HPROFEQ

        qancep(13,ism)=qancep(13,ism)*HPROFTROP
        qancep(14,ism)=qancep(14,ism)*HPROFTROP
        qancep(15,ism)=qancep(15,ism)*HPROFTROP
        qancep(16,ism)=qancep(16,ism)*HPROFTROP
        qancep(17,ism)=qancep(17,ism)*(1.0+(HPROFTROP-1.0)/2.0)
        qancep(18,ism)=qancep(18,ism)*(1.0+(HPROFTROP-1.0)/2.0)
        qancep(19,ism)=qancep(19,ism)*HPROFW
        qancep(20,ism)=qancep(20,ism)*HPROFW
        qancep(21,ism)=qancep(21,ism)*HPROFW
        qancep(22,ism)=qancep(22,ism)*HPROFW
        qancep(23,ism)=qancep(23,ism)*HPROFAN
        qancep(24,ism)=qancep(24,ism)*HPROFW
        qancep(25,ism)=qancep(25,ism)*HPROFW
        qancep(26,ism)=qancep(26,ism)*HPROFW
        qancep(27,ism)=qancep(27,ism)*HPROFW
      enddo
#if ( IMSK == 1 )
      CALL ec_masq(.TRUE.)
#endif

!***  read GHG concentrations
      if (iscenghg.eq.1) then
!* 0-nyscenmax GHG concentrations
        read(GHG_dat_id,*) nyscenmax,y1scenghg
        write(info_id,*) 'scen GHG',nyscenmax,y1scenghg
        do i=1,nyscenmax
#if ( F_PALAEO >= 1 )
         read(GHG_dat_id,*) ii, (ghgscen(j,i),j=1,19)
#else
         read(GHG_dat_id,100)(ghgscen(j,i),j=1,19)
#endif
        enddo
      else
!*  "standard" GHG concentrations, 1st line of forcing file
!    -> iscenghg = 0, 2, 3 &
        read(GHG_dat_id,*) nyscenmax,y1scenghg
        read(GHG_dat_id,100)(ghgscen(j,1),j=1,19)
      endif

!***  read 2 times CO2 concentrations
!     -> other GHGs at their "standard" concentrations (1st line of file)
      if (iscenghg.eq.2) then
        nyscenmax=2100
        do i=1,nyscenmax
         read(scenario2Xco2_dat_id,*) ii,ghgscen(1,i)
         do j=2,19
          ghgscen(j,i)=ghgscen(j,1)
         enddo
        enddo
      endif

      if (iscenghg.eq.3) then
!* get GHG concentrations corresponding to year "isghgstrt"
!    -> use: equilibrium runs
       do i=2,isghgstrt-y1scenghg+1
        read(GHG_dat_id,100)(ghgscen(j,1),j=1,19)
       enddo
      endif

!***  read O3 concentrations
      ghgscen(20,:)=25.0
      if (isceno3.eq.1) then
        do i=1,241
          read(OZONE_dat_id,*)j,ghgscen(20,i)
        enddo
      endif
!***  read 1850-2100 sulfates optical depths
      sulopt(:,:,:)=0.0
      if (iscensul.eq.1) then
        mfin=12*250
        do m=1,mfin
          do i=1,nlat
           read(SUL_dat_id) (sulopt(m,i,j),j=1,nlon)
          enddo
        enddo
      endif

!***  read 0-2000 TSI anomalies

      solartsi(:)=0.0
      if (iscentsi.eq.1) then
        do i=1,2001
         read(TSI_RM_dat_id,*) jyear,solartsi(i)
        enddo
      endif

!***  read 0-2000 anomalies associated with volcanos

      solarvol(:,:,:)=0.
      if (iscenvol.eq.1) then
        do i=1,2001
          do k=1,12
            read(VOLC_dat_id,'(I5,I3,1x,4(F7.3))') jyear,jmonth,
     &          solarvol(i,k,1),solarvol(i,k,2),
     &          solarvol(i,k,3),solarvol(i,k,4)
          enddo
        enddo
      endif

100   FORMAT(4X,3(X,F7.2),9(2X,F6.2),2X,F7.2,6(2X,F6.2))


      i=1
      call ec_ghgupdate(i)

      do ireg=1,27
        do im=1,12
          beta=(tncep(11,ireg,im)-tncep(10,ireg,im))/alog(400./300.)
          tncep12(1,ireg,im)=tncep(10,ireg,im)+beta*alog(350./300.)
          beta=(tncep(14,ireg,im)-tncep(13,ireg,im))/alog(700./600.)
          tncep12(2,ireg,im)=tncep(13,ireg,im)+beta*alog(650./600.)
        enddo
      enddo

      do k=1,17
         rlogtl(k)=log(pncep(k)/65000.)
      enddo

      do ireg=1,27
        rlogts(ireg)=log(pisccp(ireg)/65000.)
      enddo


!
! *** shortwave radiation parameters
!


      read(swrref_dat_id) costref_silp,salbref_silp
      read(swrref_dat_id) swrref_silp


      read(swrcoef_dat_id) swrcost_silp,swrsalb_silp

      costref = dble(costref_silp)
      salbref = dble(salbref_silp)
      swrref  = dble(swrref_silp)
      swrcost = dble(swrcost_silp)
      swrsalb = dble(swrsalb_silp)

      call ec_detqmtabel


310   format(9f7.2)
330   format(2f5.2)
340   format(f5.2)


! *** computation of the effective turning angle


      draglar=3.141592654/180.0*dragla
      draganr=3.141592654/180.0*dragan
      do i=1,nlat
        if (phi(i).lt.(-1.0*draglar)) then
           dragane(i)=-1.0*draganr
        else
           if (phi(i).gt.draglar) then
               dragane(i)=draganr
           else
               dragane(i)=draganr*phi(i)/draglar
           endif
        endif
      enddo

!
! *** evaporation factor
!
!--- moved to initial0      evfac=1d0

      ksw=0.042

      do i=1,32
        solarcl(i)=solarm
      enddo
      indxtsi=1
      indxvol=1
!     if (iscentsi.eq.0.) solartsi(:)=0
!     if (irunlabelf+iyear.lt.istsistrt) solartsi(:)=0.
      if ((iscentsi.eq.1).AND.(irunlabelf+iyear.ge.istsistrt) ) then
        indxtsi=irunlabelf+1-1
        if (indxtsi.gt.2001) indxtsi=2001
        do i=1,32
          solarcl(i)=solarcl(i)+solartsi(indxtsi)*facttsi
        enddo
      endif
      if ((iscenvol.eq.1).AND.(irunlabelf+iyear.ge.isvolstrt) ) then
        indxvol=irunlabelf+1-1
        if (indxvol.gt.2001) indxvol=2001
        do i=1,10
          solarcl(i)=solarcl(i)+solarvol(indxvol,imonth,1)
        enddo
        do i=11,16
          solarcl(i)=solarcl(i)+solarvol(indxvol,imonth,2)
        enddo
        do i=17,22
          solarcl(i)=solarcl(i)+solarvol(indxvol,imonth,3)
        enddo
        do i=23,32
          solarcl(i)=solarcl(i)+solarvol(indxvol,imonth,4)
        enddo
      endif

      if ((isceno3.eq.1).AND.(irunlabelf+iyear.ge.iso3strt) ) then
        do i=1,16
        solarcl(i)=solarcl(i)+(4*0.247*ksw*(ghg(20)-25.))
        enddo
        do i=17,32
        solarcl(i)=solarcl(i)+(4*0.324*ksw*(ghg(20)-25.))
        enddo
      endif

      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_ghgupdate(istep)
!-----------------------------------------------------------------------
! *** updates ghg concentrations: indxghg 1 corresponds to y1scenghg AD
!-----------------------------------------------------------------------

      use comphys, only: ghg,iscenghg,isghgstrt,nyscenmax,y1scenghg,
     &                   ghgscen,isceno3,iso3strt,ghgipcc,mag_alpha,
     &                   lwrref,lwrghg,lwrflux
      use comemic_mod, only: iyear, nstpyear
      use comrunlabel_mod, only: irunlabelf
      use atmos_composition_mod, only: get_PGA_CO2, set_PGA_CO2
      use newunit_mod, only: info_id 
      
!pbakker
#if ( INTERACT_CYCC == 1 && CYCC == 2 )
      use carbone_co2, only: PA_C
#endif
!pbakker

#if ( F_PALAEO >= 1 )
      use palaeo_timer_mod, only: indx_ghg
#endif

      implicit none


      integer i,istep,indxghg,s,r,k,l,m,indxo3,h
      real*8  logco2,sqrch4,sqrn2o
      real*8 alpho3lw(2)
      logical :: result_function_call


      if (mod(istep,nstpyear).eq.1) then


       indxghg=1

       if (iscenghg.eq.1) then
#if ( F_PALAEO >= 1 )
!         indxghg = indx_ghg
! dmr --- To avoid indxghg being 0 at first run at init
         indxghg = max(indx_ghg,1)
#else
         if (((irunlabelf+iyear).ge.isghgstrt) )then
              indxghg=irunlabelf+iyear-y1scenghg+1
            if (iyear.eq.0) indxghg=indxghg+1
            if (indxghg.gt.nyscenmax) indxghg=nyscenmax
         endif
#endif
         WRITE(*,*) "COUCOU GHG !!", irunlabelf, isghgstrt, y1scenghg,
     >      iyear, indxghg, nyscenmax
       endif

       write(info_id,*) 'iscenghg',iscenghg
       if (iscenghg.eq.2) then
            indxghg=irunlabelf+iyear-isghgstrt
            write(info_id,*) 'indxghg',irunlabelf,iyear,isghgstrt
            if (iyear.eq.0) indxghg=indxghg+1
            if (indxghg.gt.nyscenmax) indxghg=nyscenmax
            write(info_id,*) 'indxghg',indxghg
       endif

       do i=1,19
          ghg(i)=ghgscen(i,indxghg)
       enddo

#if ( INTERACT_CYCC == 1 )
       if (PA_C.eq.0) then
        result_function_call = set_PGA_CO2(ghg(1))
        write(info_id,*) 'PA_C eq 0, using reference CO2 (ghg(1))'
       else
        result_function_call = set_PGA_CO2(PA_C)
        ghg(1)=PA_C
       endif
       write(info_id,10) 'ghg forcing with interactive CO2',
     &          indxghg,get_PGA_CO2(),ghg(2),ghg(3),ghg(20)
#else
       result_function_call=set_PGA_CO2(ghg(1))
       write(info_id,11) 'ghg forcing ',indxghg,get_PGA_CO2(),
     &          ghg(2),ghg(3),ghg(20)
#endif

       indxo3=iso3strt-1859
       if ((isceno3.eq.1).and.((irunlabelf+iyear).ge.iso3strt)
     &    ) then
          indxo3=irunlabelf+iyear-1859
          if (indxo3.gt.241)indxo3=241
       endif
        ghg(20)=ghgscen(20,indxo3)
!pbakker        write(info_id,10) 'ghg forcing ',iyear+irunlabelf,ghg(1),
!pbakker     &           ghg(2),ghg(3),ghg(20)
10      format(A12,1i6,4f12.3)
11      format(A12,1i6,4f12.3)
        call flush(info_id)
!*** Update LW reference radiation fluxes using new GHG concentrations
        logco2=mag_alpha*log(ghg(1)/280.0) + log(280.0/ghgipcc(1))
        sqrch4=sqrt(ghg(2))-sqrt(ghgipcc(2))
        sqrn2o=sqrt(ghg(3))-sqrt(ghgipcc(3))
        alpho3lw(1)=153.6
        alpho3lw(2)=201.2
        do h=1,2
        do l=0,1
         do s=1,4
          do r=1,27
           do k=1,7
            lwrflux(k,r,s,l,h)=lwrref(k,r,s,l)+
     *            lwrghg(k,1,r,s,l)*logco2+
     *            lwrghg(k,2,r,s,l)*sqrch4+
     *            lwrghg(k,3,r,s,l)*sqrn2o
            do m=4,19
             lwrflux(k,r,s,l,h)=lwrflux(k,r,s,l,h)+
     *            lwrghg(k,m,r,s,l)*(ghg(m)-ghgipcc(m))
            enddo
              lwrflux(k,r,s,l,h)=lwrflux(k,r,s,l,h)+
     *              lwrghg(k,4,r,s,l)*alpho3lw(h)*(ghg(20)-25.)
           enddo
          enddo
         enddo
        enddo
        enddo
      endif
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_atmphyszero
!-----------------------------------------------------------------------
! *** initializes data arrays to zero for physics routines
!-----------------------------------------------------------------------


      use comphys, only: nlon,nlat,torain,tosnow,cormois,
     &                   dyrain,corain,dysnow,cosnow,
     &                   thforg0,thforg1,thforg2,
     &                   vhforg0,vhforg1,vhforg2

c~ #if ( ISOATM >= 1 )
c~       use iso_param_mod, only : ieau
c~ #endif

#if ( DOWNSCALING == 2 )
      use input_subgrid2L, only: reset_rain_snow_sub_grid
#endif

#if ( DOWNSTS == 1 )
      use vertDownsc_mod, only: dyrain_d, dysnow_d, corain_d, cosnow_d
     & , torain_d, tosnow_d
#endif


      implicit none

      integer i,j


      do j=1,nlon
        do i=1,nlat

          cormois(i,j,:)=0.d0
          torain(i,j,:)=0.d0
          tosnow(i,j,:)=0.d0
          dyrain(i,j,:)=0.d0
          corain(i,j,:)=0.d0
          dysnow(i,j,:)=0.d0
          cosnow(i,j,:)=0.d0

#if ( DOWNSTS == 1 )
          dyrain_d(i,j,:) = 0.0d0
          dysnow_d(i,j,:) = 0.0d0
          corain_d(i,j,:) = 0.0d0
          cosnow_d(i,j,:) = 0.0d0
          torain_d(i,j,:) = 0.0d0
          tosnow_d(i,j,:) = 0.0d0
#endif

          thforg0(i,j)=0.d0
          thforg1(i,j)=0.d0
          thforg2(i,j)=0.d0
          vhforg0(i,j)=0.d0
          vhforg1(i,j)=0.d0
          vhforg2(i,j)=0.d0
        enddo
      enddo

#if ( DOWNSCALING == 2 )
          call reset_rain_snow_sub_grid()
#endif
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_sensrad
!-----------------------------------------------------------------------
! *** computes atmospheric forcing due to sensible heat and radiation
! *** the forcing terms are computed in dimensional units (K/s)
! *** input
! ***       hflux : sensible heat flux between atmosphere and earth
! ***       hesws : short wave solar radiation in layer 1 or 2
! ***       ulrad : upward long wave radiation in layer 1 or 2
! *** output
! ***       thforg: temperature forcing in layer 1 or 2 in K/s
! ***       vhforg: diabatic forcing in layer 1 or 2 in K/s
!-----------------------------------------------------------------------


      use comatm, only: grav
      use comphys, only: nlon,nlat,cpair,dp0,dp1,dp2,thforg0,thforg1,
     &                   thforg2,vhforg0,vhforg1,vhforg2
      use comsurf_mod, only: hflux,ulrads,dlrads,hesw0,hesw1,hesw2,
     &                       ulrad0,ulrad1,ulrad2


      implicit none


      integer i,j
      real*8  halpha,halpha1,halpha2,sum1,sum2,sum0,halpha0

      halpha=grav/cpair


      halpha0 =halpha/dp0
      halpha1 =halpha/dp1
      halpha2 =halpha/dp2



!
! *** summation of forcing terms
!
      do j=1,nlon
        do i=1,nlat


          sum0 = (hesw0(i,j) - ulrad0(i,j) + ulrad1(i,j))*halpha0
          sum1 = (hesw1(i,j) - ulrad1(i,j) + ulrad2(i,j))*halpha1
          sum2 = (hflux(i,j) + hesw2(i,j) - ulrad2(i,j) +
     *           ulrads(i,j) - dlrads(i,j))*halpha2


          thforg0(i,j) = thforg0(i,j) + sum0
          thforg1(i,j) = thforg1(i,j) + sum1
          thforg2(i,j) = thforg2(i,j) + sum2


          vhforg0(i,j) = vhforg0(i,j) + sum0
          vhforg1(i,j) = vhforg1(i,j) + sum1
          vhforg2(i,j) = vhforg2(i,j) + sum2


        enddo
      enddo

      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_fluxes(nn)
!-----------------------------------------------------------------------
! *** computes energy fluxes above ocean surface
! *** short wave radiation, long wave radiation, sensible heat flux
! *** latent heat flux, evaporation
!-----------------------------------------------------------------------

      use comdiag, only: irad
      use comsurf_mod, only: noc, swrad
      use atmmois_mod, only: ec_surfmois

      implicit none

      integer nn
      logical :: success

      swrad = 0.0  !mohr

      if (irad.eq.1) call ec_swaverad2(nn)
      call ec_swaverad(nn)
      call ec_lwaverad(nn)
      if (irad.eq.1) call ec_lwaverad2(nn)
      call ec_dragcoef(nn)
      success = ec_surfmois(nn)
      call ec_sensibheat(nn)
      call ec_latentheat(nn)
      if (nn.eq.noc) call ec_momentflux(nn)

      return
      end



!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_swaverad(nn)
!-----------------------------------------------------------------------
! *** computes short wave radiation
! *** linearization of RCM with ISCCP D2 1990 clouds
!-----------------------------------------------------------------------


      use comphys, only: irn,bup,dayfr,solarf,tcc,tas1,kosz,costref,
     &                   dso4,sulopt,issulstrt,iscensul,salbref,
     &                   swrref,swrcost,swrsalb,fswusfc,fswdsfc
      use comdiag, only: nlon,nlat,irad
      use comemic_mod, only: imonth, iyear
      use comrunlabel_mod, only: irunlabelf
      use comsurf_mod, only: nse,noc,nld,fractn,albesn,alb2esn,heswsn,
     &                       hesw0n,hesw1n,hesw2n

      use commons_mod, only: fswutoa0,fswutoaGA,fswdtoa0,fswdtoaGA

      implicit none


#if ( CLAQUIN == 1 )

!dmr --- Pour le forcage a la Claquin et al., 2003
#include "radforc.h"
!dmr --- Pour le forcage a la Claquin et al., 2003

#endif

      integer i,j,k,l,ireg
      integer m, d, r, nn , nol


      real*8 f0,f1,ftot(8),fn(8,0:1)
      real*8 drs, drs2, drs3
      real*8 dcost, df,sk,sr,x,y,dfs,smsc
      real*8 fswdtoa,fswutoa!,fswdsfc(nlat,nlon),fswusfc
      real*8 fswutoa2(nlat,nlon)
c~ ,fswutoa0(nlat,nlon)
c~       real*8 fswdtoa0(nlat,nlon)
      real*8 ec_globalmean
      real*8 fswutoaG0
c~ ,fswutoaGA
      real*8 fswdtoa2,fswdtoaG0
c~ ,fswdtoaGA
      real*8 fswutoa_diff,fswutoaG,df_test,fswdtoa_diff,fswdtoaG

      integer nreg(2),indxsul
      real*8 zac(2),asup
!     real*8 zac(2),asup,bup
      common /rad_sul0 /fswutoaG,df_test,fswdtoaG
c~       common /rad_sul1 /fswutoa0,fswutoaGA,fswdtoa0,fswdtoaGA

#if ( CLAQUIN == 1 )
!dmr --- Pour le forcage a la Claquin et al., 2003
      REAL*8 dfprime
!dmr --- Pour le forcage a la Claquin et al., 2003
#endif

! *** aerosol scattering included as a correction on the upward
! *** clear sky fluxes
! *** sk,sr: empirical coefficients Dorland et al, J. Geophys. Res.,102,
! *** 28079-28100, 1997.
! *** smsc: mass scattering coefficient [m2/g]
! *** dso4: change in sulfate aerosol column integrated concentration since
! *** pre-industrial times [g/m2]


      sk=0.058d0*1370d0
      sr=0.05d0
      smsc=8.0
!     write(info_id,*) 'bup=',bup

      if (nn.eq.noc.or.nn.eq.nse) nol=1
      if (nn.eq.nld) nol=2


      do j=1,nlon
       do i=1,nlat
          indxsul=issulstrt+1-1850
          if(indxsul.lt.1) indxsul=1
          indxsul=12*(indxsul-1)+imonth
          if ((iscensul.eq.1).AND.(irunlabelf+iyear.ge.issulstrt))then
            if (iyear.eq.0) then
              indxsul=irunlabelf+1-1850
              if (indxsul.gt.250) indxsul=250
              indxsul=12*(indxsul-1)+1
            else
              indxsul=irunlabelf+iyear-1850
              if (indxsul.gt.250) indxsul=250
              indxsul=12*(indxsul-1)+imonth
            endif
          endif
            tas1(i,j) = sulopt(indxsul,i,j)
            if (irunlabelf+iyear.le.issulstrt) tas1(i,j) =0.
            if ((issulstrt.le.1850).and.(iscensul.eq.0)) tas1(i,j) =0.

!           bup=0.13
            nreg(nol)=irn(i,j,nol)
            zac(nol)=dble(costref(nreg(nol),imonth))
            if (zac(nol).GT.0.) then
             asup = (bup*tas1(i,j)*(1-albesn(i,j,nn))**2)/zac(nol)
            else
             asup =0.
            endif

            alb2esn(i,j,nn) = albesn(i,j,nn)+ asup
!        else
!          alb2esn(i,j,nn) = albesn(i,j,nn)
!        endif
          if (alb2esn(i,j,nn).ge.1.) then
!           write(*,*)alb2esn(i,j,nn),i,j,iyear,imonth,iday,zac(nol),
!    &       asup,albesn(i,j,nn),iscensul,tas1(i,j)
            alb2esn(i,j,nn)=1.
          endif
          df=dayfr(i)*solarf(i)
          ireg=irn(i,j,nol)
          dcost=kosz(i)-costref(ireg,imonth)
          do l=1,8
            do k=0,1
              fn(l,k) =
     &               swrref(l,ireg,imonth,k)
     &               +  swrcost(l,ireg,imonth,k)*dcost
            enddo
          enddo

          x=sqrt(kosz(i))
          y=sqrt(1-alb2esn(i,j,nn))
          dfs=sk*(4d0*x*y*(y-x)-sr)*dso4(i,j)*smsc
!         WRITE(*,*)dso4(i,j)
          if (dfs.gt.0d0.and.kosz(i).lt.0.05) dfs=0d0
          drs=alb2esn(i,j,nn)-salbref(ireg,imonth)
          drs2=drs*drs
          drs3=drs2*drs


          do l=1,4
            f0=fn(l,0)+swrsalb(l,ireg,imonth,0)*drs+dfs
            f1=fn(l,1)+swrsalb(l,ireg,imonth,1)*drs
     &                     +swrsalb(l,ireg,imonth,2)*drs2
     &                     +swrsalb(l,ireg,imonth,3)*drs3
            ftot(l) = (1.-tcc(i,j))*f0 + tcc(i,j)*f1
!           write(*,*)f0,f1
!           if (l.eq.1)ftot_test=f0
          enddo
          do l=5,8
            f0=fn(l,0)+swrsalb(l,ireg,imonth,0)*drs
            f1=fn(l,1)+swrsalb(l,ireg,imonth,1)*drs
     &                     +swrsalb(l,ireg,imonth,2)*drs2
     &                     +swrsalb(l,ireg,imonth,3)*drs3
            ftot(l) = (1.-tcc(i,j))*f0 + tcc(i,j)*f1
          enddo

#if ( CLAQUIN == 1 )
!dmr --- Ajout de l utilisation du forcage des poussieres
!dmr ---
!dmr --- Distinction des cas forcage transitoire ou non

!dmr --- Si transitoire ...
          IF (transitforce.NE.0) THEN

!dmr --- Si l annee est superieure au forcage total
!dmr ---   ---> on conserve la derniere annee !
            IF (iyear.GT.transitforce) THEN
              dfprime=df*(solarc+radforc(i,j)
     >     *facteurpoussieres(transitforce))/solarc

              if
     >        ((imonth.eq.1).and.(iday.eq.1).and.(i.eq.1).and.(j.eq.1))
     >        write(info_id,*) 'facteur poussieres'
     >     , facteurpoussieres(transitforce)

            ELSE
              dfprime=df*(solarc+radforc(i,j)
     >     *facteurpoussieres(iyear))/solarc

              if
     >        ((imonth.eq.1).and.(iday.eq.1).and.(i.eq.1).and.(j.eq.1))
     >        write(info_id,*) 'facteur poussieres'
     >     , facteurpoussieres(iyear)

            ENDIF

!dmr --- Si pas transitoire ...
          ELSE
            dfprime=df*(solarc+radforc(i,j))/solarc
          ENDIF
!dmr ---
#endif

! alternative calculation of upward flux at ground:
! in parameterisation no cross terms are accounted for, which are important for
! upward shortwave radiation at surface and therefore also for net flux
! heswsn(i,j)

          ftot(4)=-alb2esn(i,j,nn)*ftot(8)

#if ( CLAQUIN ==1 )
!dmr --- Ajout de dfprime

          hesw0n(i,j,nn)=(-ftot(1)-ftot(5)+ftot(2)+ftot(6))*dfprime
          hesw1n(i,j,nn)=(-ftot(2)-ftot(6)+ftot(3)+ftot(7))*dfprime
          hesw2n(i,j,nn)=(-ftot(3)-ftot(7)+ftot(4)+ftot(8))*dfprime
          heswsn(i,j,nn)=(-ftot(4)-ftot(8))*dfprime
#else
          hesw0n(i,j,nn)=(-ftot(1)-ftot(5)+ftot(2)+ftot(6))*df
          hesw1n(i,j,nn)=(-ftot(2)-ftot(6)+ftot(3)+ftot(7))*df
          hesw2n(i,j,nn)=(-ftot(3)-ftot(7)+ftot(4)+ftot(8))*df
          heswsn(i,j,nn)=(-ftot(4)-ftot(8))*df
#endif

          if (irad.eq.1) then
! for diagnostic purposes:
! (1) downward shortwave radiation at TOA
             if (nn.eq.1) then 
               fswdtoa0(i,j)= 0.
             endif
#if ( CLAQUIN ==1 )
             fswdtoa2=-ftot(5)*dfprime
#else
             fswdtoa2=-ftot(5)*df
#endif
             fswdtoa0(i,j)=fswdtoa0(i,j)+(fractn(i,j,nn)*fswdtoa2)
! (2) upward shortwave radiation at TOA
             if (nn.eq.1) then 
               fswutoa0(i,j)=0.
               fswusfc(i,j) = 0.
              endif
#if ( CLAQUIN ==1 )
             fswutoa2(i,j)=ftot(1)*dfprime
#else
             fswutoa2(i,j)=ftot(1)*df
#endif
             fswutoa0(i,j)=fswutoa0(i,j)+(fractn(i,j,nn)*fswutoa2(i,j))
! (3) downward shortwave radiation at SURFACE
             fswdsfc(i,j)=-ftot(8)*df
! (4) upward shortwave radiation at SURFACE
! incorrect, ftot(8) does not depend on surfaces ...  fswusfc(i,j)=fswusfc(i,j) + (heswsn(i,j,nn)+ftot(8)*df)

         endif
        enddo
      enddo
      if (irad.eq.1) then
       if (nn.eq.3) then
        fswutoaG0=ec_globalmean(fswutoa0)
        fswdtoaG0=ec_globalmean(fswdtoa0)
        fswutoa_diff=fswutoaG0-fswutoaG
        fswdtoa_diff=fswdtoaG0-fswdtoaG
        if (iyear.eq.0) fswutoaGA=0.
        if (iyear.eq.0) fswdtoaGA=0.
        fswutoaGA=fswutoaGA+(fswutoa_diff/(360.*6.))
        fswdtoaGA=fswdtoaGA+(fswdtoa_diff/(360.*6.))
       endif
      else
       fswutoaGA=0.
       fswdtoaGA=0.
      endif



      return
      end



!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_lwaverad(nn)
!-----------------------------------------------------------------------
! *** computes long wave radiation according to the parameterization of
! *** Chao Chou and Neelin and substantially adapted and extended
! *** for global scale and more
! *** specific ECBILT application by the one and only Michiel Schaeffer
! ***
! *** parameters: nlat   = number of gridpoints in the meridional
! ***                      direction (32)
! ***             nlon   = number of gridpoints in the zonal
! ***                      direction (64)
! ***
! *** input : dtemp(19,nlat,nlon): temperature anomalies [K] wrt ncep
! ***                              climatology tncep in common lwrscheme
! ***         dqa(nlat,nlon) : anomalies of total prec. water cont. below
! ***                          500 hPa wrt ncep climatology
! ***         tcc(nlat,nlon)  : total cloud cover
! ***         ghg(19) : concentrations of well mixed ghg s (see comphys.h)
! ***
! *** output : ulrad1(nlat,nlon): upward longwave radiation [Wm-2] at
! ***                             toa
! ***          ulrad2(nlat,nlon): net longwave radiation [Wm-2] at
! ***                             (500 hPa)
! ***          dlrads(nlat,nlon): downward longwave radiation [Wm-2] at
! ***                             the surface
! ***          ulrads(nlat,nlon): upward longwave radiation [Wm-2] at
! ***                             the surface
!-----------------------------------------------------------------------


      use comphys, only: irn,tcc,ipl,dtemp,sboltz,expir,lwrflux,lwrqa,
     &                   lwrt,lwrts,lwrqts,qancep,tncep
      use comdiag, only: nlon,nlat,irad
      use comemic_mod, only: imonth
      use comsurf_mod, only: nse,noc,nld,tsurfn,fractn,emisn,lwrmois,
     &                       ulradsn,ulrad0n,ulrad1n,ulrad2n,dlradsn
      use comunit, only:

!!    USE OMP_LIB

      implicit none

      integer i,j,l,k,m,is,ism,nol,nn,ireg,h
      real*8  lwr(7,0:1),dumts
      real*8  dqa,dqreg(27)
      real*8  ulrad0nm,ulrad1nm,ulrad2nm,ulradsnm,dlradsnm
      real*8  ec_globalmean
      real*8  ulrad0nU,ulrad1nU,ulrad2nU,ulradsnU,dlradsnU
      real*8  ulrad0nz(nlat,nlon),ulrad1nz(nlat,nlon)
      real*8  ulrad2nz(nlat,nlon),ulradsnz(nlat,nlon)
      real*8  dlradsnz(nlat,nlon), ulrad0nUz(nlat,nlon)
      real*8  ulrad1nUz(nlat,nlon)
!dmr [OMP]
      real*8  mompd,mompd2
!dmr [OMP]
      common / radO3 / ulrad0nU,ulrad1nU,ulrad2nU,ulradsnU,dlradsnU
      common / radO32 / ulrad0nUz,ulrad1nUz



      is=imonth/3+1
      if (is.gt.4) is=1
      ism=(is-1)*3+1

      do i=1,27
!dqa    dqreg(i)=qancep(i,ism)**0.3333
        dqreg(i)=qancep(i,ism)
      enddo

      if (nn.eq.noc.or.nn.eq.nse) nol=1
      if (nn.eq.nld) nol=2

c~ !$OMP PARALLEL
c~ !$OMP DO PRIVATE (i,j,l,k,m,ireg,h,dqa,lwr,dumts,mompd,mompd2) SCHEDULE(static)
      do j=1,nlon
        do i=1,nlat
          ireg=irn(i,j,nol)

!-hemispheric dependence of tropospheric ozone forcing
          if (i.le.16) then
           h=1
          else
           h=2
          endif

!dqa      dqa=lwrmois(i,j)-dqreg(ireg)
!dqa      q**1/3-qm**1/3=qm**(1/3-n)*(q**n-qm**n)
          dqa=dqreg(ireg)**(0.3333-EXPIR)*
     *          (lwrmois(i,j)**EXPIR-dqreg(ireg)**EXPIR)
!dqa      write(info_id,'(i3,3F12.5)') i,lwrmois(i,j)-
!dqa &      dqreg(ireg),dqa,lwrmois(i,j)**0.333-dqreg(ireg)**0.33333
          do l=0,1
            do k=1,7
              lwr(k,l)=lwrflux(k,ireg,is,l,h)+lwrqa(k,ireg,is,l)*dqa
              do m=1,ipl(ireg)-1
                lwr(k,l)=lwr(k,l)+
     *                   lwrt(k,m,ireg,is,l)*dtemp(m,i,j,nol)
              enddo
              lwr(k,l)=lwr(k,l)+
     *              lwrt(k,18,ireg,is,l)*dtemp(18,i,j,nol)
            enddo

! afq --    dumts=tsurfn(i,j,nn)-tncep(19,ireg,ism) ! moved above !
            dumts=tsurfn(i,j,nn)-tncep(19,ireg,ism) !crash if not here?
            do m=1,4
              do k=1,3
                lwr(k,l)=lwr(k,l)+
     *          (lwrts(k,m,ireg,is,l)+lwrqts(k,m,ireg,is,l)*dqa)
     *          *dumts
              enddo
              lwr(7,l)=lwr(7,l)+
     *        (lwrts(7,m,ireg,is,l)+lwrqts(7,m,ireg,is,l)*dqa)
     *        *dumts
              dumts=dumts*(tsurfn(i,j,nn)-tncep(19,ireg,ism))
            enddo


          enddo

          mompd=(lwr(1,0)*(1-tcc(i,j))+lwr(1,1)*tcc(i,j))
!dmr [OMP]          ulrad0n(i,j,nn)=(lwr(1,0)*(1-tcc(i,j))+lwr(1,1)*tcc(i,j))
          ulrad0n(i,j,nn)=mompd
          mompd2=(lwr(2,0)+lwr(5,0))*(1-tcc(i,j)) +
     *             (lwr(2,1)+lwr(5,1))*tcc(i,j)
!dmr [OMP]          ulrad1n(i,j,nn)=(lwr(2,0)+lwr(5,0))*(1-tcc(i,j)) +
!dmr [OMP]     *             (lwr(2,1)+lwr(5,1))*tcc(i,j)
          ulrad1n(i,j,nn)=mompd2
          ulrad2n(i,j,nn)=(lwr(3,0)+lwr(6,0))*(1-tcc(i,j)) +
     *             (lwr(3,1)+lwr(6,1))*tcc(i,j)

          ulradsn(i,j,nn)=emisn(nn)*sboltz*tsurfn(i,j,nn)**4
          dlradsn(i,j,nn)=-lwr(7,0)*(1-tcc(i,j))-lwr(7,1)*tcc(i,j)


         if (irad.eq.1) then
           if(nn.eq.1) then
!dmr [OMP]            ulrad0nUz(i,j)=0.
!dmr [OMP]            ulrad1nUz(i,j)=0.
            ulrad0nUz(i,j)=mompd*fractn(i,j,nn)
            ulrad1nUz(i,j)=mompd2*fractn(i,j,nn)
           else
            ulrad0nUz(i,j)=ulrad0nUz(i,j)+(mompd*fractn(i,j,nn))
            ulrad1nUz(i,j)=ulrad1nUz(i,j)+(mompd2*fractn(i,j,nn))
           endif

          ulrad0nUz(i,j)=ulrad0nUz(i,j)+(ulrad0n(i,j,nn)*fractn(i,j,nn))
          ulrad1nUz(i,j)=ulrad1nUz(i,j)+(ulrad1n(i,j,nn)*fractn(i,j,nn))

         endif
       enddo
      enddo
c~ !$OMP END DO
c~ !$OMP END PARALLEL

      if (irad.eq.1) then
       ulrad0nU=ec_globalmean(ulrad0nUz)
       ulrad1nU=ec_globalmean(ulrad1nUz)
      endif

!     ulrad0nm=ec_globalmean(ulrad0n)
!     ulrad1nm=ec_globalmean(ulrad1n)
!     ulrad2nm=ec_globalmean(ulrad2n)
!     ulradsnm=ec_globalmean(ulradsn)
!     dlradsnm=ec_globalmean(dlradsn)

!     write(iuo+37,*)ulrad0nm,ulrad1nm,ulrad2nm,ulradsnm,dlradsnm

! *** that s all folks


      return
      end


!2345678901234567890123456789012345678901234567890123456789012345678901
      SUBROUTINE ec_dragcoef(nn)
!----------------------------------------------------------------------
! *** drag coefficient
! *** depending on stability
!----------------------------------------------------------------------


      USE comphys, only: nlon,nlat,cdrag,cdragw,cwdrag
      use comsurf_mod, only: tsurfn,tempsgn,cdragvn


      implicit none


      integer i,j,nn


      real*8 cdrags,cdragl,tdif,cred,cdum

      cdrags=cdrag
      cdragl=cdrag
      cred=0.2
      cdum=(1-cred)*0.5
      do j=1,nlon
        do i=1,nlat
          cdragw(i,j)=cwdrag

          tdif=tsurfn(i,j,nn)-tempsgn(i,j,nn)
          cdragvn(i,j,nn)=cdrags*
     &                max(cred,cred+min(cdum+cdum*tdif,1-cred))

        enddo
      enddo


      return
      end


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
! dmr   subroutine ec_sensibheat
! ***    Computes the sensible heatflux exchanged between surface and atmosphere
! dmr& aurel
!      Update to include the vertical downscaling at nb_ipoints resolution
!      Last updated: 2016-01-25
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

      subroutine ec_sensibheat(nn)

      use comphys, only: nlon,nlat,uv10,alphad
      use comsurf_mod, only: nld,tsurfn,tempsgn,hfluxn,cdragvn,hficof


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
! dmr   Downscaling input variables for the sub-grid computations
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
#if ( DOWNSTS == 1 )
      use vertDownsc_mod       , only: tempsg_d, tsurfn_d, hfluxn_d
      use ecbilt_topography    , only: nb_levls
#endif
      implicit none


      integer :: i,j,nn  ! common indices: nlat, nlon and ntyps
      integer :: nb_down ! loop indice over nb_ipoints

! *** sensible heatflux (watt/m**2)

      if (nn.eq.2) then ! [TOBEDONE] => immediate indice number 2 is not a clean way to reference
                        !               the surface type nse = sea-ice
       do j=1,nlon
        do i=1,nlat
          hficof(i,j)=alphad*cdragvn(i,j,2)*uv10(i,j)
        enddo
       enddo
       do j=1,nlon
        do i=1,nlat
          hfluxn(i,j,2)=hficof(i,j)*(tsurfn(i,j,2)-tempsgn(i,j,2))
        enddo
       enddo
      else
       do j=1,nlon
        do i=1,nlat
           hfluxn(i,j,nn)=alphad*cdragvn(i,j,nn)*uv10(i,j)*
     &                      (tsurfn(i,j,nn)-tempsgn(i,j,nn))

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

#if ( DOWNSTS == 1 )

          if (nn.eq.nld) then

!     --- The other two surface types (sea-ice and ocean) do not admit sub-grid orography.

                do nb_down = 1, nb_levls
                   hfluxn_d(i,j,nn,nb_down)=alphad*cdragvn(i,j,nn)*
     &       uv10(i,j)* (tsurfn_d(i,j,nn,nb_down)-tempsg_d(i,j,nb_down))
                enddo ! on nb_ipoints


          endif ! on ntyps

#endif
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

        enddo   ! on nlat
       enddo    ! on nlon
      endif     ! on nn .eq. nse

      return

      end subroutine ec_sensibheat
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
! dmr  ec_sensibheat <END>
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
! dmr   subroutine ec_latentheat
! ***    Computes the latent heatflux exchange between surface and atmosphere
! dmr& aurel
!      Update to include the vertical downscaling at nb_ipoints resolution
!      Last updated: 2016-01-25
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

      subroutine ec_latentheat(nn)


      USE comatm, only: dtime,iwater
      USE comphys, only: nlon,nlat,rowat,uv10,rlatsub,evfac,alphav,alphas,
     &                   rlatvap
      use comsurf_mod, only: epss,nse,noc,nld,q10n,qsurfn,fractn,efluxn,evapn,
     &                       evfacan,adsnow,abmoisg,abmoism,cdragvn


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
! dmr   Downscaling input variables for the sub-grid computations
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

#if ( DOWNSTS == 1 )
      use vertDownsc_mod,   only: efluxn_d, q10n_d, qsurfn_d
      use ecbilt_topography,only: nb_levls
#endif

#if ( WISOATM == 1 ) 
      use comatm, only: iwat17, iwat18, iwat2h
      use iso_param_mod,    only: ratio_oceanatm, ratio_evap_ocean, deltaR
     >                          , dexcess
      use iso_alphas,       only: alpha_lv17, alpha_lv18, alpha_lv2h
      use comsurf_mod, only: fractoc, tsurfn
#endif

c~ #if (ISOATM >= 1 )
c~       USE iso_param_mod, ONLY: ieau, neauiso, ieau18,ieau17,ieaud,
c~      & alpha_diff_oc, rsmow, dexcess, delta
c~       USE isoatm_mod, ONLY: ratio_oceanatm, ratio_evap,datmini
c~       USE iso_alphas, ONLY: alpha_lv
c~       use iso_funcs,  only: mois2ratiosmj
c~ #endif
c~ #if ( ISOATM >= 2 )
c~       USE isoatm_mod, ONLY: ratio_evap_snow, ratio_evap_land
c~ #endif

#if ( IMSK == 1 )
      USE input_icemask, ONLY: icemask
#endif

      implicit none

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
! dmr   Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

      integer :: i,j,nn,n,NSTAT

#if ( DOWNSTS == 1 )
      integer :: nb_down     ! loop indice over nb_ipoints
#endif


#if ( WISOATM == 1 )
      integer                            :: quick_fix
#endif

c~ #if ( ISOATM >= 1 )
c~       integer                            :: k, quick_fix
c~       real, dimension(nlat,nlon,neauiso) :: ration_qatm
c~       real, dimension(nlat,nlon)         :: relhumtilde
c~       real, dimension(neauiso)           :: mean_evapiso
c~       real                               :: mean_evap, nb_cazes
c~       real, dimension(neauiso)           :: variso
c~       character(len=10)                  :: varisonm
c~       real                               :: esnowiso, emoisiso
c~ #endif

      double precision :: ec_qsat,db,emois,esubf,evapf,esnow,sfrac,edum
     &  ,psilai

#if ( IMSK == 1 )
      double precision :: esice
#endif

c~ #if ( ISOATM >= 2 )
c~       double precision    rainf(nlat,nlon,5),snowf(nlat,nlon,5)
c~ #else
c~       double precision    rainf(nlat,nlon),snowf(nlat,nlon)
c~ #endif
!      double precision    fswdsfc(nlat,nlon)
      double precision    dc,eflux_t,eflux_g,eflux_bare

#if ( EVAPTRS == 0 )
      double precision :: lai(2), resist(3), k0(2), rs
      double precision    fswdsfcM
#endif

      double precision st,sg,sd,snlt,anup,blai,pnpp,b12,b34,b1,b2,b3,b4,
     &       anup_moy,stock,st_moy
      double precision stR,sgR,sdR,snltR

!      common /pr_evap /fswdsfc
      common /BIOTA/
     &   ST(nlat,nlon), SG(nlat,nlon), SD(nlat,nlon), SNLT(nlat,nlon),
     &   BLAI(nlat,nlon,2), PNPP(nlat,nlon),
     &   B12(nlat,nlon),   B34(nlat,nlon),
     &   B1(nlat,nlon), B2(nlat,nlon), B3(nlat,nlon), B4(nlat,nlon),
     &   ANUP_MOY(nlat,nlon),ANUP(nlat,nlon), STOCK(nlat,nlon),
     &   st_moy(nlat,nlon),
     &   NSTAT(nlat,nlon),STR(nlat,nlon), SGR(nlat,nlon),SDR(nlat,nlon),
     &   SNLTR(nlat,nlon)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
! dmr   Start of main computation
! *** latent heatflux due to evaporation from surface (watt/m**2)
! *** and evaporation rate (m/s)
! *** limit evaporation to available snow or bottom moisture
! *** evaporation factor =1 over snow, over wet land maximal 1
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

#if ( EVAPSI == 1 )
! NO EVAPORATION ON SEA-ICE !!
      if (nn.EQ.noc) then
#else
      if (nn.ne.nld) then
#endif
!dmr --- Cas où l on est un ocean ou de la banquise ...
#if ( ISOATM == 2 )

!dmr --- dummy variables to compute debugging diagnostics

       mean_evapiso = 0.0d0
       mean_evap = 0.0d0
       nb_cazes = 0.0d0
#endif
        do j=1,nlon
          do i=1,nlat
            evfacan(i,j,nn)=evfac

            efluxn(i,j,nn)=alphav*cdragvn(i,j,nn)*uv10(i,j)*
     &                 (qsurfn(i,j,nn)-q10n(i,j,nn))
            efluxn(i,j,nn)=evfacan(i,j,nn)*max(0.d0,efluxn(i,j,nn))
c~ #if ( ISOATM >= 1 )
c~             evapn(i,j,nn,ieau)=efluxn(i,j,nn)/(rowat*rlatvap)
c~ #else
            evapn(i,j,nn,iwater)=efluxn(i,j,nn)/(rowat*rlatvap)
c~ #endif
          enddo
        enddo

#if ( WISOATM == 1 )
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
!      Computation of evaporation for the water isotopes over oceans
!     == Version without isotopic fractionnation for now: only equilibrium 
!          with surface oceanic values ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
       do j=1, nlon
        do i=1, nlat

!dmr --- Quick and dirty fix to have consistency between fractoc and
!dmr ---   nse/noc
         quick_fix = CEILING(fractoc(i,j))

         ratio_evap_ocean(i,j,iwat17) = ratio_oceanatm(i,j,iwat17)
     >                                / alpha_lv17(tsurfn(i,j,nn))

         ratio_evap_ocean(i,j,iwat18) = ratio_oceanatm(i,j,iwat18)
     >                                / alpha_lv18(tsurfn(i,j,nn))
         ratio_evap_ocean(i,j,iwat2h) = ratio_oceanatm(i,j,iwat2h)
     >                                / alpha_lv2h(tsurfn(i,j,nn))         
         evapn(i,j,nn,iwat17:iwat2h) = evapn(i,j,nn,iwater) 
     >                    *ratio_evap_ocean(i,j,iwat17:iwat2h)*quick_fix


         if (( quick_fix.gt.0 ) .and.( tsurfn(i,j,nn).gt.280.0d0 )) then
         WRITE(*,*) "Test iso evap_ocean", 
     >      deltaR(ratio_evap_ocean(i,j,iwat18),iwat18),
     >      deltaR(ratio_evap_ocean(i,j,iwat2h),iwat2h),
     >      deltaR(ratio_oceanatm(i,j,iwat18),iwat18),
     >      deltaR(ratio_oceanatm(i,j,iwat2h),iwat2h),     
     >      dexcess(ratio_evap_ocean(i,j,iwat2h),ratio_evap_ocean(i,j,iwat18)),
     >      dexcess(ratio_oceanatm(i,j,iwat2h),ratio_oceanatm(i,j,iwat18)),
     >      tsurfn(i,j,nn)
         endif

        enddo ! on spatial loop - nlon
       enddo ! on spatial loop - nlat

#endif

      else

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
!       Evaporating from the continent surface
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

        do j=1,nlon
          do i=1,nlat

! *** evaporation factor =1 over snow, over wet land maximal 1
! *** care has to be taken in case of a snowcover in order to
! *** conserve heat: the evaporation is constraint to the amount
! *** of snow and moisture available; it can happen that in
! *** one timestep, the remaining snow is sublimated and part
! *** of the latent heat flux is also used to evaporate bottom
! *** moisture: Eflux = (1-sfrac)*Esub + sfrac*Evap


            if (fractn(i,j,nld).gt.epss) then

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Evaporating from a snowy surface
!-----|--1--------2---------3---------4---------5---------6---------7-|
c~ #if ( ISOATM >= 2 )
c~               if (adsnow(i,j,ieau).gt.0.) then
c~ #else
              if (adsnow(i,j,iwater).gt.0.) then
c~ #endif

                evfacan(i,j,nld)=evfac
                edum=cdragvn(i,j,nld)*uv10(i,j)*
     &                   (qsurfn(i,j,nld)-q10n(i,j,nld))
                edum=evfacan(i,j,nld)*max(edum,0d0)
!-----|--1--------2---------3---------4---------5---------6---------7-|
!       esubf == energy available for sublimation
!-----|--1--------2---------3---------4---------5---------6---------7-|
                esubf=alphas*edum
                evapf=alphav*edum
c~ #if ( ISOATM >= 2 )
c~                 esnow=min(rowat*adsnow(i,j,ieau)*rlatsub/dtime,
c~      &                  esubf)
c~ #else
                esnow=min(rowat*adsnow(i,j,iwater)*rlatsub/dtime,
     &                  esubf)
c~ #endif
!-----|--1--------2---------3---------4---------5---------6---------7-|
!       esnow == energy needed to evaporate all snow in adsnow
!       emois == same, from moisture layer ...
!       Here, case where not enough snow for accomodating energy flux
!         ==> part of the flux is used to evaporate the moisture in
!          abmoisg, but limited to its maximum amount
!-----|--1--------2---------3---------4---------5---------6---------7-|
                if (esnow.lt.esubf) then ! all possible snow is sublimated!
                  sfrac=(esubf-esnow)/esubf
c~ #if ( ISOATM >= 2 )
c~ ! [TODO] == update this !!!!
c~ ! ## TAG FIX AURELIEN ##
c~                 if (abmoisg(i,j,ieau).LT.0d0) then
c~ !                   WRITE(*,*) "fixed abmoisg, i,j, in l.1930 atmphys0.f"
c~ !     & ,i,j
c~                   abmoisg(i,j,ieau) = 0d0
c~                 endif
c~ ! ## TAG FIX AURELIEN ##
c~ ! [TODO] == update this !!!!
c~                 emois=min(rowat*abmoisg(i,j,ieau)*rlatvap/dtime
c~      &                        ,sfrac*evapf)
c~ #else
                emois=min(rowat*abmoisg(i,j,iwater)*rlatvap/dtime,sfrac*evapf)
c~ #endif

#if ( IMSK == 1 )

                esice = 0.0d0

!-- dmr Case where there is an ice-sheet. So whatever the depth of the snow layer
!         we can use of the energy to melt either snow, or ice below.
!       Hence we do not have any energy left for moisture. Set emois to zero.
!       Store the energy left in esice to provide for efluxn
!       Do not change the amount of evaporation (potentially sublimation of ice?).
! -- dmr
                if (icemask(i,j).gt.0.9) then
                    esice = esubf-esnow
                    emois = 0.0
                endif
#endif

                  efluxn(i,j,nld)=esnow+emois
#if ( IMSK == 1 )
                  efluxn(i,j,nld)=efluxn(i,j,nld)+esice
#endif

c~ #if ( ISOATM >= 1 )
c~                   evapn(i,j,nld,ieau)=esnow/(rowat*rlatsub)+
c~      &                                           emois/(rowat*rlatvap)

c~ #else
                  evapn(i,j,nld,iwater)=esnow/(rowat*rlatsub)+
     &                                           emois/(rowat*rlatvap)
c~ #endif
#if ( ISOATM == 1 )
!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Computation of evaporation for the water isotopes
!        == Version without isotopic fractionnation ==
!      (Evaporation from the land / snow surface)
!-----|--1--------2---------3---------4---------5---------6---------7-|
!dmr Compute the evaporative ratio
         do k=ieau+1,neauiso
           ratio_evap(i,j,k) = (datmini(k)+1.0d0) * rsmow(k)
         enddo

         do k=ieau+1, neauiso
!dmr apply the evaporative ratios to the evaporative flux
           evapn(i,j,nld,k) = evapn(i,j,nld,ieau) * ratio_evap(i,j,k)
         enddo
#elif ( ISOATM == 2 )
!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Computation of sublimation for the water isotopes
!     == Version with isotopic fractionnation as in MJ79==
!       ... first version with equilibrium fractionnation, no kinetic
!      This is the case where all snow is "sublimated"
!      And possibly part of the land moisture evaporated
!-----|--1--------2---------3---------4---------5---------6---------7-|
       DO k=ieau+2, neauiso !  beware in ISOATM == 2, alpha_lv is undefined for ieau16

!dmr --- Here all snow is evaporated, no complicated scheme for the
!         ratio of isotopes in snow: take all isotopes!
!         ratio_evap_snow(i,j,k) = adsnow(i,j,k)/adsnow(i,j,ieau)
!
!         esnowiso = min(
!     &  esnow/(rowat*rlatsub) * ratio_evap_snow(i,j,k)
!     & ,(adsnow(i,j,k)/dtime)
!     &                 )
         esnowiso = (adsnow(i,j,k)/dtime)
!!
         if (abmoisg(i,j,k).LT.0.0d0) then
           emoisiso = 0.0d0
!           WRITE(*,*) "WARN @@@ : ", abmoisg(i,j,k), i,j,k
         else
         if (emois.GT.0.d0) then
          ratio_evap_land(i,j,k) = abmoisg(i,j,k)/abmoisg(i,j,ieau)
! temporarily no knietic fractionnation for moisture evap from land
!     &           /alpha_lv(tsurfn(i,j,nld),k)

         emoisiso = min(
     &  emois/(rowat*rlatvap)*ratio_evap_land(i,j,k)
     & ,(abmoisg(i,j,k)/dtime)
     &                 )
         else
           emoisiso = 0.0d0
         endif

        endif ! on abmoisg < 0

!      dmr here we take all possible snow and part of the moisture ...
!!         evapn(i,j,nld,k) =
!!     &      esnow/(rowat*rlatsub) * ratio_evap_snow(i,j,k)
!!     &    +
!!     &      emois/(rowat*rlatvap)*ratio_evap_land(i,j,k)

!dmr --- Beware, I did not compute esnowiso as esnow in unit, hence the
!        different formulation hereafter
           evapn(i,j,nld,k) = esnowiso + emoisiso

       ENDDO ! on isotopes

#elif ( ISOATM == 3 )
        WRITE(*,*) "Option non implementee !!! , atmphys0.f"
#endif
!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Case where there is some snow left after sublimation
!-----|--1--------2---------3---------4---------5---------6---------7-|
                else ! case where (esnow.lt.esubf) is not true

                  efluxn(i,j,nld)=esubf

c~ #if ( ISOATM >= 1 )
c~                   evapn(i,j,nld,ieau)=esubf/(rowat*rlatsub)
c~ #else
                  evapn(i,j,nld,iwater)=esubf/(rowat*rlatsub)
c~ #endif

#if ( ISOATM == 1 )

!dmr Compute the evaporative ratio
       DO k=ieau+1,neauiso
         ratio_evap(i,j,k) = (datmini(k)+1.0d0) * rsmow(k)
       ENDDO

       DO k=ieau+1, neauiso
!dmr Apply the evaporative ratios to the evaporative flux
         evapn(i,j,nld,k) = evapn(i,j,nld,ieau) * ratio_evap(i,j,k)
       ENDDO

#elif ( ISOATM == 2 )
!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Computation of sublimation for the water isotopes
!     == Version with isotopic fractionnation as in MJ79==
!       ... first version with equilibrium fractionnation, no kinetic
!      This is the case where there is some snow left after sublimation
!-----|--1--------2---------3---------4---------5---------6---------7-|
        DO k=ieau+2, neauiso !  beware in ISOATM == 2, alpha_lv is undefined for ieau16
         ratio_evap_snow(i,j,k) = adsnow(i,j,k)/adsnow(i,j,ieau)
! temporarily no knietic fractionnation for snow evap from land
! should it be some
!     &           /alpha_sv(tsurfn(i,j,nld),k)

        if ((evapn(i,j,nld,ieau) * ratio_evap_snow(i,j,k)).GT.
     > (adsnow(i,j,k)/dtime)) then
          WRITE(*,*) 'Error: non conservative issue in evapn from dsnow'
          WRITE(*,*) 'atmphys0.f, line 2095'
          WRITE(*,*) '(1)', i,j,k
          WRITE(*,*) '(2)', evapn(i,j,nld,ieau), adsnow(i,j,ieau)/dtime
          WRITE(*,*) '(3)', adsnow(i,j,k)
!          READ(*,*)
        endif

        evapn(i,j,nld,k) = evapn(i,j,nld,ieau) * ratio_evap_snow(i,j,k)
#if ( 0 )
        if ((k.EQ.ieau18).and.(adsnow(i,j,ieau).GT.0.0d0)) THEN
          WRITE(*,*) "Test evap snow: some snow left!"
     &    ,(ratio_evap_snow(i,j,k)/rsmow(k)-1.0d0)*1000.0d0
        ENDIF
#endif

       ENDDO

#elif ( ISOATM == 3 )
        WRITE(*,*) "Option non implementee !!! , atmphys0.f"
#endif
                endif

!-----|--1--------2---------3---------4---------5---------6---------7-|
!        Evaporating from a non - snowy surface
!-----|--1--------2---------3---------4---------5---------6---------7-|
              else
#if ( EVAPTRS == 1 )

c~ #if ( ISOATM >= 2 )
c~                 evfacan(i,j,nld)=evfac*min(1d0,abmoisg(i,j,ieau)
c~      > /max(1.0E-10,abmoism(i,j)))
c~ #else
                evfacan(i,j,nld)=evfac*min(1d0,abmoisg(i,j,iwater)
     > /max(1.0E-10,abmoism(i,j)))
c~ #endif
                efluxn(i,j,nld)=alphav*cdragvn(i,j,nld)*uv10(i,j)*
     &                   (qsurfn(i,j,nld)-q10n(i,j,nld))
#elif ( EVAPTRS == 0 )
                evfacan(i,j,nld)=evfac*min(1d0,abmoisg(i,j,iwater)
     > /max(1.0E-10,abmoism(i,j)))
                dc=ec_qsat(pgroundn(i,j,nld),tempsgn(i,j,nld))
     &                         -q10n(i,j,nld)
                psilai=0.2+(0.08*min(max
     &               (tempsgn(i,j,nld)-tzero,0.),10.))
                lai(1)=6*psilai
                lai(2)=2*psilai
                k0(1)=30.0E-5
                k0(2)=25.0E-5
                fswdsfcM=max(fswdsfc(i,j),1.0d0)
                do n=1,2
                rs=((fswdsfcM+125)/fswdsfcM)*
     &                  (23E-3+(1.5*dc))*(1/(lai(n)*k0(n)))
                rs=max(0d0,rs)
                resist(n)=rs+(1/(cdragvn(i,j,nld)*uv10(i,j)))
                resist(n)=1/resist(n)
                enddo
                resist(3)=1/(cdragvn(i,j,nld)*uv10(i,j))
                resist(3)=1/resist(3)
!	write(*,*)cdragvn(i,j,nld)*uv10(i,j), rs(i,j)

!               efluxn(i,j,nld)=alphav*(qsurfn(i,j,nld)-q10n(i,j,nld))*
!    &              ((st(i,j)*resist(1))+(sg(i,j)*resist(2))+
!    &              (sd(i,j)*resist(3)))
                eflux_bare=alphav*(qsurfn(i,j,nld)-q10n(i,j,nld))
     &           *cdragvn(i,j,nld)*uv10(i,j)*(10./30.)
                eflux_g=alphav*(qsurfn(i,j,nld)-q10n(i,j,nld))
     &           *resist(2)*(15./30.)
                eflux_t=alphav*(qsurfn(i,j,nld)-q10n(i,j,nld))
     &           *resist(1)*(30./30.)
                efluxn(i,j,nld)=eflux_bare+(sg(i,j)*eflux_g)
     &          + (st(i,j)*eflux_t)
#endif

                efluxn(i,j,nld)=evfacan(i,j,nld)*
     &                          max(efluxn(i,j,nld),0d0)
c~ #if ( ISOATM >= 2 )
c~                efluxn(i,j,nld)=min(abmoisg(i,j,ieau)*rowat*rlatvap/dtime
c~      &                       ,efluxn(i,j,nld))
c~ #else
            efluxn(i,j,nld)=min(abmoisg(i,j,iwater)*rowat*rlatvap/dtime
     &                       ,efluxn(i,j,nld))
c~ #endif
c~ #if ( ISOATM >= 1 )
c~                 evapn(i,j,nld,ieau)=efluxn(i,j,nld)/(rowat*rlatvap)
c~ #else
            evapn(i,j,nld,iwater)=efluxn(i,j,nld)/(rowat*rlatvap)
c~ #endif

#if ( ISOATM == 1 )

!dmr Compute the evaporative ratio
       DO k=ieau+1, neauiso
         ratio_evap(i,j,k) = (datmini(k)+1.0d0) * rsmow(k)
       ENDDO

       DO k=ieau+1, neauiso
!dmr Apply the evaporative ratios to the evaporative flux
         evapn(i,j,nld,k) = evapn(i,j,nld,ieau) * ratio_evap(i,j,k)
       ENDDO

#elif ( ISOATM == 2 )
!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Computation of evaporation for the water isotopes
!       ... case where all water is evaporated from "land water"
!     == Version with isotopic fractionnation as in MJ79==
!       ... first version with equilibrium fractionnation, no kinetic
!-----|--1--------2---------3---------4---------5---------6---------7-|
            DO k=ieau+2, neauiso !  beware in ISOATM == 2, alpha_lv is undefined for ieau16


            if (abmoisg(i,j,k).LT.0.0d0) then
              emoisiso = 0.0d0
            else
            if (efluxn(i,j,nld).GT.0.d0) then
            ratio_evap_land(i,j,k) = abmoisg(i,j,k)/abmoisg(i,j,ieau)
! temporarily no kinetic fractionnation for moisture evap from land
!     &           /alpha_lv(tsurfn(i,j,nld),k)

            emoisiso = min(
     &     efluxn(i,j,nld)/(rowat*rlatvap)*ratio_evap_land(i,j,k)
     &    ,(abmoisg(i,j,k)/dtime)
     &                 )
            else
              emoisiso = 0.0d0
            endif
           endif ! on abmoisg < 0

           evapn(i,j,nld,k) = emoisiso
           ENDDO

#elif ( ISOATM == 3 )
        WRITE(*,*) "Option non implementee !!! , atmphys0.f"
#endif
              endif

            endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
! dmr   Downscaling on nb_ipoints level
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

#if ( DOWNSTS == 1 )

!     --- The other two surface types (sea-ice and ocean) do not admit sub-grid orography.
! dmr --- Moreover, flag_down should only exist if ntyp == nland ...
! dmr --- This section of code is called only if nn == nld anyhow

            if (nn.eq.nld) then ! in principle the "nn.eq.nld" is useless here

             do nb_down = 1, nb_levls

c~ #if ( ISOATM >= 2 )
c~               if (adsnow(i,j,ieau).gt.0.) then ! [TODO] adsnow could be taken on the high-res
c~                                                !            land surface if to make sense with
c~                                                !            the rest of the downscaling!!
c~                                                ! 2016-01-25
c~ #else
              if (adsnow(i,j,iwater).gt.0.) then
c~ #endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

!                evfacan(i,j,nld)=evfac

                edum=cdragvn(i,j,nn)*uv10(i,j)*
     &                 (qsurfn_d(i,j,nn,nb_down)-q10n_d(i,j,nn,nb_down))
                edum=evfacan(i,j,nn)*max(edum,0d0)

                esubf=alphas*edum
                evapf=alphav*edum

c~ #if ( ISOATM >= 2 )
c~                 esnow=min(rowat*adsnow(i,j,ieau)*rlatsub/dtime,
c~      &                  esubf)
c~ #else
                esnow=min(rowat*adsnow(i,j,iwater)*rlatsub/dtime,
     &                  esubf)
c~ #endif
                if (esnow.lt.esubf) then
                  sfrac=(esubf-esnow)/esubf
c~ #if ( ISOATM >= 2 )
c~                   emois=min(rowat*abmoisg(i,j,ieau)*rlatvap/dtime,
c~      &                   sfrac*evapf)
c~ #else
                  emois=min(rowat*abmoisg(i,j,iwater)*rlatvap/dtime,
     &                   sfrac*evapf)
c~ #endif
                  efluxn_d(i,j,nn,nb_down)=esnow+emois
!                  evapn(i,j,nld)=esnow/(rowat*rlatsub)+
!     &                                           emois/(rowat*rlatvap)
                else
                  efluxn_d(i,j,nn,nb_down)=esubf
!                  evapn(i,j,nld)=esubf/(rowat*rlatsub)
                endif

              else ! adsnow < 0.0d0

!                evfacan(i,j,nld)=evfac*min(1d0,abmoisg(i,j)/abmoism(i,j))
                efluxn_d(i,j,nn,nb_down)=alphav*cdragvn(i,j,nn)
     &      *uv10(i,j)*(qsurfn_d(i,j,nn,nb_down)-q10n_d(i,j,nn,nb_down))


                efluxn_d(i,j,nn,nb_down)=evfacan(i,j,nn)*
     &                          max(efluxn_d(i,j,nn,nb_down),0d0)
c~ #if ( ISOATM >= 2 )
c~                 efluxn_d(i,j,nn,nb_down)=min(abmoisg(i,j,ieau)
c~      &                   *rowat*rlatvap/dtime
c~      &                       ,efluxn_d(i,j,nn,nb_down))
c~ #else
                efluxn_d(i,j,nn,nb_down)=min(abmoisg(i,j,iwater)
     &                   *rowat*rlatvap/dtime
     &                       ,efluxn_d(i,j,nn,nb_down))
c~ #endif
!     evapn(i,j,nld)=efluxn(i,j,nld)/(rowat*rlatvap)

             endif              ! adsnow

             enddo    ! on nb_levls

! [DELETE]          enddo

           endif ! nn == nld
#endif
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!


          enddo
        enddo
      endif

#if ( EVAPSI == 1 )

! NO EVAPORATION ON SEA-ICE !!

c~ #if ( ISOATM >= 1 )
c~       if (nn.EQ.nse) then
c~        evapn(:,:,nse,:) = 0.0d0
c~       endif
c~ #else
      if (nn.EQ.nse) then
       evapn(:,:,nse,:) = 0.0d0
      endif
c~ #endif

#endif

      return
      end subroutine ec_latentheat
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
! dmr  ec_latentheat <END>
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!


!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_momentflux(nn)
!-----------------------------------------------------------------------
! *** computation of windstress
!-----------------------------------------------------------------------



      USE comdyn, only: utot,vtot
      USE comphys, only: nlon,nlat,roair,uv10rws,uvw10,dragane,cdragw
      use comsurf_mod, only: winstu,winstv


      implicit none


      integer i,j,nn
      real*8  uv,costt,sintt,facstr

      facstr=roair * uv10rws

      do j=1,nlon
        do i=1,nlat
          costt=cos(dragane(i))
          sintt=sin(dragane(i))
          winstu(i,j)=cdragw(i,j)*facstr*uvw10(i,j)*
     *                  (utot(i,j,3)*costt-vtot(i,j,3)*sintt)
          winstv(i,j)=cdragw(i,j)*facstr*uvw10(i,j)*
     *                  (utot(i,j,3)*sintt+vtot(i,j,3)*costt)
        enddo
      enddo


      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_totwind10
!-----------------------------------------------------------------------
! *** computation of strength of 10-meter winds (uv10r* 800 hPa wind)
! *** input  u800, v800 , udivg, vdivg
! *** output uv10 strength of 10 m wind at gaussian grid with minimum
! ***        uvw10 strength of 10 m wind at gaussian grid for windstress
!-----------------------------------------------------------------------



      USE comdyn, only: utot,vtot
      USE comphys, only: nlat,nlon,utot10,vtot10,uv10,uvw10,uv10m,
     &                   uv10rfx,uv10rws,uv10rwv


      implicit none


      integer i,j,k
      real*8  uv

! *** bug fix 27 march 97: uv was declared integer

      do j=1,nlon
        do i=1,nlat
          uv=sqrt((utot(i,j,3))**2 + (vtot(i,j,3))**2)
          uv10(i,j)=uv10rfx*uv
          uvw10(i,j)=uv10rws*uv
! *** minimum value of uv10 set to uv10m m/s
          if (uv10(i,j).lt.uv10m) uv10(i,j)=uv10m

!dmr @-@ iceb0
! JONO iceberg wind
          utot10(i,j) = utot(i,j,3) * uv10rwv
          vtot10(i,j) = vtot(i,j,3) * uv10rwv
!dmr @-@ iceb0

        enddo
      enddo


      return
      end
!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_tracer
!-----------------------------------------------------------------------
! *** advection of tracer field
!-----------------------------------------------------------------------


      USE comatm, only: dtime
      USE comdyn, only: pp
      USE comphys, only: nlon,nlat,nsh2,co2


      implicit none


      integer i,j
      real*8  hdivmg(nlat,nlon)
      real*8  co2sp(nsh2)



      call ec_rggtosp(co2,co2sp)
      call ec_sptogg(co2sp,co2,pp)

! *** horizontal divergence of tracer

      call ec_trafluxdiv(hdivmg,co2sp,co2)

!
! *** time stepping forward time stepping
!
      do j=1,nlon
        do i=1,nlat
          co2(i,j)=co2(i,j)-dtime*(hdivmg(i,j))
          if (co2(i,j).lt.0d0) co2(i,j)=0d0
        enddo
      enddo


      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_cloud
!-----------------------------------------------------------------------
! *** calculates total cloud cover tcc, based on a global climatology
! *** from isccpd2. Michiel Schaeffer may 1998.
! *** diagnostic cloud scheme based on relative humidity and omega
! *** and stability. Jules Beersma
! *** depending on iradcloud, cloud climatology of diagnostic clouds are
! *** used in the calculation of radiative fluxes.
!-----------------------------------------------------------------------



      use comdyn, only: omegg
      use comphys, only: nlon,nlat,relhum,tcc,tccd,relhcrit,relhfac,
     &                   iradcloud,ccisccp
      use comemic_mod, only: imonth
      use comsurf_mod, only: epss,fractoc

      implicit none

      integer i,j
      real*8 rhc, rhfac
      real*8 cc


      do j=1,nlon
        do i=1,nlat


           rhc = relhcrit
           rhfac=relhfac


! enhance clouds in case of vertical motion
           if (omegg(i,j,2) .lt.  0.0 )  rhfac=0.95d0
           if (omegg(i,j,2) .lt. -0.04)  rhfac=0.90d0
! enhance clouds in areas of subsidence inversions
           if ((fractoc(i,j) .gt. epss) .and.
     &         (omegg(i,j,2) .gt.  0.03))
     &         rhfac=0.7d0


           cc=(relhum(i,j)/rhfac - rhc)/(1.0d0 - rhc)


           cc=max(cc,0.0d0)
           cc=min(cc,1.0d0)


           tccd(i,j)=cc
           if (iradcloud.eq.1) then
             tcc(i,j)=tccd(i,j)
           else
             tcc(i,j)=ccisccp(i,j,imonth)
             tccd(i,j) = tcc(i,j)
           endif
        enddo
      enddo


      return
      end



!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_detqmtabel
!-----------------------------------------------------------------------
!***  calculate tabel of maximum content of water in the layer below
!***  500 hPa using the Clausius-Clapeyron relation and assuming a
!***  temperature profile linear in log(p): T(p)= Tr + alpha*log(p/pr)
!***  where alpha is given by (T350-T650)/log(350/650)
!***  for given groundtemperature and 650 and 350 hPa temperatures
!***  the maximum water content in [m] is given by qmtabel
!-----------------------------------------------------------------------



      USE comatm, only: grav,rlogtl12
      USE comphys, only: tzero,rowat,cc1,cc2,cc3,tqmi,tqmj,tqmk,
     &                   tqmimin,tqmjmin,tqmkmin,dtqmi,dtqmj,dtqmk,
     &                   qmtabel,iqmtab,jqmtab,kqmtab


      implicit none

      integer i,j,k
      real*8  tmount,t500,b,qmax,ec_expint,z1,z2,bz1,bz2,hulpx
      real*8  t350,t650,rlogp500,alpha
      real*8  ec_detqmaxexact
      real*4  hulp(0:iqmtab,0:jqmtab,0:kqmtab)


      rlogp500=log(500d0/650d0)
      b=cc2*cc3-cc2*tzero


!      call system('rm outputdata/atmos/qmtabel.dat')


!      open(88,file='qmtabel.dat')
!      open(89,file='qmtabel.test')
!     *  form='unformatted',access='direct',recl=51*21*21)


      do i=0,iqmtab
        tqmi(i)=tqmimin + i*dtqmi
      enddo
      do j=0,jqmtab
        tqmj(j)=tqmjmin + j*dtqmj
      enddo
      do k=0,kqmtab
        tqmk(k)=tqmkmin + k*dtqmk
      enddo


      hulpx=cc1*exp(cc2)/(rowat*grav)


      do i=0,iqmtab
        t650=tqmi(i)


        do j=0,jqmtab
          tmount=tqmj(j)+t650


          do k=0,kqmtab
            t350=t650-tqmk(k)


            alpha=(t350-t650)*rlogtl12


            t500=t650+alpha*rlogp500


            z1=1/(tmount-cc3)
            z2=1/(t500-cc3)


            bz1=b*z1
            bz2=b*z2


            qmax=(exp(bz1)+ec_expint(1,-bz1)*bz1)/z1 -
     #           (exp(bz2)+ec_expint(1,-bz2)*bz2)/z2


            qmax=qmax*hulpx/alpha


            if (qmax.lt.0d0) qmax=0d0
            qmtabel(i,j,k)=qmax
!            write(88,111) tqmi(i),tqmj(j),tqmk(k),qmtabel(i,j,k)
          enddo
        enddo
      enddo


!      write(89) (((qmtabel(i,j,k),i=0,iqmtab),j=0,jqmtab),k=0,kqmtab)



!      do i=0,iqmtab
!        temp4g(1,1)=tqmi(i)
!        do j=0,jqmtab
!          tmount=tqmj(j)+temp4g(1,1)
!          do k=0,kqmtab
!            temp2g(1,1)=temp4g(1,1)-tqmk(k)
!            hulp(i,j,k)=ec_detqmaxexact(tmount,1,1)
!            write(89,111) tqmi(i),tqmj(j),tqmk(k),hulp(i,j,k)
!          enddo
!        enddo
!      enddo



 111  format(4F14.7)
!      close(88)
!      close(89)
!      stop


      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      FUNCTION ec_detqmaxexact(tmount,i,j)
!-----------------------------------------------------------------------
! *** determines the maximum water content in latlon point
! *** i,j for given ground- and 650 and 350 hPa temperature
! *** by linear interpolation in qmtabel
!-----------------------------------------------------------------------



      USE comatm,  only: rlogtl12,grav
      USE comphys, only: tzero,rowat,cc1,cc2,cc3,temp2g,temp4g
      use newunit_mod, only: error_id

      implicit none


      integer i,j
      real*8  tmount,alpha,t500,z1,z2,bz1,bz2,b,hulpx
      real*8  qmax,ec_detqmaxexact,ec_expint


      b=cc2*cc3-cc2*tzero


      alpha=(temp2g(i,j)-temp4g(i,j))*rlogtl12


      t500=temp4g(i,j)+alpha*log(500d0/650d0)


      z1=1/(tmount-cc3)
      z2=1/(t500-cc3)


      bz1=b*z1
      bz2=b*z2


      hulpx=cc1*exp(cc2)/(rowat*grav*alpha)


      qmax=hulpx*(exp(bz1)+ec_expint(1,-bz1)*bz1)/z1 -
     #     hulpx*(exp(bz2)+ec_expint(1,-bz2)*bz2)/z2


      if (qmax.lt.0d0) qmax=0d0


      if (qmax.gt.0.2) then
!dmr --- Tentative fix to try and stop the crash of model
        qmax = 0.2
        write(error_id,*) 'in latlon ',i,j,' qmax ',qmax,'tentative fix'
        write(error_id,*) 'qmax tentative fix, tmount:', tmount
!dmr --- Tentative fix to try and stop the crash of model
        call ec_error(121)
      endif


      ec_detqmaxexact=qmax


      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      FUNCTION ec_expint(n,x)
      implicit none
      integer n,maxit
      real*8 ec_expint,x,eps,fpmin,euler
      parameter (maxit=100,eps=1.e-10,fpmin=1.e-30,euler=.5772156649)
      integer i,ii,nm1
      real*8 a,b,c,d,del,fact,h,psi
      nm1=n-1
      if(n.lt.0.or.x.lt.0..or.(x.eq.0..and.(n.eq.0.or.n.eq.1)))then
        write(*,*) 'bad arguments in ec_expint'
        read(*,*)
      else if(n.eq.0)then
        ec_expint=exp(-x)/x
      else if(x.eq.0.)then
        ec_expint=1./nm1
      else if(x.gt.1.)then
        b=x+n
        c=1./fpmin
        d=1./b
        h=d
        do 11 i=1,maxit
          a=-i*(nm1+i)
          b=b+2.
          d=1./(a*d+b)
          c=b+a/c
          del=c*d
          h=h*del
          if(abs(del-1.).lt.eps)then
            ec_expint=h*exp(-x)
            return
          endif
11      continue
!        pause 'continued fraction failed in ec_expint'
        call ec_error(20)
      else
        if(nm1.ne.0)then
          ec_expint=1./nm1
        else
          ec_expint=-log(x)-euler
        endif
        fact=1.
        do 13 i=1,maxit
          fact=-fact*x/i
          if(i.ne.nm1)then
            del=-fact/(i-nm1)
          else
            psi=-euler
            do 12 ii=1,nm1
              psi=psi+1./ii
12          continue
            del=fact*(-log(x)+psi)
          endif
          ec_expint=ec_expint+del
          if(abs(del).lt.abs(ec_expint)*eps) return
13      continue
!        pause 'series failed in ec_expint'
        call ec_error(20)
      endif
      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_trafluxdiv(tfdiv,ctrasp,ctra)
!-----------------------------------------------------------------------
! *** computes horizontal divergence of tracer flux
!-----------------------------------------------------------------------



      USE comatm, only: cosfi,radius
      USE comdyn, only: pp,pd,u800,v800,divg,udivg,vdivg
      USE comphys, only: nlon,nlat,nsh2,umoisr


      implicit none


      integer i,j,k
      real*8  ctrasp(nsh2),vv(nsh2),ww(nsh2)
      real*8  dcdl(nlat,nlon),dcdm(nlat,nlon)
      real*8  tfdiv(nlat,nlon),ctra(nlat,nlon)


! *** 800 hPa winds are reduced with umoisr in the advection of the
! *** tracer field

! *** spatial derivatives of tracer


      call ec_ddl (ctrasp,vv)
      call ec_sptogg (vv,dcdl,pp)
      call ec_sptogg (ctrasp,dcdm,pd)


! *** advection of tracer by total wind + convergence of tracer


      do j=1,nlon
        do i=1,nlat
          tfdiv(i,j)=dcdl(i,j)*(u800(i,j) + udivg(i,j,3))/
     *               (radius*cosfi(i)) +
     *               dcdm(i,j)*(v800(i,j) + vdivg(i,j,3))/
     *               (radius/cosfi(i)) +
     *               ctra(i,j)*divg(i,j,3)
          tfdiv(i,j)=tfdiv(i,j)*umoisr
        enddo
      enddo


      call ec_rggtosp(tfdiv,vv)

      vv(1)=0d0
      call ec_sptogg (vv,tfdiv,pp)

      return
      END SUBROUTINE ec_trafluxdiv


!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_hdivspec(hduvg,ug,vg)
!-----------------------------------------------------------------------
! *** computes horizontal divergence
!-----------------------------------------------------------------------


      USE comatm, only: cosfi, radius
      USE comdyn, only: pp, pd
      USE comphys, only: nlon,nlat,nsh2


      implicit none

      integer i,j
      real*8  hduvg(nlat,nlon),ug(nlat,nlon),vg(nlat,nlon)
      real*8  dugdl(nlat,nlon),dvgdm(nlat,nlon)
      real*8  usp(nsh2),vv(nsh2),vsp(nsh2)
      real*8  dx,dy(nlat)


      do j=1,nlon
        do i=1,nlat
          vg(i,j)=vg(i,j)*cosfi(i)
        enddo
      enddo


      call ec_rggtosp(ug,usp)
      call ec_rggtosp(vg,vsp)
      call ec_ddl (usp,vv)
      call ec_sptogg (vv,dugdl,pp)
      call ec_sptogg (vsp,dvgdm,pd)


      do j=1,nlon
        do i=1,nlat
          hduvg(i,j)= (dugdl(i,j)/cosfi(i)+dvgdm(i,j))/radius
        enddo
      enddo

      return
      end



!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_hdiff(rmoisslok,hdmg)
!-----------------------------------------------------------------------
! *** horizontal diffusion of moisture
!-----------------------------------------------------------------------


      USE comdyn,  only: pp
      USE comphys, only: nlon,nlat,nsh2,tdifq


      implicit none


      integer idifq,k
      real*8  hdmoiss(nsh2),hdmg(nlat,nlon)
!dmr --- added 27th may, 2011 to avoid rmoiss to be taken from common
      REAL*8  rmoisslok(nsh2)
!dmr ---
      real*8  difq,rll

      difq=max(0.d0,1.d0/(tdifq*24d0*3600d0))


      call ec_lap(rmoisslok,hdmoiss)

      hdmoiss(1)=0d0


      do k=2,nsh2
        hdmoiss(k)=difq*hdmoiss(k)
      enddo


      call ec_sptogg(hdmoiss,hdmg,pp)


      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_moisbalance(ka)
!-----------------------------------------------------------------------
! *** fix moisture balance in the atmosphere: due to advection and
! *** spectral truncation, the moisture at a specific point can become
! *** negative. All moisture additions to reset the moisture at zero
! *** have been accumulated in cormois. In this routine these additions
! *** are subtracted from rmoisg. This is done for each latitude
! *** separately to prevent an artificial meridional transport of
! *** moisture. If at a given latitude not enough moisture is available
! *** to accomodate the subtraction, the moisture balance is fixed at
! *** the neighbouring equatorward latitude. The weights pw(nlat,1) are
! *** used to correct for changes in the gridsize with latitude.
!-----------------------------------------------------------------------


      USE comdyn,  only: pw
      USE comphys, only: nlon,nlat,cormois,rmoisg
      !use comrunlabel_mod, only: irunlabelf

c~ #if (ISOATM >= 1 )
c~        USE iso_param_mod, ONLY : ieau, neauiso
c~ #endif

      implicit none

      integer, intent(in) :: ka

      integer             :: i,j,nn ! loop indicies
      double precision    :: gmc,gmm,gfac

c~ #if ( ISOATM >= 1 )
c~       integer :: k
c~       double precision, dimension(ieau+2:neauiso) :: gmc_iso, gmm_iso
c~       double precision, dimension(ieau+2:neauiso) :: coriso_done
c~ #endif

      ! From Antarctica to Equator ...
      do i=1,nlat/2
        ! For each latitudinal band ...
        gmc=0d0
        gmm=0d0
c~ #if ( ISOATM >= 1 )
c~         gmc_iso = 0.0d0
c~         gmm_iso = 0.0d0
c~ #endif
        do j=1,nlon
c~ #if ( ISOATM >= 1 )
c~           ! gmc = sum of corrections to be substracted
c~           gmc=gmc+cormois(i,j,ka)
c~           ! gmm = sum of available moisture into that lat. band
c~           gmm=gmm+rmoisg(i,j,ka)

c~           ! ditto, for water isotopes
c~           gmc_iso(:) = gmc_iso(:) + cormois(i,j,ieau+2:neauiso)
c~           gmm_iso(:) = gmm_iso(:) +  rmoisg(i,j,ieau+2:neauiso)
c~ #else
          gmc=gmc+cormois(i,j,ka)
          gmm=gmm+rmoisg(i,j,ka)
c~ #endif
        enddo

        ! gmm > 0 (there IS moisture!) and gmc > 0 (there IS something to correct for)
        if (gmm.gt.0d0.and.gmc.gt.0d0) then
          ! gmm > gmc = case where there is enough humidity to account for the correction
          if (gmm.gt.gmc) then
            gfac=gmc/gmm ! reduction factor for the humidity (plain water) < 1
c~ #if ( ISOATM >= 1 )
c~             ! same proportions for isotopes
c~             coriso_done(:) = gfac*gmm_iso(:)
c~ #endif
            do j=1,nlon
c~ #if ( ISOATM >= 1 )
c~            ! there is nothing left in cormois for plain water after this
c~            rmoisg(i,j,ka)=rmoisg(i,j,ka)-gfac*rmoisg(i,j,ka)

c~            ! same proportion for isotopes ...
c~            rmoisg(i,j,ieau+2:neauiso)=rmoisg(i,j,ieau+2:neauiso)
c~      &                               -gfac*rmoisg(i,j,ieau+2:neauiso)
c~ #else
           rmoisg(i,j,ka)=rmoisg(i,j,ka)-gfac*rmoisg(i,j,ka)
c~ #endif
            enddo ! on nlon

c~ #if ( ISOATM >= 1 )
c~             ! For water isotopes ... did we substracted enough?
c~             do k = ieau+2, neauiso
c~               ! there is still a part to be corrected ...
c~               if (coriso_done(k).lt.gmc_iso(k)) then
c~                 ! [TODO]
c~                 ! here we could handle the issue, or not !
c~               else ! nothing left already, too much done in fact ...
c~                 ! [TODO]
c~                 ! here we could handle the issue, or not !
c~               endif

c~             enddo ! on water isotopic types
c~ #endif

          else ! not enough water in that lat. band for the corrections
            ! gmm<gmc, so gfac < 1
            gfac=gmm/gmc
            gmc=(gmc-gmm)*pw(i,1)/(pw(i+1,1)*dble(nlon))

            do j=1,nlon ! on that lat. band ...
c~ #if ( ISOATM >= 1 )
c~               rmoisg(i,j,ka)=0d0 ! suppress all water in that lat. band
c~               cormois(i,j,ka)=gfac*cormois(i,j,ka)    ! correction that is dealt with at that lat.
c~               cormois(i+1,j,ka)=cormois(i+1,j,ka)+gmc ! remaining correction, for next latitude band

c~               ! for isotopes ...
c~               rmoisg(i,j,ieau+2:neauiso)=0d0
c~               ! potentially there is a part of the correction to be sent to next lat. band
c~               ! [TODO]
c~ #else
              rmoisg(i,j,ka)=0d0
              cormois(i,j,ka)=gfac*cormois(i,j,ka)
              cormois(i+1,j,ka)=cormois(i+1,j,ka)+gmc
c~ #endif
            enddo
          endif ! on water to be accomodated
        endif   ! on something to be done ...
      enddo     ! on southern hemisphere

      ! same procedure from N. Pole to equator ...
      do i=nlat,1+nlat/2,-1
        gmc=0d0
        gmm=0d0
c~ #if ( ISOATM >= 1 )
c~         gmc_iso = 0.0d0
c~         gmm_iso = 0.0d0
c~ #endif
        do j=1,nlon
c~ #if ( ISOATM >= 1 )
c~           gmc=gmc+cormois(i,j,ka)
c~           gmm=gmm+rmoisg(i,j,ka)

c~           ! ditto, for water isotopes
c~           gmc_iso(:) = gmc_iso(:) + cormois(i,j,ieau+2:neauiso)
c~           gmm_iso(:) = gmm_iso(:) +  rmoisg(i,j,ieau+2:neauiso)
c~ #else
          gmc=gmc+cormois(i,j,ka)
          gmm=gmm+rmoisg(i,j,ka)
c~ #endif
        enddo
        if (gmm.gt.0d0.and.gmc.gt.0d0) then
          if (gmm.gt.gmc) then
            gfac=1d0-gmc/gmm
c~ #if ( ISOATM >= 1 )
c~             ! same proportions for isotopes
c~             coriso_done(:) = gfac*gmm_iso(:)
c~ #endif
            do j=1,nlon
c~ #if ( ISOATM >= 1 )
c~               rmoisg(i,j,ka)=rmoisg(i,j,ka)*gfac

c~            ! same proportion for isotopes ...
c~            rmoisg(i,j,ieau+2:neauiso)=rmoisg(i,j,ieau+2:neauiso)*gfac

c~ #else
              rmoisg(i,j,ka)=rmoisg(i,j,ka)*gfac
c~ #endif
            enddo
c~ #if ( ISOATM >= 1 )
c~             ! For water isotopes ... did we substracted enough?
c~             do k = ieau+2, neauiso
c~               ! there is still a part to be corrected ...
c~               if (coriso_done(k).lt.gmc_iso(k)) then
c~                 ! [TODO]
c~                 ! here we could handle the issue, or not !
c~               else ! nothing left already, too much done in fact ...
c~                 ! [TODO]
c~                 ! here we could handle the issue, or not !
c~               endif

c~             enddo ! on water isotopic types
c~ #endif
          else
            gfac=gmm/gmc
            gmc=(gmc-gmm)*pw(i,1)/(pw(i-1,1)*dble(nlon))
            do j=1,nlon
c~ #if ( ISOATM >= 1 )
c~               rmoisg(i,j,ka)=0d0
c~               cormois(i,j,ka)=gfac*cormois(i,j,ka)
c~               cormois(i-1,j,ka)=cormois(i-1,j,ka)+gmc
c~               ! for isotopes ...
c~               rmoisg(i,j,ieau+2:neauiso)=0d0
c~               ! potentially there is a part of the correction to be sent to next lat. band
c~               ! [TODO]
c~ #else
              rmoisg(i,j,ka)=0d0
              cormois(i,j,ka)=gfac*cormois(i,j,ka)
              cormois(i-1,j,ka)=cormois(i-1,j,ka)+gmc
c~ #endif
            enddo
          endif
        endif
      enddo

      return
      end SUBROUTINE ec_moisbalance

#if ( ISOATM >= 1 )
!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_moisbalance_globiso
!-----------------------------------------------------------------------
! *** fix moisture balance in the atmosphere: due to advection and
! *** spectral truncation, the moisture at a specific point can become
! *** negative. All moisture additions to reset the moisture at zero
! *** have been accumulated in cormois. In this routine these additions
! *** are subtracted from rmoisg. This is done for each latitude
! *** separately to prevent an artificial meridional transport of
! *** moisture. If at a given latitude not enough moisture is available
! *** to accomodate the subtraction, the moisture balance is fixed at
! *** the neighbouring equatorward latitude. The weights pw(nlat,1) are
! *** used to correct for changes in the gridsize with latitude.
!-----------------------------------------------------------------------


c~ #if ( COMATM == 1 )
      USE comatm
      USE comdyn
      USE comphys
      use comrunlabel_mod, only: irunlabelf
c~ #endif

      USE iso_param_mod, ONLY : ieau, neauiso, ieau18, rsmow

      implicit none

c~ [DEPRECATED]
c~ #if ( COMATM == 0 )
c~ #include "comatm.h"
c~ #include "comdyn.h"
c~ #include "comphys.h"
c~ #include "comrunlabel.h"
c~ #endif
c~ [DEPRECATED]

      INTEGER :: i,j,k
      INTEGER, PARAMETER :: trops = 8, tropn = 25 ! only tropics 10 , 23
      REAL(KIND=8), EXTERNAL :: ec_globalmean
      REAL(KIND=8) :: gmc,gmm
      REAL(KIND=8), PARAMETER :: globfac=1d0/dsqrt(dble(nlon))
      REAL(KIND=8), DIMENSION(neauiso) :: gfac

      gfac(:)=0d0

!      WRITE(*,*) "Mean d18O avt ::",((
!     &ec_globalmean(rmoisg(:,:,ieau18))/ec_globalmean(rmoisg(:,:,ieau))
!     & )/rsmow(ieau18) - 1d0 ) *1000d0

      DO k=ieau,neauiso

!      WRITE(*,*) 'globmeancor avt', ec_globalmean(cormois(:,:,k)), k
!      WRITE(*,*) 'globmeanwat avt', ec_globalmean(rmoisg(:,:,k)), k

      gmc = ec_globalmean(cormois(:,:,k))/globfac ! global mean of correction
      gmm = 0d0

      do i=trops, tropn
      do j=1, nlon
        gmm=gmm+rmoisg(i,j,k)*pw(i,1)
      enddo
      enddo

      IF (gmm.NE.0d0) THEN

        gfac(k) = (ec_globalmean(cormois(:,:,k))/globfac) / gmm
!        WRITE(*,*) "gfac ::", 1d0 - gfac(k)
        do i=trops, tropn
        do j=1, nlon
          rmoisg(i,j,k) = rmoisg(i,j,k) * (1d0 - gfac(k))
        enddo
        enddo

        cormois(:,:,k) = 0d0
      ENDIF

!      WRITE(*,*) 'globmeancor apr', ec_globalmean(cormois(:,:,k)), k
!      WRITE(*,*) 'globmeanwat apr', ec_globalmean(rmoisg(:,:,k)), k

      ENDDO ! on k, iso

!      WRITE(*,*) "Mean d18O apr ::",((
!     &ec_globalmean(rmoisg(:,:,ieau18))/ec_globalmean(rmoisg(:,:,ieau))
!     & )/rsmow(ieau18) - 1d0 ) *1000d0

!      WRITE(*,*) "in ec_moisbalance_globiso <return>"
!      READ(*,*)

      return

      END SUBROUTINE ec_moisbalance_globiso
#endif

!23456789012345678901234567890123456789012345678901234567890123456789012
      FUNCTION ec_qsat(press,temp)
!-----------------------------------------------------------------------
! *** saturation mixing ratio
! *** input press in [Pa], temp in K
! *** output ec_qsat: saturation mixing ratio
!-----------------------------------------------------------------------


      USE comphys, only: tzero,cc1,cc2,cc3

      implicit none

      real*8  press,temp,ec_qsat


!      print*, "/", temp, tzero, cc2, cc1, cc3, press, "/"
      ec_qsat=cc1*exp(cc2*(temp-tzero)/(temp-cc3))
     &     /press

      end



!23456789012345678901234567890123456789012345678901234567890123456789012
       SUBROUTINE ec_levtemp(tlev,plev)
!-----------------------------------------------------------------------
! *** computation temperatures at level p [Pa] assuming a constant
! *** temperature lapse rate : dt/dlnp = constant
! *** input: plev
! *** ouput: tlev
!-----------------------------------------------------------------------


      USE comatm,  only: rlogtl12
      USE comphys, only: nlon,nlat,temp2g,temp4g

      implicit none

      integer i,j
      real*8  tlev(nlat,nlon)
      real*8  r,plev


      r=log(plev/65000.d0)*rlogtl12
      do j=1,nlon
        do i=1,nlat
          tlev(i,j)=temp4g(i,j)+r*(temp2g(i,j)-temp4g(i,j))
        enddo
      enddo

      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      FUNCTION ec_levtempgp(plev,i,j)
!-----------------------------------------------------------------------
! *** computation temperatures at level p [Pa] assuming a constant
! *** temperature lapse rate : dt/dlnp = constant
! *** input: plev [Pa],i,j
! *** output: ec_levtempgp [K]
!-----------------------------------------------------------------------


      USE comatm,  only: rlogtl12
      USE comphys, only: temp2g,temp4g

      implicit none

      integer i,j
      real*8  ec_levtempgp
      real*8  r,plev


      r=log(plev/65000.d0)*rlogtl12
      ec_levtempgp=temp4g(i,j)+r*(temp2g(i,j)-temp4g(i,j))

      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_fortemp
!-----------------------------------------------------------------------
! *** computes  new atmospheric temperatures due to heatforcing
! *** apply diffusion in stratosphere: timestep smaller than:
! *** 2*diffusiontimescale/462(=21*22) fe 10 days, 1 hour timestep
!-----------------------------------------------------------------------


      USE comatm,  only: dtime
      USE comdyn,  only: rinhel,pp
      USE comphys, only: nlon,nlat,nsh2,temp0g,temp2g,temp4g,thforg0,
     &                   thforg1,thforg2


      implicit none

      integer i,j,k,it,ipd
      real*8  temp0sp(nsh2),tdifc,ec_globalmean,tstep,tdifday

      tdifday=100d0
      tdifc=1.0d0/(tdifday*24.*3600.)

      tstep=2d0*tdifday*24.*3600./462.


      ipd=1+int(dtime)/int(tstep)


      tstep=dtime/ipd

      call ec_ggtosp(temp0g,temp0sp)

      do it=1,ipd

        do k=1,nsh2
          temp0sp(k)=temp0sp(k) + tstep*rinhel(k,0)*temp0sp(k)*tdifc
        enddo
      enddo

      call ec_sptogg(temp0sp,temp0g,pp)


      do j=1,nlon
        do i=1,nlat
          temp0g(i,j)=temp0g(i,j) + dtime*thforg0(i,j)
          temp2g(i,j)=temp2g(i,j) + dtime*thforg1(i,j)
          temp4g(i,j)=temp4g(i,j) + dtime*thforg2(i,j)
        enddo
      enddo


      return
      end


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
! dmr   subroutine ec_tempprofile
! *** computation of vertical temperature profile
! *** based on reference profiles from NCEP reanalysis data
! *** also used in the LWR parameterisation
! *** assuming temperature anomalies wrt these temperature profiles
! *** vary linearly with log of the pressure
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

       subroutine ec_tempprofile
       
      USE comatm,  only: rlogtl12,grav,rgas
      USE comphys, only: nlon,nlat,irn,temp0g,temp2g,temp4g,dtemp,rlogtl,tzero,
     &                   pncep,ipl,tncep,tncep12,z500ncep
      use comemic_mod, only: imonth
      use comsurf_mod, only: noc,nse,nld,ntyps,pground,pgroundn,tempsg,tempsgn,
     & rmountn,fractn


#if ( IMSK == 1 )
      use input_icemask, ONLY: icemask
#endif

#if ( DOWNSTS == 1 )
      use vertDownsc_mod, only: pground_d, tempsg_d, rmount_ps_d
      use ecbilt_topography,only: nb_levls, rmount_virt
#endif
#if ( DOWN_T2M == 1 )
      use input_subgrid2L, only: nbpointssg, weights_low_sg, index_low_sg
     &                         , tempsg_sg  
#endif

      implicit none

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
! dmr  Local variables to the subroutine
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

      integer          :: i,j,k,l,ireg(2),is,ism,nn,k1,k2
      double precision :: ro,ro1,ro2,z,z0,dt350,dt650(2),beta(2),tsref,
     &                    dt100,z1,z2
      double precision :: dtemp_tmp, pground_tmp,tsref_tmp,z_tmp,ro_tmp,
     &                    z0_tmp
      integer          :: k2_tmp
      double precision :: beta_tmp

#if ( DOWNSTS == 1 )
      double precision :: beta_d, dtemp_d, tsref_d
      integer          :: nb_down, k2_d
#endif
#if ( DOWN_T2M == 1 )
      integer :: n_point
      double precision :: ind_low, weight_low
#endif


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
! dmr   Start of the main code of the subroutine
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

! *** Example reference profile for
! *** zonal band between 15s and 15n
! *** Month  4 : tncep(19,27,12)
! ***  nr.   pfl     tfl
! ***        (mb)    (K)
! ***  1   10.00  235.113
! ***  2   20.00  223.086
! ***  3   30.00  216.125
! ***  4   50.00  206.119
! ***  5   70.00  198.420
! ***  6  100.00  196.311
! ***  7  150.00  208.144
! ***  8  200.00  221.461
! ***  9  250.00  232.611
! *** 10  300.00  242.385
! *** 11  400.00  257.836
! *** 12  500.00  268.271
! *** 13  600.00  276.160
! *** 14  700.00  282.859
! *** 15  850.00  290.298
! *** 16  925.00  294.405
! *** 17 1000.00  299.345
! *** 18 1011.99  300.156
! *** Ps 1013.00  301.265 Ts

#if ( DOWNSTS == 1 )

! dmr Update to the new global variables ... 2016-01-26

         pground_d   = 0.0
         tempsg_d    = 0.0
         rmount_ps_d = 0.0
#endif


      is=imonth/3+1
      if (is.gt.4) is=1
      ism=(is-1)*3+1

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
! dmr  Start of the spatial loop ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

      do j=1,nlon
        do i=1,nlat

! dmr&afq ireg(1) if for the ocean surface (or sea-ice)
!         ireg(2) is for the land surface  (nn = nld)

          ireg(1)=irn(i,j,1)
          ireg(2)=irn(i,j,2)

! *** logarithmic interpolation of temperature anomalies; 200 hPa or
! *** higher the anomalies approach T100 temperature within 3
! *** pressure levels


          do nn=1,2
            dt100=temp0g(i,j)-tncep(6,ireg(nn),imonth)
            dt350=temp2g(i,j)-tncep12(1,ireg(nn),imonth)
            dt650(nn)=temp4g(i,j)-tncep12(2,ireg(nn),imonth)
            beta(nn)=(dt350-dt650(nn))*rlogtl12
            do k=1,6
              dtemp(k,i,j,nn)=dt100
            enddo
            do k=7,9
              dtemp(k,i,j,nn)=((10-k)*dt100+(k-6)*(dt650(nn)
     &            + beta(nn)*rlogtl(k)))*0.25
            enddo
            do k=10,17
              dtemp(k,i,j,nn)=dt650(nn) + beta(nn)*rlogtl(k)
            enddo
          enddo
! *** from a mean height of 500 hPa from NCEP reanalysis data, the
! *** surface pressure is found using hydrostatic equilibrium
! *** and ideal gas law.


          z=z500ncep(ireg(1),imonth)
          k1=12
          DOWHILE (z.GT.rmountn(i,j,noc).AND.k1.LT.17)
            z0=z
            ro1=pncep(k1)/(rgas*(dtemp(k1,i,j,1)+
     &                        tncep(k1,ireg(1),imonth)))
            ro2=pncep(k1+1)/(rgas*(dtemp(k1+1,i,j,1)+
     &         tncep(k1+1,ireg(1),imonth)))
            ro=(ro1+ro2)*0.5
            z=z-(pncep(k1+1)-pncep(k1))/(grav*ro)
            k1=k1+1
          ENDDO
          z1=z
          pgroundn(i,j,noc)=ro*grav*(z0-rmountn(i,j,noc))+pncep(k1-1)
          pgroundn(i,j,nse)=pgroundn(i,j,noc)


#if ( ISM == 1 )
          if (flgism) then
!         rmount_ref(i,j)=rmountn(i,j,nld)
          if(rmount_ref(i,j).lt.-900.) rmount_ref(i,j)=rmountn(i,j,nld)
          if(rmount_ref(i,j).lt.0.) rmount_ref(i,j)=0.
          z_tmp=z500ncep(ireg(2),imonth)
          k2_tmp=12

          DOWHILE ((z_tmp.GT.rmount_ref(i,j)).AND.(k2_tmp.LT.17))
            z0_tmp=z_tmp
            ro1=pncep(k2_tmp)/(rgas*(dtemp(k2_tmp,i,j,2)+
     *                        tncep(k2_tmp,ireg(2),imonth)))
            ro2=pncep(k2_tmp+1)/(rgas*(dtemp(k2_tmp+1,i,j,2)+
     *         tncep(k2_tmp+1,ireg(2),imonth)))
            ro_tmp=(ro1+ro2)*0.5
            z_tmp=z_tmp-(pncep(k2_tmp+1)-pncep(k2_tmp))/(grav*ro_tmp)
            k2_tmp=k2_tmp+1
          ENDDO
          pground_tmp=ro_tmp*grav*(z0_tmp-rmount_ref(i,j))
     &                 +pncep(k2_tmp-1)
          endif
#endif

          z=z500ncep(ireg(2),imonth)
          k2=12

          if (z.lt.rmountn(i,j,nld)) then
            rmountn(i,j,nld) = z-10.0d0
          endif

          dowhile ((z.GT.rmountn(i,j,nld)).AND.(k2.LT.17))
            z0=z
            ro1=pncep(k2)/(rgas*(dtemp(k2,i,j,2)+
     &                        tncep(k2,ireg(2),imonth)))
            ro2=pncep(k2+1)/(rgas*(dtemp(k2+1,i,j,2)+
     &         tncep(k2+1,ireg(2),imonth)))
            ro=(ro1+ro2)*0.5
            z=z-(pncep(k2+1)-pncep(k2))/(grav*ro)
            k2=k2+1
          enddo

          z2=z
          pgroundn(i,j,nld)=ro*grav*(z0-rmountn(i,j,nld))+pncep(k2-1)
          pground(i,j)=0.0

          do nn=1,ntyps
            pground(i,j)=pground(i,j)+fractn(i,j,nn)*pgroundn(i,j,nn)
          enddo


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
! *** Temperature of air near surface is then found by
! *** interpolating for this pressure level the temperature profile calculated
! *** above.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

          dtemp(18,i,j,1)=dt650(1)+beta(1)*log(pgroundn(i,j,noc)/65000.)
          beta(1)=(tncep(k1,ireg(1),imonth)-tncep((k1-1),ireg(1),imonth)
     &            )/log(pncep(k1)/pncep(k1-1))
          tsref=tncep((k1-1),ireg(1),imonth)+
     &          beta(1)*log(pgroundn(i,j,noc)/pncep(k1-1))
          tempsgn(i,j,noc)=tsref+dtemp(18,i,j,1)
!dmr --- I see no reason why temperature of the sea-ice should be above freezing
!        at any time ...
!dmr ---          tempsgn(i,j,nse)=tempsgn(i,j,noc)
          tempsgn(i,j,nse)=min(tempsgn(i,j,noc),tzero)

#if ( ISM == 1 )
          if (flgism) then
          dtemp_tmp=dt650(2)+beta(2)*log(pground_tmp/65000.)
          beta_tmp=(tncep(k2_tmp,ireg(2),imonth)-
     &             tncep((k2_tmp-1),ireg(2),imonth))/
     &         log(pncep(k2_tmp)/pncep(k2_tmp-1))
          endif
#endif

#if ( DOWNSTS == 1 )

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
! dmr On cherche à calculer une anomalie de température à ajouter par
! dmr pour avoir le profil sur la grille GRISLI et non sur ECBilt
!
! dmr On calcule un ensemble de niveaux verticaux, définis dans ecbilt_topgraphy
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
! dmr   1.0) Where to do the downscaling ...
!            2016-01-25: update to have a nb_ipoints number of vertical points
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

            k2_tmp=12
            z_tmp=z500ncep(ireg(2),imonth)

            do nb_down= nb_levls, 1, -1

              if (z_tmp.lt.rmount_virt(nb_down)) then
                 rmount_virt(nb_down) = z_tmp-10.0d0
              endif

              do while ((z_tmp.gt.rmount_virt(nb_down)).and.(k2_tmp.LT.17))

               z0_tmp=z_tmp
               ro1=pncep(k2_tmp)/(rgas*(dtemp(k2_tmp,i,j,2)+
     &                  tncep(k2_tmp,ireg(2),imonth)))
               ro2=pncep(k2_tmp+1)/(rgas*(dtemp(k2_tmp+1,i,j,2)+
     &                  tncep(k2_tmp+1,ireg(2),imonth)))
               ro_tmp=(ro1+ro2)*0.5
               z_tmp=z_tmp-(pncep(k2_tmp+1)-pncep(k2_tmp))/(grav*ro_tmp)
               k2_tmp=k2_tmp+1

              enddo

             k2_d = k2_tmp-1

                                    ! pground_d is an output of the routine
             pground_d(i,j,nb_down) =
     &            ro_tmp*grav*(z0_tmp-rmount_virt(nb_down))+pncep(k2_d)
                                    ! rmount_ps_d is an output of the routine
             rmount_ps_d(i,j,nb_down) = z0_tmp

             dtemp_d=dt650(2)+beta(2)*log(pground_d(i,j,nb_down)/65000.)
             beta_d=(tncep(k2_d+1,ireg(2),imonth)-tncep((k2_d),ireg(2),
     &         imonth))/log(pncep(k2_d+1)/pncep(k2_d))

             tsref_d=tncep((k2_d),ireg(2),imonth)+
     &          beta_d*log(pground_d(i,j,nb_down)/pncep(k2_d))

                                    ! tempsg_d is an output of the routine
             tempsg_d(i,j,nb_down) = tsref_d +dtemp_d

            k2_tmp = k2_d
            z_tmp = z0_tmp
            enddo ! nb_ipoints

#if ( IMSK == 1 )
            if ( icemask(i,j).gt.0.9) then
            ! dmr
            ! temperature should be limited to zero Celsius => no above freezing
            ! on the ice-sheet
              tempsg_d(i,j,:)=min(tzero,tempsg_d(i,j,:))
            endif
#endif

#endif

#if ( DOWN_T2M == 1 )
            do n_point = 1, nbpointssg(i,j)
               ind_low = index_low_sg(i,j,n_point)
               weight_low = weights_low_sg(i,j,n_point)
               tempsg_sg(i,j,n_point) = weight_low *tempsg_d(i,j,ind_low)
     &              + (1. - weight_low)*tempsg_d(i,j,ind_low+1)
            enddo
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
! dmr  End addition of vertical downscaling in that section
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!



          dtemp(18,i,j,2)=dt650(2)+beta(2)*log(pgroundn(i,j,nld)/65000.)
          beta(2)=(tncep(k2,ireg(2),imonth)-tncep((k2-1),ireg(2),imonth)
     &            )/log(pncep(k2)/pncep(k2-1))
          tsref=tncep((k2-1),ireg(2),imonth)+
     &          beta(2)*log(pgroundn(i,j,nld)/pncep(k2-1))


#if ( ISM == 1 )
          if (flgism) then
          tsref_tmp=tncep((k2_tmp-1),ireg(2),imonth)+
     &          beta_tmp*log(pground_tmp/pncep(k2_tmp-1))
          endif
#endif


          tempsgn(i,j,nld)=tsref+dtemp(18,i,j,2)
#if ( IMSK == 1 )
          if ( icemask(i,j).gt.0.9) then
            ! dmr
            ! temperature should be limited to zero Celsius => no above freezing
            ! on the ice-sheet
              tempsgn(i,j,nld)=min(tzero,tempsgn(i,j,nld))
            endif
#endif

#if ( ISM == 1 )
          if (flgism) then
          tempsg_ism(i,j)=(tsref_tmp+dtemp_tmp)*fractn(i,j,nld)+
     *    (tempsgn(i,j,noc)*fractn(i,j,noc))+(tempsgn(i,j,nse)*
     *    fractn(i,j,nse))
          endif
#endif

          tempsg(i,j)=0.0
          do nn=1,ntyps
            tempsg(i,j)=tempsg(i,j)+fractn(i,j,nn)*tempsgn(i,j,nn)
          enddo

#if ( ISM == 1 )
          tempsg_ism(i,j)=tempsg(i,j)-tempsg_ism(i,j)
#endif

! *** Temperature of pressure levels between diagnosed surface pressure and
! *** reference surface pressure are set equal to surface air temp.


          if(z1.le.rmountn(i,j,noc))then
            do l=ipl(ireg(1)),k1,-1
               dtemp(l,i,j,1)=tempsgn(i,j,noc)-tncep(l,ireg(1),imonth)
            enddo
          endif

          if(z2.le.rmountn(i,j,nld))then
            do l=ipl(ireg(2)),k2,-1
               dtemp(l,i,j,2)=tempsgn(i,j,nld)-tncep(l,ireg(2),imonth)
            enddo
          endif

! *** for LWR parameterisation temperature anomalies wrt seasonal mean are
! *** required:
          do nn=1,2
            do k=1,18
              dtemp(k,i,j,nn)=tncep(k,ireg(nn),imonth)
     &                  +dtemp(k,i,j,nn)-tncep(k,ireg(nn),ism)
            enddo
          enddo

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
! dmr  End of spatial loop below
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

         enddo
       enddo

      return

      end subroutine ec_tempprofile

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
! dmr  ec_tempprofile <END>
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

       SUBROUTINE ec_ptmoisgp(tmount,qmax,i,j,dqmdt)

! dmr --- removed pmount from list since unused ...
!       SUBROUTINE ec_ptmoisgp(pmount,tmount,qmax,i,j,dqmdt)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
! *** computation of ground pressure and temperature in order to
! *** to calculate the maximum precipitable water content in latlon i,j
! *** qmount contains the topography for this purpose
! *** assuming temperature varies linearly with log of the pressure
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!


      USE comatm,  only: nlat, nlon, nsh, nsh2, nvl, grav, nm, ntl, rgas
     &                , rlogtl12, tlevel, plevel
      USE comdyn,  only: geopg
      USE comphys, only: temp4g, hmoisr, gpm500, qmount, temp2g
      use ec_detqmax_mod, only: ec_detqmax      
      use newunit_mod, only: error_id

      implicit none


      integer i,j

! dmr --- removed ec_levtempgp from list since unused ...
!      real*8  t500,ec_levtempgp,hmount,hred,z500,dqmdt

      real*8  t500,hmount,hred,z500,dqmdt

! dmr --- removed pmount from list since unused ...
!      real*8  alpha,pfac,hfac,pmount,tmount,qmax,ec_detqmax

      real*8  alpha,pfac,hfac,tmount,qmax
      
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
! dmr   Main code of  the function starts here !!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

      z500=gpm500*grav
      hfac=2/rgas
      hred=hmoisr*grav
      pfac=log(plevel(2)/tlevel(2))


! *** calculate temperature at t500 assuming the temperature varies
! *** linearly with log(p) : T = Tref + alpha * log (p/pref)


      alpha=(temp2g(i,j) - temp4g(i,j))*rlogtl12
      t500 =temp4g(i,j) + alpha*pfac


! *** calculate reduced ground height in decameters
! *** reduction occurs in order to tune the amount of moisture which
! *** is allowed to pass a topographic barier


      hmount=qmount(i,j)*hred
      if (hmount.lt.0d0) hmount=0d0


! *** calculate the groundpressure assuming that the mean geopotential
! *** height at 500 hPa is gpm500 decameter
! *** calculate 10 mtr temperature in K


      tmount=t500**2 - hfac*alpha*(hmount-geopg(i,j,2)-z500)

      if (tmount.lt.0) then
        write(error_id,*) 'in latlon ',i,j
        write(error_id,*) tmount,hmount,t500,geopg(i,j,2)
        call ec_error(18)
      else
        tmount=sqrt(tmount)
      endif


!      pmount=plevel(2)*exp((tmount-t500)/alpha)

      qmax=ec_detqmax(tmount,i,j,dqmdt)


      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      FUNCTION ec_globalmean(gfield)
!-----------------------------------------------------------------------
! *** computes global mean of gfield
!-----------------------------------------------------------------------


c~ #if ( COMATM == 1 )
      use comdyn, only: nlon,nlat,pw
c~ #endif

      implicit none

c~ [DEPRECATED]
c~ #if ( COMATM == 0 )
c~ #include "comatm.h"
c~ #include "comdyn.h"
c~ #endif
c~ [DEPRECATED]

      integer i,j
      real*8  gfield(nlat,nlon),suml(nlat),globsum,globfac,ec_globalmean


      globfac=1d0/dsqrt(dble(nlon))


      do i=1,nlat
        suml(i)=0d0
      enddo


      do j=1,nlon
        do i=1,nlat
          suml(i)=suml(i)+gfield(i,j)
        enddo
      enddo


      globsum=0d0


      do i=1,nlat
        globsum=globsum+suml(i)*pw(i,1)
      enddo


      ec_globalmean=globsum*globfac

      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_meantemp
!-----------------------------------------------------------------------
! *** computes mean atmospheric temperatures
!-----------------------------------------------------------------------


c~ #if ( COMATM == 1 )
      use comphys, only: temp0g,temp2g,temp4g,tempm
c~ #endif

      implicit none

c~ [DEPRECATED]
c~ #if ( COMATM == 0 )
c~ #include "comatm.h"
c~ #include "comphys.h"
c~ #include "comrunlabel.h"
c~ #endif
c~ [DEPRECATED]


      real*8  ec_globalmean


! *** mean temperatures

      tempm(0)=ec_globalmean(temp0g)
      tempm(1)=ec_globalmean(temp2g)
      tempm(2)=ec_globalmean(temp4g)

      return
      end



!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_dyntemp
!-----------------------------------------------------------------------
! *** computes temperature distribution in K from geopotential
! *** the mean level is given by tempm
! *** input:  geopg,tempm
! *** output: temp2g,temp4g
!-----------------------------------------------------------------------


c~ #if ( COMATM == 1 )
      use comatm, only: rgas
      USE comdyn, only: ntl,geopg
      USE comphys, only: nlon,nlat,temp2g,temp4g,tempm
      !use comrunlabel_mod, only: irunlabelf
c~ #endif

      implicit none

c~ [DEPRECATED]
c~ #if ( COMATM == 0 )
c~ #include "comatm.h"
c~ #include "comdyn.h"
c~ #include "comphys.h"
c~ #include "comrunlabel.h"
c~ #endif
c~ [DEPRECATED]

      integer i,j,l
      real*8  tempfac(ntl)
      real*8  geogt(nlat,nlon,ntl),tempgt(nlat,nlon,ntl)


      tempfac(1)=350.d0/(rgas*300.d0)
      tempfac(2)=650.d0/(rgas*300.d0)


      do j=1,nlon
        do i=1,nlat
          geogt(i,j,1)=geopg(i,j,1)-geopg(i,j,2)
          geogt(i,j,2)=geopg(i,j,2)-geopg(i,j,3)
        enddo
      enddo

      do l=1,ntl


! ***  calculate temperatures and add mean temperature level


        do j=1,nlon
          do i=1,nlat
            tempgt(i,j,l)=tempfac(l)*geogt(i,j,l)+tempm(l)
          enddo
        enddo


      enddo

      do j=1,nlon
        do i=1,nlat
          temp2g(i,j)=tempgt(i,j,1)
          temp4g(i,j)=tempgt(i,j,2)
        enddo
      enddo


      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_vortfor
!-----------------------------------------------------------------------
! *** computes vorticity forcing
!-----------------------------------------------------------------------


c~ #if ( COMATM == 1 )
      use comatm, only: om,cosfi
      USE comdyn, only: pp,psi,for,dfor1,dfor2,forcgg1,rinhel,iartif,
     &                  ipvf1,ipvf2,ipvf3,ipvf4,ipvf5
      USE comphys, only: nlon,nlat,nvl,nsh2,vfor1,vfor2,vfor3,vforg1,
     &                   vforg1,vforg2,vforg3,vhforg1,vhforg2
c~ #endif

      implicit none

c~ [DEPRECATED]
c~ #if ( COMATM == 0 )
c~ #include "comatm.h"
c~ #include "comdyn.h"
c~ #include "comphys.h"
c~ #include "comrunlabel.h"
c~ #endif
c~ [DEPRECATED]

      integer i,j,k,l
      real*8  vforg(nlat,nlon,nvl)
      real*8  forcgg(nlat,nlon),forcgg2(nlat,nlon)
      real*8  zetas(nsh2,3),zetag(nlat,nlon,3)
      real*8  dimfac


      call ec_rggtosp(vhforg1,dfor1)
      call ec_rggtosp(vhforg2,dfor2)


! ***  compute the relative vorticity zeta from steamFUNCTION psi
! ***  psi is dimensionless, zeta with dimension


      do l=1,nvl
        do k=1,nsh2
          zetas(k,l)=rinhel(k,0)*psi(k,l)*om
        enddo
        call ec_sptogg(zetas(1,l),zetag(1,1,l),pp)
        do j=1,nlon
          do i=1,nlat
            vforg(i,j,l)=0d0
          enddo
        enddo
      enddo


! *** potential vorticity (pv) forcing (1/(second**2)) due to
! *** the diabatic heating


      if (ipvf1.ne.0) call ec_pvf1(vforg)


! *** pv forcing due to advection of the planetary vorticy f
! *** by the divergent wind

      if (ipvf2.ne.0) call ec_pvf2(vforg)


! *** pv forcing due to d*zeta. d is the divergence.
! *** zeta is the relative vorticity


      if (ipvf3.ne.0) call ec_pvf3(vforg,zetag)


! *** pv forcing due to advection of the relative vorticy
! *** zeta by the divergent wind


      if (ipvf4.ne.0) call ec_pvf4(vforg,zetas)


! *** pv forcing due to advection of temperature by
! *** the divergent wind


      if (ipvf5.ne.0) call ec_pvf5(vforg)


! *** the total pv forcing in nondimensional units


      dimfac=1./(om**2)


      call ec_ggtosp(vforg(1,1,1),vfor1)
      call ec_ggtosp(vforg(1,1,2),vfor2)
      call ec_ggtosp(vforg(1,1,3),vfor3)


! *** adding the artificial forcing (forcgg1)

      if (iartif .eq. 1) then

         call ec_sptogg(vfor1,forcgg2,pp)
         do i=1,nlat
           do j=1,nlon
             forcgg(i,j)=forcgg1(i,j)*cosfi(i)**8/dimfac+forcgg2(i,j)
           enddo
         enddo
         call ec_rggtosp(forcgg,vfor1)
      endif


      vfor1(1)=0.d0
      vfor2(1)=0.d0
      vfor3(1)=0.d0


      do k=2,nsh2
        for(k,1)=dimfac*vfor1(k)
        for(k,2)=dimfac*vfor2(k)
        for(k,3)=dimfac*vfor3(k)
      enddo


! *** transfer to grid point


      call ec_sptogg(vfor1,vforg1,pp)
      call ec_sptogg(vfor2,vforg2,pp)
      call ec_sptogg(vfor3,vforg3,pp)


      return
      end



!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_pvf1(vforg)
!-----------------------------------------------------------------------
! *** potential vorticity (pv) forcing (1/(second**2)) due to
! *** the adiabatic heating
!-----------------------------------------------------------------------


c~ #if ( COMATM == 1 )
      USE comatm, only: fzero, radius, om, rgas, sinfi, dp
      USE comdyn, only: pp, rrdef1,rrdef2
      USE comphys, only: nlon,nlat,nvl,nsh2,vhforg1,vhforg2
c~ #endif

      implicit none

c~ [DEPRECATED]
c~ #if ( COMATM == 0 )
c~ #include "comatm.h"
c~ #include "comdyn.h"
c~ #include "comphys.h"
c~ #include "comrunlabel.h"
c~ #endif
c~ [DEPRECATED]

      integer i,j,k,l
      real*8  pvf1s(nsh2,3),vforg(nlat,nlon,nvl)
      real*8  vhfor1(nsh2),vhfor2(nsh2)
      real*8  vhforg1x(nlat,nlon),vhforg2x(nlat,nlon)
      real*8  vorfac1,vorfac2,drdef1,drdef2


      vorfac1=+rgas*dp/3.5d+4
      vorfac2=+rgas*dp/6.5d+4
      drdef1=1./( om*fzero*(radius*rrdef1)**2 )
      drdef2=1./( om*fzero*(radius*rrdef2)**2 )


      do j=1,nlon
        do i=1,nlat
          vhforg1x(i,j)=vhforg1(i,j)*sinfi(i)/fzero
          vhforg2x(i,j)=vhforg2(i,j)*sinfi(i)/fzero
        enddo
      enddo


      call ec_rggtosp(vhforg1x,vhfor1)
      call ec_rggtosp(vhforg2x,vhfor2)


      pvf1s(1,1)=0d0
      pvf1s(1,2)=0d0
      pvf1s(1,3)=0d0


      do k=2,nsh2
        pvf1s(k,1)=-drdef1*vorfac1*vhfor1(k)
        pvf1s(k,2)=drdef1*vorfac1*vhfor1(k)-drdef2*vorfac2*vhfor2(k)
        pvf1s(k,3)=drdef2*vorfac2*vhfor2(k)
      enddo


      do l=1,nvl
        call ec_sptogg(pvf1s(1,l),vforg(1,1,l),pp)
      enddo

      return
      end



!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_pvf2(vforg)
!-----------------------------------------------------------------------
! *** pv forcing due to advection of the planetary vorticy f
! *** by the divergent wind
!-----------------------------------------------------------------------


c~ #if ( COMATM == 1 )
      USE comatm, only: om,cosfi,radius
      USE comdyn, only: vdivg
      USE comphys, only: nlon,nlat,nvl
c~ #endif

      implicit none

c~ [DEPRECATED]
c~ #if ( COMATM == 0 )
c~ #include "comatm.h"
c~ #include "comdyn.h"
c~ #include "comphys.h"
c~ #include "comrunlabel.h"
c~ #endif
c~ [DEPRECATED]

      integer i,j,k,l
      real*8  vforg(nlat,nlon,nvl)

      do l=1,nvl
        do j=1,nlon
          do i=1,nlat
            vforg(i,j,l)=vforg(i,j,l)-vdivg(i,j,l)*om*cosfi(i)/radius
          enddo
        enddo
      enddo


      return
      end



!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_pvf3(vforg,zetag)
!-----------------------------------------------------------------------
! *** pv forcing due to d*zeta. d is the divergence.
! *** zeta is the relative vorticity
!-----------------------------------------------------------------------


c~ #if ( COMATM == 1 )
      USE comdyn, only: divg
      USE comphys, only: nlon,nlat,nvl
c~ #endif

      implicit none

c~ [DEPRECATED]
c~ #if ( COMATM == 0 )
c~ #include "comatm.h"
c~ #include "comdyn.h"
c~ #include "comphys.h"
c~ #include "comrunlabel.h"
c~ #endif
c~ [DEPRECATED]

      integer i,j,k,l
      real*8  vforg(nlat,nlon,nvl),zetag(nlat,nlon,nvl)


      do l=1,nvl
        do j=1,nlon
          do i=1,nlat
            vforg(i,j,l)=vforg(i,j,l)-divg(i,j,l)*zetag(i,j,l)
          enddo
        enddo
      enddo


      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_pvf4(vforg,zetas)
!-----------------------------------------------------------------------
! *** pv forcing due to advection of the relative vorticy
! *** zeta by the divergent wind
!-----------------------------------------------------------------------


c~ #if ( COMATM == 1 )
      USE comatm, only: radius,cosfi
      USE comdyn, only: pp,pd,udivg,vdivg
      USE comphys, only: nlon,nlat,nvl,nsh2
c~ #endif

      implicit none

c~ [DEPRECATED]
c~ #if ( COMATM == 0 )
c~ #include "comatm.h"
c~ #include "comdyn.h"
c~ #include "comphys.h"
c~ #include "comrunlabel.h"
c~ #endif
c~ [DEPRECATED]

      integer i,j,k,l
      real*8  vforg(nlat,nlon,nvl),zetas(nsh2,nvl)
      real*8  x(nsh2),xhelp(nsh2),dxdl(nlat,nlon),dxdm(nlat,nlon)


      do l=1,nvl
        do k=1,nsh2
          x(k)=zetas(k,l)
        enddo
        call ec_ddl(x,xhelp)
        call ec_sptogg(xhelp,dxdl,pp)
        call ec_sptogg(x,dxdm,pd)
        do i=1,nlat
          do j=1,nlon
             vforg(i,j,l)=vforg(i,j,l)
     &                  -udivg(i,j,l)*dxdl(i,j)/(radius*cosfi(i))
     &                  -vdivg(i,j,l)*cosfi(i)*dxdm(i,j)/radius
          enddo
        enddo
      enddo


      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_pvf5(vforg)
!-----------------------------------------------------------------------
! *** computes vorticity forcing due to advection of temperature by
! *** divergent wind
!-----------------------------------------------------------------------


c~ #if ( COMATM == 1 )
      USE comatm, only: radius,om,cosfi,sinfi,fzero
      USE comdyn, only: pp,pd,psi,udivg,vdivg,rrdef1,rrdef2
      USE comphys, only: nlon,nlat,nvl,nsh2
c~ #endif

      implicit none

c~ [DEPRECATED]
c~ #if ( COMATM == 0 )
c~ #include "comatm.h"
c~ #include "comdyn.h"
c~ #include "comphys.h"
c~ #include "comrunlabel.h"
c~ #endif
c~ [DEPRECATED]

      integer i,j,k,l
      real*8  vforg(nlat,nlon,nvl),pvf7g(nlat,nlon,nvl)
      real*8  ud(nlat,nlon,2),vd(nlat,nlon,2)
      real*8  x(nsh2,2),y(nlat,nlon,3)
      real*8  xhelp(nsh2),dxdl(nlat,nlon),dxdm(nlat,nlon)
      real*8  rr1,rr2,sinfact


      rr1=1.d0/(rrdef1*radius)**2
      rr2=1.d0/(rrdef2*radius)**2


      do j=1,nlon
        do i=1,nlat
          ud(i,j,1)=(udivg(i,j,1)+udivg(i,j,2))/2.d0
          ud(i,j,2)=(udivg(i,j,2)+udivg(i,j,3))/2.d0
          vd(i,j,1)=(vdivg(i,j,1)+vdivg(i,j,2))/2.d0
          vd(i,j,2)=(vdivg(i,j,2)+vdivg(i,j,3))/2.d0
        enddo
      enddo

      do k=1,nsh2
        x(k,1)=(psi(k,1)-psi(k,2))*om*radius**2
        x(k,2)=(psi(k,2)-psi(k,3))*om*radius**2
      enddo


      do l=1,2
        call ec_ddl(x(1,l),xhelp)
        call ec_sptogg(xhelp,dxdl,pp)
        call ec_sptogg(x(1,l),dxdm,pd)

        do j=1,nlon
          do i=1,nlat
             y(i,j,l)=ud(i,j,l)*dxdl(i,j)/(radius*cosfi(i))
     &               +vd(i,j,l)*cosfi(i)*dxdm(i,j)/radius
          enddo
        enddo
      enddo


      do j=1,nlon
        do i=1,nlat
          pvf7g(i,j,1)=rr1*y(i,j,1)
          pvf7g(i,j,2)=-rr1*y(i,j,1)+rr2*y(i,j,2)
          pvf7g(i,j,3)=-rr2*y(i,j,2)
        enddo
      enddo


      do i=1,nlat
        sinfact=(sinfi(i)/fzero)**2
        do j=1,nlon
          vforg(i,j,1)=vforg(i,j,1)+pvf7g(i,j,1)*sinfact
          vforg(i,j,2)=vforg(i,j,2)+pvf7g(i,j,2)*sinfact
          vforg(i,j,3)=vforg(i,j,3)+pvf7g(i,j,3)*sinfact
        enddo
      enddo

      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_lwaverad2(nn)
!-----------------------------------------------------------------------
! *** computes long wave radiation according to the parameterization of
! *** Chao Chou and Neelin and substantially adapted and extended
! *** for global scale and more
! *** specific ECBILT application by the one and only Michiel Schaeffer
! ***
! *** parameters: nlat   = number of gridpoints in the meridional
! ***                      direction (32)
! ***             nlon   = number of gridpoints in the zonal
! ***                      direction (64)
! ***
! *** input : dtemp(19,nlat,nlon): temperature anomalies [K] wrt ncep
! ***                              climatology tncep in common lwrscheme
! ***         dqa(nlat,nlon) : anomalies of total prec. water cont. below
! ***                          500 hPa wrt ncep climatology
! ***         tcc(nlat,nlon)  : total cloud cover
! ***         ghg(19) : concentrations of well mixed ghg s (see comphys.h)
! ***
! *** output : ulrad1(nlat,nlon): upward longwave radiation [Wm-2] at
! ***                             toa
! ***          ulrad2(nlat,nlon): net longwave radiation [Wm-2] at
! ***                             (500 hPa)
! ***          dlrads(nlat,nlon): downward longwave radiation [Wm-2] at
! ***                             the surface
! ***          ulrads(nlat,nlon): upward longwave radiation [Wm-2] at
! ***                             the surface
!-----------------------------------------------------------------------

c~ #if ( COMATM == 1 )
      use comphys, only: nlon,nlat,irn,ipl,tncep,dtemp,qancep,tcc,expir,
     &                   ghgscen,ghgipcc,lwrt,lwrts,lwrqa,lwrref,lwrghg,
     &                   lwrqts
      use comemic_mod, only: iyear, imonth
      use comsurf_mod, only: nse,noc,nld,tsurfn,fractn,lwrmois
c~ #endif

      use ipcc_output_mod, only: moc,tmc,tmc0,tsurfmean,cland,thex
      use commons_mod,     only: ulrad0nz,ulrad1nz,ulrad0nT,ulrad1nT

!!    USE OMP_LIB

           implicit none

c~ [DEPRECATED]
c~ #if ( COMATM == 0 )
c~ #include "comatm.h"
c~ #include "comdyn.h"
c~ #include "comphys.h"
c~ #include "comemic.h"
c~ #include "comsurf.h"
c~ #include "comunit.h"
c~ #endif




      integer i,j,l,k,m,is,ism,nol,nn,ireg,h,r,s,igas
      real*8  lwrz(7,0:1),dumts
      real*8  dqa,dqreg(27)
      real*8  ulrad0nm,ulrad1nm,ulrad2nm,ulradsnm,dlradsnm
      real*8  ulrad0nmm,ulrad1nmm
      real*8  ulrad0nU,ulrad1nU,ulrad2nU,ulradsnU,dlradsnU
c~       real*8  ulrad0nz(nlat,nlon),ulrad1nz(nlat,nlon)
      real*8  ulrad1nzz(nlat,nlon,3)
      real*8  ulrad2nz(nlat,nlon),ulradsnz(nlat,nlon)
      real*8  dlradsnz(nlat,nlon)
c~       real*8  ulrad0nT,ulrad1nT
      real*8  ulrad2nT,ulradsnT,dlradsnT
      real*8  ec_globalmean
      real*8  logco2T,sqrch4T,sqrn2oT,ghgz(20)
      real*8  alpho3lw(2)
      real*4  lwrfluxz(7,27,4,0:1,2)
c~       real*8  moc,tmc,tmc0,tsurfmean,cland,thex
!dmr [OMP]
      real*8  mompdy
!dmr [OMP]
      common / radO3 / ulrad0nU,ulrad1nU,ulrad2nU,ulradsnU,dlradsnU
c~       common /rad031/ulrad0nz,ulrad1nz,ulrad0nT,ulrad1nT
c~       common/IPCC_out2/moc,tmc,tmc0,tsurfmean,cland,thex



!dmr ###        ghgz(1)=280.
!dmr ###        do igas=2,20
        do igas=1,20
         ghgz(igas)=ghgscen(igas,1)
        enddo
        logco2T=log(ghgz(1)/ghgipcc(1))
        sqrch4T=sqrt(ghgz(2))-sqrt(ghgipcc(2))
        sqrn2oT=sqrt(ghgz(3))-sqrt(ghgipcc(3))
        alpho3lw(1)=153.6
        alpho3lw(2)=201.2
c~ !$OMP PARALLEL
c~ !$OMP DO COLLAPSE(3) PRIVATE (h,l,s,r,k,m,mompdy) SCHEDULE(static)
        do h=1,2
        do l=0,1
         do s=1,4
          do r=1,27
           do k=1,7
!dmr [OMP]            lwrfluxz(k,r,s,l,h)=lwrref(k,r,s,l)+
             mompdy=lwrref(k,r,s,l)+
     *            lwrghg(k,1,r,s,l)*logco2T+
     *            lwrghg(k,2,r,s,l)*sqrch4T+
     *            lwrghg(k,3,r,s,l)*sqrn2oT
            do m=4,19
!dmr [OMP]             lwrfluxz(k,r,s,l,h)=lwrfluxz(k,r,s,l,h)+
             mompdy=mompdy+
     *            lwrghg(k,m,r,s,l)*(ghgz(m)-ghgipcc(m))
            enddo
!dmr [OMP]              lwrfluxz(k,r,s,l,h)=lwrfluxz(k,r,s,l,h)+
              lwrfluxz(k,r,s,l,h)=mompdy+
     *              lwrghg(k,4,r,s,l)*alpho3lw(h)*(ghgz(20)-25.)
           enddo
          enddo
         enddo
        enddo
        enddo
c~ !$OMP END DO NOWAIT
c~ !$OMP SINGLE
       is=imonth/3+1
       if (is.gt.4) is=1
       ism=(is-1)*3+1

      do i=1,27
!dqa    dqreg(i)=qancep(i,ism)**0.3333
        dqreg(i)=qancep(i,ism)
      enddo

      if (nn.eq.noc.or.nn.eq.nse) nol=1
      if (nn.eq.nld) nol=2

c~ !$OMP END SINGLE
c~ !$OMP DO PRIVATE (j,i,l,k,h,m,ireg,dqa,lwrz,dumts,mompdy) SCHEDULE(static)
      do j=1,nlon
        do i=1,nlat
          ireg=irn(i,j,nol)

!-Hemispheric dependence of tropospheric ozone forcing
          if (i.le.16) then
           h=1
          else
           h=2
          endif

!dqa      dqa=lwrmois(i,j)-dqreg(ireg)
!dqa      q**1/3-qm**1/3=qm**(1/3-n)*(q**n-qm**n)
          dqa=dqreg(ireg)**(0.3333-EXPIR)*
     *          (lwrmois(i,j)**EXPIR-dqreg(ireg)**EXPIR)

          do l=0,1
            do k=1,7
              lwrz(k,l)=lwrfluxz(k,ireg,is,l,h)+lwrqa(k,ireg,is,l)*dqa
              do m=1,ipl(ireg)-1
                lwrz(k,l)=lwrz(k,l)+
     *                   lwrt(k,m,ireg,is,l)*dtemp(m,i,j,nol)
              enddo
              lwrz(k,l)=lwrz(k,l)+
     *              lwrt(k,18,ireg,is,l)*dtemp(18,i,j,nol)
            enddo

            dumts=tsurfn(i,j,nn)-tncep(19,ireg,ism)
            do m=1,4
              do k=1,3
                lwrz(k,l)=lwrz(k,l)+
     *          (lwrts(k,m,ireg,is,l)+lwrqts(k,m,ireg,is,l)*dqa)
     *          *dumts
              enddo
!             lwrz(7,l)=lwrz(7,l)+
!    *        (lwrts(7,m,ireg,is,l)+lwrqts(7,m,ireg,is,l)*dqa)
!    *        *dumts
              dumts=dumts*(tsurfn(i,j,nn)-tncep(19,ireg,ism))
            enddo


          enddo

          if (nn.eq.1) then
            ulrad1nz(i,j)=0.
            if (iyear.eq.0) then
             ulrad1nT=0.
            endif
          endif


!dmr [OMP]          ulrad1nzz(i,j,nn)=(lwrz(2,0)+lwrz(5,0))*(1-tcc(i,j)) +
          mompdy=(lwrz(2,0)+lwrz(5,0))*(1-tcc(i,j)) +
     *             (lwrz(2,1)+lwrz(5,1))*tcc(i,j)
          ulrad1nzz(i,j,nn)=mompdy
!dmr [OMP]         ulrad1nz(i,j)=ulrad1nz(i,j)+(ulrad1nzz(i,j,nn)*fractn(i,j,nn))
          ulrad1nz(i,j)=ulrad1nz(i,j)+(mompdy*fractn(i,j,nn))

        enddo
      enddo
c~ !$OMP END DO
c~ !$OMP END PARALLEL

      if (nn.eq.3) then
      ulrad1nm=ec_globalmean(ulrad1nz)

      ulrad1nmm=ulrad1nm-ulrad1nU

      ulrad1nT=ulrad1nT+(ulrad1nmm/(360.*6.))

      endif

! *** that s all folks
      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_swaverad2(nn)
!-----------------------------------------------------------------------
! *** computes short wave radiation
! *** linearization of RCM with ISCCP D2 1990 clouds
!-----------------------------------------------------------------------

c~ #if ( COMATM == 1 )
      use comphys, only: nlon,nlat,irn,dayfr,solarf,solarm,solardref,
     &                   kosz,costref,swrref,swrcost,dso4,tcc,salbref,
     &                   swrsalb,fswdsfc
      use comemic_mod, only: imonth
      use comsurf_mod, only: nse,noc,nld,fractn,albesn,alb2esn,heswsn,
     &                       hesw0n,hesw1n,hesw2n
c~ #endif

      implicit none

c~ [DEPRECATED]
c~ #if ( COMATM == 0 )
c~ #include "comatm.h"
c~ #include "comdyn.h"
c~ #include "comphys.h"
c~ #include "comemic.h"
c~ #include "comsurf.h"
c~ #include "comunit.h"
c~ #endif
c~ [DEPRECATED]


#if ( CLAQUIN == 1 )
!dmr --- Pour le forcage a la Claquin et al., 2003
#include "radforc.h"
!dmr --- Pour le forcage a la Claquin et al., 2003
#endif
      integer i,j,k,l,ireg
      integer m, d, r, nn , nol


      real*8 f0,f1,ftot(8),fn(8,0:1)
      real*8 drs, drs2, drs3
      real*8 dcost, df,sk,sr,x,y,dfs,smsc,df2
      real*8 fswutoa(nlat,nlon)! ,fswdsfc(nlat,nlon),fswusfc
      real*8 fswutoa2(nlat,nlon),fswdtoa(nlat,nlon)
      real*8 fswutoaG,fswdtoa2,fswdtoaG


      integer nreg(2),indxsul
!     real*8 zac(2),asup,bup
      real*8 zac(2),asup
      real*8  ec_globalmean
      real*8 fswutoaGA,fswutoaG0
      real*8 fswutoa_diff,df_test
      common /rad_sul2 /fswutoa,fswdtoa
      common /rad_sul0 /fswutoaG,df_test,fswdtoaG
!      common /pr_evap /fswdsfc
! *** aerosol scattering included as a correction on the upward
! *** clear sky fluxes
! *** sk,sr: empirical coefficients Dorland et al, J. Geophys. Res.,102,
! *** 28079-28100, 1997.
! *** smsc: mass scattering coefficient [m2/g]
! *** dso4: change in sulfate aerosol column integrated concentration since
! *** pre-industrial times [g/m2]

#if ( CLAQUIN == 1 )
!dmr --- Pour le forcage a la Claquin et al., 2003
      REAL*8 dfprime
!dmr --- Pour le forcage a la Claquin et al., 2003
#endif

! dmr [TODO] In this routine, it seems that the solar constant is implicitly taken as
!             being 1370, inconsistent with the global shared value for the rest of ECBilt
!            A way to change it? Or is it just meant to be so?
! dmr

      sk=0.058d0*1370d0
      sr=0.05d0
      smsc=8.0

      if (nn.eq.noc.or.nn.eq.nse) nol=1
      if (nn.eq.nld) nol=2


      do j=1,nlon
       do i=1,nlat
           alb2esn(i,j,nn)=albesn(i,j,nn)
!           alb2esn(i,j,3) = albesnR(i,j)
          if (alb2esn(i,j,nn).ge.1.) then
            alb2esn(i,j,nn)=1.
          endif
          df=dayfr(i)*solarf(i)

! dmr [NOTA] solarm == solar constant (currently 1368 ...)
!            solardref == 1/ro**2.
!                      with ro Earth-Sun distance, semi-major axis (in meters ?)
!                      so ro**2 % m^2
!                      and solardref % m^-2
!            dayfr == ha1 / pi % rad / pi hence % degree of circle
!                  with ha1 == acos(rkosha1)
!                           and rkosha1 == -tanfi(i)*tandl
!                                       where tanfi == tangent(phi) where phi is latitude in rad
!                                       and   tandl == is an arc-sine function of obliquity and  day from the VE
          df2=dayfr(i)*solarm*solardref/1370.d0

          ireg=irn(i,j,nol)
          dcost=kosz(i)-costref(ireg,imonth)
          do l=1,8
            do k=0,1
              fn(l,k) =
     &               swrref(l,ireg,imonth,k)
     &               +  swrcost(l,ireg,imonth,k)*dcost
            enddo
          enddo

          x=sqrt(kosz(i))
          y=sqrt(1-alb2esn(i,j,nn))
          dfs=sk*(4d0*x*y*(y-x)-sr)*dso4(i,j)*smsc
          if (dfs.gt.0d0.and.kosz(i).lt.0.05) dfs=0d0
          drs=alb2esn(i,j,nn)-salbref(ireg,imonth)
          drs2=drs*drs
          drs3=drs2*drs


          do l=1,4
           f0=fn(l,0)+swrsalb(l,ireg,imonth,0)*drs+dfs
           f1=fn(l,1)+swrsalb(l,ireg,imonth,1)*drs
     &                     +swrsalb(l,ireg,imonth,2)*drs2
     &                     +swrsalb(l,ireg,imonth,3)*drs3
            ftot(l) = (1.-tcc(i,j))*f0 + tcc(i,j)*f1
          enddo
          do l=5,8
            f0=fn(l,0)+swrsalb(l,ireg,imonth,0)*drs
            f1=fn(l,1)+swrsalb(l,ireg,imonth,1)*drs
     &                     +swrsalb(l,ireg,imonth,2)*drs2
     &                     +swrsalb(l,ireg,imonth,3)*drs3
            ftot(l) = (1.-tcc(i,j))*f0 + tcc(i,j)*f1
          enddo
#if ( CLAQUIN == 1 )
!dmr --- Ajout de l utilisation du forcage des poussieres
!dmr ---
!dmr --- Distinction des cas forcage transitoire ou non

!dmr --- Si transitoire ...
          IF (transitforce.NE.0) THEN

!dmr --- Si l annee est superieure au forcage total
!dmr ---   ---> on conserve la derniere annee !
            IF (iyear.GT.transitforce) THEN
              dfprime=df*(solarc+radforc(i,j)
     >     *facteurpoussieres(transitforce))/solarc

              if
     >        ((imonth.eq.1).and.(iday.eq.1).and.(i.eq.1).and.(j.eq.1))
     >        write(info_id,*) 'facteur poussieres'
     >     , facteurpoussieres(transitforce)

            ELSE
              dfprime=df*(solarc+radforc(i,j)
     >     *facteurpoussieres(iyear))/solarc

              if
     >        ((imonth.eq.1).and.(iday.eq.1).and.(i.eq.1).and.(j.eq.1))
     >        write(info_id,*) 'facteur poussieres'
     >     , facteurpoussieres(iyear)

            ENDIF

!dmr --- Si pas transitoire ...
          ELSE
            dfprime=df*(solarc+radforc(i,j))/solarc
          ENDIF
!dmr ---
#endif

! alternative calculation of upward flux at ground:
! in parameterisation no cross terms are accounted for, which are important for
! upward shortwave radiation at surface and therefore also for net flux
! heswsn(i,j)

          ftot(4)=-alb2esn(i,j,nn)*ftot(8)
#if ( CLAQUIN == 1 )
!dmr --- Ajout de dfprime
          hesw0n(i,j,nn)=(-ftot(1)-ftot(5)+ftot(2)+ftot(6))*dfprime
          hesw1n(i,j,nn)=(-ftot(2)-ftot(6)+ftot(3)+ftot(7))*dfprime
          hesw2n(i,j,nn)=(-ftot(3)-ftot(7)+ftot(4)+ftot(8))*dfprime
          heswsn(i,j,nn)=(-ftot(4)-ftot(8))*dfprime
#else
          hesw0n(i,j,nn)=(-ftot(1)-ftot(5)+ftot(2)+ftot(6))*df
          hesw1n(i,j,nn)=(-ftot(2)-ftot(6)+ftot(3)+ftot(7))*df
          hesw2n(i,j,nn)=(-ftot(3)-ftot(7)+ftot(4)+ftot(8))*df
          heswsn(i,j,nn)=(-ftot(4)-ftot(8))*df
#endif

! for diagnostic purposes:
! (1) downward shortwave radiation at TOA
             if (nn.eq.1) then 
                fswdtoa(i,j)=0.
             endif
             fswdtoa2=-ftot(5)*df2
             fswdtoa(i,j)=fswdtoa(i,j)+(fractn(i,j,nn)*fswdtoa2)
! (2) upward shortwave radiation at TOA
            if (nn.eq.1) then 
              fswutoa(i,j)=0.
!              fswusfc(i,j) = 0.
            endif
            fswutoa2(i,j)=ftot(1)*df2
            fswutoa(i,j)=fswutoa(i,j)+(fractn(i,j,nn)*fswutoa2(i,j))
! (3) downward shortwave radiation at SURFACE
            fswdsfc(i,j)=-ftot(8)*df
! (4) upward shortwave radiation at SURFACE
! incorrect, ftot(8) does not depend on surfaces ...            fswusfc(i,j)=(heswsn(i,j,nn)+ftot(8)*df)
        enddo
      enddo
      if (nn.eq.3) then
        fswutoaG=ec_globalmean(fswutoa)
        fswdtoaG=ec_globalmean(fswdtoa)
      endif


      return
      end

#if ( CLAQUIN == 1 )
!dmr %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!dmr --- Ci-dessous ajoute pour le forcage radiatif des poussieres, 21 juin 2006

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE rad_forc_add(loc_year)

c~ #if ( COMATM == 1 )
      use global_constants_mod, only: dblp=>dp, ip     
c~ #endif

c~ [DEPRECATED]
c~ #if ( COMATM == 0 )
c~ #include "comatm.h"
c~ #include "comunit.h"
c~ #endif
c~ [DEPRECATED]

c~ [TODO] Update radforc to a Fortran90 module or deprecate!!
#include "radforc.h"

      INTEGER loc_year, i, j
      REAL videtotal
      integer(kind=ip):: Forcage_Claquin_dat_id,Dust_scale_factor_txt_id

!dmr --- ---------------------------------------------------------------
! *** Load an additional radiative forcing
! *** to take into account, e.g. the dust at the LGM
!dmr --- Created : Didier M. Roche, Cedric Van Meerbeeck
!dmr --- Last updated : July, 12th, 2007
!dmr --- ---------------------------------------------------------------

       print*, "Trying to read radiative forcing", transitforce

       open(newunit=Forcage_Claquin_dat_id,file='Forcage-Claquin.dat',
     >     form='unformatted',access='direct', recl=nlat*nlon*4)
       READ(Forcage_Claquin_dat_id, REC=loc_year) 
     &         ((radforc(i,j),j=1,nlon),i=1,nlat)
       close(Forcage_Claquin_dat_id)

       print*, "read additional radiative forcing", transitforce

       if (transitforce.ne.0) then
         open(newunit=Dust_scale_factor_txt_id,
     &          file='Dust-scale-factor.txt')
         DO i=1,transitforce
            READ(Dust_scale_factor_txt_id,*) videtotal, 
     &         facteurpoussieres(i)
         ENDDO
         close(Dust_scale_factor_txt_id)
       endif

       print*, "Read fact poussieres ! "

       END
#endif

