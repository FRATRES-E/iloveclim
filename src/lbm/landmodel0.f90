!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:21 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:21 CET 2009

      module landmodel_mod

      use global_constants_mod, only: dblp=>dp, ip

      implicit none

      public

        ! Fate of the snow above dsnm in lbm_landprecip (wiso:3a:0000)
      integer, parameter :: snow_excess_bucket = 1 ! into the soil water (bmoisg), the heat counted in heatsnown/s
      integer, parameter :: snow_excess_out    = 2 ! removed from the snow only (icebergs, ICEBERG == 2)
      integer, parameter :: snow_excess_keep   = 0 ! no cap (ICEBERG == 2 with ISM == 2)

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      SUBROUTINE ec_initlbm
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!! initialises and sets parameters of the land model
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      use global_constants_mod, only: dblp=>dp, ip

      use comatm, only: nlat,nlon, iwater
      use comemic_mod, only: fini, iatm, ilan, irunlabel, fracto, dareafac
      use comsurf_mod, only: nld, fractn, tempsgn
      use comunit, only: iuo
      use newunit_mod, only: parameterschk_id

! add use statements for call of subroutine ec_landcoverupdate, ec_landalbedo
      use comland_mod, only: albland, forestfr, fractl, albsnow, dsnow,alblbm

      use input_icemask, only: icemask


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Downscaling input variables for the sub-grid computations
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( DOWNSTS == 1 )
      use vertDownsc_mod,    only: tland_d
      use ecbilt_topography, only: nb_levls
#endif

#if ( ROUTEAU == 1 )
      USE routageEAU_mod
#endif

#if ( WISOATM == 1 )
      use comatm,      only: iwat16, nwisos
      use comland_mod, only: lbm_wiso
      use iso_nml_mod, only: land_bucket, r_soil_init, r_snow_init
#endif
      use comland_mod, only: fractl, dareas, dsnow, bmoisg, runofo, tland, albsnow, bmoismfix, dsnm, dtland, epsl, heatsnown &
                     , heatsnows, iscenland, islndstrt, lhcap, pi, radius, rdtland, rlatfus, rlatsub, rlatvap, rlhcap, rowat &
                     , tareas, tzero, alb_dat_id, forfr_dat_id

#if ( CARAIB > 1 )
      use ec_ca2lbm, only: pixarea
#endif

      implicit none

#if ( CARAIB > 1 )
#include "/home/climwork/nbouttes/caraib_dj/sources/com_18/parameter.common"
      real, parameter      :: pi_car = 2*acos(0.)
      real                 :: resg, parea,parea_tot
      real                 :: iggr,xlg,xlt,isu,clay,silt,sand,elvpix,
     &                        colour,fland
#endif

      integer(kind=ip):: inlanddat_id


      integer i,j,ija,is
#if ( WISOATM == 1 )
      integer k
#endif
      real*8  ds,db,dumwei
      real*8  phi(nlat)
      real*8  cosfi(nlat),sinfi(nlat),tanfi(nlat)
      real*8  dt,dtime,dtt,rdtime
#if ( DOWNSTS == 1 )
      INTEGER :: nb_down
#endif


      integer(kind=ip):: albsnow_dat_id, namelistland_id

      common /ec_ctstep/ dt,dtime,dtt,rdtime

      NAMELIST /landpar/ bmoismfix,dsnm,lhcap,iscenland,islndstrt

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! *** open statements of land input files:
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      open(newunit=albsnow_dat_id,file='inputdata/albsnow.dat')
      open(newunit=alb_dat_id,file='inputdata/alb.dat')
      open(newunit=forfr_dat_id,file='inputdata/forfr.dat')

      open(newunit=namelistland_id,file='namelistlbm',status='old',form='formatted')


! *** land fraction

      do j=1,nlon
        do i=1,nlat
          fractl(i,j)=1.-fracto(i,j)
        enddo
      enddo


! *** set land parameters

      bmoismfix = 0.15
      dsnm      = 1000.
      lhcap     = 2e6
      iscenland = 0
      islndstrt = 1990

!      NAMELIST /landpar/ bmoismfix,dsnm,lhcap,iscenland,islndstrt

      read(namelistland_id, NML = landpar)
      close(namelistland_id)

      write(parameterschk_id, 910) 'bmoism    =', bmoismfix
      write(parameterschk_id, 910) 'dsnm      =', dsnm
      write(parameterschk_id, 910) 'lhcap     =', lhcap
      write(parameterschk_id, 900) 'iscenland =', iscenland
      write(parameterschk_id, 900) 'islndstrt =', islndstrt

      close(parameterschk_id)



      rlhcap=1d0/lhcap

! *** grid of the surface is gaussian, read from dareafac
! *** which is initialised in initemic from file darea.dat

      ! Areas along the latitude (Gaussian grid varies only by latitude for the area)
      do i=1,nlat
        dareas(i)= dareafac(i)
      enddo

      tareas=0d0

      ! tareas seems to be the total area of the grid
      do i=1,nlat
        tareas=tareas+nlon*dareas(i)
      enddo


! *** initialisation of evaporation, soil moisture, snow cover, & land temp


      if (irunlabel.eq.0) then
        do j=1,nlon
          do i=1,nlat
            bmoisg(i,j,iwater)=bmoismfix
            dsnow(i,j,iwater)=0d0
            runofo(i,j,iwater)=0d0
            tland(i,j)=tzero+10.
          enddo
        enddo
        heatsnown=0.0
        heatsnows=0.0
      else
        open(newunit=inlanddat_id,file='startdata/inland'//fini//'.dat', form='unformatted')
        read(inlanddat_id) bmoisg(:,:,iwater),runofo(:,:,iwater), dsnow(:,:,iwater),tland
         close(inlanddat_id)
        do j=1,nlon
          do i=1,nlat
            dsnow(i,j,iwater)=min(dsnow(i,j,iwater),dsnm)
          enddo
        enddo

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|



        heatsnown=0.0
        heatsnows=0.0

      endif

#if ( DOWNSTS == 1 )
      do nb_down=1,nb_levls
         tland_d(:,:,nb_down) = tland(:,:)
      enddo
#endif

      do j=1,nlon
        do i=1,nlat
          if (fractl(i,j).lt.epsl) then
            bmoisg(i,j,:)=0d0
          endif
        enddo
      enddo

#if ( WISOATM == 1 )
! dmr&clo   wiso:3a:0002 land isotopes in the stores (land_iso_mode = 'bucket'): cold start at the namelist compositions
!           (no land isotope restart yet); non-fractionating land processes (lbm_landprecip, ec_runoff, ec_runoff_2)
      lbm_wiso = land_bucket
      if (lbm_wiso) then
        if (ilan.ne.1) then
          write(*,*) "ec_initlbm: land_iso_mode = 'bucket' needs ilan = 1 (the snow / soil split of the land evaporation", &
                     " is made with the atmospheric time step), ilan = ", ilan
          stop 1
        endif
#if ( ICEBERG == 2 || CARAIB > 0 )
        write(*,*) "ec_initlbm: land_iso_mode = 'bucket' is not available with ICEBERG = 2 or CARAIB > 0 (ROADMAP 3a)"
        stop 1
#endif
        do k = iwat16, nwisos
          bmoisg(:,:,k) = r_soil_init(k)*bmoisg(:,:,iwater)
          runofo(:,:,k) = r_soil_init(k)*runofo(:,:,iwater)
          dsnow(:,:,k)  = r_snow_init(k)*dsnow(:,:,iwater)
        enddo
        write(*,*) "ec_initlbm: land isotopes in the soil water and the snow (land_iso_mode = 'bucket'), cold start"
      endif
#endif

! *** read snow albedos

      read(albsnow_dat_id,*)
      do i=1,nlat
        read(albsnow_dat_id,*) albsnow(i)
      enddo
      close(albsnow_dat_id)

      call ec_landcoverupdate(1,1,fractl,albland,forestfr,albsnow)

![TODO] Need to check te use of this call ... what if no VECODE???
!     if (flgveg) then

! Note: the value of patmCO2 has no impact during this call.
       call veget(1,1,dtime,1d-10,280.0,fractn(1,1,nld),dareas,tempsgn(1,1,nld))

!     endif

      call ec_landalbedo(1, albland, albsnow, forestfr, dsnow,icemask, alblbm)
#if ( ROUTEAU == 0 )
      call ec_inirunoff
#else

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       WHERE(fractl.GT.epsi_lon)
         mask_lnd = 1
       ELSEWHERE
         mask_lnd = 0
       ENDWHERE

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      call ec_inirunoff_2
#endif

! *** land time step

      dtland=3600.*24./(iatm*ilan)
      rdtland=1d0/dtland

#if ( CARAIB > 1 )
!#include "/home/climwork/nbouttes/caraib_dj/sources/com_18/parameter.common"

      open(2192, file='/home/climwork/CLIM-DATA/CARAIB-DATA/PI/ecotxt.dat')

      resg=5.625
      parea_tot=0.0
      print*, 'pix_car', pix_car
      do i = 1, pix_car
         read(2192,*) iggr,xlg,xlt,isu,clay,silt,sand,elvpix,colour,fland
         !print*,  iggr,xlg,xlt,isu,clay,silt,sand,elvpix,colour,fland
         ! surface en m2
         parea=1.23604306e10*(resg*resg)*cos(xlt*pi_car/180.)
         pixarea(i)=parea*fland
         parea_tot=parea_tot+pixarea(i)
      enddo
      close(2192)
#endif


900   format(a12,1x,i6)
910   format(a12,1x,e12.5)

      return
      end

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      SUBROUTINE ec_inirunoff
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! *** initialises the runof basins
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      use global_constants_mod, only: dblp=>dp, ip

      use comatm, only: nlat, nlon

      use comland_mod, only: iocbas, ilabas, arocbas, dareas, epsl, mbasins, nbasins, fractl
      use newunit_mod, only: error_id

      use comunit, only: iuo
      use error0_mod, only: ec_error

      implicit none

      integer     i,j,k,ias,ibas,iac
      character*1 ch(nlon),space
      integer(kind=ip):: labas_dat_id

! *** Open file
      open(newunit=labas_dat_id,file='inputdata/labas.dat')

! *** computation of runoff masks

! *** asci number of small letter a
      do i=1,nlat
       do j=1,nlon
         iocbas(i,j)=0
         ilabas(i,j)=0
       enddo
      enddo

      nbasins=1
      ias=ichar('a') - 1
      iac=ichar('A') - 1
      do i=nlat,1,-1
        read (labas_dat_id,100) k,space,(ch(j),j=1,nlon)
        do j=1,nlon
          if (fractl(i,j).gt.epsl) then
            ilabas(i,j)=ichar(ch(j)) - ias
            if (ilabas(i,j).lt.1.or.ilabas(i,j).gt.26) then
              write(error_id,*) 'in lat-lon point ',i,j
              write(error_id,*) ch(j),ilabas(i,j),fractl(i,j)
              call ec_error(16)
            endif
            if (nbasins.lt.ilabas(i,j)) nbasins=ilabas(i,j)
          endif
        enddo
      enddo

      if (nbasins.gt.mbasins) then
        write(error_id,*) 'Error inirunoff: number of land basins greater than mbasins'
        STOP
      endif

      do i=nlat,1,-1
        read (labas_dat_id,100) k,space,(ch(j),j=1,nlon)
        do j=1,nlon
          if ((1d0-fractl(i,j)).gt.epsl) then
            iocbas(i,j)=ichar(ch(j)) - ias
            if (iocbas(i,j).lt.1.or.iocbas(i,j).gt.26) then
              iocbas(i,j)=0
            endif
          endif
        enddo
      enddo

      do i=nlat,1,-1
        do j=1,nlon
          if (fractl(i,j).gt.epsl) then
            ch(j)=char(ilabas(i,j)+ias)
          else
            if (iocbas(i,j).eq.0) iocbas(i,j)=ichar('0')-iac
            ch(j)=char(iocbas(i,j)+iac)
          endif
        enddo
        write(labas_dat_id,100) i-1,space,(ch(j),j=1,nlon)
      enddo

      close(labas_dat_id)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! *** computation of area of ocean runoff basins
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      do ibas=1,nbasins
        arocbas(ibas)=0.
      enddo
      do i=1,nlat
        do j=1,nlon
          if (iocbas(i,j).gt.0) then
            arocbas(iocbas(i,j))=arocbas(iocbas(i,j)) + dareas(i)*(1.-fractl(i,j))
          endif
          if (iocbas(i,j).gt.nbasins) then
            write(error_id,*) 'Error inirunof: iocbas out of range'
            STOP
          endif
        enddo
      enddo
      do ibas=1,nbasins
        if (arocbas(ibas).eq.0.) then
          write(error_id,*) 'Error inirunof: ocean basin empty ',ibas
!          STOP
        endif
      enddo


 100  format(i4,65A1)
      return
      end


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      SUBROUTINE ec_lbm(ist,jst,kst)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! *** landmodel
! *** computes: bottom moisture, snow coverage, runoff, landtemperature
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      use comatm, only: nlat, nlon
      use comemic_mod, only: iatm
      use comland_mod, only: iscenland

! add use statements for call of subroutine ec_landcoverupdate, ec_landalbedo
      use comland_mod, only: albland, forestfr, albsnow, dsnow, alblbm, fractl

      use input_icemask, only: icemask

#if ( F_PALAEO_FWF ==2 )
      use comsurf_mod, only: thi_chge
      use routageEAU_mod, only: eni, enj, mask_lnd
      use comland_mod, only: wf_ice_sheet
#endif

#if ( ISM == 1 )
#include "ismecv.com"
#endif

      integer ist,jst,kst,istep,longit

      integer :: i,j

      istep=(ist-1)*iatm+jst

#if ( ISM == 1 )
      if (flgism) then
        if (mod(istep,nstpyear).eq.1.) then
             cumprecn=0.
             cumprecs=0.
             cumprecismntemp=0.
             cumprecismstemp=0.
             bilann=0.
             bilans=0.
             dsnomeln=0.
             dsnomels=0.
        endif
      endif
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


      if (iscenland.eq.1) call ec_landcoverupdate(istep,kst,fractl,albland,forestfr,albsnow)

      call lbm_landsurface(istep)

#if ( ROUTEAU == 0 )
      call ec_runoff(istep)
#else
      call ec_runoff_2(istep)
#endif

#if ( F_PALAEO_FWF == 2 )
!nb routage de l eau venant de la fonte de la calotte  partir du
!changement d epaisseur thi_chge
! flux en m/s
      wf_ice_sheet(:,:)=0.0
          do j =1, ubound(thi_chge, dim=2)
           do i =1, ubound(thi_chge, dim=1)
             if (mask_lnd(i,j).eq.0) then
                        wf_ice_sheet(i,j)=wf_ice_sheet(i,j)+thi_chge(i,j)
             else
                        wf_ice_sheet(eni(i,j),enj(i,j))=wf_ice_sheet(eni(i,j),enj(i,j))+thi_chge(i,j)
             endif
           enddo
          enddo

#endif

      call ec_landalbedo(istep,albland, albsnow, forestfr, dsnow,icemask, alblbm)
      call ec_landcheck(istep)

      return

      end

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      SUBROUTINE lbm_landsurface(istep)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
        !! Land surface temperature and water of every grid point, one point at a time (lbm_landtemp, lbm_landprecip),
        !!   and the grid-wide sums that go with them (heat of the snow removed above dsnm, snow sent to icebergs).
        !! Replaces the whole-field ec_landtemp and ec_landprecip (wiso:3a:0000, dmr&clo).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      use comatm,        only: nlat, nlon
      use comland_mod,   only: fractl, meltheat, tland, nethfxland, landheat, dsnow, evapl, bmoisg, snowf, rainf, dareas  &
                             , betam, heatsnown, heatsnows
      use input_icemask, only: icemask

#if ( DOWNSTS == 1 )
      use vertDownsc_mod, only: tland_d, nethfxland_d, meltheat_d
#endif

#if ( ICEBERG == 2 && ISM != 2 )
      use comland_mod,   only: iceb_totl, sntoicebl
      use comemic_mod,   only: nstpyear
      use iceberg_mod,   only: iceberg_info_out_id
#endif

      implicit none

        !> Time step number in the year
      integer, intent(in) :: istep

      integer         :: i, j
      integer         :: snow_excess
      real(kind=dblp) :: snowdiff

      heatsnown = 0.0_dblp
      heatsnows = 0.0_dblp

#if ( ICEBERG == 2 && ISM != 2 )
      if (mod(istep,nstpyear).eq.1) then ! first timestep of the year
         sntoicebl = 0.
         iceb_totl = 0.
      endif
#endif

      do j=1,nlon
        do i=1,nlat

#if ( DOWNSTS == 1 )
          call lbm_landtemp(fractl(i,j),meltheat(i,j),icemask(i,j),tland(i,j),nethfxland(i,j),landheat(i,j),dsnow(i,j,:)     &
                          , tland_d(i,j,:), nethfxland_d(i,j,:), meltheat_d(i,j,:))
#else
          call lbm_landtemp(fractl(i,j),meltheat(i,j),icemask(i,j),tland(i,j),nethfxland(i,j),landheat(i,j),dsnow(i,j,:))
#endif

            ! fate of the snow above dsnm
#if ( ICEBERG < 2 )
          snow_excess = snow_excess_bucket
#elif ( ICEBERG == 2 && ISM != 2 )
          if (icemask(i,j).gt.0.9) then ! only grid cells with ice cover of more than 90%
            snow_excess = snow_excess_out
          else
            snow_excess = snow_excess_bucket
          endif
#else
          snow_excess = snow_excess_keep
#endif

          call lbm_landprecip(landheat(i,j),evapl(i,j,:),bmoisg(i,j,:),dsnow(i,j,:),snowf(i,j,:),rainf(i,j,:),meltheat(i,j)  &
                            , fractl(i,j), snow_excess, snowdiff)

          if (snowdiff.gt.0.0_dblp) then
            if (snow_excess.eq.snow_excess_bucket) then
              if (i.gt.15) then
                heatsnown=heatsnown+snowdiff*dareas(i)*fractl(i,j)/betam
              else
                heatsnows=heatsnows+snowdiff*dareas(i)*fractl(i,j)/betam
              endif
#if ( ICEBERG == 2 && ISM != 2 )
            else
              sntoicebl(i,j)=sntoicebl(i,j) + snowdiff ! accumulate excess snow over ice sheets [m]
#endif
            endif
          endif

        enddo
      enddo

#if ( ICEBERG == 2 && ISM != 2 )
      if (mod(istep,nstpyear).eq.0) then ! last timestep of the year
        do i=1,nlat
          do j=1,nlon
            iceb_totl = iceb_totl + sntoicebl(i,j)*dareas(i)*fractl(i,j)
          end do
        end do
        WRITE(iceberg_info_out_id,*) "yearly sum iceberg flux land-grid [m3]",iceb_totl,istep
      endif
#endif

      return
      end subroutine lbm_landsurface

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      SUBROUTINE lbm_landtemp(fractl,meltheat,icemask,tland,nethfxland,landheat,dsnow,tland_d, nethfxland_d, meltheat_d)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
        !! computes surface land temperature given the heat fluxes provided
        !! Revamped version to be location independant

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Downscaling input variables for the sub-grid computations
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( DOWNSTS == 1 )
!      use vertDownsc_mod, only: tland_d, nethfxland_d, meltheat_d
      use ecbilt_topography, only: rmount_virt
#endif

      use comatm, only: iwater, nwisos
      use comland_mod, only: dtland, epsl, lhcap, rdtland, rlatfus, rlhcap, rowat, tzero, betam


      implicit none

        !> Land fraction [1]
      real(kind=dblp), intent(in)    :: fractl
        !> Net heat flux over land [?]
      real(kind=dblp), intent(in)    :: nethfxland
        !> Is this an icesheet (0/1)? [1]
      real(kind=dblp), intent(in)    :: icemask

        !> Landheat flux [?]
      real(kind=dblp), intent(inout) :: landheat

        !> Heat availability to melt snow [?]
      real(kind=dblp)                 , intent(out)   :: meltheat
        !> Computed temperature over the land surface [K]
      real(kind=dblp)                 , intent(inout)   :: tland
        !> Computed temperature over the land surface [K] in downscaling
      real(kind=dblp), dimension(:), optional, intent(inout)   :: tland_d

      real(kind=dblp), dimension(:), optional, intent(inout)   :: nethfxland_d
      real(kind=dblp), dimension(:), optional, intent(inout)   :: meltheat_d


        !> Snow pack [m water equivalent], all species; read only (intent(in): with intent(out) and a strided actual
        !>   argument, gfortran -O3 does not copy the values in, wiso:3a:0000)
      real(kind=dblp), dimension(nwisos), intent(in) :: dsnow

        !> Heat availability to melt ice [?]
      real(kind=dblp) :: iceheat

        meltheat = 0._dblp

        IF (PRESENT(meltheat_d)) then
            meltheat_d(:) = 0.0d0
        ENDIF

          ! If there is a significant land surface
        if (fractl.gt.epsl) then

            !> New temperature over land is netheatflux times landheatcapacity times timestep
          tland = tland + dtland*rlhcap*(nethfxland + landheat)

          IF (PRESENT(tland_d)) then
            tland_d(:) = tland_d(:) + dtland * rlhcap * (nethfxland_d(:)+landheat)
          ENDIF

          if ((dsnow(iwater).gt.0._dblp).and.(tland.gt.tzero)) then
            meltheat = min(dsnow(iwater)/betam,lhcap*(tland-tzero))
            tland    = tland - meltheat/lhcap
            meltheat = meltheat*rdtland

            if (PRESENT(meltheat_d)) then
              meltheat_d(:) = min(dsnow(iwater)/betam,lhcap*(tland_d(:)-tzero))
              tland_d(:)    = tland_d(:) - meltheat_d(:)/lhcap
              meltheat_d(:) = meltheat_d(:)*rdtland
            endif

#if ( IMSK == 1 )
              ! dmr&afq
              ! in the case where there is an ice-sheet below the snow ...
              ! dmr&afq
            if ((icemask.gt.0.9_dblp).and.(tland.gt.tzero)) then
              iceheat    = lhcap*(tland-tzero)
              tland      = tzero
            endif
#endif


          endif ! on snow

        endif ! on fractl > epsl

      return
      end subroutine lbm_landtemp

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      SUBROUTINE lbm_landprecip(landheat,evapl,bmoisg,dsnow,snowf,rainf,meltheat,fractl,snow_excess,snowdiff)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       !! Revamped routine from ec_landprecip, for one grid point
       !! Computes the fate of precipitation on land: evaporation, snow fall and rain, snow melt, snow above dsnm
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      use comatm, only: iwater, nwisos

      use comland_mod, only:  dsnm, betam, dtland, epsl

      use newunit_mod, only: info_id
#if ( WISOATM == 1 )
      use comatm,      only: iwat16
      use comland_mod, only: lbm_wiso
#endif

      implicit none

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        !> Land heat [?]
      real(kind=dblp),                   intent(out) :: landheat
        !> Height of snow removed above dsnm [m]
      real(kind=dblp),                   intent(out) :: snowdiff
        !> Fate of the snow above dsnm: snow_excess_bucket (into bmoisg), snow_excess_out (removed, returned in snowdiff
        !>   only) or snow_excess_keep (no cap)
      integer,                           intent(in)  :: snow_excess
        !> Evaporation over land [m/s]
      real(kind=dblp), dimension(nwisos),intent(in)  :: evapl
        !> Precipitation over land, liquid and solid [m/s]
      real(kind=dblp), dimension(nwisos),intent(in)  :: rainf, snowf
        !> Heat from melting snow/ice [?]
      real(kind=dblp),                   intent(in)  :: meltheat
        !> Fraction of land in the grid cell [1]
      real(kind=dblp),                   intent(in)  :: fractl

        !> "bottom moisture" ergo water on land [m]
      real(kind=dblp), dimension(nwisos),intent(inout) :: bmoisg
        !> Thickness of snow on land [m water equivalent]
      real(kind=dblp), dimension(nwisos),intent(inout) :: dsnow

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=dblp) ::  dsnomel,dsnovap
#if ( WISOATM == 1 )
        ! ratios of the snow and of the soil water (wiso:3a:0002)
      real(kind=dblp), dimension(nwisos) :: rsn, rso
#endif

      snowdiff = 0.0_dblp
      dsnovap  = 0.0_dblp

      landheat = 0d0

          ! if on land then ...
      if (fractl.gt.epsl) then

#if ( WISOATM == 1 )
! dmr&clo   wiso:3a:0002 land isotopes (lbm_wiso): each store keeps a ratio per species, unchanged by what leaves it
!           (evaporation, melt, snow above dsnm: non-fractionating), mixed with what enters it (snow fall, rain, melt
!           water, snow above dsnm); the isotopes are the ratio times the water at the end. The evaporation leaves each
!           store at the ratio it had when the atmosphere computed the flux (ec_la2co, same step, ilan = 1).
        if (lbm_wiso) then
          rsn = lbm_wiso_ratio(dsnow)
          rso = lbm_wiso_ratio(bmoisg)
        endif
#endif

! dmr --- First deal with evaporation from the landsurface

        if (dsnow(iwater).le.0._dblp) then ! there is no snow here
          bmoisg(iwater)=bmoisg(iwater)-dtland*evapl(iwater) ! remove water that goes into evaporation from the land surface
                                                             ![NOTA] no check that there is enough water
        else                               ! there is snow
          dsnovap=dtland*evapl(iwater)           ! how much water needs to be removed in total?
          if (dsnovap.gt.dsnow(iwater)) then     ! if too much evaporation for the content of the snow layer
            bmoisg(iwater) = bmoisg(iwater)-(dsnovap-dsnow(iwater))           ! get the rest from the liquid water
            dsnow(iwater)=0._dblp                                             ! no snow left at all
          else
            dsnow(iwater) = dsnow(iwater)-dsnovap                             ! remove part of the snow layer to evaporation
          endif
        endif

#if ( WISOATM == 1 )
        if (lbm_wiso) then
          rsn = lbm_wiso_mix(rsn, dsnow(iwater), dtland*snowf)
          rso = lbm_wiso_mix(rso, bmoisg(iwater), dtland*rainf)
        endif
#endif
        dsnow(iwater) = dsnow(iwater) + dtland*snowf(iwater)                  ! add snow fall
        bmoisg(iwater)= bmoisg(iwater)+ dtland*rainf(iwater)                  ! add rainfall

! dmr the deal with the melting of snow if necessary

        if (meltheat.gt.0._dblp) then           ! there is energy to melt the snow

          dsnomel = dtland*betam*meltheat       ! potential melting

          if (dsnomel.gt.dsnow(iwater)) then    ! if too much melting for the existing snow
            landheat = landheat + (dsnomel-dsnow(iwater))/(dtland*betam) ! add the energy left to the land temperature
            dsnomel  = dsnow(iwater)                                     ! potential melting the whole layer
          endif

#if ( WISOATM == 1 )
          if (lbm_wiso) rso = lbm_wiso_mix(rso, bmoisg(iwater), rsn*dsnomel)
#endif
          dsnow(iwater) = dsnow(iwater)-dsnomel ! actually melting the snow
          bmoisg(iwater)=bmoisg(iwater)+dsnomel ! acutally adding the snow melt to the water layer

        endif ! on meltheat > 0.0

! *** if snowdepth above a thresshold, remove excessive snow
! *** artificially through the bottom moisture (or to icebergs: snow_excess, set by the caller).
! *** the heat involved is accounted for by the caller (heatsnown, heatsnows)

        if ((snow_excess.ne.snow_excess_keep).and.(dsnow(iwater).gt.dsnm)) then
          snowdiff=dsnow(iwater)-dsnm
          if (snow_excess.eq.snow_excess_bucket) then
#if ( WISOATM == 1 )
            if (lbm_wiso) rso = lbm_wiso_mix(rso, bmoisg(iwater), rsn*snowdiff)
#endif
            bmoisg(iwater)=bmoisg(iwater)+ dsnow(iwater)-dsnm
          endif
          dsnow(iwater) = dsnm
        endif

! [WARNING] This bit is not water conservative. This need to be fixed!!! dmr 2025-07-02
        if (bmoisg(iwater).lt.0d0) then
          if (bmoisg(iwater).lt.-1e-10) then
            write(info_id,*) 'bottom moisture less than zero '
            write(info_id,*) bmoisg(iwater)
            write(info_id,*) dtland*betam*meltheat
            write(info_id,*) dtland*snowf(iwater)
            write(info_id,*) dtland*rainf(iwater)
            write(info_id,*) dtland*evapl(iwater)
            write(info_id,*) dsnow(iwater),dsnovap
          endif
          bmoisg(iwater)=0d0
        endif

#if ( WISOATM == 1 )
        if (lbm_wiso) then
          dsnow(iwat16:nwisos)  = rsn(iwat16:nwisos)*dsnow(iwater)
          bmoisg(iwat16:nwisos) = rso(iwat16:nwisos)*bmoisg(iwater)
        endif
#endif

      endif ! on land ...

      return
      end subroutine lbm_landprecip

#if ( WISOATM == 1 )
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      pure function lbm_wiso_ratio(store) result(r)
        !! Isotope ratios of a land store (amounts / water); 0 for an empty store (its ratio does not matter: no weight)
        !!   wiso:3a:0002
      use comatm, only: iwater, nwisos
      real(kind=dblp), dimension(nwisos), intent(in) :: store
      real(kind=dblp), dimension(nwisos)             :: r
      if (store(iwater).gt.0.0_dblp) then
        r(:) = store(:)/store(iwater)
      else
        r(:) = 0.0_dblp
      endif
      r(iwater) = 1.0_dblp
      end function lbm_wiso_ratio

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      pure function lbm_wiso_mix(r_old, w_old, a_in) result(r)
        !! Ratios of a land store of water w_old (taken as 0 if negative) and ratios r_old after the inflow a_in (amounts
        !!   of all species, a_in(iwater) the water); unchanged without inflow. wiso:3a:0002
      use comatm, only: iwater, nwisos
      real(kind=dblp), dimension(nwisos), intent(in) :: r_old, a_in
      real(kind=dblp),                    intent(in) :: w_old
      real(kind=dblp), dimension(nwisos)             :: r
      real(kind=dblp) :: wp
      wp = max(w_old, 0.0_dblp)
      if (wp + a_in(iwater).gt.0.0_dblp) then
        r(:) = (r_old(:)*wp + a_in(:))/(wp + a_in(iwater))
      else
        r(:) = r_old(:)
      endif
      r(iwater) = 1.0_dblp
      end function lbm_wiso_mix
#endif


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      SUBROUTINE ec_runoff(istep)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! *** computes runoff from rivers and distributes it over the ocean
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      use comatm, only: nlat, nlon, nwisos, iwater
      use comemic_mod, only:
      use comunit, only: iuo


      use comland_mod, only: runofl, bmoisg, runofo, dtland, epsl, nbasins, fractl, bmoism, ilabas, dareas, iocbas, rlatfus &
                     , arocbas
#if ( WISOATM == 1 )
      use comatm,      only: iwat16
      use comland_mod, only: lbm_wiso
#endif
      use newunit_mod, only: runoff_id

      implicit none



#if ( ISM == 1 )
#include "ismecv.com"
#endif

      real*8    runo(nbasins,nwisos)
      integer   i,j,ibas,istep
#if ( WISOATM == 1 )
      real(kind=dblp) :: rso(nwisos)
#endif



! ***  total runoff for each land basin

      do ibas=1,nbasins
        runo(ibas,:)=0d0
      enddo

      do i=1,nlat
        do j=1,nlon

          runofl(i,j,:)=0.

          if (fractl(i,j).gt.epsl) then




!     MB   bmoism : Max value of bottom moisture [m]
            if (bmoisg(i,j,iwater).gt.bmoism(i,j)) then
#if ( WISOATM == 1 )
              if (lbm_wiso) rso = lbm_wiso_ratio(bmoisg(i,j,:))
#endif
!     MB dtland = 14400 [s/step]
!     MB   [m step/s]  = ( [m] - [m] ) / [s/step]
              runofl(i,j,iwater)=(bmoisg(i,j,iwater)-bmoism(i,j))/dtland
              bmoisg(i,j,iwater)=bmoism(i,j)
#if ( WISOATM == 1 )
! dmr&clo   wiso:3a:0002 overflow at the ratio of the soil water
              if (lbm_wiso) then
                runofl(i,j,iwat16:nwisos) = rso(iwat16:nwisos)*runofl(i,j,iwater)
                bmoisg(i,j,iwat16:nwisos) = rso(iwat16:nwisos)*bmoisg(i,j,iwater)
              endif
#endif
            endif ! on (bmoisg(i,j).gt.bmoism(i,j))

#if  ( ISM == 1 )
            if (flgism) then
              if(rmount_ism(i,j).LT.0.) then
              runo(ilabas(i,j))=runo(ilabas(i,j)) + dareas(i)*runofl(i,j)*fractl(i,j)
              endif
            else
#endif
!     MB land runoff to runoff, meridional convergence, fraction of land per box
!     MB  [m3 step/s] = [m3 step/s] + [m2] [m step/s] [1]
              runo(ilabas(i,j),iwater)=runo(ilabas(i,j),iwater) + dareas(i)*runofl(i,j,iwater)*fractl(i,j)
#if ( WISOATM == 1 )
              if (lbm_wiso) runo(ilabas(i,j),iwat16:nwisos)=runo(ilabas(i,j),iwat16:nwisos)                            &
                                                          + dareas(i)*runofl(i,j,iwat16:nwisos)*fractl(i,j)
#endif

#if  ( ISM == 1 )
            endif
#endif
          endif
        enddo
      enddo

#if ( ISM == 1 )
      if (flgisma) then
       if(runocorra.ne.0.) runo(23)=runocorra+runo(23)
      endif
      if (flgismg) then
       if(runocorrg.ne.0.) runo(5)=runocorrg+runo(5)
      endif
#endif
! ***  distribution of land runoff over the ocean
! correction test Greenland, Antartcica
!     runo(5)=runo(5)*0.80
!     runo(23)=runo(23)*0.80
!     MB runo is in [m3 step/s]... means m3/s per timestep!

      do i=1,nlat
        do j=1,nlon
          runofo(i,j,iwater)=0.
          if (iocbas(i,j).gt.0) then
!     MB   [m step/s] = [m3 step/s]/[m2]
             runofo(i,j,iwater)= runo(iocbas(i,j),iwater)/arocbas(iocbas(i,j))
          endif
#if ( WISOATM == 1 )
          if (lbm_wiso) then
            runofo(i,j,iwat16:nwisos)=0.
            if (iocbas(i,j).gt.0) runofo(i,j,iwat16:nwisos)= runo(iocbas(i,j),iwat16:nwisos)/arocbas(iocbas(i,j))
          endif
#endif
        enddo
      enddo
! outputs for Greenland and Antarctica
!     write(runoff_id,*) runo(5),runo(23)
! MBEDIT
! Set to 1 to use runoff output,
#if ( 0 )
      if (mod(istep,nstpyear).eq.1) then
! MBEDIT Every first atmospheric year set it zero
         runo_yearly=0.
      endif
! MBEDIT Sum it up and multiply with dtland (14400) Land timestep
! MBEDIT [m3] = [m3] + [m3 step/s]*[s/step]
      do i=1,mbasins
         runo_yearly(i) = runo_yearly(i) + (runo(i)*dtland)/1e9
      enddo
! MBEDIT

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


      if (istep.eq.1) then
         write(runoff_id,'(A)') '# Landmodel Runoff(km**3/y)'
         write(runoff_id,'(A4,26A16)') 'Year','Sibiria','Gulf.of.Alaska','Beaufort.Sea','Hudson.Bay','Greenland','Europe', &
                                       'Mediterrean','India','Okhotsk','Philippines','West.NA','Mississippi','St.Lawrence',&
                                       'West.Africa','Gulf.of.Guinea','SouthEastAfrica','SouthAfrica','Indonesia',         &
                                       'S.Australia','West.SA','Amazonas','Rio.Plata','Antarctica','Svalbard','Iceland',    &
                                       'New.Zealand'
      endif
      if (mod(istep,nstpyear).eq.0.) then
         write(runoff_id,'(I4,26E16.6)') iyear,(runo_yearly(i),i=1,mbasins)
      endif
! MBEDIT
#endif
!      if (mod(istep,nstpyear).eq.0.) then
!        write(runoff_id,'(A25,I3,2E15.5)')'runoff land (m**3/y)',
!    &   iyear,runo_yn,runo_ys
!      endif

      return
      end

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      SUBROUTINE ec_landcoverupdate(istep,kst,fractl,albland,forestfr,albsnow)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! *** updates land surface albedo and forest fraction
! *** once every 5 years
! *** starts when simulation year is equivalent to 1970, remains
! *** unchanged if simulation year exceeds 2100
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      use comatm, only: nlat, nlon
      use comemic_mod, only: iyear, nstpyear
      use comunit, only: iuo

#if ( IMSK == 1 )
      USE input_icemask
#endif

      use comland_mod, only: epsl, islndstrt, alb_dat_id, forfr_dat_id

#if ( ISM >= 2 )
      USE comland_mod, only: albkeep
#endif

      use newunit_mod, only: info_id
      implicit none

      integer i,j,is,scenyr,yr,iy
      real*8  d
      integer, intent(in) :: istep, kst
      double precision, dimension(nlat), intent(in) :: albsnow
      double precision, dimension(nlat,nlon), intent(in) :: fractl
      double precision, dimension(nlat,nlon), intent(inout) :: forestfr
      double precision, dimension(nlat,nlon,4), intent(inout)  :: albland


! *** update once every 5 years

      if (mod(istep,5*nstpyear).eq.1.and.kst.eq.1) then

! *** scenario starts in 1970, after 2100 no updates
        if (iyear.eq.0) then
          scenyr=islndstrt
        else
          scenyr=islndstrt+iyear-1
        endif
        if (scenyr.ge.1970.and.scenyr.le.2100) then

! *** read seasonal albedo data
          read(alb_dat_id,90)yr
          read(alb_dat_id,*)
          if(yr.ne.scenyr) then
            rewind(alb_dat_id)
            read(alb_dat_id,90)yr
            read(alb_dat_id,*)
          endif
          if(yr.ne.scenyr) then
            call ec_forwardfile(alb_dat_id,scenyr,yr)
          endif
          do i=1,nlat
            do j=1,nlon
              read(alb_dat_id,100)(albland(i,j,is),is=1,4)
            enddo
          enddo
#if ( LGMSWITCH == 1 || BATHY > 0 )
!! LGM or bathy
          do i=1,nlat
            do j=1,nlon
             do is=1,4
               if (albland(i,j,is).eq.0.0) albland(i,j,is)=0.18
             enddo
            enddo
          enddo
! LGM
#endif

! *** read yearly forest fraction data
          read(forfr_dat_id,90)yr
          if(yr.ne.scenyr) then
            rewind(forfr_dat_id)
            read(forfr_dat_id,90)yr
          endif
          if(yr.ne.scenyr) then
            call ec_forwardfile(forfr_dat_id,scenyr,yr)
          endif
          do i=1,nlat
            do j=1,nlon
              read(forfr_dat_id,110)forestfr(i,j)
            enddo
          enddo
        endif

#if ( IMSK == 1 )
        do i=1,nlat
          do j=1,nlon
#if ( ISM >= 2 )
            albkeep(i,j,:) = albland(i,j,:) ! afq, used for ice sheet surface mass balance
#endif
            if (icemask(i,j).gt.0.9) then
               forestfr(i,j)=0.
               albland(i,j,:) = albsnow(i) !albland(30,58,:) ! Greenland centre is at lat 30, lon 58
!              do is=1,4
!               if(albland(i,j,1).le.albsnow(i)) albland(i,j,is)=albsnow(i)
!              enddo
            endif
          enddo
        enddo
#endif

! *** addition for LGM runs ***

        d=0d0
        do j=2,25
          d=d+albland(27,j,1)
        enddo
        d=d/24.
        write(info_id,120) 'landcover update year: ',iyear,scenyr,yr,d
        do i=1,nlat
          do j=1,nlon
            if (fractl(i,j).gt.epsl) then
              do is=1,4
                if ((albland(i,j,is).lt.0.01).or.(albland(i,j,is).gt.0.99)) then
                  write(info_id,130) 'Albedo of land out of range ',i,j,is,albland(i,j,is)
                endif
              enddo
            endif
          enddo
        enddo

      endif
      call flush(info_id)
! *** FORMATS:

90    FORMAT(I5)
100   FORMAT(4(2X,F8.4))
110   FORMAT(F8.4)
120   FORMAT(A23,3I8,F12.5)
130   FORMAT(A28,3I4,F12.5)


      return
      end

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      SUBROUTINE ec_forwardfile(filenr,year,yr)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! *** forwards scenario ASCII file to find scenario year in
! *** file with number filenr
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


      use comatm, only: nlat, nlon
      use comunit, only: iuo

      use comland_mod, only: alb_dat_id
      implicit none

      integer year,filenr,yr,i,j

      do i=1,nlat
        do j=1,nlon
          read(filenr,*)
        enddo
      enddo
      read(filenr,90)yr
      if (filenr.eq.alb_dat_id) read(filenr,*)
      do while((yr.ne.year).and.(yr.lt.2100))
        do i=1,nlat
          do j=1,nlon
            read(filenr,*)
          enddo
        enddo
        read(filenr,90)yr
        if (filenr.eq.alb_dat_id) read(filenr,*)
      enddo

90    FORMAT(I5)

      return
      end

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      SUBROUTINE ec_landalbedo(istep, albland, albsnow, forestfr, dsnow,icemask, alblbm)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! *** calculates albedos as a FUNCTION of the time of the year, linearly
! *** interpolating between seasonal mean values
! *** albes  is the albedo of the earth surface, its value depends on
! ***      whether it is a land or sea point, and on whether the grid
! ***      point is snow or ice covered or not
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      use comatm, only: nlat, nlon, iwater, nwisos
      use comemic_mod, only: iday, imonth
      use comunit, only: iuo
#if ( CARAIB > 0 )
      use comsurf_mod, only: nld, fractn
      use comland_mod, only: bmoism
#endif

#if ( CARAIB > 0 )
      use ec_co2ca,  only: j_antarctique, frac_land
      use ec_ca2lbm, only: alb_transit, est_calculee, alb_for_lbm
      use ec_ca2lbm, only: tree_frac, grass_frac, desert_frac
      use comemic_mod, only: iyear
#endif

#if ( ISM >= 2 )
      USE comland_mod, only: albkeep
#endif

      implicit none

      double precision, dimension(nlat),        intent(in) :: albsnow
      double precision, dimension(nlat,nlon),   intent(in) :: forestfr
      double precision, dimension(nlat,nlon,4), intent(in) :: albland
      double precision, dimension(nlat,nlon,nwisos), intent(in) :: dsnow
      ! [TBD] dmr --- 32 64 hard coded !!!!
      real,             dimension(32,64),       intent(in) :: icemask

      double precision, dimension(nlat,nlon),   intent(out) :: alblbm

      integer i,j,id1,is1,is2,istep
      real*8  sfrac,albforfac,asnow,anosnow,snm,rsnm,snd

      integer imonth2,ktype,spv
      real*8  fracm,snowfrac,wcontent,forestfrac,hsnow
      real*8  albvec,albvecnosnow,albvecsnow,sndism

! *** reduction factor of snow albedo over forests

      albforfac=0.5

! *** if snowdepth exceeds snm, the albedo of snow
! *** is taken; for smaller snowdepths the albedo is
! *** a linear interpolation between the landalbedo
! *** and the albedo of snow.

      snm=0.05d0
      rsnm=1d0/snm

! *** interpolate between seasonal means: first shift day
! *** one to the middle of the first season (15 januari) and
! *** then interpolate

      id1=(imonth-1)*30+iday-14
      if (id1.lt.1) id1=id1+360

      is1=(id1+89)/90
      is2=is1+1
      if (is2.eq.5) is2=1

      sfrac=(id1-((is1-1)*90.+1.))/90.

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( CARAIB > 1 )
      if (iyear.gt.1) then

      do j=1,nlon
       do i=1,nlat

        if((fractn(i,j,nld).GE.frac_land).AND.(i.GE.j_antarctique))then
             forestfr(i,j) = tree_frac(i,j)

             bmoism(i,j)= (tree_frac(i,j)*0.25) + (desert_frac(i,j)*0.1)+ (grass_frac(i,j)*0.15)

        endif

       enddo
      enddo
      endif
#endif

      do j=1,nlon
        do i=1,nlat
          anosnow=albland(i,j,is1)+(albland(i,j,is2)-albland(i,j,is1))*sfrac
          asnow=albsnow(i)*(1.-forestfr(i,j))+forestfr(i,j)*albsnow(i)*albforfac
          snd=min(dsnow(i,j,iwater),snm)

          alblbm(i,j)=rsnm*((snm-snd)*anosnow+snd*asnow)

#if ( ISM >= 2 )
          albkeep(i,j,:) = alblbm(i,j) !only albkeep(:,:,nld) is used for SMB
#endif

#if ( IMSK == 1 )
          if (icemask(i,j).GT.0.9) then
            alblbm(i,j) = albsnow(30)
          endif
#endif

        enddo
      enddo

      return
      end


      SUBROUTINE ec_landcheck(istep)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! *** this routine checks the conservation of the moisture and
! *** heat budget over land
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


      use comatm, only: nlat, nlon, iwater
      use comemic_mod, only:
      use comunit, only: iuo
      use newunit_mod, only: info_id



      use comland_mod, only: dtland, epsl, lhcap, tareas, fractl, evapl, snowf, runofl, rainf, dsnow, dareas, bmoisg, runofo, tland
#if ( ICEBERG == 2 && ISM != 2 )
      use comland_mod, only: sntoicebl
#endif


      implicit none

#if ( ISM == 1 )
#include "ismecv.com"
#endif

      integer i,j,istep,nn
      real*8 iflux,oflux,labufch,labufcht,moisbdif
      real*8 totrunl,totruno,difrun,dum1,dum2,dum3
      real*8 ifluxt,ofluxt,abufchn,abufchnt,facmois,totrain,totevap
      real*8 difatmois,obufchnt,ocbud,totcor,totmois,globalmean

      real*8 bmoisgo(nlat,nlon),dsnowo(nlat,nlon),tlando(nlat,nlon)

      common /rlhch/bmoisgo,dsnowo,tlando

      if (istep.gt.2) then

         do j=1,nlon
            do i=1,nlat
               if (fractl(i,j).gt.epsl) then
!~ #if ( ISOATM >= 2 )
!~                   if (evapl(i,j,ieau)*dtland-bmoisgo(i,j)
!~      &                 -dsnowo(i,j).gt.1E-15)
!~ #else
                  if (evapl(i,j,iwater)*dtland-bmoisgo(i,j)-dsnowo(i,j).gt.1E-15) then
!~ #endif
                  write(info_id,*) 'error in evaporation over land '
                  write(info_id,*) 'latlon ',i,j,bmoisgo(i,j)+ dsnowo(i,j),evapl(i,j,iwater)*dtland

                   endif            ! evapl
                endif               ! fractl
             enddo                  ! nlat
          enddo                     ! nlon


! *** local control for surface moisture budget

!       if (.not.flgism) then
          labufcht=0d0
          do j=1,nlon
             do i=1,nlat
                if (fractl(i,j).gt.epsl) then
!~ #if ( ISOATM >= 2 )
!~                    iflux=(rainf(i,j,ieau)+snowf(i,j,ieau)
!~      &                  -evapl(i,j,ieau))
!~                    oflux=runofl(i,j,ieau)
!~                    labufch=(bmoisg(i,j,ieau)+dsnow(i,j,ieau))
!~      &                  -(bmoisgo(i,j)+dsnowo(i,j))
!~ #else
                   iflux=(rainf(i,j,iwater)+snowf(i,j,iwater)-evapl(i,j,iwater))
                   oflux=runofl(i,j,iwater)
                   labufch=(bmoisg(i,j,iwater)+dsnow(i,j,iwater))-(bmoisgo(i,j)+dsnowo(i,j))

!~ #endif
                   labufcht=labufcht + labufch*dareas(i)*fractl(i,j)
                   moisbdif=labufch - dtland*(iflux-oflux)
                   if (abs(labufch).gt.1d-6) then
                      if (abs(moisbdif/labufch).gt.1d-5) then
!dmr @-@ iceb0
!dmr --- Comment the following BLOCK IF you don't put the excess snow in runoff !
!dmr --- (to limit the output)
!dmr @-@ iceb0

                         write(info_id,*) 'error in landmoisture budget',fractl(i,j)
                         write(info_id,*) 'latlon',i,j,moisbdif,labufch,dtland*(iflux-oflux)
!~ #if ( ISOATM >= 2 )
!~                          write(info_id,*) rainf(i,j,ieau),
!~      &                        evapl(i,j,ieau),runofl(i,j,ieau)
!~                          write(info_id,*) moisbdif/dtland,
!~      &                        bmoisg(i,j,ieau)
!~ #else
                         write(info_id,*) rainf(i,j,iwater),evapl(i,j,iwater),runofl(i,j,iwater)
                         write(info_id,*) moisbdif/dtland,bmoisg(i,j,iwater)
!~ #endif
                         write(info_id,*)
                      endif
                   else
                      if (abs(moisbdif).gt.1d-5) then
                         write(info_id,*) 'error in landmoisture budget'
                         write(info_id,*) 'latlon',i,j,moisbdif,labufch
                      endif
                   endif
                endif
             enddo
          enddo
!     endif

! *** global control for runoff budget

!- Not ok in the ice sheet-coupled version because
#if ( ISM == 1 )
          if (.not.flgism) then
#endif
             totrunl=0d0
             totruno=0d0
             do j=1,nlon
                do i=1,nlat
                   if (fractl(i,j).gt.epsl) then
!~ #if ( ISOATM >= 2 )
!~                       totrunl=totrunl + runofl(i,j,ieau)*dareas(i)
!~      &                     *fractl(i,j)
!~ #else
                      totrunl=totrunl + runofl(i,j,iwater)*dareas(i)*fractl(i,j)
!~ #endif
                   endif
                   if ((1d0-fractl(i,j)).gt.epsl) then
!~ #if ( ISOATM >= 2 )
!~             totruno=totruno + runofo(i,j,ieau)*dareas(i)*(1-fractl(i,j))
!~ #else
            totruno=totruno + runofo(i,j,iwater)*dareas(i)*(1-fractl(i,j))
!~ #endif
                   endif
                enddo
             enddo

             totrunl=totrunl/tareas
             totruno=totruno/tareas
             difrun=totrunl-totruno

             if (abs(totruno).gt.0d0) then
                if (abs(difrun/totruno).gt.1d-5) then
                   write(info_id,*) 'error in runoff budget'
                   write(info_id,*) difrun,totrunl,totruno
                endif
             endif

#if ( ISM == 1 )
          endif
#endif
       endif                    ! istep.gt.2

! *** water conservation in the ISM coupled version (see forism.f)

#if ( ISM == 1 )
       if(flgism) then
          bilann=bilann+((runocorrg+cumprecismn)*dtland)
          bilans=bilans+((runocorra+cumprecisms)*dtland)
          cumprecismntemp=cumprecismntemp+(cumprecismn*dtland)
          cumprecismstemp=cumprecismstemp+(cumprecisms*dtland)
       endif
#endif

! *** local control for surface heat budget

       labufcht=0d0
       do j=1,nlon
          do i=1,nlat
             if (fractl(i,j).gt.epsl) then
!~ #if ( ISOATM >= 2 )
!~               iflux=(rainf(i,j,ieau)+snowf(i,j,ieau)-evapl(i,j,ieau))
!~               oflux=runofl(i,j,ieau)
!~ #else
              iflux=(rainf(i,j,iwater)+snowf(i,j,iwater)-evapl(i,j,iwater))
              oflux=runofl(i,j,iwater)
!~ #endif
                labufch=(tland(i,j)-tlando(i,j))*lhcap
                labufcht=labufcht + labufch*dareas(i)*fractl(i,j)
                moisbdif=labufch - dtland*(iflux-oflux)
                if (abs(labufch).gt.1d-6) then
                   if (abs(moisbdif/labufch).gt.1d-5) then
!                  write(info_id,*) 'error in landmoisture budget',
!     &                        fractl(i,j)
!                  write(info_id,*) 'latlon',i,j,moisbdif,labufch,
!     &                     dtland*(iflux-oflux)
!                    write(info_id,*) rainf(i,j),evapl(i,j),runofl(i,j)
!                  write(info_id,*) moisbdif/dtland,bmoisg(i,j)
!                  write(info_id,*)
                   endif
                else
                   if (abs(moisbdif).gt.1d-5) then
!                  write(info_id,*) 'error in landmoisture budget'
!                  write(info_id,*) 'latlon',i,j,moisbdif,labufch
                   endif
                endif
             endif
          enddo
       enddo

! *** storing of moisture variables

!~ #if ( ISOATM >= 2 )
!~        do j=1,nlon
!~           do i=1,nlat
!~              bmoisgo(i,j)=bmoisg(i,j,ieau)
!~              dsnowo(i,j)=dsnow(i,j,ieau)
!~              tlando(i,j)=tland(i,j)
!~           enddo
!~        enddo
!~ #else
       do j=1,nlon
          do i=1,nlat
             bmoisgo(i,j)=bmoisg(i,j,iwater)
             dsnowo(i,j)=  dsnow(i,j,iwater)
             tlando(i,j)=tland(i,j)
          enddo
       enddo
!~ #endif

       return
       end

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      SUBROUTINE ec_co2la
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! *** communicate between coupler and land
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      use comatm, only: nlat, nlon, iwater
      use comemic_mod, only:
      use comcoup_mod
      use comsurf_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Downscaling input variables for the sub-grid computations
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( DOWNSTS == 1 )
      use vertDownsc_mod, only: nethfxland_d, hfluxn_d, efluxn_d
#endif


      use comland_mod, only: evapoc, nethfxland, evapl, rainf, snowf, rlatfus, rowat

      implicit none

      integer i,j, k
      real*8 rsnrai
      real*8 temp_pr,temp_ev,tot_pr,pr_ism,tot_prnet

      evapoc(:,:)=evapn(:,:,noc,iwater)

      do j=1,nlon
        do i=1,nlat
          nethfxland(i,j)=(heswsn(i,j,nld)+dlradsn(i,j,nld)-ulradsn(i,j,nld)-efluxn(i,j,nld)-hfluxn(i,j,nld))

#if ( DOWNSTS == 1 )
          do k = LBOUND(nethfxland_d,dim=3), UBOUND(nethfxland_d,dim=3)
          nethfxland_d(i,j,k)=(heswsn(i,j,nld)+dlradsn(i,j,nld)-ulradsn(i,j,nld)-efluxn_d(i,j,nld,k)-hfluxn_d(i,j,nld,k))
          enddo
#endif

          evapl(i,j,:)=evapn(i,j,nld,:)
          rainf(i,j,:)=couprf(i,j,:)
          snowf(i,j,:)=coupsf(i,j,:)

        enddo
      enddo

      end


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      SUBROUTINE ec_la2co
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! *** communicate land data to coupler
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


      use comcoup_mod
      use comsurf_mod

#if ( DOWNSTS == 1 )
      use vertDownsc_mod,    only: tland_d, tsurfn_d
#endif

      use comatm, only: nlat, nlon, iwater
      use comland_mod, only: alblbm, tland, bmoisg, bmoism, dsnow
#if ( WISOATM == 1 )
      use comatm, only: iwat16, nwisos
      use comland_mod, only: lbm_wiso
#endif


      implicit none

      integer i,j

      do j=1,nlon
        do i=1,nlat
          albesn(i,j,nld)=alblbm(i,j)

          tsurfn(i,j,nld)=tland(i,j)

#if ( DOWNSTS == 1 )
          tsurfn_d(i,j,nld,:) = tland_d(i,j,:)
#endif

          abmoisg(i,j,iwater)= bmoisg(i,j,iwater)
          abmoism(i,j)       = bmoism(i,j)
          adsnow(i,j,iwater) = dsnow(i,j,iwater)
#if ( WISOATM == 1 )
! dmr&clo   wiso:3a:0002 isotopes of the stores, for the composition of the land evaporation (land_iso_mode = 'bucket')
          if (lbm_wiso) then
            abmoisg(i,j,iwat16:nwisos) = bmoisg(i,j,iwat16:nwisos)
            adsnow(i,j,iwat16:nwisos)  = dsnow(i,j,iwat16:nwisos)
          endif
#endif
!~ #endif
        enddo
      enddo
      end

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      SUBROUTINE ec_lae2co
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! *** communicate land data to coupler
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      use comcoup_mod
      use comsurf_mod

      use comatm, only: nlat, nlon, iwater
      use comland_mod, only: heatsnown, heatsnows, runofo, runofl
#if (ICEBERG == 2 && ISM != 2 )
      use comland_mod, only: sntoicebo,sntoicebl
#endif

      implicit none

      integer i,j

!~ #if ( ISOATM >= 2 && ISOOCN >= 1 )
!~       do j=1,nlon
!~         do i=1,nlat
!~           coupruno(i,j,:)=runofo(i,j,:)
!~           couprunl(i,j,:)=runofl(i,j,:)
!~         enddo
!~       enddo
!~ #elif ( ISOATM >= 2 )
!~       do j=1,nlon
!~         do i=1,nlat
!~           coupruno(i,j)=runofo(i,j,ieau)
!~           couprunl(i,j)=runofl(i,j,ieau)
!~         enddo
!~       enddo
!~ #else
      do j=1,nlon
        do i=1,nlat
          coupruno(i,j,:)=runofo(i,j,:)
          couprunl(i,j,:)=runofl(i,j,:)

#if ( ICEBERG == 2 && ISM != 2 )
          coupiceb(i,j)=sntoicebo(i,j) ! [m] excess snow summed since jan. 1st

#endif

        enddo
      enddo

      couphsnn=heatsnown
      couphsns=heatsnows
      return
      end

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      SUBROUTINE ec_wrendland
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! *** output land for start new run
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      use comland_mod, only: bmoisg, dsnow, runofo, tland
      use newunit_mod, only: newunit_id

      use comatm,      only: iwater
      use comunit,     only: iuo

      implicit none

!~ #if ( ISOATM >= 2 )
        write(newunit_id) bmoisg(:,:,iwater),runofo(:,:,iwater),dsnow(:,:,iwater),tland
!~ #else
!~       write(newunit_id) bmoisg,runofo,dsnow,tland
!~ #endif

      return
      end


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       SUBROUTINE send_caraib2lbm
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! *** communicate CARAIB data to coupler
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if (CARAIB > 0)

      use taillesGrilles,       only : nlat => iEcb, nlon => jEcb
      use global_constants_mod, only : nbdays => days_year360d_i
      use global_constants_mod, only : nbdays365 => days_year365d_i
      use ec_co2ca,             only : shift_indexes, frac_land, j_antarctique, n_pix_caraib, npft
      use ec_ca2lbm,            only : alb_transit, shift_albedo,pixarea
#if ( CYCC > 1 )
      use ec_ca2lbm,            only : stock_carbon_caraib
#endif
      use comland_mod,          only : alblbm ! alblbm is dimension(nlat,nlon)
      use comsurf_mod,          only: nld, fractn

      implicit none

!#include "comsurf.h" /* nld, fractn */
!#include "comsurf.h"

#if ( CYCC > 1 )
#include "/home/climwork/nbouttes/caraib_dj/sources/com_18/parameter.common"
#include "/home/climwork/nbouttes/caraib_dj/sources/com_18/inidata.common"
#include "/home/climwork/nbouttes/caraib_dj/sources/com_18/biomasse.common"
#include "/home/climwork/nbouttes/caraib_dj/sources/com_18/ecoin.common"


      integer(kind=4)      :: i,j,ij
      real, parameter  :: pi = 2*acos(0.)
      !real, parameter  :: pi = 3.141592654
      real                 :: pixcveg, glcveg
!     &                        parea_tot
      real                 :: fcveg, pixcsoil, glcsoil
!     &                       iggr,xlg,xlt,isu,clay,silt,sand,elvpix,colour

      ij = 0

      !open(2192, file='/home/climwork/CLIM-DATA/CARAIB-DATA/PI/ecotxt.dat')

      glcveg=0.0
      glcsoil=0.0
      !fland=1
      !resg=5.625
      !parea_tot=0.0
      print*, 'pix_car', pix_car
      do i = 1, pix_car
         !read(2192,*)  iggr,xlg,xlt,isu,clay,silt,sand,elvpix,colour,fland
         !print*,  iggr,xlg,xlt,isu,clay,silt,sand,elvpix,colour,fland
         ! surface en m2
         !parea=1.23604306e10*(resg*resg)*cos(xlt*pi/180.)
         !pixarea=parea*fland
         !parea_tot=parea_tot+pixarea

         !carbon in vegetation from caraib
         pixcveg=0.0
         do j = 1, npft
           fcveg=frac_cpl(i,j)*ybiomf_tot(i,j) !ybiomf(j) !dans biom.f xcveg(j)->ybiomf dans biomass
           pixcveg=pixcveg+fcveg
         end do
         glcveg=glcveg+pixcveg*pixarea(i)

         !carbon in soil from caraib (soil+litter) carbon (xcsoil)
         !print*, 'ysoilr_tot', ysoilr_tot(i)
         pixcsoil=ysoilr_tot(i) !(0) !xcsoil=ysoilr dans biomass
         glcsoil=glcsoil+pixcsoil*pixarea(i)
      enddo

      !carbon total
      !print*, 'surface total caraib', parea_tot
      print*,'carbon caraib dans landmodel (GtC) soil, veg', glcsoil*1e-15, glcveg*1e-15
      stock_carbon_caraib = glcsoil+glcveg


      !close(2192)


        do i=1,nlat
         do j=1,nlon
           if(fractn(i,j,nld).GE.frac_land.AND.i.GE.j_antarctique) then
             ij = ij + 1
               !stock_carbon_caraib(i,j) = sum(ycar_ini(1:26,1:3,ij))
               !write(*,*) "carbon caraib", i,j, stock_carbon_caraib(i,j)
!~             alb_transit(i,j,:) = albmth_t(1:365,ij)
!~             write(*,*) "i ==",i, "j ==", j
!~             write(*,*) "alb_transit == ", alb_transit(i,j,365)
           endif
        enddo
      enddo
      !write(*,*) "stop in landmodel0.f"
      !read(*,*)

#endif

!~       call shift_albedo()

#endif

      return
      end subroutine send_caraib2lbm


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      SUBROUTINE frac_acc
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! *** communicate CARAIB data to coupler
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if (CARAIB > 1)

      use taillesGrilles,  only : nlat => iEcb, nlon => jEcb
      use ec_co2ca,        only : frac_land
      use ec_co2ca,        only : j_antarctique
      use ec_ca2lbm,       only : tree_frac, grass_frac, desert_frac, veget_frac
      use comland_mod,     only : forestfr
#if ( COMATM == 1 )
      use comsurf_mod,     only : nld, fractn
#endif
      implicit none

#if ( COMATM == 0 )
#include "comsurf.h" /* nld, fractn */
#endif


#include "parameter.common"
#include "vegfr.common"

      integer(kind=4) :: i,j,ij

      ij = 0

      do i=1,nlat
       do j=1,nlon
         if(fractn(i,j,nld).GE.frac_land.AND.i.GE.j_antarctique) then
          ij = ij + 1
          grass_frac(i,j) = sum(caraib_frc(1:11,ij))
          tree_frac(i,j) = sum(caraib_frc(12:26,ij))
          veget_frac(i,j)=grass_frac(i,j)+tree_frac(i,j)

          desert_frac(i,j) = 1 - (grass_frac(i,j)+tree_frac(i,j))
!~           write(*,*) "tree ==",tree_frac(i,j),"grass ==",grass_frac(i,j)
!~           write(*,*) "desert ==", desert_frac(i,j)
!~           read(*,*)
         endif
       enddo
      enddo

#endif

      return
      end subroutine frac_acc


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      SUBROUTINE d18O_iso
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ***
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( OXYISO > 0 )

      use taillesGrilles, only : nlat => iEcb, nlon => jEcb
      use ec_co2ca,       only : frac_land, j_antarctique, dhu, dom
      use ec_iso2co,      only : rhu_photo, Tnew, temp_photo, epsieq, d18Oleaf, Ci, GPPO2, epsil, PHOo2, Pmehler, O2prod     &
                               , fphoto, f_darkleaves, a_darksoil, alpha18ter, d18Oter, diffO2, m_18o2, gpp_C3, gpp_C4,      &
                                 d18Oleaf_PFT2, epsieq_PFT2, write_d18O_iso

      use comsurf_mod, only: nld, fractn

      implicit none

#include "/home/climwork/textier/caraib-git/com_18/parameter.common"
#include "/home/climwork/textier/caraib-git/com_18/d18Oatm.common"
#include "../../../caraib_dj/sources/com_18/d18Oatm_PFT.common"
#include "/home/climwork/textier/caraib-git/com_18/gpp_values.common"
#include "/home/climwork/textier/caraib-git/com_18/delta.common"

      integer(kind=4)      :: i,j,ij

      ij = 0

      do i=1,nlat
       do j=1,nlon

         if(fractn(i,j,nld).GE.frac_land.AND.i.GE.j_antarctique) then

           ij = ij + 1

           rhu_photo(i,j,:) = rhu_foto(1:365,ij)
           Tnew(i,j,:) = T_new(1:365,ij)
           temp_photo(i,j,:) = temp_foto(1:365,ij)
           epsieq(i,j,:) = epsi_eq(1:365,ij)
           Ci(i,j,:) = C_i(1:365,ij)
           GPPO2(i,j,:) = GPP_O2_g(1:365,ij)
           epsil(i,j,:) = eps_il(1:365,ij)
           PHOo2(i,j,:) = PHO_O2(1:365,ij)
           Pmehler(i,j,:) = P_mehler(1:365,ij)
           O2prod(i,j,:) = O2_prod(1:365,ij)
           fphoto(i,j,:) = f_foto(1:365,ij)
           f_darkleaves(i,j,:) = f_darkl(1:365,ij)
           a_darksoil(i,j,:) = a_darks(1:365,ij)
           alpha18ter(i,j,:) = alpha_18_ter(1:365,ij)
           d18Oleaf(i,j,:) = delta_leaf(1:365,ij)
           d18Oter(i,j,:) = d18O_ter(1:365,ij)
           diffO2(i,j,:) = diff_O2(1:365,ij)
           m_18o2(i,j,:) = m18o2(1:365,ij)
           gpp_C3(i,j,:) = gpp_fracC3(1:365,ij)
           gpp_C4(i,j,:) = gpp_fracC4(1:365,ij)
           epsieq_PFT2(i,j,:) = epsieq_PFT(1:26,ij)
           d18Oleaf_PFT2(i,j,:,:) = d18Oleaf_PFT(1:365,1:26,ij)

         endif

       enddo
      enddo

      call write_d18O_iso

#endif

      return
      end subroutine d18O_iso


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      SUBROUTINE D17_iso
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ***
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( D17ISO > 0 )

      use taillesGrilles, only : nlat => iEcb, nlon => jEcb
      use ec_co2ca,       only : frac_land, j_antarctique, dhu, dom, npft
      use ec_iso2co,      only : alpha17_2, lambda_t2, R18Olw2, R17Olw2, d17Oleaf2, d17Oter_local2, d17Oleaf_PFT2, d17Otl_PFT2,  &
                                 R18ter2, R17ter2, d17Oter, R18Ogw2, R17Ogw2, write_d17O_iso
      use comsurf_mod, only    : nld, fractn

      implicit none

#include "../../../caraib_dj/sources/com_18/parameter.common"
#include "../../../caraib_dj/sources/com_18/d17O.common"
#include "../../../caraib_dj/sources/com_18/d17O_PFT.common"

      integer(kind=4)      :: i,j,ij

      ij = 0

      do i=1,nlat
       do j=1,nlon

         if(fractn(i,j,nld).GE.frac_land.AND.i.GE.j_antarctique) then

           ij = ij + 1

           alpha17_2(i,j,:) = alpha_17_ter(1:365,ij)
           lambda_t2(i,j,:) = lambda_transpi(1:365,ij)
           R18Olw2(i,j,:) = R18Olw(1:365,ij)
           R17Olw2(i,j,:) = R17Olw(1:365,ij)
           d17Oleaf2(i,j,:) = d17Oleaf(1:365,ij)
           d17Oter_local2(i,j,:) = d17Oter_local(1:365,ij)
           d17Oleaf_PFT2(i,j,:,:) = d17Oleaf_PFT(1:365,1:26,ij)
           d17Otl_PFT2(i,j,:,:) = d17Otl_PFT(1:365,1:26,ij)
           R18ter2(i,j,:) = R18ter(1:365,ij)
           R17ter2(i,j,:) = R17ter(1:365,ij)
           d17Oter(i,j,:) = d17O_ter(1:365,ij)
           R18Ogw2(i,j,:) = R18Ogw(1:365,ij)
           R17Ogw2(i,j,:) = R17Ogw(1:365,ij)


         endif

       enddo
      enddo

      call write_d17O_iso

#endif

      return
      END SUBROUTINE D17_iso

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      SUBROUTINE dDwax_iso
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ***
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( WAXISO > 0 )

      use taillesGrilles, only : nlat => iEcb, nlon => jEcb
      use ec_co2ca,       only : frac_land, j_antarctique, dhu, dom, npft
      use ec_iso2co,      only : epsiplus_2, dDlw_2, dDblw_2, dDwax_2, write_dDwax_iso

      use comsurf_mod, only    : nld, fractn

      implicit none

#include "../../../caraib_dj/sources/com_18/parameter.common"
#include "../../../caraib_dj/sources/com_18/dDwax.common"
#include "../../../caraib_dj/sources/com_18/parameter_dDwax.common"

      integer(kind=4)      :: i,j,ij

      ij = 0

      do i=1,nlat
       do j=1,nlon

         if(fractn(i,j,nld).GE.frac_land.AND.i.GE.j_antarctique) then

           ij = ij + 1

           epsiplus_2(i,j,:) = epsiplus(1:365,ij)
           dDlw_2(i,j,:) = dDlw(1:365,ij)
           dDblw_2(i,j,:) = dDblw(1:365,ij)
           dDwax_2(i,j,:,:) = dDwax(1:365,1:26,ij)

         endif

       enddo
      enddo

      call write_dDwax_iso

#endif

      return
      end subroutine dDwax_iso

      end module landmodel_mod

