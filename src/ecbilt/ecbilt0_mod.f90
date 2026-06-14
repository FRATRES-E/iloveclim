!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:39 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:39 CET 2009


!==============================================================================
! MODULE ecbilt0_mod
! Main time-stepping routines for ECBilt QG3L atmospheric model.
! Converted from ecbilt0.f (F77) to F90 module.
!==============================================================================
module ecbilt0_mod

  use global_constants_mod, only: dblp=>dp, ip

  implicit none
  private

  public :: ec_ecbilt, ec_update, ec_atmstate, ec_wrendphy, ec_wrenddyn

contains


!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_ecbilt(ist,jst)
!-----------------------------------------------------------------------
! *** this routine performs one timestep of ECBilt which is a
! *** QG3L atmospheric model with simple physical parameterizations
! ***
! *** 6 april 1999 KNMI, De Bilt
! ***
! *** joint project Hugues Goosse
! ***               Rein Haarsma
! ***               Theo Opsteegh
! ***               Thijs van Reenen
! ***               Michiel Schaeffer
! ***               Frank Selten
! ***               Xueli Wang
! ***               Nanne Weber
!-----------------------------------------------------------------------

      USE atmdyn_mod, only: ec_forcdaily, ec_forward

      use comatm, only: iadyn, iaphys
      use comdyn, only: iartif
      use comemic_mod, only: iyear, iatm, day
      use atmphys_mod, only: ec_atmphyszero, ec_sensrad, ec_fortemp, ec_meantemp

#if ( ISOATM == 0 )      
      use ec_convec_mod, only: ec_convec_two
#endif

#if ( ISOATM >= 1 )
      use atmmois_mod, only: ec_moisture=>ec_moisture_r332, &
                            ec_convec=>ec_convec_r332
#else
      use atmmois_mod, only: ec_moisture
#endif

      implicit none

      integer ist,jst,istep
      logical :: success

! *** atmospheric physics (in file atmphys.f)

      istep=(ist-1)*iatm+jst
#if ( ISM == 1 )
      if (flgism.AND.(mod(istep,nstpyear).eq.1)) call ec_topo
#endif
      call ec_atmout(istep)
      call ec_checks(istep)

      if (iaphys.eq.1) then

        call ec_atmphyszero
        call ec_sensrad
!        call ec_tracer
        success = ec_moisture()
#if ( ISOATM >= 1 )
        success=ec_convec()
#else
        call ec_convec_two
#endif
        call ec_fortemp
        call ec_meantemp

      endif

! *** atmospheric dynamics (in file atmdyn.f)

      if (iadyn.eq.1) then
        if (iartif.eq.1) call ec_forcdaily
        call ec_forward
      endif

      return
      end subroutine ec_ecbilt


!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_update(ist,jst)
!--------------------------------------------------------------------------------
! ***
! *** This routine updates the model date, incoming solar radiation and
! *** the atmospheric state and does outputting and checking
! ***
!-------------------------------------------------------------------------------

      USE comatm,          only:
      use comphys,         only:
      use comemic_mod,     only: iyear, iatm
      use ECBiltTimer_mod, only: ec_mdldate
      use atmrad_mod,      only: ec_solar
      use atmphys_mod,      only: ec_ghgupdate, ec_vortfor

      implicit none

      logical :: success

!~ [DEPRECATED]
!~ [DEPRECATED]

      integer ist,jst,istep

      istep=(ist-1)*iatm+jst

      call ec_mdldate(istep)
      call ec_ghgupdate(istep)
      success = ec_solar(istep)
      call ec_atmstate
      call ec_vortfor

      end subroutine ec_update

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_checks(istep)
!-----------------------------------------------------------------------
! *** this routine performs checks on the model formulation
!-----------------------------------------------------------------------

      use comemic_mod, only: iatm, day

      implicit none

!~ [DEPRECATED]
!~ [DEPRECATED]

      integer i,j,istep

      if ( mod(istep,iatm) .eq. 0) then
        call ec_testecbilt(istep)
      endif

!      call ec_moischeck(istep)
!      call ec_heatcheck(istep)

      call flush(99)
      return
      end subroutine ec_checks

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_testecbilt(istep)
!-----------------------------------------------------------------------
! *** testing if model variables are inside prescribed ranges
!-----------------------------------------------------------------------

      use comatm,          only: dt, nlat, nlon
      use comphys,         only: uv10, tzero
      use comemic_mod,     only: iyear, iatm, day
      use comsurf_mod,     only: eflux, hflux, tsurf
      use newunit_mod,     only: book_id, error_id, info_id
      use ipcc_output_mod, only: moc, tmc, tmc0, tsurfmean, cland, thex
      use atmphys_mod,     only: ec_globalmean
      use error0_mod,      only: ec_error

      implicit none

      integer i,j,istep
      real(kind=dblp)  dum1,dum2
      character*12 chts
!~      real(kind=dblp)  moc,tmc,tmc0,cland,thex
!~      common/IPCC_out2/moc,tmc,tmc0,tsurfmean,cland,thex

!      dum1=ec_globalmean(ulrads)
!      dum2=ec_globalmean(dlrads)
!      write(*,*) istep,dum1,dum2

      do j=1,nlon
        do i=1,nlat

          if (uv10(i,j).gt.60) then
            write(info_id,*) 'uv10 out of range'
            write(info_id,*) i,j,uv10(i,j)
          endif
          if (tsurf(i,j).gt.400.or.tsurf(i,j).lt.150) then
            write(info_id,*) 'tsurf out of range in test'
            write(info_id,*) i,j,tsurf(i,j)
          endif
          if (eflux(i,j).gt.2000.or.hflux(i,j).gt.2000) then
            write(info_id,*) 'surface flux out of range in test'
            write(info_id,*) i,j,eflux(i,j),hflux(i,j),tsurf(i,j)
          endif

        enddo
      enddo

      tsurfmean=ec_globalmean(tsurf)-tzero

      write(chts,900) tsurfmean
 900  format(E12.5)
      if (chts(1:3).eq.'nan') call ec_error(99)

      write(book_id,110) iyear,int((day+0.5*dt)/(iatm*dt))+1,tsurfmean
      call flush(book_id)

      if (tsurfmean.gt.40.0 .or. tsurfmean.lt.-10.0) then
        write(error_id,*) 'mean surface temperature ',tsurfmean
        call ec_error(3)
      endif

  110 format(i8,i8,f7.2)
  100 format(f7.2)

      return
      end subroutine ec_testecbilt


!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_atmstate
!-----------------------------------------------------------------------
! *** calculates atmospheric fields from potential vorticity
! *** field, mean atmospheric temperatures and the moisture field
!-----------------------------------------------------------------------

      use atmdyn_mod, only: ec_ddl, ec_diver, ec_divwin, &
                           ec_forward, ec_geowin, ec_ggtosp, ec_omega3, &
                           ec_psitogeo, ec_qtopsi, ec_sptogg, ec_totwind

      use atmmois_mod, only: ec_moisfields
      use atmphys_mod, only: ec_dyntemp, ec_tempprofile, ec_totwind10, &
                            ec_cloud

      implicit none

      logical :: success

      call ec_qtopsi
      call ec_psitogeo
      call ec_dyntemp
      call ec_tempprofile
      call ec_geowin
      call ec_omega3
      call ec_diver
      call ec_divwin
      call ec_totwind
      call ec_totwind10
      success = ec_moisfields()
      call ec_cloud

      end subroutine ec_atmstate

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_wrenddyn
!-----------------------------------------------------------------------
! *** output atmosphere for start new run
!-----------------------------------------------------------------------

      use comdyn,      only: for, qprime

      use comatm, only: iadyn, iaphys
      use comdyn, only: iartif
      use newunit_mod, only: newunit_id
      
      implicit none

      write(newunit_id) qprime,for

      return
      end subroutine ec_wrenddyn

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_wrendphy
!-----------------------------------------------------------------------
! *** output atmosphere for start new run
!-----------------------------------------------------------------------

      use comatm, only: iadyn, iaphys
      use comdyn, only: iartif
      use comemic_mod, only: iyear
      use comphys,     only: rmoisg, temp0g, tempm, torain, tosnow
      use comsurf_mod, only: tsurfn
      use comatm,      only: iwater
      use comemic_mod, only: iyear
      use newunit_mod, only: newunit_id
#if ( CYCC >= 2 )
      use carbone_co2, only: PA_C
#endif
#if ( ISOATM >= 1 )
      use write_isoatm_mod, only: write_isoatm
#endif

      implicit none

      real(kind=dblp) :: PGACO2, PCO2ref
      common /lo2atm/ PGACO2, PCO2ref

      write(newunit_id) tsurfn, tempm, temp0g

!~ #if ( ISOATM >= 1 )
      write(newunit_id) rmoisg(:,:,iwater), &
                        torain(:,:,iwater), &
                        tosnow(:,:,iwater)
!~ #else
!~       write(newunit_id) rmoisg,torain,tosnow
!~ #endif
#if ( CYCC == 1 )
      if (flgloch)  write(newunit_id) PGACO2,PCO2ref
#elif ( CYCC >= 2 )
      write(newunit_id) PA_C
#endif
#if ( ISOATM >= 1 )
      use write_isoatm_mod, only: write_isoatm
#endif
#if ( ISOATM >= 1 )
      call write_isoatm
#endif
      return
      end subroutine ec_wrendphy

#if ( ISM == 1 )
!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_topo
!-----------------------------------------------------------------------
! *** replace topography by ISM topography
!-----------------------------------------------------------------------

      use comatm, only: iadyn, iaphys

      use comdyn, only: iartif
      use c06_fftw_mod, only: C06FPF, C06FQF                

      use comemic_mod, only: iyear
      use comcoup_mod                ! external coupler — only: not available
      use comsurf_mod, only: epss, fractn, nld, rmountn, thi_chge
      use global_constants_mod, only: dblp=>dp, ip
      use newunit_mod, only: berg_dat_id
      
#if ( NC_BERG >= 1 )      
      use ncio, only: nc_read
#endif    
#if ( NC_BERG == 2 )
      use update_clio_bathy_tools, only: la_date
#endif

      implicit none

      real(kind=dblp) asum,spv
      integer i,j,i1,j1,ii,jj

      integer k1,k2,k,l,m,n,ifail,nn
      real(kind=dblp)  pigr4,dis,dif,rll,ininag(nlat,nlon)
      real(kind=dblp)  r1,a,b,c,d,e,sqn,rsqn
      real(kind=dblp)  rnorm,rh0,dd
      real(kind=dblp)  agg(nlat,nlon), agg1(nlat,nlon), agg2(nlat,nlon)
      real(kind=dblp)  fw(nsh2),fs(nsh2),fors(nsh2,nvl), fmu(nlat,2)
      real(kind=dblp)  forw(nsh2,nvl),wsx(nsh2),areafac

#if ( NC_BERG >= 1 )      
      real(kind=dblp) aggT(nlon,nlat)
#endif
#if ( NC_BERG == 2 )
      character*30 name_file
      character(len=5) :: charI
#endif

      integer(kind=ip):: topographie_id
      
!-pour changer topo (dyn+thermo), mettre le nouveau champ dans rmount_ism(verifier les flgism):
        spv=0.
        do i=1,nlat
          do j=1,nlon
!         rmount_ism(i,j)=
          if ((flgisma.AND.(i.lt.15)).OR.(flgismg.AND.(i.gt.15))) then
             if (rmount_ism(i,j).GE.spv) then
              rmountn(i,j,nld)=rmount_ism(i,j)
              rmount(i,j)=rmount_ism(i,j)
             endif
          endif
          if (rmountn(i,j,nld).lt.0e0_dblp) rmountn(i,j,nld)=0e0_dblp
          if (fractn(i,j,nld).lt.epss) then
            rmountn(i,j,nld)=0e0_dblp
            qmount(i,j)=0e0_dblp
          else
            qmount(i,j)=rmountn(i,j,nld)
          endif
          enddo
        enddo

        do j=1,nlon
        do i=1,nlat
          asum=0e0_dblp
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
            qmount(i,j)=asum/9e0_dblp
          enddo
        enddo

         open(newunit=topographie_id &
                 ,file='outputdata/globals/topographie',status='unknown')
         do i=1,nlat
           write(topographie_id,'(1P64E10.3)') (rmountn(i,j,nld),j=1,nlon)
         enddo
         close(topographie_id)

!        read (iuo+1) nshm, ll

! *** real parameters


          pigr4=4.e0_dblp*pi
          rl1=1.0e0_dblp/rrdef1**2
          rl2=1.0e0_dblp/rrdef2**2
          relt1=max(0.0e0_dblp,rl1/(trel*pigr4))
          relt2=max(0.0e0_dblp,rl2/(trel*pigr4))
          dis=max(0.0e0_dblp,1.0e0_dblp/(tdis*pigr4))
          rll=dble(ll(nsh))
          dif=max(0.0e0_dblp,1.0e0_dblp/(tdif*pigr4*(rll*(rll+1))**idif))

! *** zonal derivative operator

         k2=0
      do m=0,nm
        k1=k2+1
        k2=k2+nshm(m)
        do k=k1,k2
          rm(k)=dble(m)
        enddo
      enddo

! *** laplace/helmholtz direct and inverse operators

      do j=0,5
        rinhel(1,j)=0.0e0_dblp
      enddo

      diss(1,1)=0.0e0_dblp
      diss(1,2)=0.0e0_dblp

      do k=2,nsh
        r1=dble(ll(k)*(ll(k)+1))
        a=-r1-3.0e0_dblp*rl1
        b=-r1-3.0e0_dblp*rl2
        c=-r1-rl1
        d=-r1-rl2
        e=a*d+b*c
        rinhel(k,0)=-r1
        rinhel(k,1)=-1.0e0_dblp/r1
        rinhel(k,2)= d/e
        rinhel(k,3)= b/e
        rinhel(k,4)=-c/e
        rinhel(k,5)= a/e
        diss(k,2)=dis*r1
        diss(k,1)=-dif*r1**idif
      enddo

      do j=0,5
        do k=1,nsh
          rinhel(k+nsh,j)=rinhel(k,j)
        enddo
      enddo

      do j=1,2
        do k=1,nsh
          diss(k+nsh,j)=diss(k,j)
        enddo
      enddo

! *** initialization of coefficients for fft


      do j=1,nlon
        do i=1,nlat
          ininag(i,j)=1.0e0_dblp
        enddo
      enddo

      ifail=0
      call c06fpf (nlat,nlon,ininag,'i',trigd,wgg,ifail)

      ifail=0
      call c06fqf (nlat,nlon,ininag,'i',trigi,wgg,ifail)

! *** orography and dissipation terms

! *** fmu(i,1): sin(phi(i))
! *** fmu(i,2): 1-sin**2(phi(i))



         rnorm=1.0e0_dblp/sqrt(3.0e0_dblp*nlon)
      do i=1,nlat
        fmu(i,1)=rnorm*pp(i,2)
        fmu(i,2)=1.e0_dblp-fmu(i,1)**2
      enddo

#if ( NC_BERG == 1 )      
         call nc_read("inputdata/berg.nc","agg1",aggT)
         agg1 = transpose(aggT)
#elif ( NC_BERG == 2 )
      write(charI,'(I5.5)'), la_date
      name_file=='inputdata/berg_'//trim(charI)//'k.nc'
      call nc_read(name_file,"agg1",aggT)
      agg1 = transpose(aggT)

#if ( F_PALAEO_FWF == 2 )
      call nc_read(name_file,"thi_chge",thi_chge)
      thi_chge = transpose(thi_chge)
      print(*,*) '!!!!!!!!!!!!!!!!!!!!!'
      print(*,*) 'thi change', thi_chge
#endif

#else
         rewind (berg_dat_id)
         read (berg_dat_id) agg1
#endif         

         rh0=max(0.0e0_dblp,0.001e0_dblp/h0)
      do j=1,nlon
        do i=1,nlat
          if ((flgisma.AND.(i.lt.15)).OR.(flgismg.AND.(i.gt.15))) then
          if (rmount_ism(i,j).GE.spv) agg1(i,j)=rmount_ism(i,j)
          endif
        enddo
      enddo

         do j=1,nlon
          do i=1,nlat
!          agg(i,j)=fmu(i,1)*agg1(i,j)*rh0
          agg(i,j) = agg1(i,j)*rh0
          rmount(i,j)=agg1(i,j)
          if (rmount(i,j).lt.0e0_dblp) rmount(i,j)=0e0_dblp
         enddo
        enddo

! *** surface dependent friction

      lgdiss=((addisl.gt.0.0).or.(addish.gt.0.0))

      call ec_ggtosp (agg,orog)
      call ec_ddl (orog,ws)
      call ec_sptogg (ws,dorodl,pp)
      call ec_sptogg (orog,dorodm,pd)

      if (lgdiss) then

#if ( NC_BERG == 1 )      
       call nc_read("inputdata/berg.nc","agg2",aggT)
       agg2 = transpose(aggT)
#elif ( NC_BERG == 2 )
      write(charI,'(I5.5)'), la_date
      name_file=='inputdata/berg_'//trim(charI)//'k.nc'
      call nc_read(name_file,"agg2",aggT)
      agg2 = transpose(aggT)
#else
       read (berg_dat_id) agg2
#endif

        do j=1,nlon
          do i=1,nlat
            agg(i,j)=1.0e0_dblp+addisl*agg2(i,j)+ &
                     addish*(1.0e0_dblp-exp(-0.001e0_dblp*agg1(i,j)))
          enddo
        enddo

        call ec_ggtosp (agg,ws)
        call ec_ddl (ws,wsx)

        call ec_sptogg (ws,rdiss,pp)
        call ec_sptogg (wsx,ddisdx,pp)
        call ec_sptogg (ws,ddisdy,pd)

        dd=0.5e0_dblp*diss(2,2)
        do j=1,nlon
          do i=1,nlat
            ddisdx(i,j)=dd*ddisdx(i,j)/fmu(i,2)
            ddisdy(i,j)=dd*ddisdy(i,j)*fmu(i,2)
          enddo
        enddo

      endif
         end
#endif

end module ecbilt0_mod
