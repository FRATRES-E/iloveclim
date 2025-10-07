!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:26 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:26 CET 2009

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_initcoup
!-----------------------------------------------------------------------
! *** initialises the coupling of ocean with atmosphere model
!-----------------------------------------------------------------------

      USE comatm
      use comdyn
      use comphys
      use comcoup_mod
      use comsurf_mod
      use comunit
      use comemic_mod, only: fini, irunlabel, dareafac
#if ( ICEBERG == 2 && ISM !=2 )
     & ,nstpyear, ntstep
#endif

#if ( WRAP_EVOL == 1 )

      use init_emic_wrapper, only: mozaic_f
#endif

#if (ISOATM >= 1 )
       USE iso_param_mod, ONLY : ieau, ieau18,
     &                           rsmow, neauiso
       USE isoatm_mod, ONLY : datmini,isolbm_restart
#endif

#if ( BATHY >= 2 )
       use update_clio_bathy_tools, only: la_date
#endif

      use global_constants_mod, only: dblp=>dp, ip

      USE OCEAN2COUPL_COM, only: ec_oc2co

      use landmodel_mod, only: ec_lae2co, ec_la2co


      implicit none

      integer i,j,g05dyf,nn,ijo,ija,k
      real*8  areafac,dumwei,ds,db

      integer(kind=ip):: incoup_dat_id, mozaic_w_id

#if ( BATHY >= 2 )
      character*30 name_file
      character(len=5) :: charI
#endif

      if (irunlabel.eq.0) then
        do j=1,nlon
          do i=1,nlat
            clhesws(i,j) = 0.
            clhesw0(i,j) = 0.
            clhesw1(i,j) = 0.
            clhesw2(i,j) = 0.
            clulrad0(i,j)= 0.
            clulrad1(i,j)= 0.
            clulrad2(i,j)= 0.
            clulrads(i,j)= 0.
            cldlrads(i,j)= 0.
            clhflux(i,j) = 0.
            cleflux(i,j) = 0.
c~ #if (ISOATM >= 1 )
c~             clevap(i,j,:)  = 0.
c~ #else
            clevap(i,j,:)  = 0.
c~ #endif
            cohesws(i,j) = 0.
            cohesw0(i,j) = 0.
            cohesw1(i,j) = 0.
            cohesw2(i,j) = 0.
            coulrad0(i,j)= 0.
            coulrad1(i,j)= 0.
            coulrad2(i,j)= 0.
            coulrads(i,j)= 0.
            codlrads(i,j)= 0.
            cohflux(i,j) = 0.
            coeflux(i,j) = 0.
c~ #if ( ISOATM >= 1 )
c~             coevap(i,j,:)  = 0.
c~ #else
            coevap(i,j,:)  = 0.
c~ #endif
            sumohfx(i,j) = 0.
            sumoswr(i,j) = 0.
            sumihfx(i,j) = 0.
            sumiswr(i,j) = 0.
c~ #if ( ISOATM >= 2 && ISOOCN >= 1 )
c~             sumofwf(i,j,:) = 0.
c~             sumisno(i,j,:) = 0.
c~ #else
            sumofwf(i,j,:) = 0.
            sumisno(i,j,:) = 0.
c~ #endif
            sumicof(i,j) = 0.
            sumuv10(i,j) = 0.
            sumpress(i,j)= 0.
            winstua(i,j) = 0.
            winstva(i,j) = 0.
            sumtx(i,j)   = 0.
            sumty(i,j)   = 0.

!dmr @-@ iceb0
! JONO_cpl_wind march 2004 introducing sum of windvector to drive icebergs
            sumu10(i,j) = 0.
            sumv10(i,j) = 0.
!dmr @-@ iceb0

#if (ICEBERG == 2 && ISM != 2 )
            if (ntstep.eq.1) then ! first timestep of the simulation
               sumiceb(i,j) = 0.
            endif
#endif


          enddo
        enddo
        sumohsn=0.
        sumohss=0.

      else

        open(newunit=incoup_dat_id,file='startdata/incoup'//fini//'.dat',
     *              form='unformatted')
c~ #if ( ISOATM >= 2 && ISOOCN >= 1 )
c~         read(incoup_dat_id) sumohfx,sumoswr,sumihfx,sumiswr,sumofwf(:,:,ieau)
c~      &       ,sumisno(:,:,ieau),
c~      &       winstua,winstva,sumtx,sumty
c~ #else
        read(incoup_dat_id) sumohfx,sumoswr,sumihfx,sumiswr
     &       ,sumofwf(:,:,iwater),sumisno(:,:,iwater),
     &       winstua,winstva,sumtx,sumty
c~ #endif
        read(incoup_dat_id) cohesws,cohesw0,cohesw1,cohesw2,coulrad0,coulrad1,
c~ #if (ISOATM >= 1 )
c~      *       coulrad2,coulrads,codlrads,cohflux,coeflux,coevap(:,:,ieau)
c~ #else
     >    coulrad2,coulrads,codlrads,cohflux,coeflux,coevap(:,:,iwater)
c~ #endif

!--- dmr [NOTA] Why is clevap not included in the restart read?

        close(incoup_dat_id)
        sumohsn=0.0
        sumohss=0.0

        do j=1,nlon
          do i=1,nlat
            clhesws(i,j) = 0.
            clhesw0(i,j) = 0.
            clhesw1(i,j) = 0.
            clhesw2(i,j) = 0.
            clulrad0(i,j)= 0.
            clulrad1(i,j)= 0.
            clulrad2(i,j)= 0.
            clulrads(i,j)= 0.
            cldlrads(i,j)= 0.
            clhflux(i,j) = 0.
            cleflux(i,j) = 0.
c~ #if (ISOATM >= 1 )
            clevap(i,j,:)  = 0.
#if ( ISOOCN >= 1 && ISOATM >= 2 )

!dmr ### Lignes suivantes fausses ... doivent être remplacées par une relecture propre des isotopes dans le coupler ...
!dmr --- Initialisation barbare ...
! dmr ###            DO k=ieau+1,neauiso
! dmr ###              coevap(:,:,k) = coevap(:,:,ieau) * ( datmini(k) + 1.0d0 )
! dmr ###     &                          * rsmow(k)
! dmr ###            ENDDO

!mab: NEED restart file for that! if no restart wanted -> don't try to open it!

        IF (isolbm_restart .EQ. 1) THEN

        open(newunit=incoup_dat_id,file='startdata/wisocpl_restart.dat'
     >              ,form='unformatted')
        do k=ieau+1,neauiso
          read(incoup_dat_id) coevap(:,:,k)
        enddo
        do k=ieau+1,neauiso
          read(incoup_dat_id) sumofwf(:,:,k)
        enddo
        do k=ieau+1,neauiso
          read(incoup_dat_id) sumisno(:,:,k)
        enddo
        close(incoup_dat_id)


        ELSE

!dmr --- Initialisation barbare ...
          DO k=ieau+1,neauiso
            coevap(:,:,k) = coevap(:,:,ieau) * ( datmini(k) + 1.0d0 )
     &                          * rsmow(k)
            sumofwf(:,:,k) = sumofwf(:,:,k) * ( datmini(k) + 1.0d0 )
     &                          * rsmow(k)
            sumisno(:,:,k) = sumisno(:,:,k) * ( datmini(k) + 1.0d0 )
     &                          * rsmow(k)

          ENDDO
        ENDIF
!dmr ---
#elif ( ISOATM == 3 )
              WRITE(*,*) "Option non impementee !! coupling0.f"
#endif

c~ #else
c~             clevap(i,j)  = 0.
c~ #endif
            sumicof(i,j) = 0.
            sumuv10(i,j) = 0.
            sumpress(i,j)= 0.

!dmr @-@ iceb0
! JONO_cpl_wind
            sumu10(i,j) = 0.
            sumv10(i,j) = 0.
!dmr @-@ iceb0

          enddo
        enddo
      endif


! *** calculate total area for atmospheric gaussian grid darea
! *** dareafac is read in initemic.f from file darea.dat

      areafac=4*pi*radius**2/dsqrt(dble(nlon))
      tarea=0d0
      tareas=0d0

      do i=1,nlat
        darea(i)=dareafac(i)
        dareas(i)=darea(i)
        tarea=tarea + nlon*darea(i)
        tareas=tareas + nlon*dareas(i)
      enddo

! *** read adresses of neighbouring points and weights to
! *** interpolate between CLIO and ECBilt grid
! *** unit 45 connected to mozaic.w
! *** interpolate from ocean to atmosphere:
! *** do ji = 1, ijatm
! ***   zsum = 0.
! ***   do jk = 1, kamax
! ***     zsum = zsum + wo2a(ji,jk) * ocean(indo2a(ji,jk))
! ***   enddo
! ***   atmos(ji) = zsum
! *** enddo

! *** interpolate from atmosphere to ocean:
! *** do ji = 1, ijocn
! ***   zsum = 0.
! ***   do jk = 1, komax
! ***     zsum = zsum + wa2o(ji,jk) * atmos(indo2a(ji,jk))
! ***   enddo
! ***   ocean(ji) = zsum
! *** enddo

! *** kamax=14 komax=17 (see comcoup.h)


#if ( WRAP_EVOL == 1 )

      read(mozaic_f%id)
      read(mozaic_f%id) ((indo2a(ija,k),k=1,kamax),ija=1,ijatm)
      read(mozaic_f%id)
      read(mozaic_f%id) ((wo2a(ija,k),k=1,kamax),ija=1,ijatm)
      read(mozaic_f%id)
      read(mozaic_f%id) ((inda2o(ijo,k),k=1,komax),ijo=1,ijocn)
      read(mozaic_f%id)
      read(mozaic_f%id) ((wa2o(ijo,k),k=1,komax),ijo=1,ijocn)

#else

#if ( BATHY <= 1 )

      open(newunit=mozaic_w_id
     &  ,file = 'inputdata/clio/mozaic.w',status='old',form='unformatted')
#else

      write(charI,'(I5.5)'), la_date
      name_file ='inputdata/clio/mozaic_'//trim(charI)//'.w'

      open(newunit=mozaic_w_id,file = name_file,status='old',
     &        form='unformatted')

!dmr --- [TODO] add the read of KAMAX HERE !!!
#endif

      allocate(indo2a(ijatm,kamax))
      allocate(wo2a(ijatm,kamax))

      read(mozaic_w_id)
      read(mozaic_w_id) ((indo2a(ija,k),k=1,kamax),ija=1,ijatm)
      read(mozaic_w_id)
      read(mozaic_w_id) ((wo2a(ija,k),k=1,kamax),ija=1,ijatm)
      read(mozaic_w_id)
      read(mozaic_w_id) ((inda2o(ijo,k),k=1,komax),ijo=1,ijocn)
      read(mozaic_w_id)
      read(mozaic_w_id) ((wa2o(ijo,k),k=1,komax),ijo=1,ijocn)

      close(mozaic_w_id)
#endif

      do k=1,kamax
        do ija=1,ijatm
          if (indo2a(ija,k).eq.0) then
            indo2a(ija,k)=1
            wo2a(ija,k)=0d0
          endif
        enddo
      enddo

      do k=1,komax
        do ijo=1,ijocn
          if (inda2o(ijo,k).eq.0) then
            inda2o(ijo,k)=1
            wa2o(ijo,k)=0d0
          endif
        enddo
      enddo

      call ec_iniocbas

      call ec_la2co
      call ec_lae2co
      call ec_oc2co(1)
      call ec_at2co

! *** initialisation of fluxes

      do nn=1,ntyps
        call ec_fluxes(nn)
      enddo


900   format(a12,1x,i6)
910   format(a12,1x,e12.5)

      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_at2co
!-----------------------------------------------------------------------
! *** communicate coupler data to the atmosphere
!-----------------------------------------------------------------------

      USE comatm
      USE comphys
      use comcoup_mod

      implicit none

      integer i,j,k

      do j=1,nlon
        do i=1,nlat
c~ #if ( ISOATM == 1 )
c~           couprf(i,j)=torain(i,j,ieau)
c~           coupsf(i,j)=tosnow(i,j,ieau)
c~ #elif ( ISOATM >= 2 )
          couprf(i,j,:)=torain(i,j,:)
          coupsf(i,j,:)=tosnow(i,j,:)
          couptcc(i,j)=tcc(i,j)

         enddo
      enddo

! *** calculate flux correction to snow and rain over oceans
! *** to deal with too much precipitation in the arctic in ecbilt

      call preccor


      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_co2at
!-----------------------------------------------------------------------
! *** communicate coupler data to the atmosphere
!-----------------------------------------------------------------------

c~ #if (ISOATM >= 1 )
c~       use iso_param_mod,  only : ieau, neauiso, ieau18, rsmow
c~ #endif
c~ #if ( COMATM == 1 )
      use comatm
      use comphys
      use comcoup_mod
      use comsurf_mod
c~ #endif
#if ( DOWNSTS == 1 )
      use vertDownsc_mod, only : tsurf_d, tsurfn_d
#endif

      implicit none

c~ #if ( COMATM == 0 )
c~ #include "comatm.h"
c~ #include "comphys.h"
c~ #include "comcoup.h"
c~ #include "comsurf.h"
c~ #endif

      integer i,j,nn
c~ #if (ISOATM >= 1 )
c~       integer k
c~ #endif

      do j=1,nlon
        do i=1,nlat
          tsurf(i,j)=0.0
#if ( DOWNSTS == 1 )
          tsurf_d(i,j,:)=0.0
#endif
          albes(i,j) =0.0
          alb2es(i,j) =0.0
          qsurf(i,j)=0.0
          cdragv(i,j)=0.0
          do nn=1,ntyps
            tsurf(i,j) =tsurf(i,j)+fractn(i,j,nn)*tsurfn(i,j,nn)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
! dmr  Addition to account for the vertical downscaling of precipitations ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
#if ( DOWNSTS == 1 )
            tsurf_d(i,j,:) =tsurf_d(i,j,:)+fractn(i,j,nn)*tsurfn_d(i,j,nn,:)
#endif
            albes(i,j) =albes(i,j)+fractn(i,j,nn)*albesn(i,j,nn)
            alb2es(i,j)=alb2es(i,j)+fractn(i,j,nn)*alb2esn(i,j,nn)
            qsurf(i,j) =qsurf(i,j)+fractn(i,j,nn)*qsurfn(i,j,nn)
            cdragv(i,j)=cdragv(i,j)+fractn(i,j,nn)*cdragvn(i,j,nn)
          enddo
#if ( DOWNSTS == 1 )
! afq -- we might want mostly land for the downscaling ->
          if (fractn(i,j,nld).gt.0.05) then
             tsurf_d(i,j,:) = tsurfn_d(i,j,nld,:)
          end if
#endif
          hesws(i,j) =cohesws(i,j)  + clhesws(i,j)
          hesw0(i,j) =cohesw0(i,j)  + clhesw0(i,j)
          hesw1(i,j) =cohesw1(i,j)  + clhesw1(i,j)
          hesw2(i,j) =cohesw2(i,j)  + clhesw2(i,j)
          ulrad0(i,j)=coulrad0(i,j) + clulrad0(i,j)
          ulrad1(i,j)=coulrad1(i,j) + clulrad1(i,j)
          ulrad2(i,j)=coulrad2(i,j) + clulrad2(i,j)
          ulrads(i,j)=coulrads(i,j) + clulrads(i,j)
          dlrads(i,j)=codlrads(i,j) + cldlrads(i,j)
          hflux(i,j) =cohflux(i,j)  + clhflux(i,j)
          eflux(i,j) =coeflux(i,j)  + cleflux(i,j)
c~ #if (ISOATM >= 1 )
c~           DO k = ieau, neauiso
c~           evap(i,j,k)  =coevap(i,j,k)   + clevap(i,j,k)
c~           ENDDO
c~ #else
          evap(i,j,iwater)  =coevap(i,j,iwater) + clevap(i,j,iwater)
c~ #endif
! *** these two variables are communicated to the atmosphere only
! *** for the output routines
c~ #if ( ISOATM >= 2 && ISOOCN >= 1 )
c~           arunofo(i,j)=sumro(i,j,ieau)
c~           arunofl(i,j)=sumrl(i,j,ieau)
c~ #else
          arunofo(i,j)=sumro(i,j,iwater)
          arunofl(i,j)=sumrl(i,j,iwater)
c~ #endif
        enddo
      enddo

      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_sumfluxland(kst)
!-----------------------------------------------------------------------
! *** accumulate surface fluxes between atmosphere and land
!-----------------------------------------------------------------------

      use comemic_mod, only: ilan
      use comcoup_mod
      use comatm, only: nlat, nlon
      use comsurf_mod
      use comatm, only: iwater

#if ( ISOATM >= 1 )
      USE iso_param_mod, ONLY : ieau, neauiso, rsmow, ieau18, delta
#endif


      implicit none

      integer kst,i,j
c~ #if ( ISOATM >= 1 )
c~       integer k
c~       real, dimension(neauiso) :: variso
c~       character*10 varisonm
c~ #endif
      real*8  rlan

      if (kst.eq.1) then
        do j=1,nlon
          do i=1,nlat
            clhesws(i,j) = 0.
            clhesw0(i,j) = 0.
            clhesw1(i,j) = 0.
            clhesw2(i,j) = 0.
            clulrad0(i,j)= 0.
            clulrad1(i,j)= 0.
            clulrad2(i,j)= 0.
            clulrads(i,j)= 0.
            cldlrads(i,j)= 0.
            clhflux(i,j) = 0.
            cleflux(i,j) = 0.
c~ #if ( ISOATM >= 1 )
            clevap(i,j,:)  = 0.
c~ #else
c~             clevap(i,j)  = 0.
c~ #endif
c~ #if ( ISOATM >= 2 && ISOOCN >= 1 )
c~             sumrl(i,j,:)   = 0.
c~             sumro(i,j,:)   = 0.
c~ #else
            sumrl(i,j,:)   = 0.
            sumro(i,j,:)   = 0.
c~ #else
c~             sumrl(i,j)   = 0.
c~             sumro(i,j)   = 0.

c~ #endif /* on ISOATM >= 2 && ISOOCN >= 1 */
          enddo
        enddo
        sumhsn=0.
        sumhss=0.
      endif
      do j=1,nlon
        do i=1,nlat
          clhesws(i,j) = clhesws(i,j) + fractn(i,j,nld)*heswsn(i,j,nld)
          clhesw0(i,j) = clhesw0(i,j) + fractn(i,j,nld)*hesw0n(i,j,nld)
          clhesw1(i,j) = clhesw1(i,j) + fractn(i,j,nld)*hesw1n(i,j,nld)
          clhesw2(i,j) = clhesw2(i,j) + fractn(i,j,nld)*hesw2n(i,j,nld)
          clulrad0(i,j)= clulrad0(i,j)+ fractn(i,j,nld)*ulrad0n(i,j,nld)
          clulrad1(i,j)= clulrad1(i,j)+ fractn(i,j,nld)*ulrad1n(i,j,nld)
          clulrad2(i,j)= clulrad2(i,j)+ fractn(i,j,nld)*ulrad2n(i,j,nld)
          clulrads(i,j)= clulrads(i,j)+ fractn(i,j,nld)*ulradsn(i,j,nld)
          cldlrads(i,j)= cldlrads(i,j)+ fractn(i,j,nld)*dlradsn(i,j,nld)
          clhflux(i,j) = clhflux(i,j) + fractn(i,j,nld)*hfluxn(i,j,nld)
          cleflux(i,j) = cleflux(i,j) + fractn(i,j,nld)*efluxn(i,j,nld)
c~ #if ( ISOATM >= 1 )
c~        DO k = ieau, neauiso
c~        clevap(i,j,k)= clevap(i,j,k)+fractn(i,j,nld)
c~      &                  *evapn(i,j,nld,k)

c~
c~        ENDDO
c~ #else
          clevap(i,j,iwater)  = clevap(i,j,iwater)
     >                        + fractn(i,j,nld)*evapn(i,j,nld,iwater)
c~ #endif
c~ #if ( ISOATM >= 2 && ISOOCN >= 1 )
          sumrl(i,j,:)   = sumrl(i,j,:)   + couprunl(i,j,:)
          sumro(i,j,:)   = sumro(i,j,:)   + coupruno(i,j,:)
c~ #if ( 0 )
c~ !dmr --- Background checking of fluxes ...
c~           k = ieau18
c~           variso(:) = sumro(i,j,:)
c~           varisonm = "sumro"
c~
c~           if ((variso(k)*variso(ieau)).LT.0.0d0) then
c~            WRITE(*,*) "Inconsistent sign in variso: ", varisonm
c~            WRITE(*,*) variso(k),variso(ieau)
c~           endif
c~
c~           if (abs(variso(k)).gt.abs(variso(ieau)*1.0E-1)) then
c~             WRITE(*,*) "Too much isotopes in: ",varisonm, i,j
c~             WRITE(*,*) variso(k), variso(ieau)
c~             WRITE(*,*) delta(variso(:),k)
c~             READ(*,*)
c~           endif
c~ #endif
c~ #else
c~           sumrl(i,j)   = sumrl(i,j)   + couprunl(i,j)
c~           sumro(i,j)   = sumro(i,j)   + coupruno(i,j)
c~
c~ #endif /* ISOATM >= 2 && ISOOCN >= 1 */

        enddo
      enddo

      sumhsn=sumhsn+couphsnn
      sumhss=sumhss+couphsns

      if (kst.eq.ilan) then
        rlan=1d0/float(ilan)
        do j=1,nlon
          do i=1,nlat
            clhesws(i,j) = clhesws(i,j) *rlan
            clhesw0(i,j) = clhesw0(i,j) *rlan
            clhesw1(i,j) = clhesw1(i,j) *rlan
            clhesw2(i,j) = clhesw2(i,j) *rlan
            clulrad0(i,j)= clulrad0(i,j)*rlan
            clulrad1(i,j)= clulrad1(i,j)*rlan
            clulrad2(i,j)= clulrad2(i,j)*rlan
            clulrads(i,j)= clulrads(i,j)*rlan
            cldlrads(i,j)= cldlrads(i,j)*rlan
            clhflux(i,j) = clhflux(i,j) *rlan
            cleflux(i,j) = cleflux(i,j) *rlan
c~ #if ( ISOATM >= 1 )
c~        DO k = ieau, neauiso
c~             clevap(i,j,k)=clevap(i,j,k)*rlan
c~        ENDDO
c~ #else
            clevap(i,j,iwater)  = clevap(i,j,iwater)  *rlan
c~ #endif
c~ #if ( ISOATM >= 2 && ISOOCN >= 1 )
c~             sumrl(i,j,:)   = sumrl(i,j,:)   *rlan
c~             sumro(i,j,:)   = sumro(i,j,:)   *rlan
c~ #else
            sumrl(i,j,:)   = sumrl(i,j,:) * rlan
            sumro(i,j,:)   = sumro(i,j,:) * rlan

c~ #endif /* on  ISOATM >= 2 && ISOOCN >= 1 */
          enddo
        enddo
        sumhsn=sumhsn *rlan
        sumhss=sumhss *rlan
!       sumhsn=sumhsn
!       sumhss=sumhss
      endif

      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_sumfluxocean(ist,jst)
!-----------------------------------------------------------------------
! *** accumulate surface fluxes for the ocean-atmosphere interface
!-----------------------------------------------------------------------

      USE comatm
      use comphys
      use comemic_mod, only: nbclins, nbtrops, fracto
      use comcoup_mod
      use comsurf_mod

#if ( ISOATM >= 1 )
      USE iso_param_mod, ONLY : ieau, neauiso, ieau18, delta
#endif
      implicit none

      integer ist,jst,i,j
#if (ISOATM >= 1 )
      integer k
      real, dimension(neauiso) :: variso
      character*10 varisonm
#endif
      real*8  ratm,facstr,rndws,dwsm1,facsnow

      facsnow=rowat*rlatfus

      do j=1,nlon
        do i=1,nlat
          cohesws(i,j) = fractn(i,j,nse)*heswsn(i,j,nse)+
     &                                  fractn(i,j,noc)*heswsn(i,j,noc)
          cohesw0(i,j) = fractn(i,j,nse)*hesw0n(i,j,nse)+
     &                                  fractn(i,j,noc)*hesw0n(i,j,noc)
          cohesw1(i,j) = fractn(i,j,nse)*hesw1n(i,j,nse)+
     &                                  fractn(i,j,noc)*hesw1n(i,j,noc)
          cohesw2(i,j) = fractn(i,j,nse)*hesw2n(i,j,nse)+
     &                                  fractn(i,j,noc)*hesw2n(i,j,noc)
          coulrad0(i,j)= fractn(i,j,nse)*ulrad0n(i,j,nse)+
     &                                  fractn(i,j,noc)*ulrad0n(i,j,noc)
          coulrad1(i,j)= fractn(i,j,nse)*ulrad1n(i,j,nse)+
     &                                  fractn(i,j,noc)*ulrad1n(i,j,noc)
          coulrad2(i,j)= fractn(i,j,nse)*ulrad2n(i,j,nse)+
     &                                  fractn(i,j,noc)*ulrad2n(i,j,noc)
          coulrads(i,j)= fractn(i,j,nse)*ulradsn(i,j,nse)+
     &                                  fractn(i,j,noc)*ulradsn(i,j,noc)
          codlrads(i,j)= fractn(i,j,nse)*dlradsn(i,j,nse)+
     &                                  fractn(i,j,noc)*dlradsn(i,j,noc)
          cohflux(i,j) = fractn(i,j,nse)*hfluxn(i,j,nse)+
     &                                  fractn(i,j,noc)*hfluxn(i,j,noc)
          coeflux(i,j) = fractn(i,j,nse)*efluxn(i,j,nse)+
     &                                  fractn(i,j,noc)*efluxn(i,j,noc)
c~ #if ( ISOATM >= 1 )
c~           DO k = ieau, neauiso
c~           coevap(i,j,k)  = fractn(i,j,nse)*evapn(i,j,nse,k)+
c~      &                               fractn(i,j,noc)*evapn(i,j,noc,k)
c~           ENDDO
c~
c~ #else
c~           coevap(i,j)  = fractn(i,j,nse)*evapn(i,j,nse)+
c~      &                                  fractn(i,j,noc)*evapn(i,j,noc)
          coevap(i,j,iwater)  = fractn(i,j,nse)*evapn(i,j,nse,iwater)+
     &                          fractn(i,j,noc)*evapn(i,j,noc,iwater)
c~ #endif
        enddo
      enddo

      if ((jst.eq.1.and.ist.eq.1).or.iobclint.eq.0) then
        do j=1,nlon
          do i=1,nlat
            sumohfx(i,j) = 0.
            sumoswr(i,j) = 0.
            sumihfx(i,j) = 0.
            sumiswr(i,j) = 0.
c~ #if ( ISOATM >= 2 && ISOOCN >= 1 )
            sumofwf(i,j,:) = 0.
            sumisno(i,j,:) = 0.
c~ #else
c~             sumofwf(i,j) = 0.
c~             sumisno(i,j) = 0.

#if ( ICEBERG == 2 && ISM != 2  )
            sumiceb(i,j) = 0.
#endif

c~ #endif
            sumicof(i,j) = 0.
            sumuv10(i,j) = 0.
            sumpress(i,j)= 0.

!dmr @-@ iceb0
! JONO_cpl_wind
            sumu10(i,j) = 0.
            sumv10(i,j) = 0.
!dmr @-@ iceb0

          enddo
        enddo
        sumohsn=0.0
        sumohss=0.0
        iobclint=0
      endif

      if ((jst.eq.1.and.ist.eq.1).or.iobtropt.eq.0) then
        do j=1,nlon
          do i=1,nlat
            winstua(i,j) = 0.
            winstva(i,j) = 0.
          enddo
        enddo
        iobtropt=0
      endif

      iobtropt=iobtropt+1
      iobclint=iobclint+1

      do j=1,nlon
        do i=1,nlat
          if (fracto(i,j).gt.epss) then
c~ #if ( ISOATM >= 2 )
c~             sumohfx(i,j) = sumohfx(i,j)+
c~      &                         dlradsn(i,j,noc)-ulradsn(i,j,noc)
c~      &                        -efluxn(i,j,noc)-hfluxn(i,j,noc)-
c~      &                        facsnow*coupsf(i,j,ieau)
c~ #else
            sumohfx(i,j) = sumohfx(i,j)+
     &                         dlradsn(i,j,noc)-ulradsn(i,j,noc)
     &                        -efluxn(i,j,noc)-hfluxn(i,j,noc)-
     &                        facsnow*coupsf(i,j,iwater)
c~ #endif
            sumoswr(i,j) = sumoswr(i,j)+ heswsn(i,j,noc)
            sumihfx(i,j) = sumihfx(i,j)+
     &                         dlradsn(i,j,nse)
     &                        -efluxn(i,j,nse)-hfluxn(i,j,nse)
            sumiswr(i,j) = sumiswr(i,j)+ heswsn(i,j,nse)
c~ #if ( ISOATM == 1 )
c~             sumofwf(i,j) = sumofwf(i,j)+
c~      &                      (evapn(i,j,noc,ieau)*fractn(i,j,noc)+
c~      &                evapn(i,j,nse,ieau)*fractn(i,j,nse))/fracto(i,j)-
c~      &                      couprf(i,j)-coupsf(i,j)-sumro(i,j)
c~             sumisno(i,j) = sumisno(i,j)+ coupsf(i,j)-
c~      &                   evapn(i,j,nse,ieau)*fractn(i,j,nse)/fracto(i,j)
c~ #elif ( ISOATM >= 2 && ISOOCN == 0 )
c~             sumofwf(i,j) = sumofwf(i,j)+
c~      &                      (evapn(i,j,noc,ieau)*fractn(i,j,noc)+
c~      &                evapn(i,j,nse,ieau)*fractn(i,j,nse))/fracto(i,j)-
c~      &                      couprf(i,j,ieau)-coupsf(i,j,ieau)-sumro(i,j)
c~             sumisno(i,j) = sumisno(i,j)+ coupsf(i,j,ieau)-
c~      &                   evapn(i,j,nse,ieau)*fractn(i,j,nse)/fracto(i,j)
c~ #elif ( ISOATM >= 2 && ISOOCN >= 1 )
            sumofwf(i,j,:) = sumofwf(i,j,:)+
     &             (evapn(i,j,noc,:)*fractn(i,j,noc)+
     &              evapn(i,j,nse,:)*fractn(i,j,nse))/fracto(i,j)-
     &              couprf(i,j,:)-coupsf(i,j,:)-sumro(i,j,:)
            sumisno(i,j,:) = sumisno(i,j,:)+ coupsf(i,j,:)-
     &              evapn(i,j,nse,:)*fractn(i,j,nse)/fracto(i,j)

! Calculate yearly sum of excess snow to iceberg
#if ( ICEBERG == 2 && ISM != 2 )
            sumiceb(i,j) = coupiceb(i,j)! [m] excess snow summed since jan. 1st
#endif

c~ #endif
            sumicof(i,j) = sumicof(i,j)+hficof(i,j)
          endif
!-For the iceberg coupling
          sumuv10(i,j) = sumuv10(i,j)+uv10(i,j)
          sumpress(i,j)= sumpress(i,j)+pground(i,j)
!-icb
          winstua(i,j) = winstua(i,j) + winstu(i,j)
          winstva(i,j) = winstva(i,j) + winstv(i,j)

!dmr @-@ iceb0
! JONO_cpl_wind (don't think i should put this in the if fracto<epss-loop...
! because some land-surface wind also needed for interpolation at iceberg position
! land-wind is taken the same as ocean but
! height and surface could have a significant influence... parameterisation??)
          sumu10(i,j) = sumu10(i,j) + utot10(i,j)
          sumv10(i,j) = sumv10(i,j) + vtot10(i,j)
!dmr @-@ iceb0

        enddo
      enddo

      sumohsn      = sumohsn+sumhsn
      sumohss      = sumohss+sumhss

      if (iobclint.eq.nbclins) then
        ratm=1d0/float(nbclins)
        do j=1,nlon
          do i=1,nlat
            sumohfx(i,j) = sumohfx(i,j)*ratm
            sumoswr(i,j) = sumoswr(i,j)*ratm
            sumihfx(i,j) = sumihfx(i,j)*ratm
            sumiswr(i,j) = sumiswr(i,j)*ratm
c~ #if ( ISOATM >= 2 && ISOOCN >= 1 )
c~             sumofwf(i,j,:) = sumofwf(i,j,:)*ratm
c~             sumisno(i,j,:) = sumisno(i,j,:)*ratm
c~ #else
            sumofwf(i,j,:) = sumofwf(i,j,:)*ratm
            sumisno(i,j,:) = sumisno(i,j,:)*ratm

c~ #endif
            sumicof(i,j) = sumicof(i,j)*ratm
            sumuv10(i,j) = sumuv10(i,j)*ratm
            sumpress(i,j)=sumpress(i,j)*ratm

!dmr @-@ iceb0
! JONO_cpl_wind not the sum but the average (per second i guess) is exported?
            sumu10(i,j) = sumu10(i,j)*ratm
            sumv10(i,j) = sumv10(i,j)*ratm
!dmr @-@ iceb0

          enddo
        enddo
!       sumohsn=sumohsn*ratm
!       sumohss=sumohss*ratm
        sumohsn=sumohsn
        sumohss=sumohss
        iobclint=0
        if (iclimflux.eq.1) call ec_climflux(ist,jst)
      endif

      if (iobtropt.eq.nbtrops) then
        rndws=1d0/dble(ndayws)
        dwsm1=ndayws-1d0
        ratm=1d0/float(nbtrops)
        do j=1,nlon
          do i=1,nlat
            winstua(i,j) = winstua(i,j)*ratm
            winstva(i,j) = winstva(i,j)*ratm
            sumtx(i,j)=(dwsm1*sumtx(i,j)+winstua(i,j))*rndws
            sumty(i,j)=(dwsm1*sumty(i,j)+winstva(i,j))*rndws
          enddo
        enddo
        iobtropt=0
      endif

      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_iniocbas
!-----------------------------------------------------------------------

! *** initialises the oceanic bassin : based on iniocout
!-----------------------------------------------------------------------

      USE comatm,      only: nlon, nlat, darea
      use comdiag,     only: iocbasa, arocbasa, nbasa
      use comsurf_mod, only: fractoc
      use comunit,     only: iuo

      use global_constants_mod, only: dblp=>dp, ip

      implicit none


      integer     i,j,k,ias,ibas
      character*1 ch(nlon),space
      real*8 sum(nlat,nbasa)
      integer(kind=ip):: ocbasdat_id

!     (1)='ATL N  '
!     (2)='PAC N  '
!     (3)='ARCTIC '
!     (4)='INDIAN '
!     (5)='ATL S  '
!     (6)='PAC S  '
!     (7)='ANTAR  '


! *** computation of basin masks

! *** asci number of small letter a
      open(newunit=ocbasdat_id,file='inputdata/ocbas.dat',
     &                   status='old',form='formatted')

      ias=ichar('a') - 1
      do j=nlat,1,-1
        read (ocbasdat_id,100) k,space,(ch(i),i=1,nlon)
        do i=1,nlon
          iocbasa(j,i)=ichar(ch(i)) - ias
          if (iocbasa(j,i).lt.1.or.iocbasa(j,i).gt.nbasa) then
            if (fractoc(j,i).gt.0.01) then
              write(iuo+29,*) 'no ocean basin defined in',i,j
            endif
            iocbasa(j,i)=0
          endif
        enddo
      enddo

      close(ocbasdat_id)

! *** computation of area of ocean basins
      do ibas=1,nbasa
        arocbasa(ibas)=0.
        do i=1,nlat
          sum(i,ibas)=0d0
        enddo
      enddo

      do j=1,nlon
       do i=1,nlat
         if (iocbasa(i,j).gt.0) then
          sum(i,iocbasa(i,j))=sum(i,iocbasa(i,j))+fractoc(i,j)
         endif
         sum(i,nbasa)=sum(i,nbasa)+fractoc(i,j)
        enddo
      enddo

      do i=1,nlat
       do ibas=1,nbasa
        arocbasa(ibas)=arocbasa(ibas)+sum(i,ibas)*darea(i)
       enddo
      enddo

 100  format(i4,65A1)
      return
      end subroutine ec_iniocbas


!23456789012345678901234567890123456789012345678901234567890123456789012

      SUBROUTINE preccor
!-----------------------------------------------------------------------

! *** computes mean values of ocean basins
!-----------------------------------------------------------------------


      USE comatm     , only: nwisos, nlon, nlat, darea
      use comphys    , only: coran, corpn, corac, corid, coras, corps
     >                     , coraa
      use comdiag    , only: nbasa, iocbasa
      use comemic_mod, only: fracto
      use comcoup_mod, only: couprf, coupsf
      use comsurf_mod, only:


c~ #if ( ISOATM >= 2 )
c~       USE iso_param_mod, ONLY : ieau, neauiso
c~ #endif
      implicit none

      integer i,j
      real*8 amcor(nbasa)
c~ #if ( ISOATM >= 2 )
c~       REAL, DIMENSION(neauiso) :: sum1rf,sum1sf
c~       REAL :: sum2rf,sum2sf
c~ #else
      real*8, DIMENSION(nwisos):: sum1rf,sum1sf
      real*8                   :: sum2rf,sum2sf
c~ #endif

!     (1)='ATL N  '
      amcor(1)=corAN
!     (2)='PAC N  '
      amcor(2)=corPN
!     (3)='ARCTIC '
      amcor(3)=corAC
!     (4)='INDIAN '
      amcor(4)=corID
!     (5)='ATL S  '
      amcor(5)=corAS
!     (6)='PAC S  '
      amcor(6)=corPS
!     (7)='ANTAR  '
      amcor(7)=corAA

! *** apply correction to rainfall (couprf)
! *** sum volume of water that is subtracted at higher latitudes over
! *** gridpoints that cover ocean

      sum1rf(:)=0.0
      do j=1,nlon
       do i=1,nlat
!        if (iocbasa(i,j).gt.0.and.fracto(i,j).gt.epss) then
         if (iocbasa(i,j).gt.0.and.fracto(i,j).gt.0.1) then
           if (amcor(iocbasa(i,j)).lt.0 ) then
c~ #if ( ISOATM >= 2 )
c~              sum1rf(:)=sum1rf(:)+couprf(i,j,:)*darea(i)*amcor(iocbasa(i,j))
c~              couprf(i,j,:)=couprf(i,j,:)*(1.+amcor(iocbasa(i,j)))
c~ #else
            sum1rf(:)=sum1rf(:)+couprf(i,j,:)*darea(i)*amcor(iocbasa(i,j))
            couprf(i,j,:)=couprf(i,j,:)*(1.+amcor(iocbasa(i,j)))
c~ #endif
           endif
         endif
       enddo
      enddo

! *** sum area over which this water is to be spread, over full ocean
! *** grid points only
      sum2rf=0.0

      do j=1,nlon
       do i=1,nlat
         if (iocbasa(i,j).gt.0.and.fracto(i,j).eq.1d0) then
           if (amcor(iocbasa(i,j)).gt.0 ) then
             sum2rf=sum2rf+darea(i)*amcor(iocbasa(i,j))
           endif
         endif
       enddo
      enddo
! *** calculate water flux per meter squared

      sum1rf(:)=sum1rf(:)/sum2rf

! *** reduce rainfall accordingly

      do j=1,nlon
       do i=1,nlat
         if (iocbasa(i,j).gt.0.and.fracto(i,j).eq.1d0) then
           if (amcor(iocbasa(i,j)).gt.0 ) then
c~ #if ( ISOATM >= 2 )
c~              couprf(i,j,:)=couprf(i,j,:)-sum1rf(:)*amcor(iocbasa(i,j))
c~ #else
             couprf(i,j,:)=couprf(i,j,:)-sum1rf(:)*amcor(iocbasa(i,j))
c~ #endif
           endif
         endif
       enddo
      enddo

! *** apply correction to snowfall (coupsf)

      sum1sf(:)=0.0
      do j=1,nlon
       do i=1,nlat
         if (iocbasa(i,j).gt.0.and.fracto(i,j).gt.0.1) then
           if (amcor(iocbasa(i,j)).lt.0 ) then
c~ #if ( ISOATM >= 2 )
c~              sum1sf(:)=sum1sf(:)+coupsf(i,j,:)*darea(i)*amcor(iocbasa(i,j))
c~              coupsf(i,j,:)=coupsf(i,j,:)*(1.+amcor(iocbasa(i,j)))
c~ #else
             sum1sf(:)=sum1sf(:)+coupsf(i,j,:)*darea(i)*amcor(iocbasa(i,j))
             coupsf(i,j,:)=coupsf(i,j,:)*(1.+amcor(iocbasa(i,j)))
c~ #endif
           endif
         endif
       enddo
      enddo

      sum2sf=0.0

      do j=1,nlon
       do i=1,nlat
         if (iocbasa(i,j).gt.0.and.fracto(i,j).eq.1d0) then
           if (amcor(iocbasa(i,j)).gt.0 ) then
             sum2sf=sum2sf+darea(i)*amcor(iocbasa(i,j))
           endif
         endif
       enddo
      enddo

      sum1sf(:)=sum1sf(:)/sum2sf

      do j=1,nlon
       do i=1,nlat
         if (iocbasa(i,j).gt.0.and.fracto(i,j).eq.1d0) then
           if (amcor(iocbasa(i,j)).gt.0 ) then
c~ #if ( ISOATM >= 2 )
c~              coupsf(i,j,:)=coupsf(i,j,:)-sum1sf(:)*amcor(iocbasa(i,j))
c~ #else
             coupsf(i,j,:)=coupsf(i,j,:)-sum1sf(:)*amcor(iocbasa(i,j))
c~ #endif
           endif
         endif
       enddo
      enddo

      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_climflux(ist,jst)
!-----------------------------------------------------------------------
! *** this routine calculates climatogical SST's and heatfluxes
! *** from year 1 onwards and outputs these to a file
!-----------------------------------------------------------------------


      USE comatm
      use comphys
      use comemic_mod, only: iday, imonth, iobclin, iyear, nyears, undef,fracto
      use comcoup_mod
      use comsurf_mod

      implicit none


      integer i,j,k,l,index,istep,ist,jst

      real*4 hulph(nlat,nlon),hulpt(nlat,nlon)
      real*8 hefxcl(nlat,nlon,360),tmixcl(nlat,nlon,360),ryear,riatm
      real*8 hefx(nlat,nlon),hefxo(nlat,nlon)

      common /mixlayer/ hefxcl,tmixcl,hefx,hefxo


      if ( ist.eq.iobclin) then
        do k=1,360
          do j=1,nlon
            do i=1,nlat
              hefxcl(i,j,k)=0d0
              tmixcl(i,j,k)=0d0
            enddo
          enddo
        enddo
      endif

      index=(imonth-1)*30+iday

      do i=1,nlat
        do j=1,nlon
          hefxcl(i,j,index)=hefxcl(i,j,index)+sumohfx(i,j)+sumoswr(i,j)
          tmixcl(i,j,index)=tmixcl(i,j,index)+tsurfn(i,j,noc)
        enddo
      enddo

      if (index.eq.360.and.iyear.eq.nyears) then
        ryear=1/real(nyears)
        do k=1,360
          do i=1,nlat
            do j=1,nlon
              if (fracto(i,j).gt.epss) then
                hefxcl(i,j,k)=hefxcl(i,j,k)*ryear
                tmixcl(i,j,k)=tmixcl(i,j,k)*ryear
              else
                hefxcl(i,j,k)=undef
                tmixcl(i,j,k)=undef
              endif
            enddo
          enddo
        enddo
        index=0
        do k=1,360/iobclin
          do l=1,iobclin
            index=index+1
            write(51) ((real(hefxcl(i,j,index)),j=1,nlon),i=1,nlat)
            write(52) ((real(tmixcl(i,j,index)),j=1,nlon),i=1,nlat)
          enddo
        enddo
      endif

      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_wrendcoup
!-----------------------------------------------------------------------
! *** output fluxes for start new run
!-----------------------------------------------------------------------
#if ( ISOATM >= 1 )
      USE iso_param_mod, ONLY : ieau, neauiso
#endif

      use comcoup_mod
      use comatm, only: nlat, nlon, iwater
      use comunit

      use newunit_mod, only: newunit_id, wisocpl_restart_id

      implicit none

      integer i,j
      real*8  dum(nlat,nlon,12)

#if ( ISOATM >= 2 && ISOOCN >= 1 )
      integer k
#endif

      do i=1,nlat
        do j=1,nlon
           dum(i,j,1)   =cohesws(i,j)  + clhesws(i,j)
           dum(i,j,2)   =cohesw0(i,j)  + clhesw0(i,j)
           dum(i,j,3)   =cohesw1(i,j)  + clhesw1(i,j)
           dum(i,j,4)   =cohesw2(i,j)  + clhesw2(i,j)
           dum(i,j,5)   =coulrad0(i,j) + clulrad0(i,j)
           dum(i,j,6)   =coulrad1(i,j) + clulrad1(i,j)
           dum(i,j,7)   =coulrad2(i,j) + clulrad2(i,j)
           dum(i,j,8)   =coulrads(i,j) + clulrads(i,j)
           dum(i,j,9)   =codlrads(i,j) + cldlrads(i,j)
           dum(i,j,10)  =cohflux(i,j)  + clhflux(i,j)
           dum(i,j,11)  =coeflux(i,j)  + cleflux(i,j)
c~ #if ( ISOATM >= 1 )
c~            dum(i,j,12)  =coevap(i,j,ieau)   + clevap(i,j,ieau)
c~ #else
           dum(i,j,12)  =coevap(i,j,iwater)   + clevap(i,j,iwater)
c~ #endif
        enddo
      enddo

#if ( ISOATM >= 2 && ISOOCN >= 1 )
c~       write(newunit_id) sumohfx,sumoswr,sumihfx,sumiswr,sumofwf(:,:,ieau),
c~      &sumisno(:,:,ieau),
c~      &winstua,winstva,sumtx,sumty

      do k=ieau+1,neauiso
        write(wisocpl_restart_id) coevap(:,:,k)
      enddo
      do k=ieau+1,neauiso
        write(wisocpl_restart_id) sumofwf(:,:,k)
      enddo
      do k=ieau+1,neauiso
        write(wisocpl_restart_id) sumisno(:,:,k)
      enddo
#endif

      write(newunit_id) sumohfx,sumoswr,sumihfx,sumiswr
     >                 ,sumofwf(:,:,iwater),sumisno(:,:,iwater)
     >                 ,winstua,winstva,sumtx,sumty

      write(newunit_id) dum
      return
      end

