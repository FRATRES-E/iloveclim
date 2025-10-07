!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:29 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:29 CET 2009

      SUBROUTINE veget_wr(kveget,nyvegt, nwrveg,iyr0vg, iyear, nyears,
     &             temveg,gd0veg,prcveg,tpsdry,
     &             prcmin,bmtdry, epss,fracgr,darea,gd5veg,
     &             titveg)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

! Write vegetation output on ASCI formated files "veget.zav" & "veget.outp" :
!       zonal+global mean: every "nyvegt" (=every call veget)
!       2.D maps: every "nwrveg" year(s).
!----
!  iyear = 0 : (=1rst call) open & write Header (2 lines) of output files.
!            and if kveget < 0 write Equilibr. Vegetation on output files.
!  WARNING : put Special Values in array:
!                          st,sg,sd,temveg,gd0veg,prcveg,tpsdry,snlt,blai
!  use file_unit vegetzav_id & vegetoutp_id (remains opened during the run).
!---------
!  modif : 01/10/0

      use comatm, only: nlat, nlon

      use veget_mod, only: numvegvar, st, sg, sd, snlt, blai, blai
     &             , anup, stock, b1, b2, b3, b4, pnpp, st, sg, sd, snlt
     &             , blai, blai, anup, stock, b1, b2, b3, b4, pnpp
     &             , newvegvar, st_moy

      USE Vegetation_Output, only: yearly_means, output, open, close
     &                           , write, vegetzav_id, vegetoutp_id
#if ( IMSK == 1 )
      USE input_icemask
#endif
#if ( CYCC == 2 )
      USE veget_iso
#endif

      use comland_mod, only: albland,forestfr
! --- BdB 05-2019: added variables for rewriting files
      use comemic_mod, only: nwrskip, new_year_veg, current_int_veg

      use ipcc_output_mod, only: cland
      
#if ( CLM_INDICES >= 2 )
       USE global_constants_mod, ONLY: sip
       USE CLIMATE_INDICES_MOD, only: SET_VEG_VARS
#endif             

      implicit none

!--dummy variables :
!- input :
      integer kveget, nyvegt, nwrveg, iyr0vg, iyear, nyears
      real*8 temveg(nlat,nlon), gd0veg(nlat,nlon),
     &     prcveg(nlat,nlon,2), tpsdry(nlat,nlon)
     &    ,gd5veg(nlat,nlon)
      real*8 epss, prcmin, bmtdry
      real*8 fracgr(nlat,nlon), darea(nlat)
      character(len=15) :: titveg
!- output :

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!--local variables saved (in common) from one call to the other :
      integer, save :: nnctr
c~       common / cm0wvg / nnctr

      integer nvgmax, kwradd
#if ( CYCC == 2 )
      parameter ( nvgmax = 20 , kwradd = 0 )
#else
      parameter ( nvgmax = 20 , kwradd = 0 )
#endif
      real*8 vegsum, vegmap
      common / cm1wvg / vegsum, vegmap(nlat,nlon,nvgmax),vegsum2

!--data + local variables :
      integer i,j,k,n, kk,ns, nbwr1,nbwr2
      real*8 zero, one, var(nlon)
      real*8 dzero, zavsum, totsum,vegsum2
      real*8 vegzav(0:nlat,nvgmax), cfmzav(nvgmax), vegspv, ttscun, tts1
      character(len=6) ccny
      character(len=8) titzav(nvgmax)
      character(len=27) titmap(nvgmax)
      character(len=30) fmtveg, fmtzav
      real*8 soiltype2(nlat,nlon)
      
!dmr for newunit file ID
      integer :: vegetccny_id

!--For output :
      data vegspv / 999.0 /
      data fmtveg / '(65E18.5E3)                      ' /
      data fmtzav / '(65E18.5E3)                      ' /
#if ( CYCC == 2 )
      data cfmzav / 4*100. , 2*1. , 25. , 1. , 3*0.001, 7*1E-12,  100, 1./
#else
      data cfmzav / 4*100. , 2*1. , 25. , 1. , 3*0.001, 7*1E-12,  100, 0.001/
#endif
#if ( CYCC == 2 )
      data titzav / 'tree o/o', 'grass   ', 'desert  ', 't.needle',
     &              't_lai   ', 'g_lai   ', 'albe o/o',
     &              'temp oC ', 'gdd0 e-3', 'prc m/y ', 'prc_veg ',
     &              'uptake GtC/y ', 'stock GtC' , 'b1 GtC' ,
     &              ' b2 GtC' , 'b3 GtC' ,' b4 GtC' ,'NPP GtC/y ',
     &               'icesheet o/o', 'B1T13 o/oo ?'  /
#else
      data titzav / 'tree o/o', 'grass   ', 'desert  ', 't.needle',
     &              't_lai   ', 'g_lai   ', 'albe o/o',
     &              'temp oC ', 'gdd0 e-3', 'prc m/y ', 'prc_veg ',
     &              'uptake GtC/y ', 'stock GtC' , 'b1 GtC' ,
     &              ' b2 GtC' , 'b3 GtC' ,' b4 GtC' ,'NPP GtC/y ',
     &               'icesheet o/o', 'gdd5 e-3' /
#endif
#if ( CYCC == 2 )
      data titmap /
     & 'Fraction of Tree   (o/o) ; ' , 'Fraction of Grass  (o/o) ; ' ,
     & 'Fraction of Desert (o/o) ; ' , 'Fract. of Needle L.(o/o) ; ' ,
     & ' L.A.I. for Tree    (1)  ; ' , ' L.A.I. for Grass   (1)  ; ' ,
     & 'An. Albedo, snow=0 (o/o) ; ' ,
     & ' Annual Temperature (oC) ; ' , ' Annual G.D.D.0  (/1000) ; ' ,
     & ' Annual Precipit.  (m/y) ; ' , 'Reduced Prc. (Min&Dry_S) ; ' ,
     & ' Annual C Uptake (kgC/y/m2) ; ','Carbon Stock (kgC/m2) ; ' ,
     & ' B1 (kgC/m2) ; ' , ' B2 (kgC/m2) ; ' , 'B3 (kgC/m2) ; ',
     & 'B4 (kgC/m2) ; ' ,
     & ' NPP (kgC/y/m2) ','icesheet fraction o/o' , 'B1T13 o/oo ?'/
#else
      data titmap /
     & 'Fraction of Tree   (o/o) ; ' , 'Fraction of Grass  (o/o) ; ' ,
     & 'Fraction of Desert (o/o) ; ' , 'Fract. of Needle L.(o/o) ; ' ,
     & ' L.A.I. for Tree    (1)  ; ' , ' L.A.I. for Grass   (1)  ; ' ,
     & 'An. Albedo, snow=0 (o/o) ; ' ,
     & ' Annual Temperature (oC) ; ' , ' Annual G.D.D.0  (/1000) ; ' ,
     & ' Annual Precipit.  (m/y) ; ' , 'Reduced Prc. (Min&Dry_S) ; ' ,
     & ' Annual C Uptake (kgC/y/m2) ; ','Carbon Stock (kgC/m2) ; ' ,
     & ' B1 (kgC/m2) ; ' , ' B2 (kgC/m2) ; ' , 'B3 (kgC/m2) ; ',
     & 'B4 (kgC/m2) ; ' ,
     & ' NPP (kgC/y/m2) ','icesheet fraction o/o' ,
     & 'Annual G.D.D.0  (/1000) ; '/
#endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

 1000 format(A30,1x,E13.6,1x,I6,1x,I6,1x,I6,1x,I6)
 1111 format(3(F7.1,1X,F7.3,1X),I3,A)
 1112 format(3(F7.2,1X,F7.3,1X),I3,A)

      zero = 0.
      one  = 1.
      dzero = 0.d0
!-----

      if (iyear.eq.0) then
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1 ) 1rst call : Open & Write Header of veget. output files :
!-----------------------------------------------------------------------

        nnctr = 0
        do n=1,len(titveg)
          if (titveg(n:n).ne.' ') nnctr = n
        enddo
        vegsum = 0.
        vegsum2 = 0.
        vegmap(:,:,:)= 0.
        st_moy(:,:)= 0.

!- time scale unit of "veget.zav" Header = kyr
        ttscun = one/1000.
        tts1 = (nyvegt+iyr0vg)*ttscun
        nbwr1 = nyears / nyvegt
        nbwr2 = nyears / nwrveg
        if (kveget.eq.0) then
!--No vegetation : write Imposed Tree Fraction & Albedo (<- from EcBilt) :
          titmap(1) = 'EcBilt Origin. Tree Frac.; '
          nbwr2 = (5-nvgmax)*nbwr2
        elseif (kveget.lt.0) then
          tts1 = iyr0vg*ttscun
          nbwr1 = nbwr1 + 1
          nbwr2 = -nvgmax*(1+nbwr2)
        else
          nbwr2 = -nvgmax*nbwr2
        endif
#if ( FAST_OUTPUT == 0 )

        open(newunit=vegetzav_id,file='outputdata/vegetation/veget.zav'
     &   ,status='unknown')
        write(vegetzav_id,1000) fmtzav, vegspv, 1+nlat, nvgmax, nbwr1,65
        write(vegetzav_id,1112) -87.19-5.625, 5.625, 8., 0.,
     &                    tts1, nyvegt*ttscun, -5

        open(newunit=vegetoutp_id,
     &   file='outputdata/vegetation/veget.outp', status='unknown')
        write(vegetoutp_id,1000) fmtveg, vegspv, nlon, nlat, nbwr2, 65
        write(vegetoutp_id,1111) 0., 5.625, -87.19, 5.625, 1., 1., 0

#endif
        if (kveget.ge.0) return
!-------
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      if (kveget.eq.0) then
!--No vegetation : fill in "st" with EcBilt "Original" Tree Fraction  :
        do j=1,nlon
         do i=1,nlat
          if (fracgr(i,j).gt.epss) st(i,j) = forestfr(i,j)
         enddo
        enddo
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) Compute & Write Global + Zonal mean values.
!-----------------------------------------------------------------------

!-- compute Global & Zonal mean (veget+albedo+Veg_input_var) :

!-- Initialize :
      totsum = 0.
        vegzav(:,:)=0.

!-- Compute sum :
      iloop: do i=1,nlat
        zavsum = 0.
        jloop: do j=1,nlon
          if (fracgr(i,j).gt.epss) then
            zavsum = zavsum + fracgr(i,j)

#if ( ISM == 1 )
            if (flgism) then
             soiltype2(i,j)=soiltype(i,j,ca)
             if(i.le.6) soiltype2(i,j)=1.
             if(soiltype(i,j,ca).lt.0.) soiltype2(i,j)=0.
            else
#endif
             soiltype2(i,j)=0.
#if ( IMSK == 1 )
            if (icemask(i,j).gt.0.9) then
              soiltype2(i,j)=1.
            endif
#endif
#if ( ISM == 1 )
            endif
#endif
            vegzav(i,1) = vegzav(i,1) + fracgr(i,j)*st(i,j)
            vegzav(i,2) = vegzav(i,2) + fracgr(i,j)*sg(i,j)
            vegzav(i,3) = vegzav(i,3) + fracgr(i,j)*sd(i,j)
            vegzav(i,4) = vegzav(i,4) + fracgr(i,j)*st(i,j)*snlt(i,j)
            vegzav(i,5) = vegzav(i,5) + fracgr(i,j)*st(i,j)*blai(i,j,1)
            vegzav(i,6) = vegzav(i,6) + fracgr(i,j)*sg(i,j)*blai(i,j,2)
            vegzav(i,7) = vegzav(i,7) + fracgr(i,j)*
     &  (albland(i,j,1)+albland(i,j,2)+albland(i,j,3)+albland(i,j,4))
            vegzav(i,8) = vegzav(i,8) + fracgr(i,j)*temveg(i,j)
            vegzav(i,9) = vegzav(i,9) + fracgr(i,j)*gd0veg(i,j)
            vegzav(i,10)= vegzav(i,10)+ fracgr(i,j)*prcveg(i,j,1)
            vegzav(i,11)= vegzav(i,11)+ fracgr(i,j)*prcveg(i,j,2)
            vegzav(i,12)= vegzav(i,12)+ fracgr(i,j)*(anup(i,j)*darea(i))
            vegzav(i,13)= vegzav(i,13)+ fracgr(i,j)*(stock(i,j)*darea(i))
            vegzav(i,14)= vegzav(i,14)+ fracgr(i,j)*(b1(i,j)*darea(i))
            vegzav(i,15)= vegzav(i,15)+ fracgr(i,j)*(b2(i,j)*darea(i))
            vegzav(i,16)= vegzav(i,16)+ fracgr(i,j)*(b3(i,j)*darea(i))
            vegzav(i,17)= vegzav(i,17)+ fracgr(i,j)*(b4(i,j)*darea(i))
            vegzav(i,18)= vegzav(i,18)+ fracgr(i,j)*(pnpp(i,j)*darea(i))
            vegzav(i,19)= vegzav(i,19)+ fracgr(i,j)*soiltype2(i,j)
            vegzav(i,20) = vegzav(i,20) + fracgr(i,j)*gd5veg(i,j)
! #if ( CYCC == 2 )
!             vegzav(i,20)= vegzav(i,20)+ fracgr(i,j)*(B1T13(i,j)*darea(i))
! #endif
          endif
        enddo jloop

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!-- Compute Zonal Mean Value :
        totsum = totsum + zavsum*darea(i)
        kloop: do k=nvgmax,1,-1
#if ( CYCC == 2 )
          if (((k.ge.12).AND.(k.le.18)).OR.(k.eq.20)) then
#else
          if ((k.ge.12).AND.(k.le.18)) then
#endif
            vegzav(i,k) = vegzav(i,k)*cfmzav(k)
            vegzav(0,k) = vegzav(0,k) + vegzav(i,k)
            cycle kloop
          endif
          if (k.ge.4 .and. k.le.6) then
            kk = (k-2)/2
            if (vegzav(i,kk).gt.epss) then
              vegzav(0,k) = vegzav(0,k) + vegzav(i,k)*darea(i)
              vegzav(i,k) = cfmzav(k)*vegzav(i,k)/vegzav(i,kk)
            else
              vegzav(i,k) = vegspv
            endif
          elseif (zavsum.gt.epss) then
            vegzav(0,k) = vegzav(0,k) + vegzav(i,k)*darea(i)
            vegzav(i,k) = cfmzav(k)*vegzav(i,k)/zavsum
          else
            vegzav(i,k) = vegspv
          endif
       enddo kloop
      enddo iloop

!     WRITE(*,*) 'NPP tot:',vegzav(0,18)
!     WRITE(*,'(E16.6)') (vegzav(i,18),i=1,nlat)

!-- Compute Global Mean Value :
      do k=nvgmax,1,-1
#if ( CYCC == 2 )
        if (((k.ge.12).AND.(k.le.18)).OR.(k.eq.20)) cycle
#else
        if ((k.ge.12).AND.(k.le.18)) cycle
#endif
        if (k.ge.4 .and. k.le.6) then
          kk = (k-2)/2
          if (vegzav(0,kk).gt.epss) then
            vegzav(0,k) = cfmzav(k)*vegzav(0,k)/vegzav(0,kk)
          else
            vegzav(0,k) = vegspv
          endif
        elseif (totsum.gt.epss) then
          vegzav(0,k) = cfmzav(k)*vegzav(0,k)/totsum
        else
          vegzav(0,k) = vegspv
        endif
      enddo

#if ( FAST_OUTPUT == 0 )
!-- Write Global & Zonal_mean values :
      write(vegetzav_id,'(3A,I6)') 'Veget. Zonal_mean Output ; ',titveg(:nnctr)
     &                   , ' ; year=', iyear+iyr0vg
      write(vegetzav_id,'(99A)') (titzav(k),k=1,nvgmax)
      do k=1,nvgmax
        write(vegetzav_id,fmtzav) (vegzav(i,k),i=0,nlat)
      enddo
#endif
      cland=vegzav(0,13) ! used in ipcc_out ...

#if ( FAST_OUTPUT == 0 )
      write(vegetzav_id,*)
#endif

!-- compute & write Global + Zonal mean : End. ----------
!-----------------------------------------------------------------------
#if ( FAST_OUTPUT == 0 )
      if (iyear+nyvegt.gt.nyears) then
        close(vegetzav_id)
      elseif (mod(iyear,nwrveg).eq.0) then
        call flush(vegetzav_id)
      endif
#endif
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  3 ) Sum 2.D vegetation map ; Prepare output .
!-----------------------------------------------------------------------

      vegsum = vegsum + 1.
      do j=1,nlon
       do i=1,nlat
         if (fracgr(i,j).gt.epss) then
           vegmap(i,j,1) = vegmap(i,j,1) + st(i,j)
           vegmap(i,j,2) = vegmap(i,j,2) + sg(i,j)
           vegmap(i,j,3) = vegmap(i,j,3) + sd(i,j)
           vegmap(i,j,4) = vegmap(i,j,4) + snlt(i,j)
           vegmap(i,j,5) = vegmap(i,j,5) + blai(i,j,1)
           vegmap(i,j,6) = vegmap(i,j,6) + blai(i,j,2)
           vegmap(i,j,7) = vegmap(i,j,7) +
     &  (albland(i,j,1)+albland(i,j,2)+albland(i,j,3)+albland(i,j,4))
           vegmap(i,j,8) = vegmap(i,j,8) + temveg(i,j)
           vegmap(i,j,9) = vegmap(i,j,9) + gd0veg(i,j)
           vegmap(i,j,10)= vegmap(i,j,10)+ prcveg(i,j,1)
           vegmap(i,j,11)= vegmap(i,j,11)+ prcveg(i,j,2)
           vegmap(i,j,12)= vegmap(i,j,12)+ (anup(i,j)*1E+12)
           vegmap(i,j,13)= vegmap(i,j,13)+ (stock(i,j)*1E+12)
           vegmap(i,j,14)= vegmap(i,j,14)+ (b1(i,j)*1E+12)
           vegmap(i,j,15)= vegmap(i,j,15)+ (b2(i,j)*1E+12)
           vegmap(i,j,16)= vegmap(i,j,16)+ (b3(i,j)*1E+12)
           vegmap(i,j,17)= vegmap(i,j,17)+ (b4(i,j)*1E+12)
           vegmap(i,j,18)= vegmap(i,j,18)+ (pnpp(i,j)*1E+12)
           vegmap(i,j,19)= vegmap(i,j,19)+ soiltype2(i,j)
           vegmap(i,j,20) = vegmap(i,j,20) + gd5veg(i,j)
! #if ( CYCC == 2 )
!            vegmap(i,j,20)= vegmap(i,j,20)+ (B1T13(i,j)*1E+12)
! #endif
         else
!          vegmap(i,j,19)= vegmap(i,j,19)+ soiltype2(i,j)
           gd0veg(i,j) = -900.
         endif
       enddo
      enddo
#if ( CLM_INDICES >= 2 )
       CALL SET_VEG_VARS(st,snlt,sg,sd,INT(iyear,kind=sip))
#endif       




      if (mod(iyear,nwrveg).eq.0) then
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! ==> decide to write 2.D map on file ; Prepare output :

!- x scaling Factor (if chg units) ; put SPecial_Value over non-land grid point
      vegsum = 1. / vegsum
      st_moy(:,:)=st_moy(:,:)+vegmap(:,:,1)*vegsum
      vegsum2=vegsum2+1
      if (iyear+nwrveg.gt.nyears) then
       vegsum2= 1. / vegsum2
       st_moy(:,:)=st_moy(:,:)*vegsum2
      endif
      do k=1,nvgmax
       do j=1,nlon
        do i=1,nlat
          if (fracgr(i,j).gt.epss) then
            vegmap(i,j,k) = cfmzav(k)*vegmap(i,j,k)*vegsum
          else
            st_moy(i,j) = 0.
            vegmap(i,j,k) = vegspv
          endif
        enddo
       enddo
      enddo
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  4 ) write 2.D vegetation map on ascii file "veget.outp" :
!-----------------------------------------------------------------------

! --- BdB 05-2019: calculate start time
      new_year_veg = new_year_veg + 1

!---- creating and writting netcdf file, !mohr
      CALL open(Yearly_Means)

      do k=1,numvegvar
         IF (output(newvegvar(k,4)))
     &   CALL write(k,vegmap(1:nlat,1:nlon,k))
      enddo

#if ( FAST_OUTPUT == 0 )
      do k=1,nvgmax
        if (kveget.eq.0 .and. k.ne.1 .and. k.le.6) cycle
        if (kveget.eq.0 .and. k.eq.19) cycle
!--Write 2 titles :
        if (k.eq.11) then
          write(vegetoutp_id,'(3A,F6.3,A)') titmap(k), titveg(:nnctr),
     &                            ' (>', prcmin*1000., ' mm/day)'
        else
          write(vegetoutp_id,'(2A)') titmap(k), titveg(:nnctr)
        endif
        write(vegetoutp_id,*) 'year=', iyear+iyr0vg
!--Write vegetation distribution :
        do i=1,nlat
          write(vegetoutp_id,fmtveg) (vegmap(i,j,k),j=1,nlon)
        enddo
        write(vegetoutp_id,*)
      enddo
#endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
#if ( FAST_OUTPUT == 0 )
!- End of writing :
      if (iyear+nwrveg.gt.nyears) then
        close(vegetoutp_id)
      else
        call flush(vegetoutp_id)
      endif
#endif

      CALL close() !closing netcdf file, !mohr

! --- BdB 05-2019: reset year counter for vegetation
      if (mod(iyear,nwrskip).eq.0 .or. iyear.eq.nyears) then
       new_year_veg = 0
       current_int_veg = current_int_veg + nwrskip
      end if

!- Reset to zero (for next Sum):
      vegsum = 0.
      do k=1,nvgmax
       do j=1,nlon
        do i=1,nlat
          vegmap(i,j,k) = 0.
        enddo
       enddo
      enddo
      
!! if ( FAST_OUTPUT == 0 )
#if ( 0 )
      if (kwradd.ne.0) then
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  5 ) write additional 2.D map on ascii file "veget.+[no_yr]" :
!-----------------------------------------------------------------------
      if (mod(iyear,kwradd).eq.0) then

!--open & write 2.l Header :
        write(ccny,'(I6)') iyear+iyr0vg
        do k=1,len(ccny)
          if (ccny(k:k).eq.' ') ccny(k:k) = '0'
        enddo

        open(newunit=vegetccny_id
     >     , file='veget.+'//ccny, status='unknown')
        write(vegetccny_id,1000) fmtveg, vegspv, nlon, nlat, 5, 65
        write(vegetccny_id,1111) 0., 5.625, -87.19, 5.625, 1., 1., 0

!- write Number of days/yr with Dry Soil Condition, i.e. Soil.W < bmtdry
        write(vegetccny_id,'(3A,F5.3,A)') 
     >   'Nb_Day/yr with Dry soil cond.; ',
     >   titveg(:nnctr), ' (Soil.W <', bmtdry,'m)'
        write(vegetccny_id,*) 'year=', iyear+iyr0vg
        do i=1,nlat
          do j=1,nlon
            var(j) = vegspv
            if (fracgr(i,j).gt.epss) var(j) = tpsdry(i,j)
          enddo
          write(vegetccny_id,fmtveg) (var(j),j=1,nlon)
        enddo
  write(vegetccny_id,*)

!--
!- write 4 seasonal albedo map :
        do ns=1,4
          write(vegetccny_id,'(2A)') 
     >          'Season. Albedo (o/o) (no snow) ; ',
     >          titveg(:nnctr)
          write(vegetccny_id,*) 'year=', iyear+iyr0vg, '  saison=', ns
          do i=1,nlat
            do j=1,nlon
              var(j) = vegspv
              if (fracgr(i,j).gt.epss) var(j) = 100.0*albland(i,j,ns)
            enddo
            write(vegetccny_id,fmtveg) (var(j),j=1,nlon)
           enddo
     write(vegetccny_id,*)
          enddo
    
        close(vegetccny_id)

      endif
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!-- write additional 2.D map : End. -----
      endif

#endif
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!-- write 2.D map : End. ----------------
      endif

      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine veget_wr -
      end
