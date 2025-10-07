!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:48 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:48 CET 2009

      SUBROUTINE iso_aver(ns,ilat,name,phix,phiy,xmult)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- Adapted form streamv.f (modif : 26/05/96)
!- Compute components of the PHT (overturning, diffusion, eddies..)
!- ppmodif : flux phi 3D , ntrsm=>ilat,
!  modif : 30/03/99

!! START_OF_USE_SECTION
      use global_constants_mod, only: dblp=>dp, ip           

      use const_mod

      use para0_mod
      use para_mod
      use vareq_mod, only: spv
      use bloc0_mod
      use bloc_mod
      use reper_mod
      use varno_mod

      use newunit_clio_mod, only: clio3_out_id      
!! END_OF_USE_SECTION

      real(kind=dblp), dimension(jmax,0:nbsmax) :: trsm
      real(kind=dblp), dimension(imax,jmax,kmax):: phiy, phix, tt

      character*30 fmt
      character*70 line
      character*40 tit
      character*10 name

!--- locales
      integer(kind=ip):: i, ii, iilat, ilat, iover, j, jcr, jcr1, jj
     >                 , jj1, jj2, jsigne, k, nb, ns, iim1, nbas      
      real(kind=dblp) :: sum, sumt, tav, tms1, tms2, tt1, tt2, ttm
     >                 , xmult, xnbr

      integer(kind=ip):: idl_overgm_id
      
      iover=0
      if (name.eq.'idl.over'.or.name.eq.'idl.overgm') iover=1
      tit=' *** iso_aver: ppmodif: 14-04-97:  ***'

      iilat = min(ilat,jmvlat)
      if (iilat.eq.0) write(clio3_out_id,*) tit//' : J-GRID   AVERAGE '//name
      if (iilat.eq.1) write(clio3_out_id,*) tit//' : True Lat AVERAGE '//name

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      do nb=0,nbsmax
      do j=1,jmax
       trsm(j,nb)=zero
      enddo
      enddo
      do k=1,kmax
      do j=1,jmax
      do i=1,imax
       tt(i,j,k)=scal(i,j,k,ns)
      enddo
      enddo
      enddo

      do nb=0,nbsmax
      do k=ku1,ku2

      if (iilat.eq.0) then
!-****************
!- MODEL Latitude
!-****************
        do j=ju1,ju2
          sum=zero
          do i=isbas(j,nb),iebas(j,ks2,nb)
            ii=icl(i)
            sum=sum+phiy(ii,j,k)
          enddo

!- Compute Tav zonal for Phiy
         if (iover.eq.1) then
          sumt=zero
          xnbr=zero
          do i=isbas(j,nb),iebas(j,ks2,nb)
            ii=icl(i)
            tms1=tms(ii,j  ,k)
            tt1 = tt(ii,j  ,k)
            tms2=tms(ii,j-1,k)
            tt2=  tt(ii,j-1,k)
            ttm=tms1*tms2
            tav=(tms1*tt1+tms2*tt2)/(tms1+tms2+epsil)
            xnbr=xnbr+ttm
            sumt=sumt+ttm*tav
           enddo
           if (xnbr.GT.epsil) sum=sum*(sumt/xnbr)
         endif

         trsm(j,nb)=trsm(j,nb)+sum*dz(k)*dy
       enddo

      elseif (iilat.eq.1) then
!-***************************************************************
!- TRUE latitude JMAXL (streamv.f)
!- separate X & Y component (!! sign for x component)
!-***************************************************************

!-// Model-grid until Equator //
        do j=ju1,jeq-1
         sum=zero
         do i=isbas(j,nb),iebas(j,ks2,nb)
          ii=icl(i)
          sum=sum+phiy(ii,j,k)
         enddo

!- Compute Tav zonal
         if (iover.eq.1) then
          sumt=zero
          xnbr=zero
          do i=isbas(j,nb),iebas(j,ks2,nb)
            ii=icl(i)
            tms1=tms(ii,j  ,k)
            tt1=  tt(ii,j  ,k)
            tms2=tms(ii,j-1,k)
            tt2=  tt(ii,j-1,k)
            ttm=tms1*tms2
            tav=(tms1*tt1+tms2*tt2)/(tms1+tms2+epsil)
            xnbr=xnbr+ttm
            sumt=sumt+ttm*tav
           enddo
           if (xnbr.GT.epsil) sum=sum*(sumt/xnbr)
         endif

         trsm(j,nb)=trsm(j,nb)+sum*dz(k)*dy
        enddo

!-// Y-Contribution //
        do j=jeq,jmvlat
         sum=zero
         do i=isbas(j,nb),iebas(j,ks2,nb)
          ii=icl(i)
          jcr=jmaxl(ii,j)
          sum=sum+phiy(ii,jcr,k)
         enddo

!- Compute Tav zonal (Phiy)
         if (iover.eq.1) then
          sumt=zero
          xnbr=zero
          do i=isbas(j,nb),iebas(j,ks2,nb)
            ii=icl(i)
            jcr=jmaxl(ii,j)
            tms1=tms(ii,jcr  ,k)
            tt1=  tt(ii,jcr  ,k)
            tms2=tms(ii,jcr-1,k)
            tt2=  tt(ii,jcr-1,k)
            ttm=tms1*tms2
            tav=(tms1*tt1+tms2*tt2)/(tms1+tms2+epsil)
            xnbr=xnbr+ttm
            sumt=sumt+ttm*tav
           enddo
           if (xnbr.GT.epsil) sum=sum*(sumt/xnbr)
         endif

         trsm(j,nb)=trsm(j,nb)+sum*dz(k)*dy
        enddo

!-// X-Contribution //
        do j=jeq,jmvlat
         sum=zero
         do i=1+isbas(j,nb),iebas(j,k,nb)
          ii=icl(i)
          jcr=jmaxl(ii,j)
          jcr1=jmaxl(icl(i-1),j)
          jsigne=sign(1,jcr1-jcr)
          jj1=min(jcr,jcr1)
          jj2=max(jcr,jcr1)-1
          do jj=jj1,jj2
           ttm=tms(icl(i-1),jj,k)*tms(ii,jj,k)
           sum=sum+DFLOAT(jsigne)*phix(ii,jj,k)
          enddo
         enddo

!- Compute Tav zonal (Phix)
         if (iover.eq.1) then
          sumt=zero
          xnbr=zero
          do i=1+isbas(j,nb),iebas(j,ks2,nb)
            ii=icl(i)
            iim1=icl(i-1)
            jcr=jmaxl(ii,j)
            jcr1=jmaxl(icl(i-1),j)
            jj1=min(jcr,jcr1)
            jj2=max(jcr,jcr1)-1
            do jj=jj1,jj2
             tms1=tms(ii  ,jj,k)
             tt1=  tt(ii  ,jj,k)
             tms2=tms(iim1,jj,k)
             tt2=  tt(iim1,jj,k)
             ttm=tms1*tms2
             tav=(tms1*tt1+tms2*tt2)/(tms1+tms2+epsil)
             xnbr=xnbr+ttm
             sumt=sumt+ttm*tav
            enddo
          enddo
          if (xnbr.GT.epsil) sum=sum*(sumt/xnbr)
         endif

         trsm(j,nb)=trsm(j,nb)+sum*dz(k)*dx
        enddo


      endif

      enddo
      enddo

!- Insert Special Value
       do nb=0,nbsmax
       do j=1,jmax
        if (isbas(j,nb).le.iebas(j,ks2,nb)) then
          trsm(j,nb)=trsm(j,nb)*xmult
        else
          trsm(j,nb)=spv(0)
        endif
        if (j.lt.jsbas(nb).or.j.gt.jebas(nb)) trsm(j,nb)=spv(0)
       enddo
       enddo


!- Save Transport
      open(newunit=idl_overgm_id,file=name,STATUS='unknown')
      do nbas=0,nbsmax
       write(idl_overgm_id,*) 'BASSIN',nbas
       write(idl_overgm_id,'(70E11.3)') (trsm(j,nbas),j=1,jmax)
      enddo
      close(idl_overgm_id)

      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine iso_aver -
      end
