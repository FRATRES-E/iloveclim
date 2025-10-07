      subroutine iceb_moy(i,j,nit)
c mab: i,j define the position in the grid, nit=1,ntotday given by emic.f
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--ancien (23/05/00) nom = moyenne(.f)
c  modif : 10/02/00
!- modif 24/10/19: cleaned use declaration sections              

!! START_OF_USE_SECTION
      use global_constants_mod, only: dblp=>dp, ip             
      use const_mod, only: rho0
      use para0_mod, only:
      use para_mod, only:
      use bloc0_mod, only: ks2, dts, dx, dy, tms
      use bloc_mod, only: cmx, cmy
      use ice_mod, only: jmax, imax
      use iceberg_mod, only: dticeb,it_icb,nbrmoins_cum,nbrmoins_cum_n,
     >    nbrmoins_cum_s, nbrplus_cum, nbrplus_cum_n, nbrplus_cum_s,
     >    nbrmois, mois, nbrjour, fonte_icb, vol_icb, nbricb_moins,
     >    nbricb_moins_s, nbricb_plus, nbricb_plus_n, nbricb_plus_s,
     >    nitrunicb,nbricb_moins_n,roiceb,nbricb,nbricb_n,nbricb_s,lmx,
     >    fonte_zon,vol_zon, fonte_icb_cum,fonte_icb_mois,vol_icb_cum
      use reper_mod, only: xwi1,dxwi, ywj1, dywj
!! END_OF_USE_SECTION

      implicit none
      
!! START_OF_INCLUDE_SECTION

c-include
!PB#include "type.com"
! [SCRPTCOM] #include "para.com"
! [SCRPTCOM] #include "bloc.com"
! [SCRPTCOM] #include "iceberg.com"
! [SCRPTCOM] #include "reper.com"
! [SCRPTCOM] #include "ice.com"

!! END_OF_INCLUDE_SECTION

c--variables locales :
      character(len=4) :: fchnum
      integer(kind=ip) :: njcum(0:12)

!PB variables added after imposing implicit none
      integer(kind=ip) ::  j,i,ii,nit,k
      real(kind=dblp) :: deltat,fonte_nord,fonte_sud,vol_icb_nord,
     >          vol_icb_sud,surface_n,surface_s

c-formats
 1000 format(A30,1x,E13.6,1x,I6,1x,I6,1x,I6,1x,I6)
 1111 format(3(F7.1,1X,F7.3,1X),I3,A)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

C     data njcum/0,31,59,90,120,151,181,212,243,273,304,334,365/
      data njcum/0,30,60,90,120,150,180,210,240,270,300,330,360/
c mab:dticeb is defined in positb.init and read in deficeb.f and is 86400
c dts() is defined in clio/parameters/run.param and there is set to 1.at every
c layer
      deltat = dticeb
C     deltat = dts(ks2)

c mab:it_icb in iceberg.f each time modified (+1); defines the time
c step of the iceberg?!

      if (it_icb.eq.1) then
        nbrmoins_cum=0
        nbrmoins_cum_n=0
        nbrmoins_cum_s=0
        nbrplus_cum=0
        nbrplus_cum_n=0
        nbrplus_cum_s=0
        nbrmois=1
        do 215 j=1,jmax
         fonte_zon(j)=0.
         vol_zon(j)=0.
 215    continue
       mois=1
       do 205 j=1,jmax
          do 205 i=1,imax
             fonte_icb_cum(i,j)=0.
             fonte_icb_mois(i,j)=0.
             vol_icb_cum(i,j)=0.
c            fonte_icb_as(i,j)=0.
c            fonte_icb_fm(i,j)=0.
c            vol_icb_as(i,j)=0.
c            vol_icb_fm(i,j)=0.
 205   continue
       do 10 ii=1,nbrjour
          if (nbrjour.gt.njcum(mois)) mois=mois+1
 10    continue
       if (mois.eq.13) mois=1
      endif

c-calcul de la fonte totale mensuelle et du volume de glace total
c-a la fin du mois pour le nord et le sud(attention, la fonte
c-est en m**3 de glace.
c-----------------------------------------------------------------

      do 300 i=1,imax
        do 300 j=1,jmax
      fonte_icb_mois(i,j)=fonte_icb_mois(i,j)+fonte_icb(i,j)
 300  continue
      if ((nbrjour.eq.njcum(mois)).and.((it_icb*dticeb).eq.
     &         (nit*dts(ks2)))) then
      fonte_nord=0.
      fonte_sud=0.

      do 212 j=1,int(jmax/2)
        do 212 i=1,imax
           fonte_sud=fonte_sud+fonte_icb_mois(i,j)
 212  continue

      do 213 j=(int(jmax/2)+1),jmax
        do 213 i=1,imax
          fonte_nord=fonte_nord+fonte_icb_mois(i,j)
 213  continue
         vol_icb_nord=0.
         vol_icb_sud=0.
         do 210 j=1,int(jmax/2)
           do 210 i=1,imax
             vol_icb_sud=vol_icb(i,j)+vol_icb_sud
 210     continue
         do 211 j=(int(jmax/2)+1),jmax
           do 211 i=1,imax
             vol_icb_nord=vol_icb(i,j)+vol_icb_nord
 211     continue
       endif

c-b.convertion du tableau de fonte(m**3) en kg d'eau par sec par m**2
c-et de celui de volume en m
c mab: where is the convertion made?
c---------------------------------------------------------------------

c-calcul des tableaux de glace et de fonte cumules sur un an.
c-------------------------------------------------------------

      if((nbrjour.eq.njcum(mois)).AND.
     &    ((it_icb*dticeb).eq.(nit*dts(ks2)))) then
C    &    (mod(it_icb,int(86400/deltat)).eq.0)) then
      do 206 j=1,jmax
        do 206 i=1,imax
          fonte_icb_cum(i,j)=fonte_icb_cum(i,j)+(fonte_icb_mois(i,j))
          vol_icb_cum(i,j)=vol_icb_cum(i,j)+(vol_icb(i,j))
 206  continue
       do 214 j=1,jmax
         do 214 i=1,imax
           fonte_zon(j)=fonte_zon(j)+(fonte_icb_mois(i,j))
           vol_zon(j)=vol_zon(j)+vol_icb(i,j)
 214   continue
       endif

      nbrmoins_cum=nbrmoins_cum+nbricb_moins
      nbrmoins_cum_n=nbrmoins_cum_n+nbricb_moins_n
      nbrmoins_cum_s=nbrmoins_cum_s+nbricb_moins_s
      nbrplus_cum=nbrplus_cum+nbricb_plus
      nbrplus_cum_n=nbrplus_cum_n+nbricb_plus_n
      nbrplus_cum_s=nbrplus_cum_s+nbricb_plus_s

c-d.ecriture dans le fichier icbwp.out a la fin de l'annee
c--------------------------------------------------------
c-attention,nfrqicb ne peut etre plus petit qu'un mois
c- sinon on divise par zero, logiquement,nfrqicb devrait etre
c-egal a nitrun ainsi on a une moyenne sur tout le run.

      if(mod(it_icb,nitrunicb).eq.0) then
C     if(mod(nit,nitrun).eq.0) then
        do 207 j=1,jmax
          do 207 i=1,imax
            fonte_icb_cum(i,j)=(fonte_icb_cum(i,j)*roiceb/
     &                     (cmx(i,j,0)*cmy(i,j,0)*dx*dy))
 207    continue
c mab: cmx,cmy are metric coefficients (i,j,0)=surface (cmy=1,cmx=1 or defined on
c a sperhical grid
        do 208 j=1,jmax
          do 208 i=1,imax
            vol_icb_cum(i,j)=(vol_icb_cum(i,j)/(cmx(i,j,0)
     &                        *cmy(i,j,0)*dx*dy))
 208    continue
c JONO_out file was called icbwp.out, 4-2004 now icb_yrly_font.out
        open(1061,file='outputdata/icebergs/icb_yrly_font.out'
     &   ,status='unknown')
        write(1061,1000) '(1P125E8.1)', 0., imax,-jmax,2,imax
        write(1061,1111) xwi1, dxwi, ywj1, dywj, 0., 1., 0
        write(1061,'(2A)')'fonte(m/an):{if(mod(it_icb,nitrunicb)=0)',
     &';12fonte_icb_cum(k,j)[+=fi_moins]
     &/rho0*nbrmoins,k=1,imax;j=jmax,1,-1}'
        write(1061,*)
        do 1502 j=jmax,1,-1
            write(1061,'(1P125E8.1)') (((12*fonte_icb_cum(k,j))/(rho0
     &           *nbrmois)),k=1,imax)
 1502   continue
         write(1061,*)
        surface_n=0.
         surface_s=0.
        do 1510 j=1,int(jmax/2)
          do 1510 i=1,imax
            if(fonte_icb_cum(i,j).ne.0.) then
             surface_s=surface_s+(tms(i,j,ks2)*cmx(i,j,0)
     &       *cmy(i,j,0)*dx*dy)
            endif
 1510   continue
        do 1511 j=(int(jmax/2)+1),jmax
          do 1511 i=1,imax
            if(fonte_icb_cum(i,j).ne.0.) then
             surface_n=surface_n+(tms(i,j,ks2)*cmx(i,j,0)
     &       *cmy(i,j,0)*dx*dy)
            endif
 1511   continue
       write(6,*)'La surface au nord:',surface_n
       write(6,*)'La surface au sud:',surface_s
c JONO_out file was fort.61, 4-04 now icb_yrly_vol61.out
       open(61,file='outputdata/icebergs/icb_yrly_vol61.out'
     &  ,status='unknown')
       write(61,*)
       write(61,'(2A)') 'volume (m):{if(mod(it_icb,nitrunicb)=0)',
     &';12*vol_icb_cum(k,j)[+=vol_icb()]/nbrmoins,k=1,imax;j=jmax,1,-1}'
       write(61,*)
       do 1505 j=jmax,1,-1
         write(61,'(1P125E8.1)') (((12*vol_icb_cum(k,j))/nbrmois)
     &               ,k=1,imax)
 1505  continue
       write(61,*)
       close(61)
      endif

      if((nbrjour.eq.njcum(mois)).and.((it_icb*dticeb).eq.
     &          (nit*dts(ks2)))) then

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|
c-ecriture chaque mois dans le fichier "nbreiceb.out"
c-COL1:       "mois" est le mois courant (1->12)
c-COL2-3-4:   "nbricb/_n/_s" est le nombre d'icb/pour HN/HS en
c-circulation
c-COL5-6:     "vol_icb_nord/sud" est le volume total de glace pour
c-l'HN/HS
c-COL7-8:     "fonte_nord/sud" est le volume d'eau de fonte total pour
c-l'HN/HS
c-COL9-10-11: "nbrmoins_cum/_n/_s" est le nombre d'icb disparus sur le
c-mois
c-COL12-13-14:"nbrplus_cum/_n/_s" est le nombre d'icb generes sur le
c-mois
c-COL15:      "lmx" est la taille du tableau d'icb (= nombre reel d'icb
c-modelise en tenant compte du facteur de ponderation (Cfr gener_iceb.f)

         write(1055,'(4I8,1P4E10.2,7I6)')mois,nbricb,nbricb_n,nbricb_s,
     &    vol_icb_nord,vol_icb_sud,fonte_nord,fonte_sud,
     &    nbrmoins_cum,nbrmoins_cum_n,nbrmoins_cum_s,
     &    nbrplus_cum,nbrplus_cum_n,nbrplus_cum_s,lmx

       call flush(1055)
      nbrmoins_cum=0
      nbrmoins_cum_n=0
      nbrmoins_cum_s=0
      nbrplus_cum=0
      nbrplus_cum_n=0
      nbrplus_cum_s=0

      mois=mois+1
      if(mois.eq.13)mois=1
      nbrmois=nbrmois+1
      open(1062,file='outputdata/icebergs/diagzon.out',status='unknown')
      do 216 j=1,jmax
        write(1062,'(I3,1P2E12.2)')j,fonte_zon(j),vol_zon(j)
 216  continue
      close(1062)
      do 217 j=1,jmax
        fonte_zon(j)=0.
        vol_zon(j)=0.
 217  continue
      do 209 i=1,imax
        do 209 j=1,jmax
          fonte_icb_mois(i,j)=0.
 209  continue
      endif

      do 204 i=1,imax
        do 204 j=1,jmax
        vol_icb(i,j)=0.
 204  continue
C     write(6,*)mois,nbrmois,nbrjour,njcum(mois)
C     write(6,'(3I8)')nbricb,nbricb_n,nbricb_s
c--calcul de la surface d'eau sur laquelle il y a des icebergs
c JONO  this surface is area(i,j), I am commenting out NICs stuff here,
c     because it might be flawed (depending on timely resetting of fonte_icb())
c     and I do this transformation in clio/sources/ec_co2oc.f
c Nic changement des unites de mesure de m**3/an en kg/m**2.an
c      do 500 j=1,jmax
c         do 500 i=1,imax
c          fonte_icb(i,j)=(fonte_icb(i,j)*roiceb/(cmx(i,j,0)*cmy(i,j,0)*dx*dy))
c 500  continue

      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine iceb_moy -
      end subroutine iceb_moy
