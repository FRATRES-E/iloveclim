!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:57 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:57 CET 2009

      SUBROUTINE vague(nflag,tab)
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  evaluation du "bruit 2dx" sur l'elevation.
!   par la formule (eta(i)-eta(i+1))^2+(eta(j)-eta(j+1))^2
!  ecaluation du mode checkboard par etachk = Moyenne [ eta(i,j)*(-1)^(i+j) ]
!  modif : 28/07/94

!! START_OF_USE_SECTION

      use global_constants_mod, only: dblp=>dp, ip

      use const_mod

      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod

      use newunit_clio_mod, only: mouchard_id      
      
!! END_OF_USE_SECTION

!--dummy argument :
      integer(kind=ip) :: nflag
      real(kind=dblp), dimension(*) :: tab
      
!--variables locale :
      real(kind=dblp), dimension(imax,jmax) :: etabar
      
      
!--- more locales
      integer(kind=ip):: i, ii, j, jj, k

      real(kind=dblp) :: bruitx, bruity, ccs, eta2av, etac, etachk, etam
     &                  , etar, etares, fact, grdx, grdy, qq, qqmi, qqmx
     &                  , sumc, summ, surf, surfx, surfy, etabam, etadyn


!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!--initialisation :
      bruitx = 0.
      bruity = 0.
      eta2av =  0.
      etachk =  0.
      etares =  0.

!--calcul direct de la somme du carre de la 1ere composantes du gradient :
      surfx = 0.
      do j=js1,js2
       do i=is1(j),is2(j)
        ccs = tms(i-1,j,ks2) * tms(i,j,ks2) * cmx(i,j,1) * cmy(i,j,1)
        grdx = smx(i,j,1) * unsdx * (eta(i-1,j)-eta(i,j))
        surfx  = surfx  + ccs
        bruitx = bruitx + ccs * grdx * grdx
       enddo
      enddo
      tab(2) = bruitx / surfx

!--calcul direct de la somme du carre de la 2nd  composantes du gradient :
      surfy = 0.
      do j=js1,js2
       do i=is1(j),is2(j)
        ccs = tms(i,j-1,ks2) * tms(i,j,ks2) * cmx(i,j,2) * cmy(i,j,2)
        grdy = smy(i,j,2) * unsdy * (eta(i,j-1)-eta(i,j))
        surfy  = surfy  + ccs
        bruity = bruity + ccs * grdy * grdy
       enddo
      enddo
      tab(3) = bruity / surfy
!     bruit = sqrt( bruitx + bruity )

!--calcul direct de la somme du carre de l'elevation :
      surf = 0.
      do j=js1,js2
       do i=is1(j),is2(j)
        ccs = tms(i,j,ks2) * cmxy(i,j,0)
        surf =  surf + ccs
        eta2av =  eta2av + ccs *  eta(i,j) * eta(i,j)
       enddo
      enddo
      tab(1) =  eta2av / surf

!--calcul direct de la somme/difference alternee de l'elevation :
      do j=js1,js2
       do i=is1(j),is2(j)
!- calcul local de etac :
        etam = 0.
        etac = 0.
        summ = 0.
        sumc = 0.
        fact = 1.
        do jj=j-1,j+1
         do ii=i-1,i+1
          etam = etam + tms(ii,jj,ks2) * eta(ii,jj)
          summ = summ + tms(ii,jj,ks2)
          etac = etac + tms(ii,jj,ks2) * fact * eta(ii,jj)
          sumc = sumc + tms(ii,jj,ks2) * fact
          fact = -fact
         enddo
        enddo
        if(summ.ne.0.) then
          etam = etam / summ
          etac = (etac - sumc * etam) / summ
        endif
        etar = eta(i,j) - etac
        ccs = tms(i,j,ks2) * cmxy(i,j,0)
        etachk = etachk + ccs * etac * etac
        etares = etares + ccs * etar * etar
       enddo
      enddo
      tab(4) = etachk / surf
      tab(5) = etares / surf

      if(nflag.eq.4.or.nflag.eq.0) then
!--calcul de la contribution barocline a l'elevation :
        qqmi =  1.E9
        qqmx = -1.E9
        do j=js1,js2
         do i=is1(j),is2(j)
           qq = 0.
           if (tms(i,j,ks2).eq.1) then
             do k=ks2,kfs(i,j),-1
               qq = qq + q(i,j,k)*dz(k)
             enddo
             qq = qq / zw(kfs(i,j)) / gpes
             qqmi = min(qq,qqmi)
             qqmx = max(qq,qqmx)
           endif
           etabar(i,j) = qq
         enddo
        enddo

        etabam = 0.
        etadyn = 0.
        do j=js1,js2
         do i=is1(j),is2(j)
          ccs = tms(i,j,ks2) * cmxy(i,j,0)
          qq = eta(i,j) - etabar(i,j)
          etabam = etabam + ccs * etabar(i,j) * etabar(i,j)
          etadyn = etadyn + ccs * qq * qq
         enddo
        enddo
        tab(6) = etabam / surf
        tab(7) = etadyn / surf
        tab(8) = qqmi
        tab(9) = qqmx
      endif

!     write(mouchard_id,*) 'surfx, surfy, surf :'
!     write(mouchard_id,*) surfx, surfy, surf
!     write(mouchard_id,*) 'qqmin, qqmax : '
!     write(mouchard_id,*) qqmi, qqmx

      return
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
      end
