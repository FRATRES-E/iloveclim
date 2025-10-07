!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:57 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:57 CET 2009

      SUBROUTINE uvi(j)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!**AVANCEMENT IMPLICITE D UN PAS DE TEMPS BAROCLINE DE U ET V**
!  avec taux advection verticale implicite : txiadu.
!    et taux diffusion verticale implicite : txidfu.
!  advection suivant la verticale :  w(k) * [ V(k-1) + V(k) ] / 2
!  modif : 09/01/95

!! START_OF_USE_SECTION

      use global_constants_mod, only: dblp=>dp, ip

      use const_mod

      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod
!! END_OF_USE_SECTION

!--variables locales :
!     dimension ccwup(kmax), ccwdn(kmax)

      real(kind=dblp), dimension(imax,jmax):: wudtd2, aa, bb, cc
      real(kind=dblp), dimension(imax):: ffcor
      
!- variables complexes :
      complex(kind=8), dimension(imax,kmax) :: acplex, bcplex, vcplex


!--- more locales

      real(kind=dblp):: ccdif, cctdif, cctimp, ccz, corfac, vt
      integer(kind=ip):: i, j, k

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  0 ) Mise en place de coeff dependant de k :
!-------------------

      cctdif = dtu * txidfu
      corfac = 2.*txifcu*dtu

!     do k=ks1+1,ks2
!       ccwup(k) =  dz(k-1) * unsdzw(k)
!       ccwdn(k) =  dz(k)   * unsdzw(k)
!     enddo

      ccz  = dtu * unsdz(ku2)

!     loop_500: do j=ju1,ju2
!**debut boucle sur la latitude** <-- Mise a l'Exterieur.

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1 ) Prise en compte des Flux fond & surface .                       |
!-----------------------------------------------------------------------

!-  gradient d'elevation  <- deplace dans uve .

!--conditions aux limites :
!--fond :
!     do i=iu1(j),iu2(j)
!       u(i,j,kfu(i,j))=u(i,j,kfu(i,j)) +dtu*unsdz(kfu(i,j))*phifu(i,j)
!       v(i,j,kfu(i,j))=v(i,j,kfu(i,j)) +dtu*unsdz(kfu(i,j))*phifv(i,j)
!     enddo

!--surface
      do i=iu1(j),iu2(j)
        u(i,j,ku2) = u(i,j,ku2) - ccz * phisu(i,j)
        v(i,j,ku2) = v(i,j,ku2) - ccz * phisv(i,j)
      enddo

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) Traitement Implicite : Calcul et Inversion matricielle .        |
!-----------------------------------------------------------------------

!**debut construction vitesses verticales (a un fact.mult. pres)
      cctimp = dtu * 0.125 * txiadu
      do k=ku1+1,ku2
       do i=iu1(j),iu2(j)
        wudtd2(i,k) = cctimp *
     &            ( (w(i-1,j-1,k)+w(i,j,k))+(w(i-1,j,k)+w(i,j-1,k)) )
!    &                (w(i,j,k)+w(i-1,j,k)+w(i,j-1,k)+w(i-1,j-1,k))
       enddo
      enddo
!**fin construction vitesses verticales**

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--Traitement qqmvt semi-implicite : coriolis, advection, diffusion.

!--construction I : partie reelle (advect.+diff.) de la matrice du systeme .
      do k=ku1,ku2
       do i=iu1(j),iu2(j)
        aa(i,k)= 1.0
       enddo
      enddo
      do k=ku1+1,ku2
       do i=iu1(j),iu2(j)
        ccdif = cctdif * avudz(i,j,k)
!       wwup  = ccwup(k) * wudtd2(i,k)
!       wwdn  = ccwdn(k) * wudtd2(i,k)
!- effet du flux phiz(k) sur V(k) : a(k)*V(k) + b(k)*V(k-1)
        vt = unsdz(k) * tmu(i,j,k-1)
        bb(i,k) = vt * (-ccdif - wudtd2(i,k))
        aa(i,k) = aa(i,k) +  vt * ( ccdif - wudtd2(i,k))
!       bb(i,k) = vt * (-ccdif - wwdw)
!       aa(i,k) = aa(i,k) +  vt * ( ccdif - wwup)
!- effet du flux phiz(k) sur S(k-1) : a(k-1)*V(k-1) + c(k-1)*V(k)
        vt = unsdz(k-1) * tmu(i,j,k-1)
        aa(i,k-1) = aa(i,k-1) + vt * (ccdif + wudtd2(i,k))
        cc(i,k-1) = vt * (-ccdif + wudtd2(i,k))
!       aa(i,k) = aa(i,k) +  vt * ( ccdif + wwdn)
!       cc(i,k) = vt * (-ccdif + wwup)
       enddo
      enddo

      if (cdbot.ne.zero) then
!- effet de le tension de fond
      do  i=iu1(j),iu2(j)
        aa(i,kfu(i,j)) = aa(i,kfu(i,j)) + avudz(i,j,1)
      enddo
      endif

!**fin construction de la matrice du systeme** --(partie reelle)

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--construction II : partie imaginaire (coriolis) + decomposition LU :

!**debut decomposition LU** --(cf: Numerical Methods, pp 165-167)
!- calcul de 1/alpha(k) dans acplex(i,k) et de beta(k) dans bcplex(i,k)
      do i=iu1(j),iu2(j)
       ffcor(i) = corfac * fs2cor(i,j)
       acplex(i,ku1) = 1.0 / DCMPLX( aa(i,ku1), ffcor(i) )
      enddo
      do k=ku1+1,ku2
       do i=iu1(j),iu2(j)
        bcplex(i,k) = bb(i,k) * acplex(i,k-1)
        acplex(i,k) = 1.0 / ( DCMPLX( aa(i,k), ffcor(i) )
     &                      - bcplex(i,k)*cc(i,k-1) )
       enddo
      enddo
!**fin decompostion LU **

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  3 ) Solution du Systeme Complex - Conversion Relle .                |
!-----------------------------------------------------------------------

!- charge f(k) dans vcplex(i,k) :
      do k=ku1,ku2
       do i=iu1(j),iu2(j)
        vcplex(i,k) = DCMPLX( u(i,j,k), v(i,j,k) )
       enddo
      enddo

!**debut substitutions avant et arriere**
!--calcul de g(k) dans vcplex(i,k) :
      do k=ku1+1,ku2
       do i=iu1(j),iu2(j)
        vcplex(i,k) = vcplex(i,k) - bcplex(i,k) * vcplex(i,k-1)
       enddo
      enddo

!--calcul de x(k) dans vcplex(i,k) :
      do i=iu1(j),iu2(j)
       vcplex(i,ku2) = vcplex(i,ku2) * acplex(i,ku2)
      enddo
      do k=ku2-1,ku1,-1
       do i=iu1(j),iu2(j)
        vcplex(i,k) = (vcplex(i,k)-cc(i,k)*vcplex(i,k+1)) * acplex(i,k)
       enddo
      enddo

!--mise en place des valeurs reelles dans u et v :
      do k=ku1,ku2
       do i=iu1(j),iu2(j)
        u(i,j,k) = tmu(i,j,k) * DREAL(vcplex(i,k))
        v(i,j,k) = tmu(i,j,k) * DIMAG(vcplex(i,k))
       enddo
      enddo
!**fin substitutions avant et arriere**

!     enddo loop_500
!**fin boucle sur la latitude**

      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine uvi -
      end
