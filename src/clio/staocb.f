!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:54 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:54 CET 2009

      SUBROUTINE staocb(nnt,ccfile)
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  redemarrage a partir de l'etat definit par le fichier binaire "rest.om"
!  et de conditons arbitraires sur la glace.
!  modif : 20/12/93


!! START_OF_USE_SECTION

      use global_constants_mod, only: dblp=>dp, ip

      use const_mod

      use para0_mod
      use para_mod
      use moment_mod, only: vicmom, copy_from_vicmon

      use bloc0_mod
      use bloc_mod
      use ice_mod
      use dynami_mod

      use newunit_clio_mod, only: clio3_out_id, mouchard_id

!! END_OF_USE_SECTION

!! START_OF_INCLUDE_SECTION


!- dummy variables :
      character*(*) ccfile
!- variables locales equivalentes :
      real(kind=dblp), dimension(imax,jmax,kmax) :: wloc
!dmr [NOEQUI]      equivalence ( w(1,1,1) ,  wloc(1,1,1) )

!- local variables :
      character*6 cc6exp

!--- loads of locales
      real(kind=dblp) ::  alinn, alins, egajc, hginn, hgins, hmajc
     &                  , hninn, hnins, ttest, zidto
      integer(kind=ip) :: i, j, k, kmaxp1, nnt

      integer(kind=ip) :: ccfile5_id, inice_param_id
      kmaxp1 = kmax + 1
      if (nnt.ge.2) then
!--initialisation de w(kmax+1) :
        do j=1,jmax
         do i=1,imax
         w(i,j,kmaxp1) = 0.
         enddo
        enddo
        endif

!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  1 ) Bloc d ouverture et de Lecture du fichier binaire "rest.om" .   |
!-----------------------------------------------------------------------

      open(newunit=ccfile5_id,
     &   file=ccfile,status='old',form='UNFORMATTED')

!--1.1 lecture de la premiere partie :
!-------------------------------------

      read (ccfile5_id) numit
      read (ccfile5_id) tpstot
      read (ccfile5_id) eta
      read (ccfile5_id) ub
      read (ccfile5_id) vb
      read (ccfile5_id) u
      read (ccfile5_id) v
      read (ccfile5_id) scal
!- option avec ou sans TKE :
!     if (kstart.eq.3) then
!       write(66,*) 'Pas de lecture de TKE '
!     else
!       read(ccfile5_id) tked
!     endif

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! correction celcius-kelvin
      do k=1,kmax
       do j=1,jmax
         do i=1,imax
!           scal(i,j,k,1) = scal(i,j,k,1) + 273.15
         enddo
       enddo
      enddo
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if (nnt.ge.2) then

!--1.2. lecture de la seconde partie :
!-------------------------------------

        read (ccfile5_id) cc6exp
        read (ccfile5_id) q
        read (ccfile5_id) bvf
        read (ccfile5_id) avsdz
        read (ccfile5_id) avudz
        read (ccfile5_id) wloc
        read (ccfile5_id) egajc
        read (ccfile5_id) hmajc
!--fin de la seconde partie .
        write(clio3_out_id,'(3A,I11,A)') 'Exp '//refexp//' , File ', ccfile,
     &      ' , Iter', numit, ' : FULL RESULTS read - OK.'

      elseif (refexp.eq.'      ') then
        read (ccfile5_id) refexp
        write(clio3_out_id,
     &   '(3A,I11,A)') 'Lecture de ', ccfile, ' terminee (Exp '
     &                    //refexp//', Iter', numit, ' ) .'
      else
        write(clio3_out_id,
     &   '(3A,I11,A)') 'Lecture de ', ccfile, ' terminee (Iter',
     &       numit, ' ) .'
      endif

      close(ccfile5_id)

!dmr [NOEQUI]
      w(1:imax,1:jmax,1:kmax) = wloc(:,:,:)
!dmr [NOEQUI]

!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  2 ) Traitement du cas d'un changement de stockage .                 |
!-----------------------------------------------------------------------

      if(nnt.lt.0) then
        do j=1,jmax
         do i=1,imax
          ub(i,j) = ub(i,j) * hu(i,j)
          vb(i,j) = vb(i,j) * hu(i,j)
         enddo
        enddo
      endif

!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  3 ) Sea ice initial conditions                                      |
!-----------------------------------------------------------------------
!
      numit=0
!     tpstot=86400*150
      tpstot=0.0
!
!  File of parameters.
!
      open(newunit=inice_param_id,file='inice.param')
      read(inice_param_id,*)
      read(inice_param_id,*)
      read(inice_param_id,*)
      read(inice_param_id,*)
      read(inice_param_id,*)

      read(inice_param_id,*)
      read(inice_param_id,*) ttest
      read(inice_param_id,*)
      read(inice_param_id,*) hninn
      read(inice_param_id,*)
      read(inice_param_id,*) hginn
      read(inice_param_id,*)
      read(inice_param_id,*) alinn
      read(inice_param_id,*)
      read(inice_param_id,*) hnins
      read(inice_param_id,*)
      read(inice_param_id,*) hgins
      read(inice_param_id,*)
      read(inice_param_id,*) alins
      close(inice_param_id)

        do j=jeq,js2
          do i=is1(j),is2(j)
!
!  Northern hemisphere.
!
!  tfu: Melting point of sea water.
!
!
            tfu(i,j)  = abs(273.15-0.0575*scal(i,j,ks2,2)+
     &                  1.710523e-03*sqrt(scal(i,j,ks2,2))**3-
     &                        2.154996e-04*scal(i,j,ks2,2)**2)
!
!          Criterion for presence (zidto=1) or absence (zidto=0) of ice.
!
            zidto     = tms(i,j,ks2)*(1.0-max(zero,
     &                  sign(one,scal(i,j,ks2,1)-tfu(i,j)-ttest)))
!
            scal(i,j,ks2,1)  = zidto*tfu(i,j)
     &                         +(1.0-zidto)*scal(i,j,ks2,1)
            scal(i,j,ks2-1,1)  = zidto*tfu(i,j)
     &                           +(1.0-zidto)*scal(i,j,ks2-1,1)
            hgbq(i,j)   = zidto*hginn
            albq(i,j)   = zidto*alinn+(1.0-zidto)*1.0
            hnbq(i,j)   = zidto*hninn
            ts(i,j)     = tfu(i,j)
            tbq(i,j,1)  = tfu(i,j)
            tbq(i,j,2)  = tfu(i,j)
            tbq(i,j,3)  = tfu(i,j)
            firg(i,j)   = 0.0
            fcsg(i,j)   = 0.0
            fleg(i,j)   = 0.0
            fsbbq(i,j)  = 0.0
            qstobq(i,j) = 0.0
            xzo(i,j)    = 0.001d0
            ug(i,j)     = 0.0
            vg(i,j)     = 0.0
!          Moments for advection.
          do k=1,35
            vicmom(i,j,k) =0.0
          enddo
!
          enddo
        enddo

!dmr [NOEQUI]
       call copy_from_vicmon()

!        write(mouchard_id,*) tfu(10,js2-10),scal(10,js2-10,ks2,1)
!    &               ,scal(10,js2-10,ks2,2)
!           write(mouchard_id,*) hginn,ttest,hgbq(10,js2-10)

!
!  Southern hemisphere.
!
        do j=js1,jeq-1
          do i=is1(j),is2(j)
!
            tfu(i,j)  = abs(273.15-0.0575*scal(i,j,ks2,2)+
     &                  1.710523e-03*sqrt(scal(i,j,ks2,2))**3-
     &                  2.154996e-04*scal(i,j,ks2,2)**2)
            zidto     = tms(i,j,ks2)*(1.0-max(zero,
     &                        sign(one,scal(i,j,ks2,1)-tfu(i,j)-ttest)))
!
            scal(i,j,ks2,1) = zidto*tfu(i,j)
     &                        +(1.0-zidto)*scal(i,j,ks2,1)
            scal(i,j,ks2-1,1) = zidto*tfu(i,j)
     &                          +(1.0-zidto)*scal(i,j,ks2-1,1)
            scal(i,j,ks2-2,1) = zidto*tfu(i,j)
     &                          +(1.0-zidto)*scal(i,j,ks2-2,1)
            hgbq(i,j)   = zidto*hgins
            albq(i,j)   = zidto*alins+(1.0-zidto)*1.0
            hnbq(i,j)   = zidto*hnins
            ts(i,j)     = tfu(i,j)
            tbq(i,j,1)  = tfu(i,j)
            tbq(i,j,2)  = tfu(i,j)
            tbq(i,j,3)  = tfu(i,j)
            firg(i,j)   = 0.0
            fcsg(i,j)   = 0.0
            fleg(i,j)   = 0.0
            fsbbq(i,j)  = 0.0
            qstobq(i,j) = 0.0
            xzo(i,j)    = 0.001d0
            ug(i,j)     = 0.0
            vg(i,j)     = 0.0
!  Moments for advection.
!
            do k=1,35
              vicmom(i,j,k) =0.0
            enddo
          enddo
        enddo

!dmr [NOEQUI]
       call copy_from_vicmon()


      if (ltest.ge.1) then
!--raccord cyclique pour hgbq,albq,hnbq,ts,tbq,firg,fcsg,fleg,
!                     fsbbq,qstobq,scal:
      call raccord(hgbq(1,1),0.0,1,8)
      call raccord(albq(1,1),0.0,1,8)
      call raccord(hnbq(1,1),0.0,1,8)
      call raccord(ts(1,1),0.0,1,8)
      call raccord(tbq(1,1,1),0.0,nkb0,8)
      call raccord(firg(1,1),0.0,1,8)
      call raccord(fcsg(1,1),0.0,1,8)
      call raccord(fleg(1,1),0.0,1,8)
      call raccord(fsbbq(1,1),0.0,1,8)
      call raccord(qstobq(1,1),0.0,1,8)
      call raccord(scal(1,1,ks2,1),0.0,1,8)
      call raccord(scal(1,1,ks2-1,1),0.0,1,8)
      call raccord(scal(1,1,ks2-2,1),0.0,1,8)
      endif

      return
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!- fin de la routine staocb -
      end
