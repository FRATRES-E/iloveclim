!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:51 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:51 CET 2009

      SUBROUTINE rahmflx(ns,nn99,scaldt)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! Calcule (en nitrah ss.iteration) et integre le flux (thermique)
!  issu de la diffusion des anomalies de T. (cf. S.Rahmstorf).
! Rappel Explic. calcule et incorpore apres Diffus.Anom(<-avec ss.iter)
!--
!  modif : 25/05/99

!! START_OF_USE_SECTION

      use const_mod, only: zero, one, epsil

      use para0_mod, only: imax, jmax, kmax, nsmax
      use para_mod,  only:

      use bloc0_mod, only: is1, js1, scalr, tms, is1, is2, scal, ahrap
     & , ibera, iberam, iberp, iberpm, ims1, ims2, jbera, jberam, jberp
     & , jberpm, jcl1, jcl2, js2, ks2, spvr, unsdx, unsdy, isf1, rappes
     & , fss, dts, isf1, isf2

      use bloc_mod,  only: nitrap, nstart, numit, smxy, cmxy
      use global_constants_mod, only: dblp=>dp, ip

      use newunit_clio_mod, only: mouchard_id, clio3_out_id

!! END_OF_USE_SECTION

      implicit none

!! START_OF_INCLUDE_SECTION

c~ #include "type.com"
! [SCRPTCOM] #include "para.com"
! [SCRPTCOM] #include "bloc.com"

!! END_OF_INCLUDE_SECTION

!--dummy variables :
      real(dblp), dimension(imax,jmax,kmax) :: scaldt

!--common local (variables a conserver d'un appel a l'autre) :
      common / rhmflx /
     &  xfrz, tfrz, tfrz1, dtfrz, xvois, xxfrz,
     &  cphix(imax,jmax), cphiy(imax,jmax), unsmvs(imax,jmax)

!--variables locales :
      real(dblp), dimension(imax,jmax) :: tanom, tanom1, tanom0, tmfrz
      real(dblp), dimension(nsmax+4)   :: buffer
!-

! --- dmr local variables arising from implicit none construct

      integer(ip) :: i, j, k, n, nbline, nit, nnloc, ns, nn99

      real(dblp)  :: ccdif, cphix, cphiy, dtfrz, fflx, flxrah, sumtm
     &   , tfrz, tfrz0, tfrz1, ttm1ob, ttm2ob, ttmfrz, unsmvs, unsrap
     &   , xfrz, xvois, xxfrz

      integer(ip) :: run_param_id
      
      k = ks2
      if (numit.eq.nstart) then
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1 ) Mise en place parametres & coeff. <= only the 1rst time step    |
!-----------------------------------------------------------------------

!--Lecture des parametres (additionel) dans le fichier "run.param" :
        open(newunit=run_param_id, file='run.param', status='old')
        nbline = 10 + 2*nsmax
        do n=1,nbline
          read(run_param_id,*)
        enddo
        read(run_param_id,*,err=900) (buffer(n),n=1,nsmax+4), nnloc,
     &             xfrz, tfrz, tfrz0, dtfrz, xvois
        close(run_param_id)
        if (nnloc.ne.nitrap) goto 900
        tfrz1 = tfrz0 + dtfrz
!---
!   xfrz = part relative (maximale) de (T*) remplacee par "tfrz"
!   tfrz = temp. de remplacement (a la place de T*)
!   tfrz0,tfrz1 = intervalle de transition (=> decide de remplacer T*)
!   xvois = poids des voisins (centre = 1) dans la moyenne lissee
!   par ex. : xfrz,tfrz,tfrz0,tfrz1,xvois = 1., -2., -2.5, -1.5,  3.
!---

        if (nn99.eq.2) then
         write(mouchard_id,'(2A,I5,1PE14.6)') 'rahmflx : Vers=c(f) ',
     &          '; Rappes apres Ah_Rap ; nitrap,ahrap=', nitrap, ahrap
         write(mouchard_id,
     &    '(A,2F10.6,3F10.5)') ' xfrz,xvois ,tfrz,tfrz0 & 1 :',
     &           xfrz, xvois, tfrz, tfrz0, tfrz1
        endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--Definit les constantes une fois pour toutes :
      xxfrz = (xfrz / dtfrz) / (1.0 + xvois)

!- precaution : cyclicite (et autre) du tab. ``T.star'' :
      call raccord(scalr(1,1,k,ns), spvr, 1, 8)

!--Mise en place des Coeff. pour  Diffus.Anom.Temp. (<-S.Rahmstorf) :
!-  cphix&y <-deja.mult par deltaT/Dz [Rq. Dz etait deja retiree ds run.param]
      unsrap = one / DFLOAT(nitrap)
      ccdif = dts(k) * ahrap * unsdx * unsrap
      do j=js1,js2
       do i=is1(j),1+is2(j)
        ttm1ob = min(tms(i-1,j,k), (scalr(i-1,j,k,1)-spvr),
     &               tms(i,j,k),   (scalr(i,j,k,1)-spvr)   )
        cphix(i,j) = ccdif * smxy(i,j,1) * ttm1ob
       enddo
      enddo
      ccdif = dts(k) * ahrap * unsdy * unsrap
      do j=js1,1+js2
       do i=isf1(j),isf2(j)
        ttm2ob = min(tms(i,j-1,k), (scalr(i,j-1,k,1)-spvr),
     &               tms(i,j,k),   (scalr(i,j,k,1)-spvr)   )
        cphiy(i,j) = ccdif * cmxy(i,j,2) * ttm2ob
       enddo
      enddo

!- Pour le lissage : masque des 4 voisins :
      do j=js1,js2
       do i=is1(j),is2(j)
        sumtm = (tms(i-1,j,k)+tms(i+1,j,k)+tms(i,j-1,k)+tms(i,j+1,k))
        unsmvs(i,j) = (xvois / dtfrz) * tms(i,j,k) / max(epsil,sumtm)
       enddo
      enddo
!-----
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) Initialisation.                                                 |
!-----------------------------------------------------------------------

!- mise en place de "tanom" et du masque "tmfrz" (temp. below freezing point)
      do j=1,jmax
       do i=1,imax
        tanom(i,j) = scalr(i,j,k,ns) - scal(i,j,k,ns)
        tmfrz(i,j) = tms(i,j,k) *
     &           max(zero, min(dtfrz, tfrz1-scal(i,j,k,ns)) )
       enddo
      enddo

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--Modifie temperature "T*" where T < tfrz ; -> remplacee par (1-X).T* + X.tfrz
!  avec : if T<tfrz, X = xfrz (1 + xvois.S_4(if)/S_4tm)/(1+xvois) ; 0 si non

      do j=js1,js2
       do i=is1(j),is2(j)
        ttmfrz = xxfrz * tmfrz(i,j) * ( 1.0 + unsmvs(i,j) *
     &       (tmfrz(i-1,j)+tmfrz(i+1,j)+tmfrz(i,j-1)+tmfrz(i,j+1)) )
        tanom(i,j)  = tanom(i,j)
     &              + ttmfrz * max(zero, tfrz - scalr(i,j,k,ns))
        tanom0(i,j) = tanom(i,j)
       enddo
      enddo
!-----

!     if (nn99.eq.2 .and. numit.eq.nlast) write(99,'(A,I9,I5,3F12.6)')
!    &   ' numit, nit, tanom =', numit, 0, tanom(icheck,jcheck),
!    &       scalr(icheck,jcheck,k,ns), scal( icheck,jcheck,k,ns)

      do nit=1,nitrap
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  3 ) Calcule et somme (ds tanom) les flux par iteration sur "nit".   |
!-----------------------------------------------------------------------

!     call raccord(tanom, zero, 1, 8)
!- "inlining explicite" :
      do j=jcl1,jcl2
        tanom(ims1-1,j) = tanom(ims2,j)
        tanom(ims2+1,j) = tanom(ims1,j)
      enddo
        tanom(iberpm,jberp) = tanom(ibera, jberam)
        tanom(iberp, jberp) = tanom(iberam,jberam)
        tanom(iberam,jbera) = tanom(iberp, jberpm)
        tanom(ibera, jbera) = tanom(iberpm,jberpm)

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--Forcage thermique : Diffus.Anom.Temp. (<-S.Rahmstorf) :
!-  cphix & cphiy <-deja.mult par deltaT/Dz et / par nitrap
!----
!- calcul du flux et l'ajoute a "tanom" dans "tanom1":
!$DIR SPP LOOP_PARALLEL
!$DIR SPP LOOP_PRIVATE(i,flxrah)
      do j=js1,js2
       do i=is1(j),is2(j)
        flxrah = unsdx*(cphix(i,j)*( tanom(i-1,j) - tanom(i,j) )
     &                 -cphix(i+1,j)*( tanom(i,j) - tanom(i+1,j) ))
     &         + unsdy*(cphiy(i,j)*( tanom(i,j-1) - tanom(i,j) )
     &                 -cphiy(i,j+1)*( tanom(i,j) - tanom(i,j+1) ))
        tanom1(i,j) = tanom(i,j) + smxy(i,j,0) * flxrah
       enddo
      enddo

!- Mise a jour de "tanom" :
      do j=js1,js2
       do i=is1(j),is2(j)
        tanom(i,j) = tanom1(i,j)
       enddo
      enddo

!     if (nn99.eq.2 .and. numit.eq.nlast) write(99,'(A,I9,I5,F12.6)')
!    &   ' numit, nit, tanom =', numit, nit, tanom(icheck,jcheck)
!----
      enddo

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  4 ) Mise a jour des tableaux fss & scaldt.                          |
!-----------------------------------------------------------------------

!  Ajoute la nouvelle evolution de T.surf calculee par nitrap iterations :
      do j=js1,js2
       do i=is1(j),is2(j)
        fflx = (tanom(i,j) - tanom0(i,j)) - rappes(i,j,ns) *
     &        ( tanom(i,j) + (scal(i,j,k,ns) - scalr(i,j,k,ns)) )
        fss(i,j,ns) = fss(i,j,ns) + fflx
        scaldt(i,j,k) = scaldt(i,j,k) - fflx
       enddo
      enddo

      return

 900  continue                      ! Error handling for read errors on file
                                    ! 'run.param' or when nnloc (read from file)
                                    ! /= nitrap
!- cas d'erreur :
      write(clio3_out_id,'(2A,I4)') 'ARRET, rahmflx :',
     &     ' Error in reading "run.param", line', nbline+1
      stop

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine rahmflx -
      end
