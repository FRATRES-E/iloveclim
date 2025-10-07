!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:53 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:53 CET 2009

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      module scale_mod
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      
      use const_mod   , only: one, yeaday
      use global_constants_mod, only: dblp=>dp, ip
      use para0_mod   , only: kmax, nsmax, nlpmax


      implicit none


      integer(kind=ip), parameter        :: k1000 = 10, j50S = 13, k2000 = 7
      real(kind=dblp)                    :: unsyr, txeadv, txedif, ccwe
      real(kind=dblp), dimension(kmax)   :: ccwabs, ccwdmy, difzmx, cczdt, cczdte, cdtsxz, cslpmx
      real(kind=dblp), dimension(nsmax)  :: ccflxt
      integer(kind=ip), dimension(nlpmax):: kslpdw
      contains
      
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|      
      SUBROUTINE scale_init
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      
      use bloc0_mod   , only: txiads, txidfs, alphaz, unsdzw, dts, dz, unsdz, ks1, ks2, unsdx, dx
!~       use isodiffu_mod, only: isodiffu_init
      use const_mod,    only: rho0, cpo
      
      integer(kind=ip) :: k
      unsyr = one / (yeaday * 86400.)
      
      !--preparation des coeffs dependant de la verticale :

      txeadv = 1.0 - txiads
      txedif = 1.0 - txidfs
      ccwe = 0.5 * txeadv

      do k=ks1+1,ks2
        ccwabs(k) = ccwe * alphaz(k)
        ccwdmy(k) = ccwe * unsdzw(k) * (dts(k-1) + dts(k)) * 0.5
        
      enddo
      
!--calcul de la diffusivite/dz explicite maximale (restant stable) :
      do k=ks1,ks2
        difzmx(k) = 0.5 * dz(k) / dts(k)
        cczdt(k)  = unsdz(k) * dts(k)
        cczdte(k) = cczdt(k) * txeadv
        cdtsxz(k) = unsdx * cczdt(k)
        cslpmx(k) = 0.125 * dx / cczdt(k)
      enddo

      ccflxt(1) = dts(ks2)*unsdz(ks2) / (rho0*cpo)
      ccflxt(2:nsmax) = 0.0_dblp
      
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      END SUBROUTINE scale_init
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      SUBROUTINE scale(nsew, nn99)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!**AVANCEMENT EXPLICITE D'UN PAS DE TEMPS SCALAIRE DES SCALAIRES**
! traitement explicite a "pas de temps fractionnaires", pour les 3 D.
!  et avec taux d'implicitete (a choisir) selon la verticale.
! ordre : 1er=Vert.expli. Si nsew=1 : 2eme=Est-West, 3eme=Nord-Sud, #0 inverse
!--
! Tableaux : scal <- Valeur conserve tel quel jusque fin de "scale" .
!          scalat <- Valeur actualise a chaque "fragment" de  pas de temps.
!          scaldt <- Variation (+correction) totalement explicite.
!--
! En entree : avsdz = Coeff.Diffus.Vert. entier ; en sortie = implicite only
!  advection suivant la verticale :  w(k) * [ S(k-1) + S(k) ] / 2
!  modif : 06/10/99
!dmr --- : 2018-11-20, 2020-06-16
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      use const_mod   , only: zero

      use para0_mod   , only: nsmax, ixjmax, kmax, imax, jmax!, nlpmax
      use para_mod    , only: 
      use bloc0_mod   , only: ks1, ks2, dts, scpme, spvr, scal, is1, is2  &
                      , phiss, phivs, rappes, fss, kfs, phifs, scalr, js1, js2, tms, alphgr, ai, rapint
                    ! , alphaz, avsdz, dx, dz, txiads, txidfs, unsdx, unsdz, unsdzw, w
      use bloc_mod    , only: numit, nstart, deriv, nrap, ijrap
      use ice_mod     , only: fcm1

      use isodiffu_mod, only: isodiffu, isodiffu_init
      
      use slopez_mod, only: slopez

#if ( GEOTHERMAL == 1 )
      use ncio        , only: nc_read
      use const_mod,    only: rho0, cpo
#endif

      use newunit_clio_mod, only: clio3_out_id

!$    use omp_lib

      implicit none

!--dmr By reference variables
      integer, intent(in) :: nsew, nn99

!--variables locales :

      real, dimension(ixjmax,kmax,2)   :: alphhk
      real, dimension(imax,jmax,kmax)  :: scalat, scaldt, scal_copie
      real, dimension(imax,kmax)       :: difzex
      real, dimension(imax,kmax+1)     :: phiz

#if ( GEOTHERMAL == 1 )
      real, dimension(imax,jmax), save :: geoheat
#endif

!~  dmr Moved module wide
!~       real, dimension(kmax)            :: ccwabs, cczdt, cczdte, ccwdmy, difzmx, cdtsxz, cslpmx
!-
!~       integer, dimension(nlpmax), save :: kslpdw


!dmr ---
      integer                          :: kfd

      real                             :: sscal, fflx
!      real                             :: phdivz
      integer                          :: i, j, k, n, ns
!$      REAL :: t_ref, t_precision, t_inter, t_final

!dmr ---

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  1 ) preparation des coeff. intervenant dans les flux d'advection .  |
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!~ !dmr ---
!~       k1000 = 10 ! 1 is all depths
!~ !-dmr      k1000 = 8 ! 8 is 1000 meters
!~       k2000 = 7 ! 5 is 2300 meters, 7 is 1225
!~       j50S = 13 ! 10 is 50°S
!~ !dmr ---

! convention :
!    CCW* coeff relatif a l'interface (= place des W).
!    CCZ* coeff relatif au centre des boites (= place de T,S)
!    CC*E : tient compte du taux Explicite (= TXE*) .
!    CCWABS(k) intervient avec |w| ;
!              CCWUP(k) avec scal(k) ; CCWDowN(k) avec scal(k-1).
!  pour le taux de Decentrement calcule a partir du Nb de courant : CCWDUP(k)
!    avec dts(k), CCWDDowN(k) avec dts(k-1), CCWDMY(k) avec la moyenne des 2.
!    CCZDT**(k) intervient pour le bilan des flux V, boite "k" ;

!~       unsyr = one / (yeaday * 86400.)

#if ( GEOTHERMAL == 1 )
      if (numit.eq.nstart) then
!      !! Read geothermal heat flux (in mW / m^2) in netcdf file
      call nc_read("inputdata/clio/geothermal_heating.nc","heat_flux",geoheat)
!      !! Convert heat flux into appropriate units
      geoheat(:,:) = geoheat(:,:) / ( 1000 * rho0 * cpo )
      endif
!      !! Add heat flux to phifs(:,:,1)
      phifs(:,:,1) = phifs(:,:,1) + geoheat(:,:)
#endif

!--initialisation (indispensable) :

      do k=1,kmax
       do i=1,imax
        difzex(i,k) = 0.0
       enddo
      enddo
      
      do k=1,kmax+1
       do i=1,imax
        phiz(i,k) = 0.0
       enddo
      enddo

!dmr --- could be moved up ? Called once per call of scale potentially ...
      if (numit.eq.nstart) then
       write(clio3_out_id,*) 'scale : Adv. Zex,(X,Y) Alterne, + Corrigee.'
       call isodiffu_init()
      endif
!dmr --- could be moved up ? Called once per call of scale potentially ...

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  2 ) Detecte les courant de Downsloping .
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!~       call slopez(nn99)
      call slopez(cslpmx, kslpdw, nn99)      

!--Debut de la boucle externe sur tous les scalaires, indices "ns" :
!- N.B : Parallelisable SAUF si txiads ou txidfs different de 1.
!- NOTE A BENETS : txiads == 1.0 et txidfs == 1.0 dans iLOVECLIM (dmr)

!dmr --- [test OMP]
!~ !$      t_ref = OMP_GET_WTIME()
!~ !$      t_precision = OMP_GET_WTICK()


!!!!$OMP DO PRIVATE(ns,i,j,k,kfd,scaldt,scalat,sscal,ccflxt,fflx,alphhk)
!$OMP PARALLEL PRIVATE(ns,i,j,k,kfd,scaldt,scalat,scal_copie,sscal,fflx,alphhk)
!$OMP DO SCHEDULE(STATIC)
      do ns=1,nsmax

!!- !$      t_inter = OMP_GET_WTIME()
!-----

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  3 ) Conditions Limites : Surface & Fond - Modification de "avsdz" . |
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!- Debut de la 1ere boucle externe sur la latitude, indice "j" :

      do j=1,jmax

!--Initialisation de scalat (<- scal) et scaldt (par le terme de "derive") :
       do k=1,kmax
!dmr orig sscal = -deriv(ns) * dts(k)
!dmr ---
        if (scpme(ns).eq.spvr) then
          sscal = -deriv(ns) * dts(k)
        else
           IF ((j.LE.j50S).AND.(k.GE.k2000)) THEN
                          ! only south of 50°S in all basins           
             sscal = -deriv(ns) * dts(k)
           ELSEIF (k.GE.k1000) THEN
             sscal = -deriv(ns) * dts(k)
           ELSE
             sscal = 0.0d0
           ENDIF
        endif
!dmr ---
        do i=1,imax
         scalat(i,j,k) = scal(i,j,k,ns)
         scaldt(i,j,k) = sscal - phivs(i,j,k,ns)
        enddo ! on i
        enddo ! on k

!--Calcul du Forcage(explicite) en Surface : (phiss <-deja.mult par deltaT/Dz)
      k = ks2
! fflx flux = impact
! rappes = restoring (should be zero) for experiments to force it with a known scalar field
! e.g. mediterrean sea example.
       do i=is1(j),is2(j)
         fflx = phiss(i,j,ns) + rappes(i,j,ns) * ( scal(i,j,k,ns) - scalr(i,j,k,ns) )

         fss(i,j,ns) = fss(i,j,ns) + fflx - ccflxt(ns) * fcm1(i,j)
         scaldt(i,j,k) = scaldt(i,j,k) - fflx
       enddo

!--Expansion par le fond : (<- Not yet multipl. by deltaT)

       do i=is1(j),is2(j)
         kfd = kfs(i,j)
         scaldt(i,j,kfd) = scaldt(i,j,kfd) + cczdt(kfd) * phifs(i,j,ns)
       enddo

!--fin de la 1ere boucle sur l'indice de latitude "j".
      enddo ! on j

      scal_copie(:,:,:) = scalat(:,:,:)


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr --- NOT THE CASE IN iLOVECLIM
#if ( NIT_RAP > 0 )
!--Forcage thermique : Diffus.Anom.Temp. (<-S.Rahmstorf) :
      if (ns.eq.1) call rahmflx(ns,nn99,scaldt)
#endif
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! dmr [NOTA] below this is only if txiads != 1.0, never the case for us ...
! dmr ---> COMMENTED


!~ !-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!~ !  4 ) Direction Verticale : Traitement Explicite .                    |
!~ !-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!~       if ( (txiads.ne.1.0).or.(txidfs.ne.1.0) ) then
!~ !-----

!~       do j=js1,js2
!~ !--debut de la 2nd boucle sur l'indice de latitude "j".

!~       if (txidfs.ne.1.0) then
      
!~ !--4.1 calcul du coefficient de diffusion verticale explicite "difzex" :
!~       do k=ks1+1,ks2
!~        do i=is1(j),is2(j)
!~          difzex(i,k) = txedif * min( avsdz(i,j,k) , difzmx(k) )
!~        enddo
!~       enddo
      
!~       if (ns.eq.nsmax) then

!~ ! --- [NOTA] PROBLEME si PARALLELISATION SUR NS ???      
!~ !- transformation du tableaux de Diffusivite Verticale "avsdz" .
!~         do k=ks1+1,ks2
!~          do i=is1(j),is2(j)
!~            avsdz(i,j,k) = avsdz(i,j,k) - difzex(i,k)
!~          enddo
!~         enddo
!~       endif
!~ !--
!~       endif

!~ !-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!~ !--4.3 Flux Verticaux advect. & diffus. - Integre les Flux Explicites. |

!~ !--calcul des flux verticaux explicite advection et diffusion.
!~       do k=ks1+1,ks2
!~        do i=is1(j),is2(j)
!~          phiz(i,k) = tms(i,j,k-1) * (  ( difzex(i,k) + ccwabs(k) * abs(w(i,j,k)) )  * ( scalat(i,j,k-1) - scalat(i,j,k) )     &
!~                    + ccwe * w(i,j,k) * (scalat(i,j,k-1) + scalat(i,j,k))    )
!~        enddo
!~       enddo
!~ !--fin du calcul des flux verticaux.

!~ !--Bilan des flux Verticaux explicites :
!~       do k=ks1,ks2
!~        do i=is1(j),is2(j)
!~         phdivz = cczdte(k) * (w(i,j,k) - w(i,j,k+1)) * scal(i,j,k,ns)
!~         scalat(i,j,k) = scalat(i,j,k) - phdivz + cczdt(k) * (phiz(i,k) - phiz(i,k+1))
!~         scaldt(i,j,k)  = scaldt(i,j,k)  + phdivz
!~        enddo
!~       enddo

!~ !--fin de la 2nd boucle sur l'indice de latitude "j".
!~       enddo ! on j
      
!~ !--raccord cyclique et autre (scalat) :
!~ !     call raccord(scalat(1,1,1),zero,kmax,0)

!~ !--fin distinction tout implicite / en partie explicite .
!~ !----------------------------
!~       endif


! dmr <---- COMMENTED


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  5 ) Calcul des Coef. pour le decentrement Horizontal .              |
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


      if (alphgr(ns).lt.-1.0) then
        call alph2dc(alphhk(1,1,1),alphhk(1,1,2),ns)
      else
        call alphdec(alphhk(1,1,1),alphhk(1,1,2),ns)
      endif

!- Dowsloping : Modif par Permutation :

      call slopes(scalat,alphhk,cdtsxz,kslpdw,ns)

!- Isopycn.Diffus.: Compute (fr scal) & Incorporate (to scalat) Diagon.Flx.
!~       if (ai(kmax).gt.zero) call isodiffu(scalat,scal_copie,ns)
      if (ai(kmax).gt.zero) call isodiffu(scalat,scal_copie)
!--raccord cyclique et autre (scalat) :
      call raccord(scalat(1,1,1),zero,kmax,0)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  6 ) Flux Horizontaux advect. & diffus - Integre les Flux Explicites |
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!--bifurcation sur l'ordre de traitement des 2 directions N-S , E-W :
      if(nsew.ne.1) then

!--Resolution des Flux advect. & diffus., E-W en 1er, N-S en 2nd :
        call scadew(scalat,scaldt,alphhk(1,1,1),alphhk(1,1,2),ns)

      else

!--Resolution des Flux advect. & diffus., N-S en 1er, E-W en 2nd :
        call scadns(scalat,scaldt,alphhk(1,1,1),alphhk(1,1,2),ns)

      endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  7 ) Incorpore la partie totalement explicite (stokee dans scaldt).  |
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      do k=ks1,ks2

!--Ajoute le terme de rappel explicite interne (= ailleurs qu'en surf.) :
       do n=1,nrap(k,ns)
         i = 1 + mod(ijrap(n,k,ns),imax)
         j = 1 + ijrap(n,k,ns)/imax
         scaldt(i,j,k) = scaldt(i,j,k) + rapint(n,k,ns) * ( scalr(i,j,k,ns) - scal(i,j,k,ns) )
       enddo
!--
       do j=js1,js2
        do i=is1(j),is2(j)
         scal(i,j,k,ns) = scalat(i,j,k) + scaldt(i,j,k) * tms(i,j,k)
        enddo
       enddo
!----
      enddo ! on k

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!--Fin de la boucle externe sur tous les scalaires, indicies "ns" .

!~ !$      t_final = OMP_GET_WTIME()
!~ !$      WRITE(*,*) "Times === ", t_final-t_inter
      enddo ! on ns
!$OMP ENDDO      
!~ !$      WRITE(*,*) "Times === ", t_final-t_ref, t_precision


!$OMP END PARALLEL

      return
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!- fin de la routine scale -
      end subroutine scale
      
      end module scale_mod
