!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:53 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:53 CET 2009

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Ce module est un portage en FORTRAN 90 du module slopez initial
!       cree en FORTRAN 77.
!      (dans l'environnement logiciel LUDUS)
!      Auteur : ??, Didier M. Roche
!      Date   : ??
!      Derniere modification : 23 juin 2020
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module slopez_mod

      implicit none

      contains
      
      SUBROUTINE slopez(cslpmx,kslpdw,nn99)
      
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  appelee par "scale" ;
! Detecte les cas de "Down-Sloping" ; Consigne Volume & Niveaux Concernes.
!  modif : 22/08/97
!dmr --- : 2018-11-14, 2020-06-23
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      use global_constants_mod, only: dblp=>dp, ip

      use const_mod, only: epsil, zero
      use para0_mod, only: imax, kmax, nlpmax, ixjmax
      use bloc0_mod, only: kfs, u, v, b, xslop, dz, unsdz
      use bloc_mod,  only: cmx, cmy, nxslp, ijslp, kslp, lslp, nslp, alpslp, nxyslp, nslpdw, uvcslp, numit, nxyslp, nstart, refexp

      use newunit_clio_mod, only: clio3_out_id
      
      implicit none

!--by reference variables :
      integer(kind=ip)       , intent(in)              :: nn99
      real(kind=dblp), dimension(kmax)  , intent(in)   :: cslpmx
      integer(kind=ip), dimension(nlpmax), intent(out) :: kslpdw


!--variables locales :
      real(kind=dblp), dimension(ixjmax)      :: kfshk, cmx1hk, cmx2hk, cmy1hk, cmy2hk
      real(kind=dblp), dimension(ixjmax,kmax) :: uhk, vhk, bhk


!--variables locales :

      integer(kind=ip) :: k, kk, nl, nldw, l, ij, iij, ijj, nnx
      real(kind=dblp)  :: u1, v2, uuslp, vvslp

!--variables locales pour le controle de premiere iteration (fin de routine)

#if ( CTRL_FIRST_ITER == 1 )
      integer, dimension(0:kmax,4) :: nnkslp
      real, dimension(0:kmax)      :: uuavr, vvavr
      integer                      :: nny, kkx, kky, jj0
      real                         :: cc00
#endif      
      
      integer(kind=ip) :: xyslope_id
!     if (numit.eq.nstart) then
!       write(clio3_out_id,*) 'slopez : Uslp = xslop * dz * Db'
!     endif
      
!dmr [NOEQUI]
      kfshk(:) = RESHAPE(kfs,(/ ixjmax /))
      cmx1hk(:)= RESHAPE(cmx(:,:,1),(/ ixjmax /))
      cmx2hk(:)= RESHAPE(cmx(:,:,2),(/ ixjmax /))
      cmy1hk(:)= RESHAPE(cmy(:,:,1),(/ ixjmax /))
      cmy2hk(:)= RESHAPE(cmy(:,:,2),(/ ixjmax /))

      do k=1,kmax
        uhk(:,k) = RESHAPE(u(:,:,k),(/ ixjmax /))
        vhk(:,k) = RESHAPE(v(:,:,k),(/ ixjmax /))
        bhk(:,k) = RESHAPE(b(:,:,k),(/ ixjmax /))
      enddo
!dmr [NOEQUI]     
      
      
      
#if ( XSLOP == 0 )
! [PREPROCESS]      if (xslop.lt.epsil) then ! dmr [NOTA] this is never the case so far for iLOVECLIM where xslop == 1
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  2 ) Detecte les cas de DownSloping pour Decentrer (upwing) .        |
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      nldw = 0
!--Direction X , Dowsloping => Decentrement :
      do nl=1,nxslp
        ij = ijslp(nl)
        k = kslp(nl)
        l = lslp(nl)
        iij = ij + max(0,l)
        u1 = DFLOAT(-l) * (uhk(iij,k) + uhk(iij+imax,k))
        if (u1.gt.zero .and. bhk(ij+l,k).gt.bhk(ij,k) ) then
          nldw = nldw + 1
          nslp(nldw) = nl
          kslpdw(nldw) = k
          alpslp(nldw) = 0.25 * cmy1hk(iij) * u1
        endif
      enddo ! sur nxslp
      nnx = nldw

!--Direction Y , Dowsloping => Decentrement :
      do nl=nxslp+1,nxyslp
        ij = ijslp(nl)
        k = kslp(nl)
        l = lslp(nl)
        ijj = ij + max(0,l)
        v2 = vhk(ijj,k) + vhk(ijj+1,k)
        if (v2*DFLOAT(l).lt.zero .and. bhk(ij+l,k).gt.bhk(ij,k) ) then
          nldw = nldw + 1
          nslp(nldw) = nl
          kslpdw(nldw) = k
          alpslp(nldw) = 0.25 * cmx2hk(ijj) * abs(v2)
        endif
      enddo ! sur nxyslp
      nslpdw = nldw

! [PREPROCESS]      else

#else /*  XSLOP > 0 */
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  3 ) Detecte cas de DownSloping pour Decentrer et Permuter (xslop)   |
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      nldw = 0
!--Direction X , Dowsloping => permutation (xslop) :
      do nl=1,nxslp
        ij = ijslp(nl)
        k = kslp(nl)
        l = lslp(nl)
        iij = ij + max(0,l)
        if ( bhk(ij+l,k).gt.bhk(ij,k) ) then
          nldw = nldw + 1
          nslp(nldw) = nl
!--Recherche du 1er Niveau kk / rho_pot(i,j,kk) < rho_pot(i+l,j,k)
          do kk=k-1,NINT(kfshk(ij)),-1
            if (bhk(ij+l,kk).le.bhk(ij,kk) ) then
              kslpdw(nldw) = kk + 1
              goto 335              ! exit loop, and skip first command after
            endif
          enddo
          kslpdw(nldw) = kfshk(ij)
 335      continue                  ! Landing point from the exit out of the kk-loop above 
!- calcul de la vitesse :
          uuslp = xslop * dz(k) * (bhk(ij+l,k) - bhk(ij,k))
          uvcslp(nldw) = cmy1hk(iij) * min( uuslp, cmx1hk(iij) * cslpmx(k) )
        endif
      enddo
      nnx = nldw

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!--Direction Y , Dowsloping => permutation (xslop) :
      do nl=nxslp+1,nxyslp
        ij = ijslp(nl)
        k = kslp(nl)
        l = lslp(nl)
        ijj = ij + max(0,l)
        if ( bhk(ij+l,k).gt.bhk(ij,k) ) then
          nldw = nldw + 1
          nslp(nldw) = nl
!--Recheche du 1er Niveau kk / rho_pot(i,j,kk) < rho_pot(i,j+l,k)
          do kk=k-1,NINT(kfshk(ij)),-1
            if (bhk(ij+l,kk).le.bhk(ij,kk) ) then
              kslpdw(nldw) = kk + 1
              goto 355              ! exit loop, and skip first command after
            endif
          enddo
          kslpdw(nldw) = kfshk(ij)
 355      continue                  ! Landing point from the exit out of the kk-loop above 
!- calcul de la vitesse :
          vvslp = xslop * dz(k) * (bhk(ij+l,k) - bhk(ij,k))
          uvcslp(nldw) = cmx2hk(iij) * min( vvslp, cmy2hk(iij) * cslpmx(k) )
        endif
      enddo
      nslpdw = nldw

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! [PREPROCESS]      endif ! on xslop ...

#endif
            
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  4 ) Verification a la 1ere iteration :                             |
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!dmr [NOTA] the following lines have no impact on the calculation, but cost
!dmr        an if statement per call.
!dmr        -> Deactivated by default


#if ( CTRL_FIRST_ITER == 1 )

      if_numit_eq_nstart: if (numit.eq.nstart .and. nxyslp.ge.1) then
        nny = nslpdw - nnx
        write(clio3_out_id,'(2A,F8.3,4I5)') &
          ' slopez(U=xslop*dz*Db) : xslop,',' nXslp,nYslp,Nx,Ny=', &
          xslop, nxslp,nxyslp-nxslp, nnx,nny
       
        if_nn99_eq_2: if (nn99.eq.2) then
!--Evaluation du Nb de "boites" impliquees : (potentiel & effectif)
        do k=0,kmax
         uuavr(k) = 0.
         vvavr(k) = 0.
         do l=1,4
          nnkslp(k,l) = 0
         enddo
        enddo

        do nl=1,nxslp
          k = kslp(nl)
          nnkslp(k,1) = nnkslp(k,1) + 1
        enddo
        
        do nl=nxslp+1,nxyslp
          k = kslp(nl)
          nnkslp(k,3) = nnkslp(k,3) + 1
        enddo
!-
        kkx = nnx
        do nldw=1,nnx
          nl = nslp(nldw)
          k = kslp(nl)
          nnkslp(k,2) = nnkslp(k,2) + 1
          kkx = kkx + k - kslpdw(nldw)
          uuavr(k) = uuavr(k) + abs(uvcslp(nldw))
        enddo
        
        kky = nny
        do nldw=nnx+1,nslpdw
          nl = nslp(nldw)
          k = kslp(nl)
          nnkslp(k,4) = nnkslp(k,4) + 1
          kky = kky + k - kslpdw(nldw)
          vvavr(k) = vvavr(k) + abs(uvcslp(nldw))
        enddo
!       write(clio3_out_id,'(2A,F10.4,4I6)') ' slopez(U=xslop*dz*Db) : ',
!    &      'xslop, Nx,Kx,Ny,Ky=', xslop, nnx, kkx, nny, kky

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!--Ecriture sur fichier "xyslope" :
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
        open(newunit=xyslope_id, file='xyslope', status='unknown' )
        write(xyslope_id,'(A,I11,A,F10.4)') &
          'Exp '//refexp//', It', numit,' , slopez : xslop=', xslop
        write(xyslope_id,'(A,6I6)') &
          ' nXslp,Nx,Kx, nYslp,Ny,Ky=',nxslp, nnx, kkx, nxyslp-nxslp, nny, kky

        write(xyslope_id,'(2A)') &
          '  k   Pot.X,Dwn.X   Pot.Y,Dwn.Y ', '    U.moy  (m/s)  V.moy'
        do k=kmax,1,-1
          uuavr(k) = uuavr(k) * unsdz(k)
          vvavr(k) = vvavr(k) * unsdz(k)
          uuavr(0) = uuavr(0) + uuavr(k)
          vvavr(0) = vvavr(0) + vvavr(k)
          if (nnkslp(k,2).ge.1) uuavr(k) = uuavr(k)/DFLOAT(nnkslp(k,2))
          if (nnkslp(k,4).ge.1) vvavr(k) = vvavr(k)/DFLOAT(nnkslp(k,4))
          write(xyslope_id,'(I3,2(2X,2I6),2X,1P2E11.3)') &
            k, (nnkslp(k,l),l=1,4), uuavr(k), vvavr(k)
         do l=1,4
          nnkslp(0,l) = nnkslp(0,l) + nnkslp(k,l)
          enddo
         enddo



          if (nnkslp(0,2).ge.1) uuavr(0) = uuavr(0)/DFLOAT(nnkslp(0,2))
          if (nnkslp(0,4).ge.1) vvavr(0) = vvavr(0)/DFLOAT(nnkslp(0,4))
          write(xyslope_id,'(A3,2(2X,2I6),2X,1P2E11.3)')'T: ', &
            (nnkslp(0,l),l=1,4), uuavr(0), vvavr(0)
        write(xyslope_id,*)
        write(xyslope_id,'(A)') ' Max V.slp (cslpmx/dz) :'
        write(xyslope_id,'(1P8E10.3)') (cslpmx(k)*unsdz(k),k=kmax,1,-1)
        write(xyslope_id,*)
!--
        write(xyslope_id,'(A)') &
          '  i   j   k   l    b(l)-b     uv.slp  '//'      u       xslop.eq'
        do nldw=1,nnx
          nl = nslp(nldw)
          ij = ijslp(nl)
          k = kslp(nl)
          l = lslp(nl)
          iij = ij + max(0,l)
          u1 = DFLOAT(-l) * 0.5 * (uhk(iij,k) + uhk(iij+imax,k))
          jj0 = (ij - 1) / imax
          if ( u1.gt.epsil ) then
            cc00 = max(epsil, bhk(ij+l,k)-bhk(ij,k) )
            cc00 = u1 / cc00
            write(xyslope_id,'(4I4,1P4E11.3)') &
              ij - imax*jj0,1+jj0, k, l, bhk(ij+l,k)-bhk(ij,k), &
              uvcslp(nldw)*unsdz(k), u1, cc00
          else
            write(xyslope_id,'(4I4,1P4E11.3)') &
              ij - imax*jj0,1+jj0, k, l, bhk(ij+l,k)-bhk(ij,k), &
              uvcslp(nldw)*unsdz(k)
          endif
        enddo
        
        write(xyslope_id,*)
!--
        write(xyslope_id,'(A)') &
          '  i   j   k   l    b(l)-b     uv.slp  '//'      v       xslop.eq'
        do nldw=nnx+1,nslpdw
          nl = nslp(nldw)
          ij = ijslp(nl)
          k = kslp(nl)
          l = lslp(nl)
          ijj = ij + max(0,l)
          v2 = 0.5 * (vhk(ijj,k) + vhk(ijj+1,k))
          v2 = sign(v2, DFLOAT(-l))
          jj0 = (ij - 1) / imax
          if ( v2.gt.epsil ) then
            cc00 = max(epsil, bhk(ij+l,k)-bhk(ij,k) )
            cc00 = v2*dz(k) / cc00
            write(xyslope_id,'(4I4,1P4E11.3)') &
              ij - imax*jj0, 1+jj0, k, l, bhk(ij+l,k)-bhk(ij,k), &
              uvcslp(nldw)*unsdz(k), v2, cc00
          else
            write(xyslope_id,'(4I4,1P4E11.3)') &
              ij - imax*jj0, 1+jj0, k, l, bhk(ij+l,k)-bhk(ij,k), &
              uvcslp(nldw)*unsdz(k)
          endif
        enddo
        write(xyslope_id,*)
!--fin d'ecriture.
        close(xyslope_id)
        endif if_nn99_eq_2
      endif if_numit_eq_nstart

#endif            
      return
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!- fin de la routine slopez -
      end subroutine slopez
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

     end module slopez_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
