!dmr -- Added optional components choice - Tue Jun 14 15:59:05 CEST 2011
#include "choixcomposantes.h"
!dmr -- Added optional components choice - Tue Jun 14 15:59:05 CEST 2011

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Ce module est un portage en FORTRAN 90 du module CLIO initial
!       cree en FORTRAN 77.
!      Il contient une routine suplementaire pour simuler les
!      equivalences du code FORTRAN 77 via des POINTER en FORTRAN 90
!
!      (dans l'environnement logiciel LUDUS)
!
!      Auteur : ??, Didier M. Roche
!      Date   : ??, 20 Aout 2014
!      Derniere modification : 19 juin 2018
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module vareq_mod

      use para0_mod, only: imax, jmax, kmax, nsmax
      use para_mod,  only: nvmax, krlmax, krlmin

      implicit none

      private:: imax, jmax, kmax, nvmax, nsmax, krlmax, krlmin

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! originally written with IMPLICIT INTEGER (I-N), REAL (A-H, O-Z)
! transferred as such in fortran90
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  fichier "vareq.com" : incorpore par instruction 'include' dans les routines :
!     CLASS, GRADP, STATIS, ncdfout, moyen, local, binout, foroutp,
!      lisstab, (unigl, option)
!  (commun a toutes les routines de sortie (output) de resultats)
! - inclus apres "type.com", "para.com", "bloc.com", complete "varno.com".
!------------------
! Les variables sont rangees par niveau, les uns a la suite des autres,
!   dans un seul tableau 3D : vrl
!   et chaque couche d un tableau 3D est reperee par un indice general (k).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  modif : 30/06/98

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr Bloc "rtowrite"
!-----
!--common regroupant les variables necessaires pour l'ecriture sur fichier :
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=8), dimension(0:nvmax) :: spv, cmult, cadit, cliss

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr Bloc "itowrite"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      integer                     :: nfrc
      integer, dimension(0:3)     :: irl1, irl2, jrl1, jrl2
      integer, dimension(imax,8)  :: irn, jrn
      integer, dimension(0:nvmax) :: nabs

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr Bloc "ctowrite"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      character(len=10), dimension(-2:nvmax) :: fmt1
      character(len=13)                      :: titex1, titex2
      character(len=14)                      :: titex3
      character(len=40)                      :: titexp

!~ !-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!~ ! dmr Definitions du bloc equivalents pre-existant
!~ !-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        real(kind=8), dimension(imax,jmax)               :: psh, hmajc, egajc
        real(kind=8), dimension(imax,jmax,kmax)          :: vaflux,vdflux, alphxo, alphyo, haterm, hdterm
        real(kind=8), dimension(imax,jmax,-1:nsmax)      :: uslpfx,vslpfx
        real(kind=8), dimension(imax,jmax,krlmin:krlmax) :: vrl

      contains

!~ !-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!~ !       Procedure vareq_init is designed to associate the pointers defined above with their target.
!~ !-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine copy_to_vrl

       use para_mod,  only: krlfw, krlps, krlet, krlub, krlvb, krls, krlt, krlu, krlv, krlb, krln2, krlas, krlau, krltke, krlw     &
                   , krlajc, krlvaf, krlvdf, krlaxt, krlayt,krlhat, krlhdt, krlhac, krlajc, krlusl, krlvsl
       use bloc0_mod, only: fss, scal, eta, ub, vb, u, v, b, bvf, avsdz, avudz, w, q2turb, fqajc

       integer n, k

       vrl(:,:,krlfw) = fss(:,:,0)
       vrl(:,:,krlps) = psh(:,:)
       vrl(:,:,krlet) = eta(:,:)
       vrl(:,:,krlub) = ub(:,:)
       vrl(:,:,krlvb) = vb(:,:)
       vrl(:,:,krlu ) = u(:,:,1)
       vrl(:,:,krlv ) = v(:,:,1)
       vrl(:,:,krlt ) = scal(:,:,1,1)
       vrl(:,:,krls ) = scal(:,:,1,2)
       vrl(:,:,krlb ) = b(:,:,1)
       vrl(:,:,krln2) = bvf(:,:,1)
       vrl(:,:,krlas) = avsdz(:,:,1)
       vrl(:,:,krlau) = avudz(:,:,1)
       vrl(:,:,krlw ) = w(:,:,1)
       vrl(:,:,krltke)= q2turb(:,:,1)
       vrl(:,:,krlajc)= fqajc(:,:,1)
       vrl(:,:,krlvaf)= vaflux(:,:,1)
       vrl(:,:,krlvdf)= vdflux(:,:,1)
       vrl(:,:,krlaxt)= alphxo(:,:,1)
       vrl(:,:,krlayt)= alphyo(:,:,1)
       vrl(:,:,krlhat)= haterm(:,:,1)
       vrl(:,:,krlhdt)= hdterm(:,:,1)
       vrl(:,:,krlhac)= hmajc(:,:)
       vrl(:,:,krlajc)= egajc(:,:)
       vrl(:,:,krlusl)= uslpfx(:,:,-1)
       vrl(:,:,krlvsl)= vslpfx(:,:,-1)

      end subroutine copy_to_vrl

!~ !-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!~ !
!~ !-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      end module vareq_mod
!--fin du fichier "vareq.com"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr --- The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
