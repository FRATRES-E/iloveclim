!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Ce sous-programme provient du couplage LOVECLIM (AG)ISM
!       de P. Huybrecht. * L'import est pirate ! *
!      Modifications pour utilisation "offline" a partir des variables
!       de LOVECLIM vers GRISLI Nord
!
!      Auteur : P. Huybrecht
!      Date   : Inconnue
!      Derniere modification : 09 Juin 2008, Didier M. Roche
!-----|--1--------2---------3---------4---------5---------6---------7-|

      SUBROUTINE cforc(ni,nj,nx,ny,nph,npv,fGCM,iE,jE,lgh,lgv,fc,vp     &
                      ,timstpl)

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree :
!       ni, nj  = dimensions du modele Atmosphere (ici ECBilt) = 64, 32
!       nx, ny  = dimension de la grille de l'ISM vers lequel on va
!       nph,npv = taille du recouvrement de grille = 4,4
!       fGCM    = variable climatique sur la grille avec recouvrement
!       iE, jE  = tableaux des longitudes, latitudes des cases ECBilt
!                 sur la grille de l'ISM considere
!       lgh,lgv = ? ? ?
!       vp      = drapeau de travail a valeurs positives
!       timstpl = nombre de pas de temps du fichier
!       Variables de sorties :
!       fc = tableau de valeurs climatiques sur la grille de l'ISM
!
!         --  initialement pour GISM : nx, ny = 165, 281
!         --  pour GRISLI Nord :       nx, ny = 241, 241
!-----|--1--------2---------3---------4---------5---------6---------7-|

        IMPLICIT NONE

!        INTEGER, PARAMETER :: mois = 12
        INTEGER, INTENT(IN) :: timstpl

        INTEGER vp

       integer ni,nj,nx,ny,nph,npv,hX,hY,i,j,k,l,m
       real iE(nx,ny),jE(nx,ny),fGCM(-nph:ni+1+nph,-npv:nj+1+npv,timstpl)
       real lgh(nx,ny,nph),lgv(nx,ny,npv)
       real rect(nph,npv)
       real lgcf(nph,npv),forc(nph,npv),help
       real fc(nx,ny,timstpl)

!-----|--1--------2---------3---------4---------5---------6---------7-|
!
!-----|--1--------2---------3---------4---------5---------6---------7-|

        do m=1,timstpl
          do i=1,nx
            do j=1,ny
              hX=FLOOR(iE(i,j))-nph/2
              hY=FLOOR(jE(i,j))-npv/2
              help=0.0 
              do k=1,nph
                do l=1,npv
                  rect(k,l)=fGCM(hX+k,hY+l,m)
                  lgcf(k,l)=lgh(i,j,k)*lgv(i,j,l)
                  forc(k,l)=rect(k,l)*lgcf(k,l)
                end do
              end do
              do k=1,nph
                do l=1,npv
                  help=help+forc(k,l)
                end do
              end do
              if ((vp.EQ.1).and.(help.lt.0.0)) help = 0.0
              fc(i,j,m)=help
            end do
          end do
        end do

      END SUBROUTINE cforc
