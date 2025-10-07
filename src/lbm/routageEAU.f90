!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Cette routine routageISM sert a router l'eau de l'ISM vers les
!       bords de la calotte, developpe pour GRISLI en vue du couplage
!       avec ECBilt-CLIO-VECODE
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche 
!      Date   : 08 Juillet 2008
!      Derniere modification : 04 décembre 2010, adaptation à ECBilt 
!                                           ... pour le mode en ligne
!-----|--1--------2---------3---------4---------5---------6---------7-|

       SUBROUTINE routageEAU(nx,ny,topo,mask,nexu,exunbi,exunbj,        &
                            exuti,exutj)

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree  : 
!        nx, ny         : taille de la grille ISM
!        topo           : altitude de la surface de l'ISM
!        mask           : mask des cases utilisees de l'ISM
!        latism, lonism : latitude, longitude des cases ISM
!       Variables de sortie : 
!        exuti, exutj   : point de sortie de la case donnee
!        exunbi, exunbj : liste unitaire des exutoires
!-----|--1--------2---------3---------4---------5---------6---------7-|

       IMPLICIT NONE

       INTEGER, INTENT(in) :: nx, ny, nexu
       REAL, DIMENSION(nx,ny), INTENT(in) :: topo
       INTEGER, DIMENSION(nx,ny), INTENT(in) :: mask
       INTEGER, DIMENSION(nx,ny), INTENT(out) :: exuti, exutj
       INTEGER, DIMENSION(nexu), INTENT(out) :: exunbi, exunbj

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Variables locales
!-----|--1--------2---------3---------4---------5---------6---------7-|

       INTEGER :: i, j, k, l, ii, jj, kk, ll, m, n
       INTEGER :: mm,nn
       INTEGER, PARAMETER :: s = 1
       INTEGER :: nb_iter, succes, creux, nb_exu
       INTEGER, DIMENSION(nx,ny) :: est_sortie
       REAL :: topo_min
       INTEGER, DIMENSION(nx,ny) :: lmask
       LOGICAL, PARAMETER :: bricole_caspienne = .TRUE.


       exuti(:,:) = 0
       exutj(:,:) = 0
       exunbi(:) = 0
       exunbj(:) = 0
       est_sortie(:,:) = 0
       nb_exu = 0

       lmask(:,:) = mask(:,:)


       IF (bricole_caspienne) THEN ! we add an output for the Caspian Sea
         j=10
         DO i=23,25
           lmask(i,j) = 0
         ENDDO
         i=24
         DO j=6,7
           lmask(i,j) = 0
         ENDDO
       ENDIF ! end on bricole_caspienne


       DO j=1, ny
         DO i=1, nx

           IF (lmask(i,j).EQ.1) THEN
 
             nb_iter = 0 
             succes = 0
             kk = i
             ll = j
             ii = i
             jj = j
             m = -s
             n = s

             DO WHILE (succes.NE.1)

               nb_iter = nb_iter + 1 
               topo_min = topo(ii,jj)
               creux = 0

               DO k=m, n, n
                 DO l=m, n, n

                     mm = ii+k
                     nn = jj+l

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Attention, il ne faut pas etre circulaire en latitude !!! 
!-----|--1--------2---------3---------4---------5---------6---------7-|
                     if (mm.LT.1) mm = 1  
                     if (mm.GT.nx) mm = nx

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       En revanche il FAUT etre circulaire en longitude !!! 
!-----|--1--------2---------3---------4---------5---------6---------7-|
                     if (nn.LT.1) nn = ny + nn
                     if (nn.GT.ny) nn = nn - ny
                    
!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Si on trouve un creux ...
!-----|--1--------2---------3---------4---------5---------6---------7-|
                   if ((topo(mm,nn).LE.topo_min).AND.                   &
                      (abs(topo(mm,nn)-topo_min).GT.1.E-3)) then


                     topo_min = topo(mm,nn)
                     kk = mm  
                     ll = nn  

                     creux = 1
                   endif

                 ENDDO
               ENDDO

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Arrive ici, soit on a trouve un creux soit non
!-----|--1--------2---------3---------4---------5---------6---------7-|

               if (creux.eq.1) then
                 m = -s
                 n = s

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      On a trouve un creux : s'il a deja une sortie, on la copie
!-----|--1--------2---------3---------4---------5---------6---------7-|

                 if ((exuti(kk,ll).NE.0).AND.(exutj(kk,ll).NE.0)) THEN
!       i,j est toujours la case de depart de notre raisonnement
                   exuti(i,j) = exuti(kk,ll)
                   exutj(i,j) = exutj(kk,ll)
                   succes = 1
                   nb_iter = 0
                   m = -s
                   n = s

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Sinon, on est peut etre arrive a une nouvelle sortie ? (i.e. au
!      bord du mask de la calotte / au bord de la mer)
!-----|--1--------2---------3---------4---------5---------6---------7-|
                 else if (lmask(kk,ll).eq.0) then
                   exuti(i,j) = kk
                   exutj(i,j) = ll
                   est_sortie(kk,ll) = 1
                   succes = 1
                   nb_iter = 0
                   m = -s
                   n = s
                   nb_exu = nb_exu + 1
                   exunbi(nb_exu) = kk
                   exunbj(nb_exu) = ll
                   
!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Il n'a pas de sortie deja existante : on continue de chercher
!        a partir des nouvelles coordonnees
!-----|--1--------2---------3---------4---------5---------6---------7-|
                 else
                   ii = kk
                   jj = ll
                 endif 
                 
!-----|--1--------2---------3---------4---------5---------6---------7-|
!      On a pas trouve de creux : on elargi la recherche d'un point 
!       A faire dans un developpement futur : traitement des bassins
!        endoreiques
!-----|--1--------2---------3---------4---------5---------6---------7-|
               else ! Aie (!) pas de creux

!-----|--1--------2---------3---------4---------5---------6---------7-|
!        On tente une sortie vers l'océan le plus proche s'il existe 
!-----|--1--------2---------3---------4---------5---------6---------7-|
               DO k=m, n, n
                 DO l=m, n, n

                     mm = ii+k
                     nn = jj+l

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Attention, il ne faut pas etre circulaire en latitude !!! 
!-----|--1--------2---------3---------4---------5---------6---------7-|
                     if (mm.LT.1) mm = 1  
                     if (mm.GT.nx) mm = nx

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       En revanche il FAUT etre circulaire en longitude !!! 
!-----|--1--------2---------3---------4---------5---------6---------7-|
                     if (nn.LT.1) nn = ny + nn
                     if (nn.GT.ny) nn = nn - ny
                    
                     if (lmask(mm,nn).eq.0) then
                       topo_min = topo(mm,nn)
                       kk = mm  
                       ll = nn  
                       exuti(i,j) = kk
                       exutj(i,j) = ll
                       est_sortie(kk,ll) = 1
                       succes = 1
                       nb_iter = 0
                       m = -s
                       n = s
                       nb_exu = nb_exu + 1
                       exunbi(nb_exu) = kk
                       exunbj(nb_exu) = ll
                     endif

                 ENDDO
               ENDDO


                 m = m - 1
                 n = n + 1
               
             endif 

             ENDDO ! sur le while ... succes
#if( 1 )
           ELSE ! sur le mask, cas ou il est zero d entree.

               ! si utilisation sans prendre garde, on a exuti(i,j) = 0
               ! ... pas glop! 
               ! solution : on impose exuti(i,j) = i et exutj(i,j) = j
               exuti(i,j) = i
               exutj(i,j) = j
               ! inconvenient : le fichier genere par out_routage n'a
               ! plus de sens !!
#endif
           ENDIF   ! sur le mask ISM


         ENDDO     ! sur nx
       ENDDO       ! sur ny

       IF (bricole_caspienne) THEN ! dump water from Caspian Sea in Med. Sea
       DO j=1, ny
         DO i=1, nx
              if (                                                      &
      (exutj(i,j).EQ.10).AND.((exuti(i,j).GE.23).AND.(exuti(i,j).LE.25))&
                 ) then
               exutj(i,j) = exutj(24,4)
               exuti(i,j) = exuti(24,4)
             endif
             if (                                                       &
      (exuti(i,j).eq.24).AND.((exutj(i,j).eq.6).OR.(exutj(i,j).eq.7))   &
                ) then
               exutj(i,j) = exutj(24,4)
               exuti(i,j) = exuti(24,4)
             endif
         ENDDO     ! sur nx
       ENDDO       ! sur ny

       j=10
       DO i=23,25
         exutj(i,j) = exutj(24,4)
         exuti(i,j) = exuti(24,4)
       ENDDO
       i=24
       DO j=6,7
         exutj(i,j) = exutj(24,4)
         exuti(i,j) = exuti(24,4)
       ENDDO

       ENDIF ! on bricole_caspienne

       END SUBROUTINE routageEAU
