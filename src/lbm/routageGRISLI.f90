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

       SUBROUTINE routageGRISLI(nx,ny,topo,mask,nexu,exunbi,exunbj,     &
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
       INTEGER, DIMENSION(nx,ny) :: mask
!       INTEGER, DIMENSION(nx,ny) :: est_sortie
       INTEGER, DIMENSION(nx,ny), INTENT(out) :: exuti, exutj
       INTEGER, DIMENSION(nexu), INTENT(out) :: exunbi, exunbj

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Variables locales
!-----|--1--------2---------3---------4---------5---------6---------7-|

       INTEGER :: i, j, k, l, ii, jj, kk, ll, m, n
       INTEGER :: mm,nn
       INTEGER, PARAMETER :: s = 1
       INTEGER :: nb_iter, succes, creux, nb_exu, nb_exo
       REAL :: topo_min
       CHARACTER(LEN=5) :: string_size
       INTEGER, DIMENSION(nx,ny) :: riverbasins
       INTEGER :: curr_rivbas

       exuti(:,:) = 0
       exutj(:,:) = 0
       exunbi(:) = 0
       exunbj(:) = 0
!       est_sortie(:,:) = 0
       curr_rivbas = 0
       riverbasins(:,:) = 0
       nb_exu = 0
       nb_exo = 0
       
       WRITE(string_size,'(i5)') nx
       open(100,file='outlet'//TRIM(ADJUSTL(string_size))//'.txt',      &
              form='formatted')
       open(101,file='mask'//TRIM(ADJUSTL(string_size))//'.txt',        &
              form='formatted')

       DO j=1, ny
         DO i=1, nx

           IF (mask(i,j).EQ.1) THEN
 
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
                   riverbasins(i,j) = riverbasins(kk,ll)
                   succes = 1
                   nb_iter = 0
                   m = -s
                   n = s

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Sinon, on est peut etre arrive a une nouvelle sortie ? (i.e. au
!      bord du mask de la calotte / au bord de la mer)
!-----|--1--------2---------3---------4---------5---------6---------7-|
                 else if (mask(kk,ll).eq.0) then

                   exuti(i,j) = kk
                   exutj(i,j) = ll

                   succes = 1
                   nb_iter = 0
                   m = -s
                   n = s
                   nb_exu = nb_exu + 1
                   exunbi(nb_exu) = kk
                   exunbj(nb_exu) = ll

                   riverbasins(i,j) = nb_exu
                   
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
                    
                     if (mask(mm,nn).eq.0) then
                       topo_min = topo(mm,nn)
                       kk = mm  
                       ll = nn  
                       exuti(i,j) = kk
                       exutj(i,j) = ll
                       succes = 1
                       nb_iter = 0
                       m = -s
                       n = s
                       nb_exu = nb_exu + 1
                       exunbi(nb_exu) = kk
                       exunbj(nb_exu) = ll
                       riverbasins(i,j) = nb_exu
                     endif

                 ENDDO
               ENDDO


                 m = m - 1
                 n = n + 1
               
             endif 

             ENDDO ! sur le while ... succes
            
           ELSE ! mask is zero
              exuti(i,j) = i
              exutj(i,j) = j
              nb_exo = nb_exo - 1
              riverbasins(i,j) = nb_exo
           ENDIF   ! sur le mask ISM
           write(100,'(2I5,1F25.10)') i,j,FLOAT(riverbasins(i,j))
           write(101,'(2I5,1F25.10)') i,j,FLOAT(mask(i,j))
!           WRITE(*,'(2I5,12F25.10)') i,j,&
!                                     FLOAT(riverbasins(i,j))
         ENDDO     ! sur nx
       ENDDO       ! sur ny
           close(100)
           close(101)

       END SUBROUTINE routageGRISLI

