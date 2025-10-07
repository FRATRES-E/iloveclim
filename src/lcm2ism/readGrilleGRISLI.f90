!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Ce sous-programme lit la grille du modele GRISLI en entree
!       et le stocke dans des tableaux de lat, lon
!
!      Auteur : Didier M. Roche
!      Date   : 15 Juin 2008
!      Derniere modification : 15 Juin 2008
!-----|--1--------2---------3---------4---------5---------6---------7-|

       SUBROUTINE                                                       &
          readGrilleGRISLI(nx,ny,klat,klon,latGRISLI,lonGRISLI,fich,num)

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree  : 
!          nx, ny : taille des tableaux a remplir par la routine
!          fich : nom du fichier a lire
!          num : numero d'unite pour l'ouverture du fichier
!       Variables de sortie : 
!          klat(nx,ny) : contient la distance en latitude  lue
!          klon(nx,ny) : contient la distance en longitude lue
!          latGRISLI(nx,ny) : contient la coordonnee latitude  lue
!          lonGRISLI(nx,ny) : contient la coordonnee longitude lue
!-----|--1--------2---------3---------4---------5---------6---------7-|

       IMPLICIT NONE

       INTEGER, INTENT(in) :: nx, ny, num
       REAL, DIMENSION(nx,ny), INTENT(out) :: klat,klon
       REAL, DIMENSION(nx,ny), INTENT(out) :: latGRISLI,lonGRISLI
       CHARACTER(*), INTENT(in) :: fich

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Variables locales
!-----|--1--------2---------3---------4---------5---------6---------7-|

       INTEGER :: i,j,k, nbpoints
       REAL :: f1, f2, f3, f4


       nbpoints=nx*ny

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Lecture du fichier externe, grille du modele GRISLI
!-----|--1--------2---------3---------4---------5---------6---------7-|

        OPEN(num,file=fich) 
     
        DO k=1, nbpoints 
           READ(num,*) i, j, f1, f2, f3, f4
           klat(i,j) = f2
           klon(i,j) = f1
           latGRISLI(i,j) = f4
           lonGRISLI(i,j) = f3
        ENDDO

        CLOSE(num)


       END SUBROUTINE readGrilleGRISLI
