!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Cette routine <squelette> sert a ???
!       et caetera, et caetera
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche 
!      Date   : 25 mai 2011 
!      Derniere modification : 25 mai 2011 
!-----|--1--------2---------3---------4---------5---------6---------7-|

       SUBROUTINE format_nc_to_ec(data_in,data_out,nx,ny)

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree  : 
!                            data_in, le champs au format nc
!                            nx,ny ; la taille de tableau d'entree
!       Variables de sortie : 
!                            data_out, le champs au format ECBilt
!-----|--1--------2---------3---------4---------5---------6---------7-|

!   INSERER ICI LES EVENTUELS "USE MODULE"

       IMPLICIT NONE

!   INSERER ICI LES DECLARATIONS DES VARIABLES D'ENTREE / SORTIE

       REAL, DIMENSION(nx,ny), INTENT(IN) :: data_in
       REAL, DIMENSION(ny,nx), INTENT(out) :: data_out
 
       INTEGER, INTENT(IN) :: nx, ny

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Variables locales
!-----|--1--------2---------3---------4---------5---------6---------7-|
       REAL, DIMENSION(ny,nx) :: var_temp
       INTEGER :: i

!-----|--1--------2---------3---------4---------5---------6---------7-|
!-----|--1--------2---------3---------4---------5---------6---------7-|

       var_temp = TRANSPOSE(data_in)
       
       FORALL(i=1:UBOUND(data_out,1))                                   &
         data_out(ny+1-i,:) = var_temp(i,:)
      
       END SUBROUTINE format_nc_to_ec
