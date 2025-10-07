!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Cette routine "moyennage" sert a moyenner les donnees en pas
!       de temps iatm vers de la moyenne mensuelle
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche
!      Date   : 17 Novembre 2008
!      Derniere modification : 17 Novembre 2008, 01 mars 2016
!-----|--1--------2---------3---------4---------5---------6---------7-|

       subroutine MOYEN_AGE(entree,sortie,step)

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree  : entree
!       Variables de sortie : sortie
!-----|--1--------2---------3---------4---------5---------6---------7-|

       implicit none

       integer, intent(in) :: step

       ! Be warned that sortie is TRANSPOSED w.r.t. entree ...
       double precision, dimension(:,:)  , intent(in)    :: entree    ! entree is ny,nx
       double precision, dimension(:,:,:), intent(inout) :: sortie    ! sortie is nx,ny,nz

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Variables locales
!-----|--1--------2---------3---------4---------5---------6---------7-|

       integer :: i,j

!-----|--1--------2---------3---------4---------5---------6---------7-|
! dmr   Main computation
!-----|--1--------2---------3---------4---------5---------6---------7-|


       write(*,*) "Testing moyen_age = "
       write(*,*) "sortie : ", UBOUND(sortie,1), UBOUND(sortie,2), UBOUND(sortie,3), SIZE(sortie)
       write(*,*) "entree : ", UBOUND(entree,1), UBOUND(entree,2), SIZE(entree)

       forall (i=LBOUND(sortie,1):UBOUND(sortie,1),                     &
               j=LBOUND(sortie,2):UBOUND(sortie,2))                     &
                   sortie(i,j,step) = sortie(i,j,step)+entree(j,i)

       return
       end subroutine MOYEN_AGE
