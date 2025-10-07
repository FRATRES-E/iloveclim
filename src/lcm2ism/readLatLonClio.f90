!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Cette routine <squelette> sert a ???
!       et caetera, et caetera
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche 
!      Date   : 15 Juin 2008
!      Derniere modification : 15 Juin 2008
!-----|--1--------2---------3---------4---------5---------6---------7-|

      SUBROUTINE readLatLonClio(nx,ny,zlatt,zlont,fichlat,fichlon)

      use global_constants_mod, only: dblp=>dp, ip

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree  : 
!         nx, ny
!       Variables de sortie : 
!         zlatt, zlont
!-----|--1--------2---------3---------4---------5---------6---------7-|

      IMPLICIT NONE

      INTEGER, INTENT(in) :: nx, ny
      REAL, DIMENSION(nx,ny), INTENT(out) :: zlatt, zlont
      CHARACTER(*), INTENT(in) :: fichlat, fichlon 

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Variables locales
!-----|--1--------2---------3---------4---------5---------6---------7-|

      INTEGER :: i, j
      integer(kind=ip):: fichlat_id, fichlon_id


      open (newunit=fichlat_id,file=fichlat, status='old')
      do j=1,ny
         read(fichlat_id,'(122( F10.5))' ) (zlatt(i,j),i=1,nx)
      enddo
      read(fichlat_id,*)
      close (fichlat_id)

      open (newunit=fichlon_id,file=fichlon,status='old')
      do j=1,ny
         read(fichlon_id,'(122( F10.5))' ) (zlont(i,j),i=1,nx)
      enddo
      write(fichlon_id,*)
      close (fichlon_id)

      END SUBROUTINE readLatLonClio
