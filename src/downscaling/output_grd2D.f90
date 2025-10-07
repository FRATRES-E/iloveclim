!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Cette routine sert a ecrire une variable 2D dans un fichier
!      compatible avec GRADS (binaire non formatte)
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche 
!      Date   : 18 Juillet 2008
!      Derniere modification : 19 Juin 2009, Didier M. Roche
!-----|--1--------2---------3---------4---------5---------6---------7-|

      SUBROUTINE output_grd2D(tableau,nx,ny,nomfich,nbrec,valundef)

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree  : 
!       Variables de sortie : 
!-----|--1--------2---------3---------4---------5---------6---------7-|

       IMPLICIT NONE

       INTEGER, INTENT(IN) :: nx, ny, nbrec
       REAL, DIMENSION(nx,ny), INTENT(IN) :: tableau
       REAL, INTENT(IN) :: valundef
       CHARACTER(*), INTENT(IN) :: nomfich
!       INTEGER, INTENT(IN), OPTIONAL :: inverse
 

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Variables locales
!-----|--1--------2---------3---------4---------5---------6---------7-|

      INTEGER i,j
      INTEGER nomfichchk_id, nomfichctl_id
      REAL*4, DIMENSION(ny,nx) :: outdata
      CHARACTER*4, PARAMETER :: extchk='.chk',extctl='.ctl'

      outdata = 0.0

       do i=1,ny
         do j=1,nx
           outdata(i,j)=tableau(j,i)
         enddo
       enddo


      open(newunit=nomfichchk_id,file=nomfich//extchk,                  &
       form='unformatted', access='direct', status='unknown',           &
       convert='little_endian', recl=Size(outdata)*Kind(outdata(1,1)))

      write(nomfichchk_id,REC=nbrec) outdata

      close(nomfichchk_id)


! dmr Ecriture du .ctl qui va avec      
      open(newunit=nomfichctl_id,file=nomfich//extctl)
      write(nomfichctl_id,fmt="('dset   ^"//nomfich//extchk//"')")
      write(nomfichctl_id,fmt="('undef ',1p,e12.4)") valundef
      write(nomfichctl_id,fmt="('title Automatic output ECBilt var1')")
      write(nomfichctl_id,                                              &
       fmt="('xdef ',i3,' linear ',2f7.3)") 64,0.00,5.625
      write(nomfichctl_id,fmt="('ydef ',i3,' levels')") 32
      write(nomfichctl_id,                                              &
       fmt="(' -85.7606 -80.2688 -74.7445 -69.2130 -63.6786')")
      write(nomfichctl_id,                                              &
       fmt="(' -58.1430 -52.6065 -47.0696 -41.5325 -35.9951')")
      write(nomfichctl_id,                                              &
       fmt="(' -30.4576 -24.9199 -19.3822 -13.8445 -8.30670')")
      write(nomfichctl_id,                                              &
       fmt="(' -2.76890 2.76890 8.30670 13.8445 19.3822')")
      write(nomfichctl_id,                                              &
       fmt="('  24.9199 30.4576 35.9951 41.5325 47.0696')")
      write(nomfichctl_id,                                              &
       fmt="('  52.6065 58.1430 63.6786 69.2130 74.7445')")
      write(nomfichctl_id,fmt="('  80.2688 85.7606')")
      write(nomfichctl_id,                                              &
       fmt="('zdef ',i3,' linear ',2f7.2)") 1,0.00,1.00
      write(nomfichctl_id,                                              &
       fmt="('tdef ',i4,' linear 1jan0001  1YR')") nbrec
      write(nomfichctl_id,fmt="('vars 1')")
      write(nomfichctl_id,                                              &
       fmt="('var1    1    99 ECBilt variable var1 check')")
      write(nomfichctl_id,fmt="('endvars')")
      close(nomfichctl_id)
      
      END SUBROUTINE output_grd2D
