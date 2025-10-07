      subroutine calclagrdenomg(gr,ld)
!dmr orig      subroutine calclagrdenomg(gr,vxtt,ld)

        IMPLICIT NONE

        INTEGER :: i,j

        integer gr,n  !dmr orig   ,vxtt(gr)
        real ld(gr)
        do i=1,gr
          n=1
          do j=1,gr
            if(j.ne.i) then
              n=n*(i-j)
            endif
          end do
          if (n.ne.0) then
            ld(i)=1./n
          else
            ld(i)=999999.
          endif  
        end do
      end 
 
      subroutine calclagrnomg(gr,x,ln)
!dmr orig      subroutine calclagrnomg(gr,vxtt,x,ln)

        IMPLICIT NONE

        INTEGER :: i,j

        integer gr  !dmr orig ,vxtt(gr)
        real x,ln(gr),t
        do i=1,gr
          t=1.
          do j=1,gr
            if(j.ne.i) then
              t=t*(x-FLOOR(x)-(j-gr/2))
            endif
          end do
          ln(i)=t
        end do
      end 

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Cette routine calclgcfg est piratee du couplage entre LOVECLIM
!       et les modele de calotte de glace (AG)ISM de P. Huybrechts
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : P. Huybrechts
!      Date   : Inconnue    
!      Derniere modification : 18 Juin 2008, Didier M. Roche
!-----|--1--------2---------3---------4---------5---------6---------7-|

      SUBROUTINE calclgcfg(nx,ny,np,iEcb,lgcf,lc,ldn)
!dmr orig      SUBROUTINE calclgcfg(nx,ny,np,iEcb,lgcf,vtt,lc,ldn)

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree  : 
!        np : taille de repli des bords de grille (orig : 4)
!       Variables de sortie : lgcf
!-----|--1--------2---------3---------4---------5---------6---------7-|

        IMPLICIT NONE

        INTEGER, INTENT(in) :: nx,ny,np
!dmr orig        INTEGER, DIMENSION(np) :: vtt
        REAL, DIMENSION(nx,ny), INTENT(in) :: iEcb
        REAL, DIMENSION(nx,ny,np), INTENT(out) :: lgcf
        REAL, DIMENSION(np), INTENT(out) :: lc,ldn

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Variables locales
!-----|--1--------2---------3---------4---------5---------6---------7-|

        INTEGER :: i,j,k

!-----|--1--------2---------3---------4---------5---------6---------7-|
!                          
!-----|--1--------2---------3---------4---------5---------6---------7-|

!dmr orig        do i=1,np
!dmr orig          vtt(i)=i-np/2
!dmr orig        end do
!dmr orig        call calclagrdenomg(np,vtt,ldn)
        CALL calclagrdenomg(np,ldn)
        do i=1,nx
          do j=1,ny
            call calclagrnomg(np,iEcb(i,j),lc)
            do k=1,np                         
              lgcf(i,j,k)=lc(k)*ldn(k)
            end do
          end do
        end do
      end  
