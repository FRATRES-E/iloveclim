!-----|--1--------2---------3---------4---------5---------6---------7-|
!
!      Auteur : Didier M. Roche
!      Date   : 17 juin 2009
!      Derniere modification : 16 septembre 2010, Didier R. & Christophe D.
!
!-----|--1--------2---------3---------4---------5---------6---------7-|

      SUBROUTINE downscale_vertical(champi,nx,ny,tmax,tmin,hmax,hmin,   &
                        ex,ey,iGRID,jGRID,sGRIS,sECB,mmmois, maskGRIS)

      USE input_icemask

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nx, ny, ex, ey

      REAL, DIMENSION(nx,ny), INTENT(INOUT)  :: champi
      REAL, DIMENSION(nx,ny), INTENT(IN)  :: sGRIS, maskGRIS
      REAL, DIMENSION(nx,ny), INTENT(IN)  :: iGRID, jGRID
      REAL, DIMENSION(nx,ny), INTENT(IN) :: sECB
      REAL, DIMENSION(ey,ex), INTENT(IN) :: tmax,tmin
      REAL, DIMENSION(ex,ey), INTENT(IN) :: hmax,hmin
      INTEGER, INTENT(IN) :: mmmois

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Variables locales
!-----|--1--------2---------3---------4---------5---------6---------7-|
      INTEGER :: i,n,ii,jj, i1, j1, NOMBRE

      REAL, DIMENSION(nx,ny) :: dtemp

      REAL, DIMENSION(ex,ey) :: A_HC, A_EC, SOMME, Q
      REAL, DIMENSION(ex,ey) :: aveE,adevE,sdevE,varE,skewE,kurtE
      REAL, DIMENSION(ex,ey,300) :: distrib_alt
      INTEGER, DIMENSION(ex,ey) :: nb_remplis

      REAL :: A
      REAl, PARAMETER :: petit = -10000000.0
      CHARACTER*2 tmois

      REAL, DIMENSION(9) :: tableaux, tableauy
      LOGICAL, DIMENSION(9) :: tableau_masque
      REAL :: coefa, coefb, siga, sigb, chi, q_l, val_max, val_min
      REAL :: coef_lin, coef_min, coef_max
      INTEGER :: nb_vals, iiiii, indx
      INTEGER, DIMENSION(1) :: indx_max, indx_min

      INTEGER, PARAMETER :: mon_choix = 0
!! attention une valeur autre que zero est impossible avec sECB(nx,ny)
!!  (reliquat de sECB(sur ECBilt)

      INTEGER :: checkgammaECBmoistxt_id

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Initialisations (tableaux) ...
!-----|--1--------2---------3---------4---------5---------6---------7-|

      nb_remplis = 0
      distrib_alt = 0.0

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Etablissement d'un tableau de distribution de la topographie
!       GRISLI dans chaque case ECBilt
!-----|--1--------2---------3---------4---------5---------6---------7-|

      DO i=1,nx
        DO n=1,ny
            jj = MIN(ex,FLOOR(jGRID(i,n)))
            ii = FLOOR(iGRID(i,n))

            nb_remplis(jj,ii) = nb_remplis(jj,ii) + 1
            distrib_alt(jj,ii,nb_remplis(jj,ii)) = sGRIS(i,n)

      END DO
      END DO

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Calcul des moments d'altitude par case ECBilt
!-----|--1--------2---------3---------4---------5---------6---------7-|

      DO jj=1,ex
      DO ii=1,ey

        IF (nb_remplis(jj,ii).GT.0) THEN
          CALL moment(distrib_alt(jj,ii,1:nb_remplis(jj,ii)),           &
                     nb_remplis(jj,ii),aveE(jj,ii),adevE(jj,ii),        &
                     sdevE(jj,ii),varE(jj,ii),skewE(jj,ii),kurtE(jj,ii) &
                     )
        ENDIF

      END DO
      END DO

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Sortie directe des moments si necessaire ...
!-----|--1--------2---------3---------4---------5---------6---------7-|

!      if (mmmois.EQ.1) then
!        open(1,file='checkmoments_sdevE.txt',form='formatted')
!        open(2,file='checkmoments_skewE.txt',form='formatted')
!        open(3,file='checkmoments_kurtE.txt',form='formatted')
!
!        DO i=1,nx
!        DO n=1,ny
!            jj = MIN(ex,FLOOR(jGRID(i,n)))
!            ii = FLOOR(iGRID(i,n))
!
!          write(1,'(2I5,1F15.4)') i,n, sdevE(jj,ii)
!          write(2,'(2I5,1F15.4)') i,n, skewE(jj,ii)
!          write(3,'(2I5,1F15.4)') i,n, kurtE(jj,ii)
!
!        END DO
!        END DO
!
!        close(1)
!        close(2)
!        close(3)
!      endif

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Calcul de la descente d'échelle proprement dite
!-----|--1--------2---------3---------4---------5---------6---------7-|


      dtemp = petit

!-----|--1--------2---------3---------4---------5---------6---------7-|

      WRITE(tmois,'(I2.2)') mmmois
      open(newunit=checkgammaECBmoistxt_id,file='outputdata/coupler/checkgammaECB_'//tmois//'.txt'    &
            ,form='formatted')

!-----|--1--------2---------3---------4---------5---------6---------7-|

      DO jj=1,ex
      DO ii=1,ey

       IF ((hmax(jj,ii)-hmin(jj,ii)).NE.0.0) THEN
       A_HC(jj,ii) = (tmax(ii,jj)-tmin(ii,jj))/(hmax(jj,ii)-hmin(jj,ii))
       ELSE
       A_HC(jj,ii) = 0.0
       ENDIF

      END DO
      END DO


!-----|--1--------2---------3---------4---------5---------6---------7-|
      IF (mon_choix.NE.0) THEN
!-----|--1--------2---------3---------4---------5---------6---------7-|


      DO jj=1,ex
        DO ii=1,ey

         tableaux = 0.0
         tableauy = 0.0

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       On considere une case ECBilt et ses plus proche voisins
!-----|--1--------2---------3---------4---------5---------6---------7-|

         DO i=-1,1,1
          DO n=-1,1,1
            i1 = ii+i
            if (i1.gt.ey) i1 = i1-ey
            if (i1.lt.1)  i1 = ey+i1

            j1 = jj+n
            if (j1.gt.ex) then
              j1 = ex - (j1-ex)
              i1 = i1 + nint(ey/2.)
            endif

            if (j1.lt.1) then
              j1 = abs(j1)
              i1 = i1 + nint(ey/2.)
            endif

             tableaux(5+3*i+n) = sECB(j1,i1) ! altitude
             tableauy(5+3*i+n) = champi(j1,i1) ! temperature
          ENDDO
         ENDDO


!-----|--1--------2---------3---------4---------5---------6---------7-|
!       On ne considere pas les niveaux inférieurs à 50 mètres ...

         tableau_masque = .TRUE.

         WHERE(tableaux.LT.50)
           tableau_masque = .FALSE.
         ENDWHERE

!-----|--1--------2---------3---------4---------5---------6---------7-|


!-----|--1--------2---------3---------4---------5---------6---------7-|
!       On ne s'intéresse qu'aux zones ou il y a plus de deux points
!        définis et on calcule la droite qui passe par ces points

         nb_vals = COUNT(tableau_masque)

         IF (nb_vals.GE.2) THEN

         indx = 0
         DO iiiii = 1, 9
           if (tableau_masque(iiiii)) then
           indx = indx + 1
           tableaux(indx) = tableaux(iiiii)
           tableauy(indx) = tableauy(iiiii)
           endif
         ENDDO

         CALL fit(tableaux,tableauy,nb_vals ,0.0 ,0 , coefa, coefb,     &
                  siga, sigb, chi, q_l)
! Attention le "fit" du Numerical Recipes inverse les coef "a" et "b"
         A_EC(jj,ii) = coefb
         Q(jj,ii) = q_l

        ELSE

          A_EC(jj,ii) = 0.0
          Q(jj,ii) = 0.0


        ENDIF
!-----|--1--------2---------3---------4---------5---------6---------7-|


       ENDDO
       ENDDO


!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Calcul de la descente d'échelle proprement dite
!       a partir du coefficient A_EC précédemment calculé
!       et du deuxième moement de la distribution (variance)

      DO i=1,nx
        DO n=1,ny
            jj = MIN(ex,FLOOR(jGRID(i,n)))
            ii = FLOOR(iGRID(i,n))


            coef_min = 300
            coef_max = 500

            coef_lin = (sdevE(jj,ii)-coef_min)/(coef_max-coef_min)
            coef_lin = MAX(coef_lin,0.0)
            coef_lin = MIN(coef_lin,1.0)

            if (sECB(jj,ii).GT.2500.) coef_lin = 0.0


            if (maskGRIS(i,n).GT.0.0) THEN
              A = coef_lin*A_HC(jj,ii) + (1-coef_lin)*A_EC(jj,ii)
            ELSE
              A = A_HC(jj,ii)
            ENDIF

            dtemp(i,n) = A*(max(sGRIS(i,n),0.0)-sECB(jj,ii))

            write(checkgammaECBmoistxt_id,'(2I5,1F15.4)') i,n, A*1000.

        END DO
      END DO


!-----|--1--------2---------3---------4---------5---------6---------7-|
       ELSE ! mon_choix = 0, downscaling brutal !!

      DO i=1,nx
        DO n=1,ny
            jj = MIN(ex,FLOOR(jGRID(i,n)))
            ii = FLOOR(iGRID(i,n))

            A = A_HC(jj,ii)

            dtemp(i,n) = A*(max(sGRIS(i,n),0.0)-sECB(i,n))

            write(checkgammaECBmoistxt_id,'(2I5,1F15.4)') i,n, A*1000.

        END DO
      END DO

!-----|--1--------2---------3---------4---------5---------6---------7-|
       ENDIF
!-----|--1--------2---------3---------4---------5---------6---------7-|

      CLOSE(checkgammaECBmoistxt_id)

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Application de la descente d'échelle

      WHERE (dtemp.GT.petit)
        champi = champi + dtemp
      END WHERE
!-----|--1--------2---------3---------4---------5---------6---------7-|

      END SUBROUTINE downscale_vertical
