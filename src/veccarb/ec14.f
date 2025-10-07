!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:21 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:21 CET 2009
!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Cette routine C14_PROD sert a calculer la production cosmogenique
!       de 14C en haute atmosphere
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Veronique Mariotti (vm)
!      Date   : 09 octobre 2012
!      Derniere modification : 10 janvier 2013
!-----|--1--------2---------3---------4---------5---------6---------7-|

#if ( KC14 == 1 )

       SUBROUTINE C14_PROD

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree  : 
!       Variables de sortie : 
!-----|--1--------2---------3---------4---------5---------6---------7-|

       use global_constants_mod, only: dblp=>dp, ip

       USE carbone_co2, ONLY: PRODC14, MPRODC14, LIMC14MASSE, N14LIM    
     >                ,C14DEC

       IMPLICIT NONE

!   INSERER ICI LES DECLARATIONS DES VARIABLES D'ENTREE / SORTIE


cvm --- be10.dat has to contains prodC14 in the PC14M column,
cvm --- where the unit is so that prodC14 for preindustrial
cvm --- is equal to 1. In the TPSC14 column, years have to be
cvm --- written with a minus, for example "21ka BP" has to be
cvm --- written "-21000"
#if ( KC14P == 1 )
      open (newunit=bel10dat_id,file='be10.dat',STATUS='unknown')
cvm      PRINT*,'lecture de be10.dat NC14max = ',NC14max
      DO n=1,NC14max
        READ(bel10dat_id,*,END=100) TPSC14(n),PC14M(n)
cvm      PRINT*, n, TPSC14(n), PC14M(n)
      END DO
 100  CONTINUE
      NC14 = n-1
cvm      PRINT*,'fichier be10.dat lu, NC14=',NC14
      close (bel10dat_id)
!-----|--1--------2---------3---------4---------5---------6---------7-|
!vm determination de l'annee courante 
!-----|--1--------2---------3---------4---------5---------6---------7-|

        if (KTIME.eq.1) then !real time
          NYR0 = NYRSR + NYR
        else !not real time : fixed
          NYR0 = NYRSR
        endif
        NYR01 = -NYR0

!-----|--1--------2---------3---------4---------5---------6---------7-|
!vm interpolation des donnees lues
!-----|--1--------2---------3---------4---------5---------6---------7-|
        IF(NYR0.le.TPSC14(1))THEN
            PC14VAR = PC14M(1)
c            print*,'NYR',NYR,NYR0
        ELSE IF (NYR0.ge.TPSC14(NC14))THEN
            PC14VAR = PC14M(NC14)
        ELSE

            DO n = 1,NC14-1
              IF((NYR0.ge.TPSC14(n)).AND.(NYR0.lt.TPSC14(n+1)))
     &            PC14VAR=PC14M(n)+(PC14M(n+1)-PC14M(n))
     &            *(NYR0-TPSC14(n))/(TPSC14(n+1)-TPSC14(n))
            END DO

        ENDIF

        PRODC14 = 2.625e-12 * 40918./48296.75 * PC14VAR
cvm      PRINT*,'PC14VAR, PRODC14: ', PC14VAR, PRODC14

#elif ( KC14P == 0 )

!-----|--1--------2---------3---------4---------5---------6---------7-|
! C14 cosmogenic production from data (Masarik and Beer 1999)
!-----|--1--------2---------3---------4---------5---------6---------7-|
c***************************************************************
c production rate from data = 2.02 atomsC14/cm2/s (+/-10% of error)
c conversion into production rate in total number of C14 atoms producted in one year 
c in the entire atmosphere: production rate = 2.02 * Earth surface (cm2) * 3600 * 24 * 365
c Earth surface = 4 * Pi * R^2 with R = 6371 km = 6.371e8 cm
c production rate = 3.25e26 atomsC14/yr
c new production rate = (prodC14/yr / Nombre d'Avogadro) * (14 / 1e15) * 0.40 ppmC14
c avec Nombre d'Avogadro = 6.0214179e23
c new production rate = 3.0225e-12 ppmC14/yr 
c production rate for Climber (in order to have DC14atm = 0) = 2.533e-12 ppmC14/yr (-16% 
c compared to data value, more than 10%)
c the ratio 40918./48296.75 corresponds to the ratios of global budget of
c carbon in iLOVECLIM vs CLIMBER in PI settings (in GtC)

!        PRODC14 = 2.625e-12                  !ppmC14/yr
!        PRODC14 = 2.625e-12 * 40918./48296.75 * 1.50 !ppmC14/yr
!        PRODC14 = 2.625e-12 * 40918./48296.75 * 1.465 !ppmC14/yr
!        PRODC14 = 2.625e-12 * 40918./48296.75 * 1.428 !ppmC14/yr
!        PRODC14 = 2.625e-12 * 40918./48296.75 * 1.405 !ppmC14/yr
!        PRODC14 = 2.625e-12 * 40918./48296.75 * 1.397 !ppmC14/yr
         PRODC14 = 3.216e-12 * 1.010          !ppmC14/yr
!         PRODC14 = 3.216e-12*1.018d0    !ppmC14/yr 03/07/2013 -- To match atmospheric 14C after change in remineralisation profile
#endif
        mPRODC14 = PRODC14 / 0.40 * 1e15     !mass of carbon in gC14/yr
        limC14masse = mPRODC14 / c14dec      !limit in gC14
        n14lim = (mPRODC14 / 14) / c14dec    !limit in moles C14       

      RETURN
      END SUBROUTINE C14_PROD

!-----|--1--------2---------3---------4---------5---------6---------7-|
c******************************************************************** 
c      C14 ATM MODIFICATION BY COSMOGENIC PRODUCTION 
c                AND RADIOACTIVE DECAY           
c******************************************************************** 
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Veronique Mariotti (vm)
!      Date   : 09 octobre 2012
!      Derniere modification : id.
!-----|--1--------2---------3---------4---------5---------6---------7-|

       SUBROUTINE C14ATM_DP

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree  : 
!       Variables de sortie : 
!-----|--1--------2---------3---------4---------5---------6---------7-|

       USE carbone_co2, ONLY: C14DEC, PRODC14, C14ATM
       USE MOD_SYNC_TIME, ONLY: KENDY

       IMPLICIT NONE

!   INSERER ICI LES DECLARATIONS DES VARIABLES D'ENTREE / SORTIE

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Variables locales
!-----|--1--------2---------3---------4---------5---------6---------7-|

c********************************************************************

cvm ajout de la production cosmogenique dans le C14ATM
#if ( KC14 == 1 ) 
      CALL C14_PROD
#endif

      if(KENDY.eq.1) then
        print*,'C14ATM avt prod', C14ATM, c14dec
        C14ATM = (1-c14dec)*C14ATM
        print*,'c14atm apdec', C14ATM  
        C14ATM = C14ATM + PRODC14 
cvm        print*,'prodc14',PRODC14
        print*,'C14ATM aps prod', C14ATM
      endif
       

      return
      END SUBROUTINE C14ATM_DP

#endif
