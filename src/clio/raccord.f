!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:51 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:51 CET 2009

      SUBROUTINE raccord(var,spv,krac,ltyp)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  Appele par "class", "?",
!  Mise en place des raccords cycliques et autres raccords
!   pour les "krac" niveaux du tableau "var" .
! tient compte de spv si ltyp > 36 et changement de signe.
!  modif : 25/05/90

!! START_OF_USE_SECTION

      use const_mod, only:

      use para0_mod, only: imax, jmax
      use para_mod,  only: 
      use bloc0_mod, only: jcl1, jcl2, ims1, ims2, ijsdl, ijudl
     &             , iberpm, jberp, ibera, jberam, iberp, iberam
     &             , jberpm, jbera
      use bloc_mod,  only:
!! END_OF_USE_SECTION

      implicit none

!-By reference variables :
      real, dimension(imax,jmax,*), intent(inout) :: var
      real,                         intent(in)    :: spv
      integer,                      intent(in)    :: krac, ltyp

!- variables locales :
      integer, dimension(0:1) :: ijdl
      integer                 :: k, j, ltp

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1 ) raccord cyclique .                                              |
!-----------------------------------------------------------------------

#if ( L_TEST >= 1 )
!dmr [CLIOpH]      if (ltest.ge.1) then
!-
        do k=1,krac
         do j=jcl1,jcl2
          var(ims1-1,j,k) = var(ims2,j,k)
          var(ims2+1,j,k) = var(ims1,j,k)
           enddo
         enddo

!dmr [CLIOpH]      endif
#endif

!dmr [NOEQUI] --- copy of the values instead of pointing to ...
      ijdl(0) = ijsdl
      ijdl(1) = ijudl
!dmr [NOEQUI] ---

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) raccord des 2 grilles .                                         |
!-----------------------------------------------------------------------

#if ( L_TEST == 2 )
!dmr [CLIOpH]      if (ltest.eq.2) then
!-
        llm = mod(ltyp,4) / 3
        llc = ltyp / 12
!-
        if (llc.eq.0) then
          do k=1,krac
           do j=jdl1,jdl2
            ii = ijdl(llm) - j
            var(ims1-1,j,k) = var(ii,jeq-1,k)
            var(ii,jeq,k) = var(ims1,j,k)
           enddo
         enddo
!-
        elseif (llc.eq.1) then
          do k=1,krac
           kk = k + krac
           do j=jdl1,jdl2
            ii = ijdl(llm) - j
            var(ims1-1,j,k) = var(ii,jeq-1,kk)
            var(ii,jeq,kk) = var(ims1,j,k)
           enddo
         enddo

!-
        elseif (llc.eq.2) then
          do k=1,krac
           kk = k + krac
           do j=jdl1,jdl2
            ii = ijdl(llm) - j
            var(ims1-1,j,kk) = -var(ii,jeq-1,k)
            var(ii,jeq,k) = -var(ims1,j,kk)
           enddo
         enddo

!-
        else
          do k=1,krac
           kk = k + krac
           do j=jdl1,jdl2
            ii = ijdl(llm) - j
            var(ims1-1,j,kk) = -var(ii,jeq-1,k)
            var(ii,jeq,k) = -var(ims1,j,kk)
            if (var(ii,jeq-1,k).eq.spv) var(ims1-1,j,kk) = spv
            if (var(ims1,j,kk) .eq.spv) var(ii,jeq,k) = spv
           enddo
         enddo


        endif ! on llc cases ...

!dmr [CLIOpH]      endif
#endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  3 ) raccord pour Bering .                                           |
!-----------------------------------------------------------------------

#if ( L_TEST == 3 )
!     if (ltest.eq.3 .and. iberp.ne.ibera) then
!dmr [CLIOpH]      if (ltest.eq.3) then
      
!-
        ltp = mod(ltyp,4)
!-
!-----position 0 :
        if (ltp.eq.0) then
          do k=1,krac
            var(iberpm,jberp,k) = var(ibera, jberam,k)
            var(iberp, jberp,k) = var(iberam,jberam,k)
            var(iberam,jbera,k) = var(iberp, jberpm,k)
            var(ibera, jbera,k) = var(iberpm,jberpm,k)
           enddo
!-----position 1 :
        elseif (ltp.eq.1) then
          if (ltyp.lt.12) then
            do k=1,krac
              var(iberp,jberp,k) = var(ibera,jberam,k)
              var(ibera,jbera,k) = var(iberp,jberpm,k)
           enddo
          elseif (ltyp.lt.36) then
            do k=1,krac
              var(iberp,jberp,k) = -var(ibera,jberam,k)
              var(ibera,jbera,k) = -var(iberp,jberpm,k)
            enddo
          else
            do k=1,krac
              if (var(iberp,jberpm,k).eq.spv) then
                var(ibera,jbera,k) = spv
              else
                var(ibera,jbera,k) = -var(iberp,jberpm,k)
              endif
              if (var(ibera,jberam,k).eq.spv) then
                var(iberp,jberp,k) = spv
              else
                var(iberp,jberp,k) = -var(ibera,jberam,k)
              endif
            enddo
          endif
!-----position 2 :
        elseif (ltp.eq.2) then
          if (ltyp.lt.12) then
            do k=1,krac
              var(ibera ,jbera,k) = var(iberpm,jberp,k)
              var(iberam,jbera,k) = var(iberp ,jberp,k)
            enddo
          elseif (ltyp.lt.36) then
            do k=1,krac
              var(ibera ,jbera,k) = -var(iberpm,jberp,k)
              var(iberam,jbera,k) = -var(iberp ,jberp,k)
            enddo
          else
            do k=1,krac
              if (var(iberpm,jberp,k).eq.spv) then
                var(ibera, jbera,k) = spv
              else
                var(ibera ,jbera,k) = -var(iberpm,jberp,k)
              endif
              if (var(iberp,jberp,k).eq.spv) then
                var(iberam,jbera,k) = spv
              else
                var(iberam,jbera,k) = -var(iberp ,jberp,k)
              endif
            enddo
          endif
!-----position 3 :
        elseif (ltyp.lt.12) then
          do k=1,krac
            var(ibera,jbera,k) = var(iberp,jberp,k)
          enddo
        elseif (ltyp.lt.36) then
          do k=1,krac
            var(ibera,jbera,k) = -var(iberp,jberp,k)
          enddo
        else
          do k=1,krac
            if (var(iberp,jberp,k).eq.spv) then
              var(ibera,jbera,k) = spv
            else
              var(ibera,jbera,k) = -var(iberp,jberp,k)
            endif
          enddo
        endif

!dmr [CLIOpH]      endif
#endif

      return

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine raccord -
      end subroutine raccord
