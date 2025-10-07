#include "choixcomposantes.h"
!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Cette routine sert a faire pointer les variables de sortie de
!       la partie couplage vers les variables internes de GRISLI en
!       entree
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche
!      Date   : 15 Juin 2008
!      Derniere modification: 27 July 2011, Marianne Bugelmayer, 4 Mars
!      2014 Didier M. Roche
!-----|--1--------2---------3---------4---------5---------6---------7-|

       SUBROUTINE faitPointer

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree  :
!       Variables de sortie :
!-----|--1--------2---------3---------4---------5---------6---------7-|

! afq deprecated 11-07-17       USE input_CLIO, ONLY: sstCLIO

       USE input_subgrid2L, ONLY: tfyearnord, pfyearnord                      &
                           , SMB_iLCMnord                                     &
                           , topoGRIS_sansnivomernord,nivo_mer                &
                           , epaisglaceGRISnord

#if ( SHELFMELT == 1 )
       USE input_subgrid2L, ONLY: bmshelf_clionord
#endif


!afq for CONSEAU       USE MODULE3D_PHY, ONLY: TANN, ACC, TJULY, S, H, BSOC, BM,        &
!afq for CONSEAU                               sealevel,calvin_GRIS,RUNOF_OC, BMELT_OC
       USE MODULE3D_PHY, ONLY: TANN, ACC, TJULY, S, H, BSOC, BM, sealevel

#if ( SHELFMELT == 1 )
       USE MODULE3D_PHY, ONLY: bmshelfCLIO
#endif

#if ( CONSEAU == 1 )
       use module3D_phy, only: trendWAC,smbWAC,bmeltWAC,calvingWAC
       use varsCONSEAU_mod, only: trendgrisnord,smbgrisnord,bmeltgrisnord,calgrisnord
#endif

       IMPLICIT NONE

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Lien vers les variables GRISLI
!-----|--1--------2---------3---------4---------5---------6---------7-|

!      Equivalent GRISLI : TANN
       tfyearnord=>TANN

!      Equivalent GRISLI : ACC
! -- afq: this is not necessary, only here for output purposes
       pfyearnord=>ACC

!      Equivalent GRISLI : BM
       SMB_iLCMnord=>BM

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Lien vers les variables CLIO
!-----|--1--------2---------3---------4---------5---------6---------7-|
#if ( SHELFMELT == 1 )
        bmshelf_clionord=>bmshelfCLIO
#endif

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Pointage coupleur => variables de GRISLI
!-----|--1--------2---------3---------4---------5---------6---------7-|

!      Equivalent GRISLI : S (!!)
        topoGRIS_sansnivomernord=>S

!      Equivalent GRISLI :
        nivo_mer=>sealevel

!      Equivalent GRISLI : H (!!)
        epaisglaceGRISnord=>H

#if ( CONSEAU == 1 )
        trendgrisnord=>trendWAC
        smbgrisnord=>smbWAC
        bmeltgrisnord=>bmeltWAC
        calgrisnord=>calvingWAC
#endif

        return
       END SUBROUTINE faitPointer
