!dmr -- Added optional components choice - Tue Jun 14 15:59:05 CEST 2011
#include "choixcomposantes.h"
!dmr -- Added optional components choice - Tue Jun 14 15:59:05 CEST 2011

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Ce module est un portage en FORTRAN 90 du module VECODE initial
!       cree en FORTRAN 77.
!      (dans l'environnement logiciel LUDUS)
!
!      Auteur : ??, Didier M. Roche
!      Date   : ??, 06 juin 2018
!      Derniere modification : 16 octobre 2019
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module veget_mod

      use global_constants_mod, only: dblp=>dp, silp=>sp, ip
      use comatm, only: nlat, nlon

      implicit none

      private:: nlat, nlon

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  modif : 19/05/00 <- size of array replaced by (nlat,nlon)
!        : may 2007 <- different enrichment factors for tree and grass (/riche/)
!        : july 2008 <- deforestation scenario (/forets/, parameter mdfor) and cleaning (suppressed array data_crop)
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      integer, parameter           ::  nm=21, nvl=3, nsh=((nm+1)*(nm+2))/2, nsh2=2*nsh, iveg = 600 , mdfor = 1500

      integer                      :: numvegvar

      real(dblp), dimension(nlat)  :: phi
      real(dblp), dimension(80,20) :: newvegvar
      real(dblp)                   :: veg_fill_value, veg_missing_value
      character*60 namevegvar(80,6)

!c---- new variables included to use netCDF output in vecode !mohr

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!~ c 2 files "bio.inc" & "buffer.inc" combined in one : "veget.h" (09/12/99)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(dblp) :: a,bet,gamm,gamm2,fmax,avecube,tmin,npp,nppmax,v1,v2,v3,c1t,c2t,c3t,c1g,c2g,c3g,d1t,d2t,d3t,d1g,d2g,d3g         &
                   ,e1t,e2t,e3t,e1g,e2g,e3g,f1t,f2t,f3t,f1g,f2g,f3g,k1t,k2t,k3t,k1g,k2g,k3g,t1t,t2t,t3t,t4t,t1g,t2g,t3g,t4g
#if ( PF_CC > 0 )
      real(dblp) :: t5t,t6t,t5g,t6g
#endif
      real(dblp) :: ps1,ps2,ps3,ps4,ps5,soilt,forshare_st,t1tn,t1td, desshare_st,nlshare_st,deng,dentd,dentn,laig,lait,ave_t       &
                   ,ave_pr,ave_pr05,desmin,desmax,ave_pr05_des,ades,acr,k0t,k0g,k4g

      integer(ip):: lon, lat
      integer(ip), dimension(nlat,nlon) :: init_flag

      real(dblp), dimension(nlat,nlon) :: b1t,b2t,b3t,b4t,b1g,b2g,b3g,b4g,carea

      real(dblp) :: gdd0,gdd0_min,gdd0_max,fr_ndx,co2ghg,acwd,acwt,acwg,acwn,zrd,zrt,zrg,zrn &
                   ,rsd,rst,rsg,rsn
#if ( PF_CC > 0 )
      real(dblp), dimension(nlat,nlon) :: pf_percent,b5t, b6t, b5g, b6g
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!~ *********************************************************************
!~ ***  BUFFER COMMON: DATA TRANSFER CLIMATE <-> TERRESTR VEG MODEL  ***
!~ *********************************************************************
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      integer(ip):: nstat,ieqveg,ieqvegwr,iscendef,ivegstrt
      real(dblp), dimension(nlat,nlon) :: st = 0.0_dblp,sg =0.0_dblp,sd= 0.0_dblp, snlt,anup,pnpp,b12,b34,b1,b2,b3,b4,anup_moy, &
                    stock,st_moy,stR,sgR,sdR,snltR
#if ( FROG_EXP > 0 )
      real(dblp), dimension(nlat,nlon) :: Fv, Fv_t, Fv_g
#endif

!~ [DEPRECATED] -- variables for the LOCH model
!~ anuploch,stockloch,

      real(dblp), dimension(nlat,nlon,2) :: blai

!~ c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!~ c-AM (2007)
!~ c CO2 enrichment :
!~ c  betat for tree, betag for grass
!~ c  nppt = tree npp, nppg = grass npp
       real(dblp) :: BETAG,BETAT,NPPG,NPPT

!~ c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!~ c-AM (2008)
!~ c Deforestation scenario (data read in VEGET.dat):
!~ c  farea = fraction of the mesh occupied by crops according to R&F (1999)
!~ c  ndfor = nbr of records in deforestation scenario (ndfor =< mdfor)
!~ c  i0dfor = year-1 of first available data in VEGET.dat (date A.D.)
!~ c Reference vegetation distribution against which deforestation is evaluated:
!~ c  sd_const(nlat,nlon) = reference desert distribution.
!~ c  st_const(nlat,nlon) = reference tree distribution.
!~ c Would CO2 flux from vegetation influence atm. CO2 or not : fco2veg

       real(dblp), dimension(nlat,nlon,mdfor) :: farea
       real(dblp), dimension(nlat,nlon)       :: sd_const,st_const
       real(dblp)                             :: fco2veg
       real(silp), dimension(mdfor)           :: VegetTime
       integer(ip)                            :: i0dfor,ndfor
!~        COMMON /FORETS/farea(nlat,nlon,mdfor),sd_const(nlat,nlon),
!~      &                st_const(nlat,nlon),fco2veg,i0dfor,ndfor,VegetTime(mdfor)
!~ c---- VegetTime variable added by !mohr
!~ c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
       end module veget_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr --- The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
