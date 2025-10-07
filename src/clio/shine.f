!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:53 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:53 CET 2009

!
      SUBROUTINE shine(ih,zmue,tfsn,tfsg,ts,hgbq,
     &                 hnbq,zalb,zalcn,zalbp,zalcnp)
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!
!  Computes albedo of snow-sea ice following SHINE &
!  HENDERSSON-SELLERS [1985].
!
!  albin: Albedo of melting ice in the arctic.
!  albis: Albedo of melting ice in the antarctic (SHINE
!         & HENDERSSON-SELLERS, 1985).
!  cgren: Correction of the snow or ice albedo to take into account
!         effects of cloudiness (GRENFELL & PEROVICH, 1984)
!

!! START_OF_USE_SECTION

      use global_constants_mod, only: dblp=>dp, ip


!! END_OF_USE_SECTION
      integer(kind=ip)::  ih
      real(kind=dblp) ::  al, albice, albin, albis, alphd, alphdi, alphs
     &                  , cgren, hgbq, hnbq, tfsg, tfsn, ts, zalb, zalbp
     &                  , zalcn, zalcnp, zmue

!
!     albin = 0.53
!     albis = 0.53
!     cgren = 0.06
!  albice: Albedo of melting ice.
!     alphd  = 0.80
!     alphdi = 0.72
!     alphs  = 0.65
      albin = 0.45
      albis = 0.45
      albice = 0.45
      alphd  = 0.72
      alphdi = 0.64
      alphs  = 0.55
      cgren = 0.04
!
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  1) Computation of surface albedo.                                   |
!-----------------------------------------------------------------------
!
!  zalbp: Albedo for clear sky.
!  zalb:  Albedo for overcast sky.
!
!
!--1.1 Case of ice or snow.
!--------------------------
!
      if (hnbq.gt.0.0) then
!
!  a) Case of ice covered by snow.
!
        if (ts.lt.tfsn) then
!
!     Freezing snow.
!
          if (hnbq.gt.0.05) then
            zalbp = alphd
          else
            if (hgbq.gt.1.5) then
              zalbp = alphdi+(hnbq*(alphd-alphdi)/0.05)
            else if (hgbq.gt.1.0.and.hgbq.le.1.5) then
                   al = 0.472+2.0*(alphdi-0.472)*(hgbq-1.0)
            else if (hgbq.gt.0.05.and.hgbq.le.1.0) then
                   al = 0.2467+(0.7049*hgbq)-(0.8608*(hgbq*hgbq))+
     &                 (0.3812*(hgbq*hgbq*hgbq))
            else
                   al = 0.1+3.6*hgbq
            endif
            if (hgbq.le.1.5) zalbp=al+(hnbq*(alphd-al)/0.05)
          endif
        else
!
!     Melting snow.
!
          if (hnbq.ge.0.1) then
            zalbp = 0.65
            zalbp = alphs
          else
            zalbp = albice+((alphs-albice)/0.1)*hnbq
          endif
        endif
      else
!
!  b) Case of ice free of snow.
!
        if (ts.lt.tfsg) then
!
!     Freezing ice.
!
          if (hgbq.gt.1.5) then
            zalbp = alphdi
          else if (hgbq.gt.1..and.hgbq.le.1.5) then
                 zalbp = 0.472+2.*(alphdi-0.472)*(hgbq-1.)
          else if (hgbq.gt.0.05.and.hgbq.le.1.) then
                 zalbp = 0.2467+
     &                   (0.7049*hgbq)-(0.8608*(hgbq*hgbq))+
     &                   (0.3812*(hgbq*hgbq*hgbq))
          else
            zalbp = 0.1+3.6*hgbq
          endif
        else
!
!     Melting ice.
!
          if (hgbq.gt.1.5) then
            zalbp = albice
          else if (hgbq.gt.1..and.hgbq.le.1.5)  then
                 zalbp = 0.472+(2.*(albice-0.472)*(hgbq-1.))
          else if (hgbq.gt.0.05.and.hgbq.le.1.) then
                 zalbp = 0.2467+0.7049*hgbq
     &                  -(0.8608*(hgbq*hgbq))
     &                  +(0.3812*(hgbq*hgbq*hgbq))
          else
            zalbp = 0.1+3.6*hgbq
          endif
        endif
      endif
      zalb=zalbp+cgren
!
!--1.2. Case of the ocean.
!-------------------------
!
      zalcnp=0.05/(1.1*zmue**1.4+0.15)
!
!  Parameterization of BRIEGLED AND RAMANATHAN (1982)
!
      zalcn=0.06
!  see KONDRATYEV (1969) AND PAYNE (1972)
      return
      end
