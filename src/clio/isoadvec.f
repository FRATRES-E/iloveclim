!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:48 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:48 CET 2009

      SUBROUTINE isoadvec
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- // BALANCE of DIAGONEL FLUXES //
! ============================================================
!   Implementation of Isopycnal thickness diffusion
!   References :
!   -Danabasoglu, G., McWilliams, J.C., P. Gent, 1994:
!    The role of mesoscale tracer transports in the global
!    ocean circulation. Science, 24, 1123-1128
!   -Gent, P.R., and J.C. McWilliams, 1990: Isopycnal mixing
!    in ocean circulation models. JPO, 20, 150-155
! ============================================================
!  modif : 25/05/99

!! START_OF_USE_SECTION
      use global_constants_mod, only: dblp=>dp, ip           

      use const_mod

      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod
      use isoslope_mod

      use newunit_clio_mod, only: clio3_out_id      
!! END_OF_USE_SECTION

!--variables locales conservees d'un appel a l'autre :
      real(kind=dblp), dimension(kmax), save :: aitds2

!--variables locales :
      real(kind=dblp), dimension(imax) :: psix
      real(kind=dblp), dimension(jmax) :: psiy

      integer(kind=ip) :: i, j, k

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      if(numit.eq.nstart) then
        write(clio3_out_id,'(2A,F5.0)') ' *** isoadvec ;  G&M.90  *** :',
     &   ' aitd=', aitd(kmax)

!- initialisation :
        do 10 k=1,kmax
          aitds2(k) = 0.5 * aitd(k)
 10     continue

!- fin du traitement 1ere iter.
      endif

!- compute slopes c4xgm & c4ygm (position 4):
!pp   call isofilter(c4xgm)
!pp   call isofilter(c4ygm)

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- Start  external loop on all Meridional Section, index "j".
!$DIR SPP LOOP_PARALLEL
!$DIR SPP LOOP_PRIVATE(i,k,psix)
      do 300 j=js1,js2

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) Merid. Section : Compute "Bolus velocity" uiso & wiso(1rst part).
!-----------------------------------------------------------------------

!- initialisation :
      do 210 i=is1(j),is2(j)+1
        uiso(i,j,ks1) = 0.
 210  continue

!- Compute Stream Function "psix" (a dy pres) - fill in uiso & wiso :
!- NB : ai*alpha2 is nul at the bottom

      do 230 k=ks1+1,ks2
        do 220 i=is1(j),is2(j)+1
          psix(i) = ttm1(i,j,k-1) * aitds2(k) *
     &            ( c4xgm(i,j,k) + c4xgm(i-1,j,k) )
          uiso(i,j,k-1) = uiso(i,j,k-1) + unsdz(k-1)*psix(i)
          uiso(i,j,k) = -unsdz(k)*psix(i)
          psix(i) = cmy(i,j,1)*psix(i)
 220    continue
        do 225 i=is1(j),is2(j)
          wiso(i,j,k) = unsdx * (psix(i+1)-psix(i))
 225    continue
 230  continue

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- End of external loop on all Meridional Section, index "j".
 300  continue

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- Start  external loop on all Zonal Section, index "i".
!$DIR SPP LOOP_PARALLEL
!$DIR SPP LOOP_PRIVATE(j,k,psiy)
      do 400 i=ims1,ims2

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  3 ) Zonal Section : Compute "Bolus velocity" viso & wiso(2nd part).
!-----------------------------------------------------------------------

!- initialisation :
      do 310 j=jsdom1(i),jsdom2(i)+1
        viso(i,j,ks1) = 0.
 310  continue

!- Compute Stream Function "psiy" (a dx pres) - fill in viso & wiso :
!- NB : ai*alpha2 is nul at the bottom

      do 330 k=ks1+1,ks2
        do 320 j=jsdom1(i),jsdom2(i)+1
          psiy(j) = ttm2(i,j,k-1) * aitds2(k) *
     &            ( c4ygm(i,j,k) + c4ygm(i,j-1,k) )
          viso(i,j,k-1) = viso(i,j,k-1) + unsdz(k-1)*psiy(j)
          viso(i,j,k) = -unsdz(k)*psiy(j)
          psiy(j) = cmx(i,j,2)*psiy(j)
 320    continue
        do 325 j=jsdom1(i),jsdom2(i)
          wiso(i,j,k) = -smxy(i,j,0) *
     &      ( wiso(i,j,k) + unsdy * (psiy(j+1)-psiy(j)) )
 325    continue
 330  continue

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- End of external loop on all Zonal Section, index "i".
 400  continue

!pp   call debugrac(uiso,viso,wiso,'GM90_RACCORD')
      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine isoadvec -
      end
