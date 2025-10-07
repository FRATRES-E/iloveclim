!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:45 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:45 CET 2009

      SUBROUTINE diftur
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! Computation of mixing length : vlturb
! Computation of the coefficients for vertical DIFFUsivities (divided by dz) :
!  avsdz, avudz,avqdz.
! DIFFUsivities :  avxdz = sxturb * vlturb * vcturb (=sqrt(q2turb) )
!
!--
!  modif : 2019-12-17, dmr
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      use global_constants_mod, only: dblp=>dp

      use const_mod, only: zero, one, epsil

      use para0_mod, only: imax, jmax, kmax
      use para_mod,  only:
      use bloc0_mod, only: ims1, ims2, jeq, js1, js2, ju1, ju2, ks1, ks2, is1, is2, iu1, iu2, z, avudz, avsdz, tms, avnub &
                         , avkb, q2turb, bvf, zw, kfs, tmu, u, v, unsdzw, unsdz, ust2s
      use bloc_mod,  only: ghmax, ghmin, nstart, numit, q2tmin, vkappa, zlotur, vlturb, avuloc, avqdz, vlmin      &
                         , sqrghm, tm2tur
      use ice_mod,  only: sdvt

      use newunit_clio_mod, only: clio3_out_id

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


      implicit none

      real(kind=dblp), dimension(imax,jmax,kmax):: vcturb

      real(kind=dblp),                     save :: ref0n2, ref0n
      real(kind=dblp), dimension(kmax),    save :: avq0ri, avs0ri, avu0ri
      real(kind=dblp) :: avcom,  db, ds, epsil2, ghturb, q2blmi, q2slon, riums1, riums2, riums3, riums4, sqbvf, sqturb, ssturb &
                       , suturb, tm4u, uz, vlneut, vz, zinbot, zintop, zintr, zrimax
      integer       :: i, j, k
      integer, save :: kk0ri, kk1ri


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  1 ) 1ere Iter, Determine 1er niv. a calculer ; Initialisation .     |
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      q2blmi = 1.d-6
      zrimax = 0.7d0
      q2slon = 6.51d0
      if (numit.le.nstart) then
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!dn2   ci1=1.0 - 1.0/2.0
!dn2   ci2=1.0 - sqrt(2.0)/2.0
        write(clio3_out_id,'(2A,F6.2,A1,I4,A1)') ' diftur : mixing length from algebric formula '

!--1er Niveau a calculer :
      kk0ri = ks1
      avu0ri(ks1) = 0.0
      avs0ri(ks1) = 0.0
      avq0ri(ks1) = unsdz(ks1) * avkb(ks1+1)
      do k=ks1+1,ks2
        avu0ri(k) = unsdzw(k) * avnub(k)
        avs0ri(k) = unsdzw(k) * avkb(k)
        avq0ri(k) = unsdz(k)  * avkb(k)
      enddo
      kk1ri = kk0ri + 1

       do j=js1,js2
        do i=ims1,ims2
           q2turb(i,j,ks2+1) = max(q2tmin,q2slon*ust2s(i,j)*(one+sdvt(i,j)))
        enddo
       enddo

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!- Fin du traitement specifique de la 1ere Iter.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      endif

      epsil2 = epsil * epsil
!- valeur typique de N2 : 1.E-5 (s-2)
      ref0n2 = 1.d-5 * epsil
      ref0n  = sqrt(ref0n2)

!--Initialisation :
      do k=ks1+1,kk0ri
       do j=js1,js2
        do i=ims1,ims2
         avudz(i,j,k) = avu0ri(k)
         avsdz(i,j,k) = avs0ri(k)
         avqdz(i,j,k) = avq0ri(k)
        enddo
       enddo
      enddo

      do k=kk0ri,ks2
       do j=js1,js2
        do i=ims1,ims2
         avqdz(i,j,k) = 0.0
        enddo
       enddo
      enddo


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  2 ) Mixing Length :                                                 |
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!

      do k=kk1ri,ks2
       do j=js1,js2
        do i=is1(j),is2(j)
          vcturb(i,j,k) = sqrt ( q2turb(i,j,k) )
          sqbvf  = sqrt (max(zero,bvf(i,j,k))) + ref0n

! Mixing length as an algebric fonction
!
          ds     = zw(ks2+1) - zw(k)
          db     = max (epsil, zw(k) - zw(kfs(i,j)))
          zintr  = vkappa*(ds*db)/(ds+db)
          vlneut = zlotur*zintr/(zintr+zlotur)
          vlturb(i,j,k) = min (vlneut, sqrghm*vcturb(i,j,k)/sqbvf) + vlmin
        enddo
       enddo
      enddo

      do j=js1,js2
       do i=is1(j),is2(j)
        vlturb(i,j,ks2) = max(vlturb(i,j,ks2),-vkappa*zw(ks2))
       enddo
      enddo


!--Debut de la boucle externe sur l'indice de niveau k :
      do k=kk1ri,ks2
!-----

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  3 ) M2turb                                                          |
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      do j=js1,js2
        do i=is1(j),is2(j)

          tm4u = tmu(i,j,k-1) + tmu(i+1,j+1,k-1) + tmu(i,j+1,k-1) + tmu(i+1,j,k-1) + epsil2

          uz = u(i+1,j+1,k) - u(i+1,j+1,k-1)
          vz = v(i+1,j+1,k) - v(i+1,j+1,k-1)
          riums1 = (uz*uz + vz*vz) * tmu(i+1,j+1,k-1)

          uz = u(i+1,j,k) - u(i+1,j,k-1)
          vz = v(i+1,j,k) - v(i+1,j,k-1)
          riums2 = (uz*uz + vz*vz) * tmu(i+1,j,k-1)

          uz = u(i,j+1,k) - u(i,j+1,k-1)
          vz = v(i,j+1,k) - v(i,j+1,k-1)
          riums3 = (uz*uz + vz*vz) * tmu(i,j+1,k-1)

          uz = u(i,j,k) - u(i,j,k-1)
          vz = v(i,j,k) - v(i,j,k-1)
          riums4 = (uz*uz + vz*vz) * tmu(i,j,k-1)

          tm2tur(i,j,k) = max ( riums1 , riums2 , riums3 , riums4 )*max(zero,sign(one,tm4u-one))* unsdzw(k) * unsdzw(k)
!    &                   * (1.0+max(zero,sign(one,dfloat(k-kajul))) )
!var &    *(1.0+varfor*sdvt(i,j)*max(zero,sign(one,dfloat(k-kajul))) )

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  4 ) Stabilty FUNCTIONs .                                            |
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


          ghturb = - (vlturb(i,j,k)*vlturb(i,j,k))*bvf(i,j,k)/max(q2turb(i,j,k),epsil2)
          ghturb = max(ghmin,min(ghmax,ghturb))
          sqturb = 0.2d0
          ssturb = (0.494d0)/(1.-34.7d0*ghturb)
          suturb = (0.393d0-3.09d0*ghturb)/(1.0 - 40.8d0*ghturb + 212.*ghturb*ghturb)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  5 ) Diffusion coefficients .                                        |
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

          avcom = vlturb(i,j,k) * vcturb(i,j,k)
          avuloc(i,j,k)= max(suturb * avcom * unsdzw(k), avu0ri(k))
          avsdz(i,j,k) = max(ssturb * avcom * unsdzw(k), avs0ri(k))

! background diffusivity as a FUNCTION of N2

! Computation of an averaged N2
!dn2      zaxmoy = max(epsil,tms(i,j,k)
!dn2 &             +ci1*(tms(i-1,j,k)+tms(i,j-1,k)
!dn2 &                  +tms(i+1,j,k)+tms(i,j+1,k))
!dn2 &             +c12*(tms(i-1,j-1,k)+tms(i-1,j+1,k)
!dn2 &                  +tms(i+1,j-1,k)+tms(i+1,j+1,k))  )
!dn2      bvfmoy = (tms(i,j,k)*bvf(i,j,k)
!dn2 &             +ci1*(tms(i-1,j,k)*bvf(i-1,j,k)
!dn2 &                  +tms(i,j-1,k)*bvf(i,j-1,k)
!dn2 &                  +tms(i+1,j,k)*bvf(i+1,j,k)
!dn2 &                  +tms(i,j+1,k)*bvf(i,j+1,k))
!dn2 &             +c12*(tms(i-1,j-1,k)*bvf(i-1,j-1,k)
!dn2 &                  +tms(i-1,j+1,k)*bvf(i-1,j+1,k)
!dn2 &                  +tms(i+1,j-1,k)*bvf(i+1,j-1,k)
!dn2 &                  +tms(i+1,j+1,k)*bvf(i+1,j+1,k))  )
!dn2 &             /zaxmoy
!dn2      avsmir = max((1.D-7/sqrt(max(4.D-6*one,bvfmoy)) )
!dn2 &               * unsdzw(k), avs0ri(k) )
!dn2      avsdz(i,j,k) = max (ssturb * avcom * unsdzw(k), avsmir)


! !! location of avqdz(k) : idem scal (k)

          zintop = (zw(k+1)-z(k))*unsdz(k)
          zinbot = (z(k-1)-zw(k-1))*unsdz(k-1)
          avqdz(i,j,k) = avqdz(i,j,k) + zintop * max (sqturb * avcom *unsdz(k), avq0ri(k))
          avqdz(i,j,k-1) = avqdz(i,j,k-1) + zinbot * max (sqturb * avcom *unsdz(k-1), avq0ri(k))

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  6 ) Diffusion in the highly sheared region below the ML             |
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!         ztest = abs (-1.*max(zero,sign(one,q2turb(i,j,k)-q2blmi))
!    &         + one-max(zero,sign(one,q2turb(i,j,k+1)-q2blmi)))
!         zri   = bvf(i,j,k)/max(tm2tur(i,j,k),epsil)
!         zri   = max(zero,min(zri,zrimax) )
!         zkhri = (5.0d-3*(1-(zri/zrimax)*(zri/zrimax) )**3 )*unsdzw(k)
!
!         avuloc(i,j,k)= avuloc(i,j,k)+ (1.0-ztest)*zkhri*tms(i,j,k-1)
!         avsdz(i,j,k) = avsdz(i,j,k) + (1.0-ztest)*zkhri*tms(i,j,k-1)

        enddo
       enddo

!-raccord cyclique pour avuloc
       do j=1,jeq
          avuloc(1,j,k)=avuloc(imax-1,j,k)
          avuloc(imax,j,k)=avuloc(2,j,k)
       enddo

       do j=ju1,ju2
        do i=iu1(j),iu2(j)
!        avudz(i,j,k)=.25* ( avuloc(i,j,k)+avuloc(i-1,j,k)+
!    &                       avuloc(i,j-1,k)+avuloc(i-1,j-1,k) )
         avudz(i,j,k)=( avuloc(i,j,k)*tms(i,j,k)                 &
                       +avuloc(i-1,j,k)*tms(i-1,j,k)             &
                       +avuloc(i,j-1,k)*tms(i,j-1,k)             &
                       +avuloc(i-1,j-1,k)*tms(i-1,j-1,k) )/      &
                     max(tms(i,j,k)+tms(i-1,j,k)+tms(i,j-1,k)    &
                        +tms(i-1,j-1,k),q2blmi)

        enddo
       enddo



      enddo

      return

      end subroutine diftur

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the subroutine here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
