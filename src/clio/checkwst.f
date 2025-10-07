!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:43 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:43 CET 2009

      SUBROUTINE checkwst(speval,iirl1,iirl2,jjrl1,jjrl2,nfrc)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! ecrit les donnees de tension de vent, masquee, presentees comme u et v.
! directement appelee par "local" appele par "class"
!  modif : 17/09/94

!! START_OF_USE_SECTION
      use global_constants_mod, only: dblp=>dp, ip

      use const_mod

      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod
      use reper_mod
      
      use newunit_clio_mod, only: clio3_out_id
!! END_OF_USE_SECTION


      real(kind=dblp), dimension(imax,jmax) :: uloc, vloc
c~       dimension uloc(imax,jmax), vloc(imax,jmax)

      real(kind=dblp) :: speval, xrl1, yrl1, zrl1, dzrl, cmult, cadit

      integer(kind=ip):: iirl1, iirl2, jjrl1, jjrl2, nfrc, iim, jjm
     >                 , ltypv, i, j
      
      character*3  fmt3
      character*10 fmt1
      character*30 fmt
      character*34 titx,tity

      integer(ip) :: wst_chekc_id

 1000 format(A30,1x,E13.6,1x,I6,1x,I6,1x,I6,1x,I6)
 1111 format(3(F7.1,1X,F7.3,1X),I3,A)

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      iim = iirl2 - iirl1 + 1
      jjm = jjrl2 - jjrl1 + 1
      ltypv = 3
      xrl1 = xlon1 + dlong * 0.5
      yrl1 = ylat1 + dlat  * 0.5
      zrl1 = -1.
      dzrl = 0.001

      cmult = -1.d+04
      cadit = 0.0
!--definition des formats :
      fmt1 = 'F8.4'
      write(fmt3,'(I3)') nfrc
      fmt = '('//fmt3//fmt1(:8)//')'

      if (speval.eq.0.0) then
        write (clio3_out_id,*) 'fichier de Tension de vent SANS valeur-speciale'
      else
        write (clio3_out_id,*) 'format , valeur-specale :', fmt, speval
      endif

!--definition des titres :
      titx = 'Zonal. Wind Stess wstx  (cm2/s2)  '
      tity = 'Merid. Wind Stess wsty  (cm2/s2)  '

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) transfert dans les 2 tableaux uloc et vloc :
!-----------------------------------------------------------------------

      if (speval.eq.zero) then
        do j=jjrl1,jjrl2
         do i=iirl1,iirl2
           uloc(i,j) = (phisu(i,j) + cadit) * cmult
           vloc(i,j) = (phisv(i,j) + cadit) * cmult
         enddo
        enddo
      else
        do j=jjrl1,jjrl2
         do i=iirl1,iirl2
           if(kniv(i,j,-1).le.ks2) then
              uloc(i,j) = (phisu(i,j) + cadit) * cmult
              vloc(i,j) = (phisv(i,j) + cadit) * cmult
           elseif(kniv(i,j,1).le.ks2) then
              uloc(i,j) = 0.
              vloc(i,j) = 0.
           else
              uloc(i,j) = speval
              vloc(i,j) = speval
           endif
         enddo
        enddo
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  3 ) ecriture sur fichier .
!-----------------------------------------------------------------------

      open(newunit=wst_chekc_id,file='wst.check',status='unknown')

!--ecriture de l'en-tete :
      write(wst_chekc_id,1000) fmt, speval, iim,-jjm, -2, nfrc
      write(wst_chekc_id,1111) xrl1, dlong, yrl1, dlat, zrl1, dzrl, ltypv

!--ecriture de la conposante X :
      write(wst_chekc_id,'(A)') titx
      write(wst_chekc_id,*)
      do j=jjrl2,jjrl1,-1
        write(wst_chekc_id,fmt) (uloc(i,j),i=iirl1,iirl2)
      enddo
      write(wst_chekc_id,*)

!--ecriture de la conposante Y :
      write(wst_chekc_id,'(A)') tity
      write(wst_chekc_id,*)
      do j=jjrl2,jjrl1,-1
        write(wst_chekc_id,fmt) (vloc(i,j),i=iirl1,iirl2)
      enddo
      write(wst_chekc_id,*)

      close(wst_chekc_id)

      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine checkwst -
      end
