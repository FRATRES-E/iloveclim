!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:44 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:44 CET 2009

      SUBROUTINE debugrac(tabx,taby,tabz,name)
!- ppmodif: 25-02-97

!! START_OF_USE_SECTION

      use global_constants_mod, only: dblp=>dp, ip

      use const_mod

      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod
      use newunit_clio_mod, only: clio3_out_id
      
!! END_OF_USE_SECTION

         IMPLICIT NONE

         real(kind=dblp), dimension(imax,jmax,kmax) :: tabx, taby, tabz
         integer(ip) :: nn
         
         CHARACTER*50 name

         integer(kind=ip) :: nncf, n, j
         integer(ip) :: name_id

        DO nn=1,len(name)
         if (name(nn:nn).ne.' ') nncf=nn
        ENDDO

        OPEN(newunit=name_id,file=name(:nncf),status='UNKNOWN')
        write(name_id,*) '****************************  '
        write(name_id,*)  name(:nncf)
        write(clio3_out_id,*)  '** DEBUG RACCORD for '//name(:nncf)//' **'
        write(name_id,*) '****************************  '
        write(name_id,*) 'RACCORD BERING Tabx,Taby,Tabz'
        write(name_id,*) 'Pacific A,B    Atlantic E,F'
        write(name_id,*) '        C,D             G,H'
        write(name_id,*) 'Raccord SCALAR, E<=D,F<=C,A<=H,B<=G'
        write(name_id,*) 'Raccord VECTOR_Y, E<=-B,F<=-A'
        write(name_id,*) '****************************  '

        write(name_id,*)' A',tabx(iberp-1,jberp,kmax),
     &            ' A',taby(iberp-1,jberp,kmax),
     &            ' A',tabz(iberp-1,jberp,kmax)

        write(name_id,*)' B',tabx(iberp,jberp,kmax),
     &            ' B',taby(iberp,jberp,kmax),
     &            ' B',tabz(iberp,jberp,kmax)

        write(name_id,*)' C',tabx(iberp-1,jberp-1,kmax),
     &            ' C',taby(iberp-1,jberp-1,kmax),
     &            ' C',tabz(iberp-1,jberp-1,kmax)

        write(name_id,*)' D',tabx(iberp,jberp-1,kmax),
     &            ' D',taby(iberp,jberp-1,kmax),
     &            ' D',tabz(iberp,jberp-1,kmax)

        write(name_id,*)' E',tabx(ibera-1,jbera,kmax),
     &            ' E',taby(ibera-1,jbera,kmax),
     &            ' E',tabz(ibera-1,jbera,kmax)

        write(name_id,*)' F',tabx(ibera,jbera,kmax),
     &            ' F',taby(ibera,jbera,kmax),
     &            ' F',tabz(ibera,jbera,kmax)

        write(name_id,*)' G',tabx(ibera-1,jbera-1,kmax),
     &            ' G',taby(ibera-1,jbera-1,kmax),
     &            ' G',tabz(ibera-1,jbera-1,kmax)

        write(name_id,*)' H',tabx(ibera,jbera-1,kmax),
     &            ' H',taby(ibera,jbera-1,kmax),
     &            ' H',tabz(ibera,jbera-1,kmax)

        write(name_id,*) '****************************  '
        write(name_id,*) 'RACCORD CYCLIC Tabx '
        write(name_id,*) '****************************  '
        write(name_id,*) ' TABx(is1)  TABx(is2+1)   '
        DO j=jcl1,jcl2
         write(name_id,*) j,tabx(is1(j),j,kmax),
     &                tabx(is2(j)+1,j,kmax)
        ENDDO
        write(name_id,*) ' TABx(is1-1)  TABx(is2)   '
        DO j=jcl1,jcl2
         write(name_id,*) j,tabx(is1(j)-1,j,kmax),
     &                tabx(is2(j),j,kmax)
        ENDDO
        CLOSE(name_id)
        RETURN
        END
