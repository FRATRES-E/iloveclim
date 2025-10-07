!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Sous-programme de calcul de correspondances lat,lon ECbilt sur
!       la grille ISM consideree 
!      Modifications pour utilisation "offline" a partir des variables
!       de LOVECLIM vers GRISLI Nord
!
!      Auteur : Didier M. Roche
!      Date   : Juillet 2008
!      Derniere modification : 13 Janvier 09, Didier M. Roche, C. Dumas
!                              20 feb 17, A. Quiquet (latitudes...)
!-----|--1--------2---------3---------4---------5---------6---------7-|

      SUBROUTINE calcijEcbg(nx,ny,latit,longit,iEcb,jEcb)

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree : 
!       nx, ny     = dimension de la grille de l'ISM vers lequel on va
!       latit      = tableaux des latitudes  des cases ISM
!       longit     = tableaux des longitudes des cases ISM
!       Variables de sorties : 
!       iEcb, jEcb = tableaux des longitudes, latitudes des cases ECBilt
!                         sur la grille de l'ISM considere
!
!         --  pour GRISLI Nord :       nx, ny = 241, 241
!-----|--1--------2---------3---------4---------5---------6---------7-|

        IMPLICIT NONE
  
        INTEGER :: i,j,l

        integer nx,ny
        integer :: longitiEcbtxt_id, latitiEcbtxt_id
        real latit(nx,ny),longit(nx,ny),iEcb(nx,ny),jEcb(nx,ny)

!-- afq -- the grid spacing is irregular in latitude...:
        real*8 latEcb(32)  !afq latitude of ECB cell centre
        real*8 bordEcb(32) !afq upper edge of ECB cell
        real*8 dlEcb(32)   !afq length of ECB cell

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      On effectue une repartition simple sur la grille ECBilt
!-----|--1--------2---------3---------4---------5---------6---------7-|

        latEcb(1:32) = (/                                        &
        -85.7606, -80.2688, -74.7445, -69.2130, -63.6786,        & 
        -58.1430, -52.6065, -47.0696, -41.5325, -35.9951,        & 
        -30.4576, -24.9199, -19.3822, -13.8445, -8.30670,        & 
        -2.76890, 2.76890, 8.30670, 13.8445, 19.3822,            & 
        24.9199, 30.4576, 35.9951, 41.5325, 47.0696,             & 
        52.6065, 58.1430, 63.6786, 69.2130, 74.7445,             & 
        80.2688, 85.7606  /)  
        
        do l=1,31
           bordEcb(l)=latEcb(l)+(latEcb(l+1)-latEcb(l))/2.
        enddo
        bordEcb(32)=90.

        dlEcb(1)=bordEcb(1)+90.
        do l=2,32
           dlEcb(l)=bordEcb(l)-bordEcb(l-1)
        enddo

     OPEN(newunit=longitiEcbtxt_id,file="outputdata/coupler/longitiEcb.txt",form="formatted")
     OPEN(newunit=latitiEcbtxt_id,file="outputdata/coupler/latitjEcb.txt",form="formatted")


        do i=1,nx
          do j=1,ny
            !jEcb(i,j)=2+(latit(i,j)+84.375)/5.625   ! afq -- this is wrong!
            iEcb(i,j)=2+(longit(i,j)-5.625/2)/5.625

            jEcb(i,j)=dble(2)+(latit(i,j)-bordEcb(1))/dlEcb(1)
            do l=2,32
               if (((latit(i,j)-bordEcb(l)).le.0.).and.(latit(i,j).gt.bordEcb(l-1))) then
                  jEcb(i,j)=dble(l+1)+(latit(i,j)-bordEcb(l))/dlEcb(l)
               endif
            end do
          end do
        end do  

        where(iEcb(:,:).ge.65.) iEcb(:,:)=iEcb(:,:)-64.

        do i=1,nx
          do j=1,ny
            WRITE(longitiEcbtxt_id,*) iEcb(i,j)
            WRITE(latitiEcbtxt_id,*) jEcb(i,j)
          end do
        end do

        CLOSE(longitiEcbtxt_id)
        CLOSE(latitiEcbtxt_id)
      END SUBROUTINE calcijEcbg
