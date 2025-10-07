      subroutine gener_iceb
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!-pas_t est le pas de temps(numero de l'iteration) en cours
!-place_libre est le numero de la premiere place libre
!-fin tavleau est la position max du tableau des iceberg (apres il n'y a que du vide)
!-
!-la variable nd_cotes= la distance minimale le long de la cotes lorsqu'
!-on lache un berg
!-
!- modif 05/06/00
!- modif 24/10/19: cleaned use declaration sections              

!! START_OF_USE_SECTION
      use global_constants_mod, only: dblp=>dp, ip           
      use const_mod, only: yeaday
      use para0_mod, only:
      use para_mod, only:
      use bloc0_mod, only: tpstot, ks2, zw
      use bloc_mod, only:
      use datadc_mod, only:
      use iceberg_mod, only: dticeb, it_icb, nbricb_plus,
     >    nbricb_plus_n, nbricb_plus_s, nbrjour, nsource,
     >    numiticb, t_lach, lmx, xn, yn, kiceb, Debk, berg_mk,
     >    pond_icb, pond, srcx, srcy, hiceb, srch, srcw, viceb,
     >    wiceb, uiceb, flag_write_iceb, nbtot_write_iceb
      use reper_mod, only: xwi1, ywj1, dxwi, dywj
      
!! END_OF_USE_SECTION

      implicit none

!--variables locales :
      character(len=4) :: fchnum

!PB variables added after imposing implicit none
      integer(kind=ip) :: min_prod,max_prod,j,kflag,l,ki,kj,k
      real(kind=dblp) :: alpha,alpha2,beta,beta2, deltat,dt_lach
      real(kind=dblp) :: dt_lach_old,yrsec,prodk,year_mk,prod_1,teta
      real(kind=dblp) :: prod_2,prod_3,p,equ1

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----| 
!--instructions "data" :
! [SCRPTCOM] #include "datadc.com"

!-definition des parametres du cycle saisonnier de production
      alpha=0.
      alpha2=0.
      beta=1.
      beta2=1.
      min_prod=30
      max_prod=200

      deltat = dticeb
!     deltat = dts(ks2)

!-initialisations
      if(it_icb.ne.1) then
      nbricb_plus=0
      nbricb_plus_n=0
      nbricb_plus_s=0
      endif
      dt_lach=0.
      dt_lach_old=0.

!-main
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!-calcul du numero du jour:
       yrsec=yeaday*86400
       nbrjour=(mod(tpstot,yrsec)/86400)
       if (nbrjour.eq.0) nbrjour=360
!      if (nbrjour.eq.0) nbrjour=365

      do 400 j=1,nsource
! mab: nsource read in in deficeb.f from the restart file/source.param file
      kflag=0
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!-1.generation
!--------------
        if (numiticb.eq.int(t_lach(j))) then

!-a.recherche de place vide dans le tableau
        do 200 l=1,lmx
!-on ne remplace que si kiceb=-2

          if ((xn(l).eq.999).AND.(kiceb(l).gt.(-1))) kiceb(l)=(-1)
          if (kiceb(l).eq.(-2)) then
               kflag=l
               goto 201
          endif
 200    continue
 201    continue


! mab: LMX IS MODIFIED

        if (kflag.eq.0) then
!-----pas d'espace vide dans le tableau on met le berg a la fin
          kflag=l
          lmx=lmx+1
!- nbricb=nbre de iceb reels
!-linit est le nombre initial de bergs
!-lmx est le numero du dernier iceb dans le tableau, il est augmente si
!-generation et que le tableau est plein
!-il ne l'est pas si un iceberg est remplace a l'interieur du tableau
        endif
! JONO_icb_dbg limiting clio3.out for comprehensive debugging
! JONO_icb_dbg          write(6,*)'generation a l''iter',it_icb,j
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!-b.generation de l'iceberg (si le nombre max n'est pas atteint)

!---limitation du nombre d'iceberg a 1/6 pour arctic et 5/6 pour antar
!         if (((j.lt.36).and.(mod(kflag,6).eq.0)).or.((j.gt.36)
!     &      .and.(mod(kflag,6).ne.0))) then

!        if (kflag.lt.(icbmax-10)) then

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!-c.calcul de l'iteration du prochain lache

!-rem: calculs pour un pas de temps de 1 jour mais rectification a la
!-fin->dt peut etre qcq.


            dt_lach_old=(t_lach(j)-int(t_lach(j)))*(deltat/86400)

!--cycle saisonier:min_prod/max_prod=jour du min/max de production
!----------------- alpha/alpha2 = pente descendante/montante
!                  beta/beta2= termes independants

!--calcul de l'integrale entre min_prod et max_prod:

           prodk=Debk(j)*(0.5*alpha2*(max_prod*max_prod-min_prod*
     &               min_prod)+beta2*(max_prod-min_prod))
!          write (6,'(A,1P1E10.2)')'prodk', prodk
!-l'int de max_prod a min_prod est debk*yeaday-prodk
!-si jamais la source relache moins d'un iceberg par an:
!-year_mk est la masse de glace qui reste a relacher ds l'annee du lache
           year_mk=berg_mk(j)
           if (berg_mk(j).gt.(Debk(j)*yeaday)) year_mk=mod(berg_mk(j)
     &        ,(Debk(j)*yeaday))
!           write (6,'(a,I4,1P2E10.2)') 'masses',j,year_mk,berg_mk(j)
!           write (6,*)'jour',nbrjour

           if (nbrjour.le.min_prod) then
!---------------------------------------

             prod_1=Debk(j)*(0.5*alpha*(min_prod*min_prod-nbrjour*
     &              nbrjour)+ beta*(min_prod-nbrjour))


!             write (6,'(A,1P1E10.2)')'prod_1',prod_1
                if (prod_1.eq.year_mk)  teta=min_prod

                if (prod_1.gt.year_mk) then
                  teta=equ1(alpha,beta,nbrjour,year_mk,Debk(j))
                endif

                if (prod_1.lt.year_mk) then
                  if ((prod_1+prodk).eq.year_mk) teta=max_prod
                  if ((prod_1+prodk).gt.year_mk) then
!                   write (6,'(a,1P2E10.2,a,1P1E10.2)')'mk_new-prod_1'
!    &                      ,year_mk-prod_1, Debk(j)
                    teta=equ1(alpha2,beta2,min_prod,(year_mk-prod_1)
     &                    ,Debk(j))
                  endif
                  if ((prod_1+prodk).lt.year_mk) then
                    teta=equ1(alpha,beta,max_prod
     &                   ,(year_mk-prodk-prod_1),Debk(j))
                  endif
               endif
               dt_lach=dt_lach_old+(teta-nbrjour)+yeaday*
     &                   (int(berg_mk(j)/((alpha2*nbrjour+beta2)
     &                   *Debk(j)*yeaday)))
            endif
!           write(6,'(A,I6,1P1E10.2)') 'teta',j,teta
!           write(6,'(A,1P1E10.2)') 'nbre ans',int(berg_mk(j)
!    &                   /((alpha2*nbrjour+beta2)*Debk(j)*yeaday))

            if ((nbrjour.gt.min_prod).and.(nbrjour.le.max_prod)) then
!--------------------------------------------------------------------
!
               prod_2=Debk(j)*(0.5*alpha2*(max_prod*max_prod-nbrjour*
     &              nbrjour)+ beta2*(max_prod-nbrjour))
               if (prod_2.eq.year_mk) teta=max_prod

               if (prod_2.gt.year_mk) then
                 teta=equ1(alpha2,beta2,nbrjour,year_mk,Debk(j))
               endif

               if (prod_2.lt.year_mk) then
                 if ((prod_2+(Debk(j)*yeaday-prodk)).eq.year_mk) then
                   teta=min_prod+yeaday
                 endif
                 if ((prod_2+(Debk(j)*yeaday-prodk)).gt.year_mk) then
                   teta=equ1(alpha,beta,max_prod,(year_mk-prod_2)
     &                   ,Debk(j))
                 endif
                 if ((prod_2+(Debk(j)*yeaday-prodk)).lt.year_mk) then
                    teta=equ1(alpha2,beta2,min_prod,
     &                  (year_mk-prod_2-(Debk(j)*yeaday-prodk)),Debk(j))
                 endif
              endif
              dt_lach=dt_lach_old+(teta-nbrjour)+yeaday*
     &         (int(berg_mk(j) /((alpha*nbrjour+beta)*Debk(j)*yeaday)))
           endif

            if (nbrjour.gt.max_prod) then
!----------------------------------------

               prod_3=Debk(j)*(0.5*alpha*((min_prod+yeaday)*
     &                (min_prod+yeaday)-nbrjour*nbrjour)
     &                 +beta*((min_prod+yeaday)-nbrjour))

               if (prod_3.eq.year_mk) teta=min_prod+yeaday

               if (prod_3.gt.year_mk) then
                  teta=equ1(alpha,beta,nbrjour,year_mk,Debk(j))
               endif

               if (prod_3.lt.year_mk) then
                 if ((prod_3+prodk).eq.year_mk) teta=max_prod+yeaday
                 if ((prod_3+prodk).gt.year_mk) then
                   teta=equ1(alpha2,beta2,nint(min_prod+yeaday)
     &                  ,(year_mk-prod_3),Debk(j))
                 endif
                 if ((prod_3+prodk).lt.year_mk) then
                   teta=equ1(alpha,beta,nint(max_prod+yeaday),
     &                  (year_mk-prod_3-prodk),Debk(j))
                 endif
              endif
              dt_lach=dt_lach_old+(teta-nbrjour)+yeaday*(int
     &             (berg_mk(j)/((alpha*nbrjour+beta)*Debk(j)*yeaday)))
              if (nbrjour.gt.yeaday) nbrjour = nbrjour-yeaday
           endif
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!-d.calcul de la ponderation

           if (dt_lach.lt.1) then
             if ((nbrjour.gt.min_prod).and.(nbrjour.lt.max_prod)) then
                p=int((alpha2*nbrjour+beta2)*Debk(j)/berg_mk(j))
                dt_lach = (p+1)*(berg_mk(j)/(alpha2*nbrjour+beta2)
     &                      *Debk(j))+dt_lach_old
             else
              p=int((alpha*nbrjour+beta)*Debk(j)/berg_mk(j))
              dt_lach = (p+1)*(berg_mk(j)/(alpha*nbrjour+beta)
     &                         *Debk(j))+dt_lach_old
             endif
              pond_icb(kflag)=(p+1.)*pond(j)
           else
              pond_icb(kflag)=1.*pond(j)
           endif
           dt_lach=dt_lach*(86400/deltat)
           t_lach(j) =dt_lach+int(t_lach(j))

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!-e.l'iceberg est genere


            xn(kflag)=srcx(j)
            yn(kflag)=srcy(j)
            hiceb(kflag)=srch(j)
            wiceb(kflag)=srcw(j)
            uiceb(kflag)=0.
            viceb(kflag)=0.
            kiceb(kflag)=0.
            nbricb_plus= nbricb_plus+(1*pond_icb(kflag))
            if(srcy(j).ge.0)nbricb_plus_n=nbricb_plus_n
     &                   +(1*pond_icb(kflag))
            if(srcy(j).lt.0)nbricb_plus_s=nbricb_plus_s
     &                   +(1*pond_icb(kflag))
!           write (6,*) 'prod:',kflag,xn(kflag),yn(kflag),hiceb(kflag),
!    &         wiceb(kflag)
!           if (pond_icb(kflag).ne.1) write(6,*)'ponderation source',j,
!    &         '=',pond_icb(kflag)

            ki=1+nint((xn(kflag)-xwi1)/dxwi)
            kj=1+nint((yn(kflag)-ywj1)/dywj)
! JONO_icb_dbg , commmenting to prevent lots of output
! JONO_icb_dbg           if (tms(ki,kj,ks2).eq.0) write(6,*)'source terre',j,'prod:'
! JONO_icb_dbg    &      ,kflag,xn(kflag),yn(kflag),ki,kj,xwi1,ywj1

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!-f. Def. du niveau dans lequel se trouve l'iceberg

            kiceb(kflag)=ks2

            do 220 k=ks2,2,-1
              if (-zw(kiceb(kflag)).le.hiceb(kflag)) then
               kiceb(kflag)=k-1

              else
               goto 300
              endif
 220        continue
 300        continue

            if (kiceb(kflag).le.-1) kiceb(kflag)=ks2

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!-2.ecriture des productions sur icbtraj.out
!--------------------------------------------
! JONO_dbg debugging writing to fort.54:
! commenting similar lines in deficeb.f (accesed just once)
! renaming and closing 1054 from icbtraj.out (unformatted) to icbtraj.header (formatted)
! naming, opening and closing this file HERE (was54) to 1054=icbtraj.out
! JONO 14-6 added format in write line
           if (flag_write_iceb) then
              open(1054,file='icbtraj.out',status='unknown')
              write(1054, '(2I7,4F11.3)') kflag,numiticb,real(xn(kflag))
     &                  ,real(yn(kflag)),
     &                  real(hiceb(kflag)),real(wiceb(kflag))
              nbtot_write_iceb = nbtot_write_iceb+1
              close(1054)
            endif

!-3.test pour la cote
!---------------------

!            ki=1+nint((xn(kflag)-xwi1)/dxwi)
!            kj=1+nint((yn(kflag)-ywj1)/dywj)
!            call dcotes(kflag,ki,kj,xn(kflag),yn(kflag))

           endif

 400  continue

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!     write (6,*) 'fin generation'
!      do 500 l=1,lmx

!        if (xn(l).eq.0) write (*,'(A5,I4,4(F7.3),2I4)') 'null',l,
!     &         xn(l),yn(l),hiceb(l),wiceb(l),lmx,nbricb

! 500  continue

      write(99,'(A,2I6)') 'it_icb nbre icb en plus', it_icb,
     &  nbricb_plus

      return
! fin de routine gener_iceb
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      end subroutine gener_iceb

      function equ1 (a,b,iday_1,amk,Debk) result(equ_un)
!--resout l'eq du second degre si l'inconnue est la borne superieure
      use global_constants_mod, only: dblp=>dp, ip           

c~ !-include
c~ #include "type.com"
      real(kind=dblp) :: equ_un, a, b, amk, Debk, c, day_2, det
      integer(kind=ip):: iday_1
!-l'equation a resoudre est 0.5*a*(day_2**2-iday_1**2)+b*(day_2-iday_1)=amk/Debk

      c=(-0.5*a*iday_1*iday_1)-b*iday_1-(amk/Debk)
      det=b*b+(2*a*c)
      if (a.eq.0.)then
         day_2=((amk/Debk/b)+iday_1)
      else
         day_2=(-b+sqrt(det))/(2*a)
      endif
      equ_un=day_2

      end function equ1
