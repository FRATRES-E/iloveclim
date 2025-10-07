!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:43 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:43 CET 2009

#if ( BATHY == 3 )

      subroutine change_grid(irunlabelclio)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  modif : 2019
!  base sur init_clio
!  pour modifier grille clio


      use comemic_mod, only: fracto, dareafac
      use comcoup_mod, only: indo2a, wo2a, inda2o, wa2o
     >         , ijatm
     >         , kamax, komax, ijocn
      use comsurf_mod, only: fractoc, fractn
     >         , nld, ntyps
      USE comatm, ONLY: nlat, nlon, darea,dareas

#if ( DOWNSTS == 1 || DOWNSCALING == 2 || ISM >= 2 )
      use sgout_mass_balance_mod, only: initakkuVars
      use transfer_subgrid_ecb_mod, only: subgrid_ecb_wrapper
#endif     

#if ( ISM == 2 || ISM == 3 || F_PALAEO >= 1 || DOWNSCALING == 2 )
      use ecbilt_topography, only: ec_topo
#endif
 

      use const_mod, only: zero, one, cstmin, cstmax

#if ( PATH >= 1 )
      use path_mod, only: path_init, scalstart_path, scalend_path
#endif

      use para0_mod, only: ijkmax, ltest, nsmax_TS

      use bloc0_mod, only: bering, iberp, tms, tmu, scal, fss, eta
     >                     , ub, vb, fqajc, b, u, v, q2turb, w, bvf
     >                     , avsdz, avudz, phifu, phifv, phifs, phiss
     >                     , rappes, phimnx, rappel, scalr, phivs, scs
     >                     , scpme, ks2

      use bloc_mod, only: q2tmin, umoy, vmoy, phihhh, fub, fvb, q
     >                     , phizzz, unsvol

      use ice_mod, only: hgbq, fsbbq, qstobq, albq, hnbq, ts, tbq, xzo
     >                     , tenagx, tenagy
     >                     , imax, jmax, kmax

      use dynami_mod, only: sber, ug, vg

      use update_clio_bathy_tools, only: get_clio_masks_fnm, read_clio_masks
     >         , update_diff_masks, update_prev
     >         , la_date
     >         ,  mean_neighbours_with_mask, mise_a_zero
     >         , mise_a_valeur
     >         ,bering_nb, bering_date, bering_value

      use global_constants_mod, only: str_len, dblp=>dp, ip

      use unix_like_libs, only: is_file

      use vareq_mod,  only: uslpfx, vslpfx

      use moment_mod, only: vicmom, copy_from_vicmon, copy_to_vicmon 

      use comland_mod, only: albland,fractl,bmoisg,bmoismfix,dsnow
     >              ,runofo, tland,dsnm, forestfr, epsl, albsnow
     >              , tzero

#if ( ISM >= 2 )
      USE comland_mod, only: albkeep
#endif

      USE comatm, ONLY: nlat, nlon

      use comland_mod, only: albland,fractl,bmoisg,bmoismfix,dsnow
     >              ,runofo, tland,dsnm

#if ( ROUTEAU == 1 )
      USE routageEAU_mod, only: epsi_lon, mask_lnd
#endif

#if (IMSK == 1)
      USE input_icemask, only: icemask
#endif

      use comunit, only: iuo
      use global_constants_mod, only: ip


      implicit none

!- nn99=2 => writing in the file "mouchard", unit=99
!      common / mchd99 / nn99
!      common /clio_control_int/ ktvar,ntrmax,mixage
!      common /clio_control_int/ mixage
!      common /clio_control_fp/ dtsd2,yrsec,daysec,unsplt


      integer(ip) :: irunlabelclio

      integer(ip) ::   ija, i, j, k, nb_var
      real(dblp), dimension (ijatm)  ::    fractocn

!      dimension irn(imax,8), jrn(jmax,8)
!      character*4 fchnum
!      character*6 chf
      character*30 name_file
      character(len=5) :: charI


      integer(ip) :: nn

      character(len=str_len) :: file_path
      logical                :: logical_elmt

      integer(ip) :: kmaxb

      real(dblp) ::  distance_before, distance_after
      real(dblp), dimension(nlat,nlon) :: fractl_prev

      integer(ip),  parameter :: nszmax = ijkmax
!- variables locales equivalentes :
       real(dblp), dimension(nszmax) :: vloc !dmr on dirait que vloc est intent(out) ...

      real(dblp), dimension(ijatm) :: fractocn_prev
      real(dblp), dimension(nlat,nlon) :: fracto_prev
      integer(ip) :: nv, ns, is, ijo
      logical :: logic_elmt

      integer(ip) :: fractoc_dat_id, mozaic_w_id

!1. Fractoc
!!!!!!!!!!!!!!!

!modify fractoc
!--------------
! a l initialisation c dans emic.f (ec_initemic)


      write(charI,'(I5.5)'), la_date

      name_file ='inputdata/fractoc_'//trim(charI)//'.dat'

      open(newunit=fractoc_dat_id,file = name_file,status='old',
     &        form='formatted')

      fractocn_prev(:) = fractocn(:)
      fracto_prev(:,:) = fracto(:,:)
      fractl_prev(:,:)=fractl(:,:)

      read(fractoc_dat_id,*)
      read(fractoc_dat_id,*) (fractocn(ija),ija=1,ijatm)
      rewind(fractoc_dat_id)
      do ija=1,ijatm
         j=int((ija-1)/nlon)+1
         i=ija-(j-1)*nlon
         fracto(j,i)=fractocn(ija)
         if (fracto(j,i).gt.0.990) fracto(j,i)=1.0d0
      enddo

      close(fractoc_dat_id)

      do j=1,nlon
        do i=1,nlat
           if (abs(fracto(i,j)-fracto_prev(i,j)).gt.0.0) then
              write(*,*) i,j,fracto(i,j),fracto_prev(i,j)
           endif
        enddo
      enddo


cnb from initial0.f
      do j=1,nlon
        do i=1,nlat

          fractoc(i,j)=fracto(i,j)
          
!dmr&nb [INUTILE] ... a priori          
c~           fractn(i,j,noc)=fractoc(i,j)
c~           fractn(i,j,nld)=1-fractoc(i,j)
          
!dmr&nb [SEAICE]
        if ((fracto(i,j).gt.0.0).and.(fracto_prev(i,j).le.0.0)) then ! new_ocean
!dmr&nb ... [KOIFAIRE???] 
c~             fractn(i,j,noc)=fractoc(i,j)
            fractn(i,j,nld)=1-fractoc(i,j)
cnb --- parts of ec_initlbm
cnb from landmodel0.f
!*** land fraction
            fractl(i,j)=fractn(i,j,nld)
c~             fractn(i,j,nse) = 0.0d0 !dmr&nb pourrait être mis à jour avec la concentration "alb" des voisins
            write(*,*) "NEW OCEAN in: ", i, j
        else if ((fracto(i,j).le.0.0).and.(fracto_prev(i,j).gt.0.0)) then ! old_ocean
c~              fractn(i,j,noc)=fractoc(i,j)
             fractn(i,j,nld)=1-fractoc(i,j)
             fractl(i,j)=fractn(i,j,nld)
c~              fractn(i,j,nse) = 0.0d0
             write(*,*) "OLD OCEAN in: ", i, j
        else if (fracto(i,j).ne.fracto_prev(i,j)) then ! fract_ocean different
             !dmr&nb cas ou changement de fraction ocean mais pas de type (land<->ocean)
c~              albq_lok = fractn(i,j,nse)/fracto_prev(i,j)
c~              fractn(i,j,nse) = albq_lok*fracto(i,j)
c~              fractn(i,j,noc) = (1.0-albq_lok)*fracto(i,j)
             fractn(i,j,nld)=1-fractoc(i,j)
             fractl(i,j)=fractn(i,j,nld)
        else
             !dmr&nb rien a priori (pas de changement de fractoc)
        endif    
!dmr&nb [SEAICE]
        
        enddo
      enddo

      do j=1,nlon
        do i=1,nlat
cc          fractl_prev=fractl(i,j)
!???          fractl(i,j)=1.-fracto(i,j)
cc          !new land cell
          
          if ((fractl(i,j).gt.epsl).and.(fractl_prev(i,j).le.epsl)) then
          
            write(*,*) 'initialize new continent cell', fractl(i,j), 
     >                fractl_prev(i,j)
            bmoisg(i,j)=bmoismfix
            dsnow(i,j)=0d0
            runofo(i,j)=0d0
            
!dmr&nb ... [KOIFAIRE???]             
            tland(i,j)=tzero+10 !+10. !dmr&nb ---plutôt prendre la moyenne des voisins !!!!!
            
c~             dsnow(i,j)=min(dsnow(i,j),dsnm)
c~             if (fractl(i,j).lt.epsl) then
c~                 bmoisg(i,j)=0d0
c~             endif

cnb --- parts of ec_landcoverupdate
ccnb albedo on land
             do is=1,4
               if (albland(i,j,is).eq.0.0) then
                 albland(i,j,is)=0.18
               endif
             enddo

          
          else if ((fractl(i,j).le.epsl).and.(fractl_prev(i,j).gt.epsl)) then
!dmr&nb ... [KOIFAIRE???] land -> ocean !!!          
          endif
          

        enddo
      enddo

c~ #if ( IMSK == 1 )
c~         do i=1,nlat
c~           do j=1,nlon
c~ #if ( ISM >= 2 )
c~             albkeep(i,j,:) = albland(i,j,:) ! afq, used for ice sheet surface mass balance
c~ #endif
c~             if (icemask(i,j).gt.0.9) then
c~                forestfr(i,j)=0.
c~                albland(i,j,:) = albsnow(i) !albland(30,58,:) ! Greenland centre is at lat 30, lon 58
c~ !              do is=1,4
c~ !               if(albland(i,j,1).le.albsnow(i)) albland(i,j,is)=albsnow(i)
c~ !              enddo
c~             endif
c~           enddo
c~         enddo
c~ #endif /* on icemask */

c~ ! *** addition for LGM runs ***

c~         d=0d0
c~         do j=2,25
c~           d=d+albland(27,j,1)
c~         enddo
c~         d=d/24.
c~         write(iuo+99,120) 'landcover update year: ',iyear,scenyr,yr,d
c~         do i=1,nlat
c~           do j=1,nlon
c~             if (fractl(i,j).gt.epsl) then
c~               do is=1,4
c~                 if (albland(i,j,is).lt.0.01.or.
c~      *            albland(i,j,is).gt.0.99) then
c~                   write(iuo+99,130)
c~      *            'Albedo of land out of range ',i,j,is,albland(i,j,is)
c~                 endif
c~               enddo
c~             endif
c~           enddo
c~         enddo

c~ 120   FORMAT(A23,3I8,F12.5)
c~ 130   FORMAT(A28,3I4,F12.5)



cnb --- fin parts of ec_landcoverupdate


!dmr&nb --- a priori pas utile ... forestfr non initialisee de toute facon
c~       call ec_landalbedo(1)


!appele dans ec_initecbilt (initial0.f):
!pas utile      call ec_inierror ! initialise erreur
!pas utile      call ec_inimdldate !initialise le temps
!pas utile      call ec_iatmpar !initialise constantes
!pas utile      call ec_iatmdyn !initialise operators
!pas utile      call ec_iatmphys
!pas utile      call ec_inioutparat

      call ec_atmstate ! useful?

!2. Bering
!!!!!!!!!!!

!Bering is initially read in defcst.f, we replace it
!nb update Bering value
       do i=2,bering_nb-1
           if ((la_date.le.bering_date(i-1))
     >        .and.(la_date.ge.bering_date(i))) then
              distance_before=bering_date(i-1)-la_date
              distance_after=la_date-bering_date(i)
              if (distance_before.le.distance_after) then
                  bering=bering_value(i-1)
              else
                  bering=bering_value(i)
              endif
           endif
       enddo
       !write(*,*) 'Bering value at year ', la_date, '  updated to ', bering

!cnb cf defcst.f
!--Bering Strait opening:
      iberp = 0
cnb le pprobleme venait d'ici
      !if ( ltest.eq.3 .and. bering.le.0. ) iberp = 1 !lt
      if ( ltest.eq.3 .and. bering.ne.0. ) iberp = 1

       sber   = 1.-max(zero,sign(one,0.1-bering))


       !write(*,*) 'IBERP 3', iberp, ltest, bering


!3. change grid: tms and tmu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

cnb update prev masks
      write(*,*) 'call update prev and diff'
      logical_elmt = update_prev(tms,tmu)

!dmr&nb [NOTA]

! nn99 is a global variable ; bad idea to change it?
! only controls some writing in mouchard, should nonetheless has no impact (?)

!dmr&nb [UNUSED??]      jflag = 0
!dmr&nb      nn99 = 0

      !call defgrid(nn99)

#if (1)
cnb instead od defgrid: parts of defgrid, to update kbath from bath_xxx.txt and modify tms and tmu
cnb note that bering is modified in this call
c~       call re_defgrid(0)
      write(*,*) 'call defgrid'
      call defgrid(3)
#endif 

!dmr --- updating the differences in tms/tmu masks
      logical_elmt = update_diff_masks(tms,tmu)
!dmr

!dmr [NOEQUI]
      call copy_to_vicmon()

!nb modifications des valeures pour nouvelles cells
!pour les autres traceurs: nsmax??

           write(*,*) 'track error scal'     
!           call track_error(scal(:,:,:,1), 270, 320) 
!           call track_error(scal(:,:,kmax-1,1), 270, 320) 
       do nv=1,nsmax_TS
cnb dans defgrid: scal(i,j,k,ns) = scal0(k,ns)
           logic_elmt = mean_neighbours_with_mask(scal(:,:,:,nv),1) ! last argument 1 = tms2D
c         do k=1,kmax
c           logic_elmt = mise_a_valeur(scal(:,:,k,nv),1,scal0(k,nv)) ! last argument 1 = tms2D
c         enddo
       enddo
c           call track_error(scal(:,:,kmax,1), 270, 320) 

      do ns=0, nsmax_TS
cnb mettre a 0 (comme dans defgrid) ?
           logic_elmt = mean_neighbours_with_mask(fss(:,:,ns),1) ! last argument 1 = tms2D
c           logic_elmt = mise_a_zero(fss(:,:,ns),1) ! last argument 1 = tms2D
      enddo

           logic_elmt = mean_neighbours_with_mask(eta(:,:),1) ! last argument 1 = tms2D
c           logic_elmt = mise_a_zero(eta(:,:),1) ! last argument 1 = tms2D

           logic_elmt = mise_a_zero(ub(:,:),2) ! last argument 1 = tms2D
!???           ub(:,:) = 0.0d0

           logic_elmt = mise_a_zero(vb(:,:),2) ! last argument 1 = tms2D
!???           vb(:,:) = 0.0d0

           logic_elmt = mise_a_zero(fqajc(:,:,:),2) ! last argument 1 = tms2D !diffusion?
c           logic_elmt = mean_neighbours_with_mask(fqajc(:,:,:),2) ! last argument 1 = tms2D

           logic_elmt = mean_neighbours_with_mask(b(:,:,:),1) ! last argument 1 = tms2D
c           logic_elmt = mise_a_zero(b(:,:,:),1) ! last argument 1 = tms2D

           logic_elmt = mise_a_zero(u(:,:,:),2) ! last argument 1 = tms2D
!???           u(:,:,:) = 0.0d0

           logic_elmt = mise_a_zero(v(:,:,:),2) ! last argument 1 = tms2D
!???           v(:,:,:) = 0.0d0

cnb dans defgrid: q2turb(i,j,k)=q2tmin
           logic_elmt = mise_a_zero(q2turb(:,:,:),1) ! last argument 1 = tms2D !!!! HUUUUUUUMMMMMM PAS SUR DU TOUT !!!!!
c           logic_elmt = mise_a_valeur(q2turb(:,:,:),1,q2tmin) ! last argument 1 = tms2D
c           logic_elmt = mean_neighbours_with_mask(q2turb(:,:,:),1) ! last argument 1 = tms2D

           logic_elmt = mise_a_zero(w(:,:,:),2) ! last argument 1 = tms2D
c~            do i=1,imax
c~              do j=1,jmax
c~                do k=1,kmax
c~                  if (k.eq.kfs(i,j))then
c~                    w(i,j,k)=0.0
c~                  endif
c~                enddo
c~              enddo
c~            enddo
!???           w(:,:,:) = 0.0d0

cnb N2 =Br.Vais.Freq
cnb mise a zero ?
           logic_elmt = mean_neighbours_with_mask(bvf(:,:,:),1) ! last argument 1 = tms2D
c           logic_elmt = mise_a_zero(bvf(:,:,:),1) ! last argument 1 = tms2D

cnb coefficients de DIFFUsion Verticale (divise par dz) : avsdz, avudz
cnb mise a zero ?
           logic_elmt = mean_neighbours_with_mask(avsdz(:,:,:),1) ! last argument 1 = tms2D
c           logic_elmt = mise_a_zero(avsdz(:,:,:),1) ! last argument 1 = tms2D

cnb mise a zero ?
           logic_elmt = mean_neighbours_with_mask(avudz(:,:,:),2) ! last argument 1 = tms2D
c           logic_elmt = mise_a_zero(avudz(:,:,:),2) ! last argument 1 = tms2D

cnb pour le downsloping
           logic_elmt = mise_a_zero(uslpfx(:,:,:),2) ! last argument 1 = tms2D

           logic_elmt = mise_a_zero(vslpfx(:,:,:),2) ! last argument 1 = tms2D

           logic_elmt = mise_a_zero(umoy(:,:),2) ! last argument 1 = tms2D
!???           umoy(:,:) = 0.0d0

           logic_elmt = mise_a_zero(vmoy(:,:),2) ! last argument 1 = tms2D
!???           vmoy(:,:) = 0.0d0

cnb ice thickness
           logic_elmt = mean_neighbours_with_mask(hgbq(:,:),1) ! last argument 1 = tms2D
c           logic_elmt = mise_a_zero(hgbq(:,:),1) ! last argument 1 = tms2D

           logic_elmt = mise_a_zero(fsbbq(:,:),1) ! last argument 1 = tms2D
c           logic_elmt = mean_neighbours_with_mask(fsbbq(:,:),1) ! last argument 1 = tms2D

cnb energy stored in brine pocket
           logic_elmt = mise_a_zero(qstobq(:,:),1) ! last argument 1 = tms2D
c           logic_elmt = mean_neighbours_with_mask(qstobq(:,:),1) ! last argument 1 = tms2D

cnb lead fraction
           logic_elmt = mean_neighbours_with_mask(albq(:,:), 1) ! last argument 1 = tms2D
c           logic_elmt = mise_a_zero(albq(:,:),1) ! last argument 1 = tms2D

cnb snow thickness
           logic_elmt = mean_neighbours_with_mask(hnbq(:,:),1) ! last argument 1 = tms2D
c           logic_elmt = mise_a_zero(hnbq(:,:),1) ! last argument 1 = tms2D

           logic_elmt = mean_neighbours_with_mask(ts(:,:),1) ! last argument 1 = tms2D
c           logic_elmt = mise_a_zero(ts(:,:),1) ! last argument 1 = tms2D

           logic_elmt = mise_a_zero(ug(:,:),2) ! last argument 1 = tms2D
!???           ug(:,:) = 0.0d0

           logic_elmt = mise_a_zero(vg(:,:),2) ! last argument 1 = tms2D
!???           vg(:,:) = 0.0d0

cnb Temperature inside the ice/snow layer
           logic_elmt = mean_neighbours_with_mask(tbq(:,:,:),1) ! last argument 1 = tms2D
c           logic_elmt = mise_a_zero(tbq(:,:,:),1) ! last argument 1 = tms2D

cnb  rugosity of the ice
           logic_elmt = mean_neighbours_with_mask(xzo(:,:),1) ! last argument 1 = tms2D
c           logic_elmt = mise_a_zero(xzo(:,:),1) ! last argument 1 = tms2D

cnb Wind stress at the ice surface (x)
           logic_elmt = mean_neighbours_with_mask(tenagx(:,:),1) ! last argument 1 = tms2D
c           logic_elmt = mise_a_zero(tenagx(:,:),1) ! last argument 1 = tms2D

           logic_elmt = mean_neighbours_with_mask(tenagy(:,:),1) ! last argument 1 = tms2D
c           logic_elmt = mise_a_zero(tenagy(:,:),1) ! last argument 1 = tms2D

c           logic_elmt = mean_neighbours_with_mask(vicmom(:,:,:),1) ! last argument 1 = tms2D
c           logic_elmt = mise_a_zero(vicmom(:,:,:),1) ! last argument 1 = tms2D
           do nb_var=1,UBOUND(vicmom,DIM=3) 
             logic_elmt = mean_neighbours_with_mask(vicmom(:,:,nb_var),1) ! last argument 1 = tms2D
           enddo

cnb added (and corresponding parts removed from defgrid dans le cas nflag 3)
        !phifu(i,j)  = 0.0
        !phifv(i,j)  = 0.0
         logic_elmt = mise_a_zero(phifu(:,:),2) ! last argument 1 =tms2D
         logic_elmt = mise_a_zero(phifv(:,:),2) ! last argument 1 =tms2D
         !phifs(i,j,ns)  = 0.0
         !phiss(i,j,ns)  = 0.0
         !rappes(i,j,ns) = 0.0
         !phimnx(i,j,0,ns) = cstmin
         !phimnx(i,j,1,ns) = cstmax
      do ns=1, nsmax_TS
         logic_elmt = mise_a_zero(phifs(:,:,ns),1) ! last argument 1 = tms2D (=w)
         logic_elmt = mise_a_zero(phiss(:,:,ns),1) ! last argument 1 = tms2D
c         logic_elmt = mean_neighbours_with_mask(phiss(:,:,ns),1) ! last argument 1 = tms2D
         logic_elmt = mise_a_zero(rappes(:,:,ns),1) ! last argument 1 = tms2D
c         logic_elmt = mean_neighbours_with_mask(phimnx(:,:,0,ns),1) ! last argument 1 = tms2D
c         logic_elmt = mean_neighbours_with_mask(phimnx(:,:,1,ns),1) ! last argument 1 = tms2D
         logic_elmt = mise_a_valeur(phimnx(:,:,0,ns),1,cstmin) ! last argument 1 = tms2D
         logic_elmt = mise_a_valeur(phimnx(:,:,1,ns),1,cstmax) ! last argument 1 = tms2D
      enddo
      do nn=1,6
cnb adv diff
        !phihhh(i,j,nn) = 0.0
        logic_elmt = mise_a_zero(phihhh(:,:,nn),1) ! last argument 1 = tms2D
      enddo
         !fub(i,j,k) = 0.0
         !fvb(i,j,k) = 0.0
         !q(i,j,k) = 0.0
         !rappel(i,j,k) = 0.0
        logic_elmt = mise_a_zero(fub(:,:,:),2) ! last argument 1 =tms2D
        logic_elmt = mise_a_zero(fvb(:,:,:),2) ! last argument 1 =tms2D
        logic_elmt = mise_a_zero(q(:,:,:),2) ! last argument 1 =tms2D
        logic_elmt = mise_a_zero(rappel(:,:,:),1) ! last argument 1 =tms2D
         !phizzz(i,j,k,1) = 0.0
         !phizzz(i,j,k,2) = 0.0
        logic_elmt = mise_a_zero(phizzz(:,:,:,1),2) ! last argument 1=tms2D
        logic_elmt = mise_a_zero(phizzz(:,:,:,2),2) ! last argument 1=tms2D
      do ns=1,nsmax_TS
cnb ???        !nrap(k,ns) = 0
        !logic_elmt = mise_a_zero(nrap(:,ns),1) ! last argument 1=tms2D
          !scalr(i,j,k,ns)=0.0
          logic_elmt = mise_a_zero(scalr(:,:,:,ns),1) ! last argument 1=tms2D
c          logic_elmt = mean_neighbours_with_mask(scalr(:,:,:,ns),1) ! last argument 1 = tms2D
          !phivs(i,j,k,ns) = 0.0
          logic_elmt = mise_a_zero(phivs(:,:,:,ns),1) ! last argument 1=tms2D
c          logic_elmt = mean_neighbours_with_mask(phivs(:,:,:,ns),1) ! last argument 1 = tms2D
      enddo
      do ns=1,nsmax_TS
!         scs(i,j,ns) = scpme(ns)
c         logic_elmt = mean_neighbours_with_mask(scs(:,:,ns),1) ! last argument 1 = tms2D
         logic_elmt = mise_a_valeur(scs(:,:,ns),1,scpme(ns)) ! last argument 1 = tms2D
      enddo

!dmr [NOEQUI]
           call copy_from_vicmon()
!dmr [NOEQUI]

c           call track_error(scal, min_v, max_v) 

c~         vloc(1) = 0.

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  Traitement des cas d'un changement de stockage / bathymetrie    |
!-----------------------------------------------------------------------

!     if(nnt.lt.0) then
!       do 710 j=1,jmax
!        do 710 i=1,imax
!         ub(i,j) = ub(i,j) * hu(i,j)
!         vb(i,j) = vb(i,j) * hu(i,j)
!710    continue
!     endif

!--Tester si la bathymetrie a ete lue :
      if (unsvol.gt.zero) then
!--precaution : vitesse nulle en dehors du domaine .
        do 730 k=1,kmax
         do 730 j=1,jmax
          do 730 i=1,imax
            u(i,j,k) = u(i,j,k) * tmu(i,j,k)
            v(i,j,k) = v(i,j,k) * tmu(i,j,k)
 730    continue
        do 740 j=1,jmax
         do 740 i=1,imax
           ub(i,j) = ub(i,j) * tmu(i,j,ks2)
           vb(i,j) = vb(i,j) * tmu(i,j,ks2)
 740    continue
      endif


      call etat(0, 0)
      call informe(0) !nb a modifier !!!


!dmr&nb [ICICESTBON]


!suite ec_initlbm

#if ( ROUTEAU == 0 )
      call ec_inirunoff
#else

       WHERE(fractl.GT.epsi_lon)
         mask_lnd = 1
       ELSEWHERE
         mask_lnd = 0
       ENDWHERE

      call ec_inirunoff_2
#endif


cnb --- end of parts of ec_initlbm


      !call ec_initcoup

!4. Mozaic
!!!!!!!!!!

cnb --- parts of ec_initcoup

! mozaic.w
!---------
!       kamaxb=kamax

      !if ((la_date.le.20000).and.(la_date.ge.10000)) then
      !    kamaxb=kamax2 !kamaxb=13
      !else
      !    kamaxb=kamax !kamaxb=14
      !endif
      !write(*,*) 'kamax effectif ', kamaxb


      name_file ='inputdata/clio/mozaic_'//trim(charI)//'.w'

      open(newunit=mozaic_w_id,file = name_file,status='old',
     &        form='unformatted')

!      if (la_date.ge.21000) then

      read(mozaic_w_id)
      read(mozaic_w_id) ((indo2a(ija,k),k=1,kamax),ija=1,ijatm)
      read(mozaic_w_id)
      read(mozaic_w_id) ((wo2a(ija,k),k=1,kamax),ija=1,ijatm)
      read(mozaic_w_id)
      read(mozaic_w_id) ((inda2o(ijo,k),k=1,komax),ijo=1,ijocn)
      read(mozaic_w_id)
      read(mozaic_w_id) ((wa2o(ijo,k),k=1,komax),ijo=1,ijocn)

      close(mozaic_w_id)
      do k=1,kamax
        do ija=1,ijatm
          if (indo2a(ija,k).eq.0) then
            indo2a(ija,k)=1
            wo2a(ija,k)=0d0
          endif
        enddo
      enddo

      do k=1,komax
        do ijo=1,ijocn
          if (inda2o(ijo,k).eq.0) then
            inda2o(ijo,k)=1
            wa2o(ijo,k)=0d0
          endif
        enddo
      enddo

!      else
!
!      read(mozaic_w_id)
!      read(mozaic_w_id) ((indo2a2(ija,k),k=1,kamax2),ija=1,ijatm)
!      read(mozaic_w_id)
!      read(mozaic_w_id) ((wo2a2(ija,k),k=1,kamax2),ija=1,ijatm)
!      read(mozaic_w_id)
!      read(mozaic_w_id) ((inda2o(ijo,k),k=1,komax),ijo=1,ijocn)
!      read(mozaic_w_id)
!      read(mozaic_w_id) ((wa2o(ijo,k),k=1,komax),ijo=1,ijocn)
!
!      close(mozaic_w_id)
!      do k=1,kamax2
!        do ija=1,ijatm
!          indo2a(ija,k)=indo2a2(ija,k)
!          wo2a(ija,k)=wo2a2(ija,k)
!          if (indo2a(ija,k).eq.0) then
!            indo2a(ija,k)=1
!            wo2a(ija,k)=0d0
!          endif
!        enddo
!      enddo
!
!      do k=1,komax
!        do ijo=1,ijocn
!          if (inda2o(ijo,k).eq.0) then
!            inda2o(ijo,k)=1
!            wa2o(ijo,k)=0d0
!          endif
!        enddo
!      enddo
!
!      endif



!ocean basins
      call ec_iniocbas

c  more  parts of ec_initcoup
!land data to coupler
      call ec_la2co
      
c~       call ec_lae2co

!oceanic data to coupler
!here it compiles again fractn(nse) et fractn(noc)
      call ec_oc2co(1)
      

!coupler data to atmosphere
c~       call ec_at2co


! *** initialisation of fluxes
      do nn=1,ntyps
        call ec_fluxes(nn)
      enddo

cnb --- end of parts of ec_initcoup
  

!5. Other things in emic.f
!!!!!!!!!!!!!!!!!!!!!!!!!! 
cnb --- dans emic.f

c~ #if (DOWNSTS == 1 || DOWNSCALING == 2)
c~       call Init_subgrid
c~ #endif

#if ( F_PALAEO == 1 )
c~          CALL load_topo_masq(.TRUE.) ! parameter given is init
c~          CALL ec_masq(.FALSE.)       ! not an init
         CALL update_masq
         CALL ec_topo                ! update topography
#elif ( F_PALAEO == 3)
c~          call load_topo_masq_hr(.TRUE.) ! afq -- we read a prescribedtopography
c~          CALL ec_topo
c~          call subgrid_ecb_wrapper
c~          call ec_masq(.FALSE.)
         CALL update_masq
         call ec_topo
         
#elif ( ISM >= 2  || DOWNSCALING == 2 )
         call ec_topo                ! update topography, afq (rmountetc.)
c~          call subgrid_ecb_wrapper
         CALL update_masq
c~          call ec_masq(.FALSE.)
         call ec_topo
#endif

!dmr&nb [ICICESTBON]

         CALL clio_masq

      end subroutine change_grid
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!endif BATHY==3
#endif

#if ( BATHY >= 2 )

!------------------------------------------------------------------------------
      subroutine read_bering
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      use comunit, only: iuo
      use update_clio_bathy_tools, only: bering_nb, bering_date,
     >       bering_value
      use global_constants_mod, only: ip

      integer(ip) :: Bering_scenario_txt_id

        !write(*,*) 'bering_nb= ', bering_nb
        open(newunit=Bering_scenario_txt_id,
     >        file='inputdata/Bering_scenario.txt')
        do i=1,bering_nb
         read(Bering_scenario_txt_id,*) bering_date(i), bering_value(i)
         !write(*,*) 'BERING ', bering_date(i), bering_value(i), i
        enddo

        close(Bering_scenario_txt_id)

      end subroutine read_bering
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!------------------------------------------------------------------------------
      subroutine read_kamax
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      use comunit, only: iuo
      use update_clio_bathy_tools, only: kamax_nb, kamax_date,
     >       kamax_value
      use global_constants_mod, only: ip

      integer(ip) :: Kamax_scenario_txt_id

        !write(*,*) 'kamax_nb= ', kamax_nb
        open(newunit=Kamax_scenario_txt_id,file='inputdata/Kamax_scenario.txt')
        do i=1,kamax_nb
         read(Kamax_scenario_txt_id,*) kamax_date(i), kamax_value(i)
         !write(*,*) 'KAMAX ', kamax_date(i), kamax_value(i), i
        enddo

        close(Kamax_scenario_txt_id)

      end subroutine read_kamax
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!------------------------------------------------------------------------------
      subroutine track_error(value_v, min_v, max_v)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      use global_constants_mod, only: dblp=>dp, ip

      real,    dimension(:,:), intent(inout)  :: value_v
!      real,    dimension(:,:,:), intent(inout)  :: value_v
      real(dblp) :: min_v, max_v

      integer(ip) :: i, j, k
      integer(ip) :: i_min, i_max, j_min, j_max
      integer(ip) :: k_min, k_max


       i_min = LBOUND(value_v, dim=1)
       i_max = UBOUND(value_v, dim=1)

       j_min = LBOUND(value_v, dim=2)
       j_max = UBOUND(value_v, dim=2)

!       k_min = LBOUND(value_v, dim=3)
!       k_max = UBOUND(value_v, dim=3)

       do i=i_min, i_max
         do j=j_min, j_max
!           do k=k_min, k_max

!            if ((value_v(i,j,k) .gt. max_v) .or. 
!     >          (value_v(i,j,k) .lt. min_v)) then
            if ((value_v(i,j) .gt. max_v) .or. 
     >          (value_v(i,j) .lt. min_v)) then
              write(*,*) 'min/max value exceeded ' 
!     >               , value_v(i,j,k), i, j, k
     >               , value_v(i,j), i, j
            endif

!           enddo
         enddo
       enddo


      end subroutine track_error
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|



!endif BATHY==2
#endif
