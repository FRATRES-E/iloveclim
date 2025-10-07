!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:21 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:21 CET 2009
c********************************************************************
      SUBROUTINE INITMBIOPAR
c********************************************************************

      use global_constants_mod, only: stderr, ip, dblp=>dp

      use declars_mod, only: jt, jx, lt, noc_cbr

      use loveclim_transfer_mod, only: mgt, zx, zz
      use mod_sync_time, only: tday, tyer

      use mbiota_mod, only: SCALE_M, SCALE_B, SCANU, C14RA
     >                   , RR, N0_M, PI_M, PAR_M, LEF_M, PMIN_M, DP_M
     >                   , ER_DOC, A_M, B_M, C_M, P0_M, G0_M, EX_DOC
     >                   , DZ_M, ZMIN_M, L0_M, D0_M, KD_M, RDOC_MS
     >              , OrgCFlxAttFactor, SUE_MCA, B_SH_FR, oc_bottom_cell
     >                   , SWR_FRAC
     >                   , O2min1, O2min2

#if ( REMIN == 1 )
      use mbiota_mod, only: SUE_3D, kremin
     >                   , Ea, Rgaz, betaPOM
     >                   , compute_remin
#endif

#if ( REMIN == 1 || REMIN_CACO3 == 1 )
      use mbiota_mod, only: dt, w_sink
#endif

#if ( REMIN_CACO3 == 1 )
      use mbiota_mod, only: kremin_ca, kremin_ar
     >                   , k_diss
     >                   , SUE_ca_3D, SUE_ar_3D
     >                   , compute_remin_ca
#endif

#if ( ARAG == 1 )
      use mbiota_mod, only: RR_ar, Kmax, SUE_MAR
#endif

      use marine_bio_mod, only: eher, zinges, ecan, sigma_m, sigma_md
     >                   , jprod
     >                   , OPO4_INI, ONO3_INI, OSI_INI, OALK_INI
     >                   , ODOCS_INI

!REFACTORING DONE: Oeta(:,1:5) replaced by OetaC*, OetaN, OetaO2*
!  Removed Oeta_floor, which has become partially obsolete, ans was
!  partially replaced OetaC_POMsedin, OetaN_POMsedin, OetaO2_POMsedin
!  for usage in MEDUSA
      use marine_bio_mod, only: 
     >   OetaC_POMoxid, OetaN_POMoxid, OetaO2_POMoxid,
     >   OetaC_DOMoxid_1D, OetaN_DOMoxid_1D, OetaO2_DOMoxid_1D,
     >   OetaC_POMrain, OetaN_POMrain

#if ( MEDUSA == 1 )
      use marine_bio_mod, only:
     >   OetaC_POMsedin, OetaN_POMsedin, OetaO2_POMsedin
#endif

#if ( REMIN_CACO3 == 1 )
      use omega_mod, only : calc_omega_ca, omega_calc3D
#endif

#if ( ARAG == 1 )
      use omega_mod, only : calc_omega_ar, omega_arag3D
#endif


      implicit none

c********************************************************************

       INTEGER(kind=dblp) :: i, j, n
       REAL(kind=dblp)    :: prom

!REFACTORING: the export depth was previously hard coded
!REFACTORING CORRECTION (08OCT2024): removed parameter attribute of zp_xp
!                                    -> needs to be adjusted to zx(jprod+1) below
       REAL(kind=dblp) :: zp_xp = 100.D+00      ! Export depth

       REAL(kind=dblp) :: z_xpratio
       REAL(kind=dblp), dimension(JX) :: OrgCFlxAtt_1D     ! 1D flux attenuation profile (auxiliary)

                                    ! Auxiliary arrays
       REAL(kind=dblp) :: OetaC_POMxp, OetaN_POMxp   ! C:P and N:P ratios in the POM export
       REAL(kind=dblp) :: OetaC_DOMxp, OetaN_DOMxp   ! C:P and N:P ratios in the DOM "export"

       REAL(kind=dblp), dimension(JT) :: OetaC_POMoxid_1D  ! C:P of oxidized POM products
       REAL(kind=dblp), dimension(JT) :: OetaN_POMoxid_1D  ! N:P of oxidized POM products
       REAL(kind=dblp), dimension(JT) :: OetaO2_POMoxid_1D ! \DeltaO2:\DeltaP during for POM oxidation
       REAL(kind=dblp), dimension(JX) :: OetaC_POMrain_1D  ! C:P in POM rain (water column flux)
       REAL(kind=dblp), dimension(JX) :: OetaN_POMrain_1D  ! C:P in POM rain (water column flux)

       REAL(kind=dblp) :: OrgCFlxAttRatio

       integer :: jrain, joxid      ! Convenience index variables


! dmr [2024-10-14] Moved parameter values into a file ... 

       CHARACTER(len=18), PARAMETER :: file_path="mbio_parameter.nml"
       INTEGER(kind=ip)             :: rc, fu
       REAL(kind=dblp)              :: dp_mfac, er_docfac, g0_mfac
     &                               , ex_docfac, dz_mfac, l0_mfac
     &                               , d0_mfac, rdoc_msfac

       NAMELIST /mbiopar/ SCALE_M, SCALE_B, SCANU, C14RA, RR
     &                 , n0_m, PI_m, PAR_m, lef_m, pmin_m
     &                 , dp_mfac, er_docfac, a_m, b_m, c_m
     &                 , p0_m, g0_mfac, eher, zinges, ex_docfac
     &                 , dz_mfac, ecan, zmin_m, l0_mfac 
     &                 , d0_mfac, kd_m, O2min1, O2min2 
     &                 , rdoc_msfac, sigma_m, sigma_md


c********************************************************************

! dmr [2024-10-14] Moved parameter values into a file ... reading that one now

        INQUIRE (file=file_path, iostat=rc)

        IF (rc /= 0) THEN
            WRITE (stderr, '("Error: input file ", a, " does not exist")') 
     &            file_path
            STOP
        ENDIF

        ! Open and read Namelist file.
        OPEN (action='read', file=file_path, iostat=rc, newunit=fu)
        READ (nml=mbiopar, iostat=rc, unit=fu)
        IF (rc /= 0) WRITE (stderr, '("Error: invalid Namelist format")')

        CLOSE (fu)

#if ( ARAG == 1 )
        RR_ar=0.35*RR
        RR=(1-RR_ar)
        Kmax=0.4
#endif

! dmr initialisation of parameters that have a time dependence: rates are per second

c phytoplankton
        dp_m=dp_mfac/TDAY
        er_doc=er_docfac/TDAY

        g0_m=g0_mfac/TDAY

        ex_doc=ex_docfac/TDAY
        dz_m=dz_mfac/TDAY

c detritus
        l0_m=l0_mfac/TYER

c DOC
        d0_m=d0_mfac/TDAY

c slow DOC
        rdoc_ms=rdoc_msfac/TYER


! Flux attenuation factors 
! ------------------------
!
!   - OrgCFlxAttFactor(i,j,n): attenuation factor for cell (i,j,n),
!     such that OrgC flux raining into (i,j,n) is equal to
!     OrgCFlxAttFactor(i,j,n) * OrgC_export_flux(i,n)
!
!   - We start with a 1D profile OrgCFlxAtt_1D(j) and replicate it
!     at each in each ocean (i,n) point, and clipping it as a function
!     of depth there

!REFACTORING DONE: OrgCFlxAttFactor was SUE_M(:)

!REFACTORING NOTE: 
! The following should go into a separate function subroutine
! allowing us to write OrgCFlxAtt_1D = fct(zx),
! and have the different flux attenuation formulations
! located there, together with their parameters 

!REFACTORING CORRECTION (08OCT2024)
! export depth zp_xp adjusted to the bottom of the photic/productive zone
      zp_xp = zx(JPROD+1)

      OrgCFlxAtt_1D(1:JPROD+1) = 1D0

      do j = JPROD + 1, JT
        ! Currently: Martin profile with b = 0.858
        ! Alternatives previously used:
        ! OrgCFlxAtt_1D(j+1) = min(max((0.8-1.0)*z_xpratio*(100./1000.)+1.,0.0),1.0)
        ! OrgCFlxAtt_1D(j+1) = z_xpratio**(-0.5)! cnb test

        z_xpratio = zx(j+1)/zp_xp
        OrgCFlxAtt_1D(j+1) = z_xpratio**(-0.858D0)
      enddo

      ! Now replicate the 1D function to the 3D array OrgCFlxAttFactor
      ! and clip with local depth as required

!REFACTORING NOTE:
! The following can be boiled down to a few lines in a loop
! over i, n if the index of the bottom cell is used:
! jbot = index of bottom cell at (i,:,n):
! #if ( MEDUSA == 1 )
! jrain = MIN(jbot+2, JX)
! #else
! jrain = jbot+1
! #endif
! OrgCFlxAttFactor(i, 1:jrain, n) = OrgCFlxAtt_1D(1:jrain)
! OrgCFlxAttFactor(i, jrain+1:JX, n) = 0D0

      do n = 1, NOC_CBR
        do i = 1, LT

          if (MGT(i,1,n).eq.1) then

            OrgCFlxAttFactor(i, 1, n) = OrgCFlxAtt_1D(1)

            do j = 1, JT

              OrgCFlxAttFactor(i, j+1, n) = OrgCFlxAtt_1D(j+1)

              if (oc_bottom_cell(i,j,n)) then

!REFACTORING NOTE:
! The separate treatment of the bottom cell better had to be removed,
! and, if no resolved sediment module is coupled to OCYCC, a virtual
! one implementing reflective boundary conditions should be adopted
! instead. In this case, only the attenuation factors for the fluxes
! leaving cells below the bottom cell would be set to zero. The flux
! processing in the water column and at or below the seafloor can then
! be cleanly separated.
#if ( MEDUSA == 1 )
                jrain = j+2         ! If a sediment module is coupled to
                                    ! OCYCC, the OrgC flux leaving the
                                    ! bottom cell is kept (i.e., flux(i,j+1,n))
                                    ! and only all those below that one
                                    ! are set to zero (i.e., flux(i,j+2:JX,n) = 0)
#else
                jrain = j+1         ! If no sediment module is coupled
                                    ! to OCYCC, the OrgC flux leaving the
                                    ! bottom cell and all those below are
                                    ! set to zero (i.e., flux(i,j+1:JX,n) = 0)
#endif
                OrgCFlxAttFactor(i, jrain:JX, n) = 0D0

                exit

              endif

            enddo

          else

!            FlxAttProf_OrgC(i,:,n) = 0D0

          endif

        enddo
      enddo


! OetaC*, OetaN*, OetaO2* function construction
! ---------------------------------------------
! For each cell (i,j,n),
!  -  OetaC_POMoxid(i,j,n) is the molar C:P in the remineralisation products
!     of organic matter oxidation in the cell;
!  -  OetaN_POMoxid(i,j,n) is the molar N:P in the remineralisation products
!     of organic matter oxidation in the cell;
!  -  OetaO2_POMoxid(i,j,n) is the molecular O2:P consumption ratio (O2:P < 0)
!     for organic matter oxidation in the cell
!  -  We start with generic 1D profiles for each of these, and
!     replicate them first in each ocean (i,n) point, and then take
!     corrections relevant for depth into account


      OetaC_POMoxid_1D(:) =         ! was column 4
     & (/ 106., 106., 106., 106., 106., 106., 106., 106., 106., 106.,
     &    106., 106., 115., 120., 125., 125., 125., 125., 125., 125.  /)

      OetaN_POMoxid_1D(:) =         ! was column 1
     & (/16., 16., 16., 16., 16., 16., 16., 16., 16., 16.,
     &   16., 16., 16., 16., 16., 16., 16., 16., 16., 16. /)

      OetaC_DOMoxid_1D(:) =         ! was column 5
     & (/ 106., 106., 106., 106., 106., 106., 106., 106., 106., 106.,
     &    106., 106., 106., 106., 106., 106., 106., 106., 106., 106.  /)

      OetaN_DOMoxid_1D(:) =         ! new column (duplicate of column 1)
     & (/16., 16., 16., 16., 16., 16., 16., 16., 16., 16.,
     &   16., 16., 16., 16., 16., 16., 16., 16., 16., 16. /)

!     Oeta1Dref_POM_Si(:) =         ! was column 3, currently not used
!    & (/ 10., 10., 10., 10., 10., 10., 10., 10., 10., 10.
!    &    10., 10., 10., 10., 10., 10., 10., 10., 10., 10. /)

!      write(*,*) 'dans initmbiopar OetaC_POMoxid_1D', OetaC_POMoxid_1D(:)

! Derived values
!  1. OetaC_POMoxid_1D(1:JPROD) and OetaN_POMoxid_1D(1:JPROD)
!     must be constant, because else the export OetaC and OetaN
!     values may vary in time, which makes it impossible to reliably
!     pre-calculate OetaC_POMrain(:,:,:) etc.
!     We therefore copy the top value into 2:JPROD

      OetaC_POMxp = OetaC_POMoxid_1D(1)
      OetaC_POMoxid_1D(2:JPROD) = OetaC_POMxp

      OetaN_POMxp = OetaN_POMoxid_1D(1)
      OetaN_POMoxid_1D(2:JPROD) = OetaN_POMxp

!  2. Similarly for the Oeta arrays of DOM
      OetaC_DOMxp = OetaC_DOMoxid_1D(1)
      OetaC_DOMoxid_1D(2:JPROD) = OetaC_DOMxp

      OetaN_DOMxp = OetaN_DOMoxid_1D(1)
      OetaN_DOMoxid_1D(2:JPROD) = OetaN_DOMxp

!  3. OetaO2_X(:) = -(OetaC_X(:) + 2 * OetaN_X(:))

      OetaO2_POMoxid_1D(:) =        ! was column 2
     &  -OetaC_POMoxid_1D(:) - 2D0*OetaN_POMoxid_1D(:)

      OetaO2_DOMoxid_1D(:) =        ! new column (duplicate of column 2)
     &  -OetaC_DOMoxid_1D(:) - 2D0*OetaN_DOMoxid_1D(:)


! Set up the auxiliary 1D profiles to derive the OetaC and OetaN values
! in the POM *rain*

      ! - export ratio used for 1:JPROD+1
      OetaC_POMrain_1D(1:JPROD+1) = OetaC_POMxp
      OetaN_POMrain_1D(1:JPROD+1) = OetaN_POMxp

      ! ratios in the rain derived recursively below
      ! the productive zone

      do j = JPROD + 1, JT

        OrgCFlxAttRatio = OrgCFlxAtt_1D(j)/OrgCFlxAtt_1D(j+1)
!        write(*,*) 'initmbiopar, OrgCFlxAttRatio ', OrgCFlxAttRatio

!nb&gm [NOTA] Refactored POMrain ...

        OetaC_POMrain_1D(j+1) =
     &     1D0 / (OrgCFlxAttRatio / OetaC_POMrain_1D(j)
     &             - (OrgCFlxAttRatio-1D0) / OetaC_POMoxid_1D(j))

!        OetaC_POMrain_1D(j+1) = 1D0 / 
!     &    (OrgCFlxAttRatio / OetaC_POMrain_1D(j) -
!     &    (OrgCFlxAttRatio-1D0) / OetaC_POMoxid_1D(j))
!
!        OetaN_POMrain_1D(j+1) = 1D0 /
!     &    (OrgCFlxAttRatio / OetaN_POMrain_1D(j) -
!     &    (OrgCFlxAttRatio-1D0) / OetaN_POMoxid_1D(j))

        OetaC_POMrain_1D(j+1) = OetaC_POMrain_1D(j) +
     &    (OetaC_POMoxid_1D(j) - OetaC_POMrain_1D(j))
     &    * (1d0 - OrgCFlxAttRatio)
     &      / ( OrgCFlxAttRatio/OetaC_POMrain_1D(j)
     &          * (OetaC_POMoxid_1D(j) - OetaC_POMrain_1D(j)) + 1d0 )

        OetaN_POMrain_1D(j+1) = OetaN_POMrain_1D(j) - 
     &    (OetaN_POMoxid_1D(j) - OetaN_POMrain_1D(j))
     &    * OrgCFlxAttRatio
     &    * (OetaC_POMrain_1D(j+1)/OetaC_POMrain_1D(j))

      enddo

!      write(*,*) 'initmbiopar, OetaC_POMrain_1D ', OetaC_POMrain_1D(:)


! Now proceed to the OetaC_POMrain, OetaN_POMrain and OetaO2_POMrain
! definitions, and adjust them as a function of the chosen seafloor
! flux processing

      do n = 1, NOC_CBR
        do i = 1, LT

          if (MGT(i,1,n).eq.1) then

            OetaC_POMoxid(i,:,n) = OetaC_POMoxid_1D(:)
            OetaN_POMoxid(i,:,n) = OetaN_POMoxid_1D(:)

            OetaC_POMrain(i, 1, n) = OetaC_POMrain_1D(1)
            OetaN_POMrain(i, 1, n) = OetaN_POMrain_1D(1)

            do j = 1, JT

              OetaC_POMrain(i, j+1, n) = OetaC_POMrain_1D(j+1)
              OetaN_POMrain(i, j+1, n) = OetaN_POMrain_1D(j+1)

              if (oc_bottom_cell(i,j,n)) then

#if ( MEDUSA == 1 )
                ! Initialise the C:P, N:P characteristics (2D)
                ! to be transmitted to MEDUSA 
                OetaC_POMsedin(n, i) = OetaC_POMrain(i, j+1, n)
                OetaN_POMsedin(n, i) = OetaN_POMrain(i, j+1, n)
                OetaO2_POMsedin(n, i) = 
     &            -OetaC_POMsedin(n, i) - 2D0*OetaN_POMsedin(n, i)

                ! Set to zero all ratios for cells deeper than
                ! j_bottom+1 and all oxidation ratios for cells deeper than
                ! j_bottom+1
                jrain = j+2
                joxid = j+1
#else
                !REFACTORING Remove this pre-compile option by adding
                ! "fake" sediment module to implement reflective
                ! boundary condition

                ! No sediment? Override standard oxidation ratios so that
                ! they match the ratios in the incoming POM rain
!                write(*,*) 'initmbiopar', j, OetaC_POMrain(i, j, n)
                OetaC_POMoxid(i, j, n) = OetaC_POMrain(i, j, n)
                OetaN_POMoxid(i, j, n) = OetaN_POMrain(i, j, n)

                jrain = j+1
                joxid = j+1
#endif
                OetaC_POMoxid(i, joxid:JT, n) = 0D0
                OetaN_POMoxid(i, joxid:JT, n) = 0D0

                OetaC_POMrain(i, jrain:JX, n) = 0D0
                OetaN_POMrain(i, jrain:JX, n) = 0D0

                exit                ! Skip the rest of this column
                                    ! as we are done

              endif

              ! Finally adapt OetaO2_POMoxid data
              OetaO2_POMoxid(i, :, n) = 
     &          -OetaC_POMoxid(i, :, n) - 2D0 * OetaN_POMoxid(i, :, n)

            enddo

          else

            OetaC_POMoxid(i,:,n) = 0D0
            OetaN_POMoxid(i,:,n) = 0D0
            OetaO2_POMoxid(i,:,n) = 0D0

            OetaC_POMrain(i,:,n) = 0D0
            OetaN_POMrain(i,:,n) = 0D0

          endif

        enddo
      enddo

!      write(*,*) 'dans initmbiopar OetaC_POMoxid ', OetaC_POMoxid(:,:,:)


! Flux attenuation factor for CaCO3
! ---------------------------------

!         do j=1,JPROD+1 
!            SUE_MCA(j)=1.
!         enddo 
!
!         do j=JPROD+2,JX 
!         prom=(ZX(j)-100.)/3000.  
!         SUE_MCA(J)=dexp(-prom)
!          enddo

      SUE_MCA(1:JPROD+1) = 1d0

      do j = JPROD+2, JX

        prom = (ZX(j)-zp_xp)/3000.d0
        SUE_MCA(J) = dexp(-prom)

      enddo

! Correction for big shell fraction:
! - b_sh_fr: fraction of CaCO3 shells that are big enough to fall to
!   sufficiently fast to the seafloor that they escape
!   dissolution in the water column

      b_sh_fr = 0.2
!      b_sh_fr = 0.0 ! test

#if ( ARAG == 1 )
      SUE_MAR(:) = SUE_MCA(:)
#endif


c-----
#if ( REMIN == 1 )
cnb   3D fixed remineralisation profile (test)
      dt=1.*TDAY ! 1 day in s
      ! Ref : Gangsto; Crichton, 2021

!test case to have same values as in original fixed profile
!      w_sink(:)=125 !7.0 !m/day vertical speed, here fixed
!      write(*,*) ' R '
!      i=1
!      n=1
!      do j=1,JT
!        kremin(i,j,n)=SUE_M(j)/dt*zz(j)/w_sink(j)
!        write(*,*), j, kremin(i,j,n)
!      enddo


!      kremin=   (/0.000000000000000E+000,
!     >            0.000000000000000E+000,
!     >            0.000000000000000E+000,
!     >            0.000000000000000E+000,
!     >            0.000000000000000E+000,
!     >            0.000000000000000E+000,
!     >            0.286720000000000,
!     >            0.286775381872049,    
!     >            0.305998655712954,    
!     >            0.337330838260894,    
!     >            0.381091640376637,    
!     >            0.432941417314840,    
!     >            0.475567768019538,     
!     >            0.481632215248118,     
!     >            0.444462296679600,     
!     >            0.386877036145944,     
!     >            0.330747248936669,     
!     >            0.284053117399321,     
!     >            0.247003862864335,     
!     >            0.217768123631704/)
!       write(*,*), kremin

      w_sink(:)=125.0! in m/day !!*1./TDAY ! m/day vertical speed -> in m/s 
      Ea=56*1e3 !55*1e3 ! activation energy in J/mol 50 to 60
      Rgaz=8.314 ! gas constant in J/K/mol
      betaPOM=1*1e12/TYER*TDAY !/TYER !9*1e11 !ou 1*1e14 rate constant for POM remineralisation as temperature approaches infinity in /yr -> in /s
      do n=1,NOC_CBR
       do j=1,JT
        do i=1,LT
           call compute_remin(i,j,n) 
        enddo
       enddo
      enddo

!test with original remin profile
!      write(*,*) ' SUE_3D and SUE_M '
!       do j=1,JT
!!           SUE_3D(1,j,1)= dt/zz(j)*w_sink(j)*kremin(1,j,1)
!           write(*,*) j, SUE_3D(1,j,1), SUE_M(j)
!       enddo


#endif


#if ( REMIN_CACO3 == 1 )
!remineralisation for CaCO3, depends on omega
! see Gehlen 2007, Gangsto 2008
      dt=1.!1 day !!*TDAY ! 1 day in s
      w_sink(:)=125.0 !in m/day !! *1./TDAY ! m/day vertical speed -> in m/s 
      k_diss=10.9 !6 ! in per day !!/TDAY ! dissolution rate constant per day- > per s
      !n_reac=1
      do n=1,NOC_CBR
       do j=1,JT
        do i=1,LT
           call calc_omega_ca(i,j,n)
           call compute_remin_ca(i,j,n,omega_calc3D(i,j,n),
     >                          kremin_ca(i,j,n), kremin_ca(i,j-1,n))
#if ( ARAG == 1 )
           call calc_omega_ar(i,j,n)
           call compute_remin_ca(i,j,n,omega_arag3D(i,j,n),
     >                          kremin_ar(i,j,n),kremin_ar(i,j-1,n))
#endif
        enddo
       enddo
      enddo


#endif
c-----

cnb - Initialisation des valeurs initiales dans l ocean
c initial concentrations; units are mumol/kg for PO4, NO3, DOC;
c                         mol/kg for DIC;
c                         eq /kg for ALK

cnb -  Phosphate
#if ( LGMSWITCH == 0 )
        OPO4_ini=2.08 ! LH
#else
        OPO4_ini=2.15 ! LGM avec sea level change 3.3 percent   
#endif

cnb - Nitrates
#if ( LGMSWITCH == 0 )
        ONO3_ini=33.6 ! LH
#else
        ONO3_ini=34.7 ! LGM avec sea level change 3.3 percent
#endif

#if ( OXNITREUX == 1 )

#if ( LGMSWITCH == 0 )
cvm - Nitrous oxide
        ON2O_ini=31.9 !LH
#else
        ON2O_ini=33.0 !LGM avec sea level change 3.3 percent
#endif

#endif 

cnb - Silice
       OSI_ini=0.0

cnb - Alkalinite

#if ( LGMSWITCH == 0 )

        OALK_ini=2373.0e-6 ! LH
#else
        OALK_ini=2451.0e-6 ! LGM avec sea level change 3.3 percent
#endif

       ODOCS_ini=4.434

!dmr --- Dimension on "JT"
!    --- Transmissivity of swr into the ocean

        do j=1,UBOUND(SWR_FRAC,DIM=1)
          SWR_FRAC(j) = PI_M*PAR_m/zz(j)/lef_m
     >                          *(exp(-lef_m*zx(j))-exp(-lef_m*zx(j+1)))
        enddo

        return
        END SUBROUTINE INITMBIOPAR
