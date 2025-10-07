!=======================================================================
      MODULE MOD_ILOVECLIM_SETUP
!=======================================================================

! This module is part of the interface with iLoveclim.
! It serves as an interface to communicate fundamental information
! about

! It is included in the main program that controls the iLoveclim
! sequence of calls. That program should fill in the geometrical
! information required for the setup into the variables and arrays, and
! upon completion, set iloveclim_setup_init to .TRUE.


      IMPLICIT NONE

      ! CLIO related information (grid)
      ! ===============================

      INTEGER :: iloveclim_nix = -1  ! Number of "longitudes" (extension of first dim of 2D arrays)
      INTEGER :: iloveclim_njy = -1  ! Number of "latitudes" (extension of second dim of 2D arrays)

      INTEGER,          DIMENSION(:,:), ALLOCATABLE, SAVE
     &   :: iloveclim_ocngrid_mask
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, SAVE
     &   :: iloveclim_ocngrid_surf
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, SAVE
     &   :: iloveclim_ocngrid_dbsl
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, SAVE
     &   :: iloveclim_ocngrid_lon
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, SAVE
     &   :: iloveclim_ocngrid_lat


      ! OCYCC related information (organic matter stoechiometry)
      ! ========================================================

      INTEGER,          DIMENSION(:,:), ALLOCATABLE, SAVE
     &   :: kfs_fond

      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, SAVE
     &   :: iloveclim_om_c, iloveclim_om_n, iloveclim_om_p,
     &      iloveclim_om_ro2


      LOGICAL :: iloveclim_setup_init = .FALSE.


      CONTAINS



!-----------------------------------------------------------------------
      SUBROUTINE setup_iloveclim_medusa_xchange
!-----------------------------------------------------------------------

!-----|--1--------2---------3---------4---------5---------6---------7-|
!     Cette subroutine contient les variables a echanger entre
!     iLOVECLIM et MEDUSA au moment de l'initialisation. Elle est
!     appelée une fois au début de simulation
!     (dans l'environnement logiciel LUDUS)
!
!     Auteur : Didier M. Roche & Guy Munhoven
!     Modified by: Santiago Moreira
!     Date   : 03 fevrier 2014
!     Derniere modification : 21 octobre 2016 (Guy Munhoven), 22 sept. 2021 (dmr)
!-----|--1--------2---------3---------4---------5---------6---------7-|


      USE mod_defines_medusa, ONLY: jp_stderr
      USE mod_execontrol_medusa, ONLY: ABORT_MEDUSA

      USE declars_mod, ONLY: LT, JT, NOC_CBR
!tbd       USE marine_bio_mod, ONLY: Oeta_floor, JPROD
      USE marine_bio_mod, ONLY: JPROD, OetaC_POMsedin, OetaN_POMsedin
     &                          , OetaO2_POMsedin

      use const_mod, only:
      use para_mod, only:
      use para0_mod, only: imax, jmax, kmax
      use bloc0_mod, only: kfs, tms, dz, xslon, yslat
      use bloc_mod, only:
      use dynami_mod, only: area

      CHARACTER(LEN=*), PARAMETER :: c_fmterr_a =
     &  '("[MOD_ILOVECLIM_SETUP/' //
     &  'setup_iloveclim_medusa_xchange]] error: ", A)'

!dmr !!!
!dmr --- Commentaire
!dmr       Pour CLIO : VarOcean(imax,jmax,kmax)   => Lon, Lat, Levels
!dmr       Dans CLIO, la surface est à k = kmax
!dmr !!!

      INTEGER i,j

      IF (iloveclim_setup_init) THEN
        WRITE(jp_stderr,c_fmterr_a) 'Setup already done -- aborting!'
        CALL ABORT_MEDUSA()
      ENDIF

!dmr --- [NOTA] imax contains Halo Points: 1 and 122 are replica of 2 and 121 swapped

      iloveclim_nix = imax ! 122
      iloveclim_njy = jmax ! 65
                                    
                                    ! Allocation of arrays on the CLIO grid
                                    
!dmr --- Le masque terre -- océan vu par le fond de CLIO
      ALLOCATE(iloveclim_ocngrid_mask(iloveclim_nix,iloveclim_njy))
!dmr --- La surface des cases CLIO
      ALLOCATE(iloveclim_ocngrid_surf(iloveclim_nix,iloveclim_njy))
!dmr --- La profondeur de l'océan en metre de la surface au fond
      ALLOCATE(iloveclim_ocngrid_dbsl(iloveclim_nix,iloveclim_njy))
!dmr --- Longitude on the CLIO grid
      ALLOCATE(iloveclim_ocngrid_lon(iloveclim_nix,iloveclim_njy))
!dmr --- Latitude on the CLIO grid      
      ALLOCATE(iloveclim_ocngrid_lat(iloveclim_nix,iloveclim_njy))
   
!dmr --- Variables in OCYCC are on grid:
! --- LT = 65, JT = 20 ,NOC_CBR = 122
!dmr --- However, MEDUSA is working on CLIO, so remove that dependencie as much as possible

! --- kfs_fond is on the Carbon Cycle (OCYCC grid)
      ALLOCATE(kfs_fond(LT,NOC_CBR)) 

! --- those variables are on the CLIO grid ...
      ALLOCATE(iloveclim_om_c(iloveclim_nix,iloveclim_njy))
      ALLOCATE(iloveclim_om_n(iloveclim_nix,iloveclim_njy))
      ALLOCATE(iloveclim_om_p(iloveclim_nix,iloveclim_njy))
      ALLOCATE(iloveclim_om_ro2(iloveclim_nix,iloveclim_njy))


                                    ! Initialise kfs_fond
                                    ! (Transpose CLIO's kfs and change the direction by using kmax+1-z
      do i=1, LT
        do j=1, NOC_CBR
          kfs_fond(i,j) = JT+1-kfs(j+1,i)
        enddo
      enddo

#ifdef DEBUG
      print*,"shape kfs_fond", shape(kfs_fond)
#endif

      FORALL(i=1:iloveclim_nix,j=1:iloveclim_njy)

        iloveclim_ocngrid_mask(i,j) = tms(i,j,kmax) ! the ocean grid is where there is water at the surface ...
        iloveclim_ocngrid_surf(i,j) = area(i,j) 
        iloveclim_ocngrid_dbsl(i,j) = SUM(dz(kfs(i,j):kmax))
        iloveclim_ocngrid_lon(i,j)  = xslon(i,j)
        iloveclim_ocngrid_lat(i,j)  = yslat(i,j)

      END FORALL

!--- dmr ignore halo points in MEDUSA (mask them in MEDUSA)
      iloveclim_ocngrid_mask(1,:)             = 0
      iloveclim_ocngrid_mask(iloveclim_nix,:) = 0

c~       print*, "values of dbsl..."
c~       print*, iloveclim_ocngrid_dbsl      
c~       print*, "end of values dbsl..."
      
      DO i = 2, iloveclim_nix-1
        DO j = 1, iloveclim_njy

                                    ! variables Oeta_floor are: dim(JT,5)
          IF (iloveclim_ocngrid_mask(i,j) > 0) THEN
                                    ! [TODO] not sure to agree with that below ... transpose for i,j on kfs_fond ... why not kfs?
            iloveclim_om_c(i,j)
!nb     &        = Oeta_floor(MAX(INT(kfs_fond(j,i-1)),JPROD+1),4)
     &        = OetaC_POMsedin(i-1,j)
            iloveclim_om_n(i,j)
!nb     &        = Oeta_floor(MAX(INT(kfs_fond(j,i-1)),JPROD+1),1)
     &        = OetaN_POMsedin(i-1,j)
            iloveclim_om_p(i,j)
     &        = 1.0D+00
            iloveclim_om_ro2(i,j)
!nb     &        = ABS(Oeta_floor(MAX(INT(kfs_fond(j,i-1)),JPROD+1),2))
     &        = ABS(OetaO2_POMsedin(i-1,j))

          ENDIF
        ENDDO
      ENDDO

! Halo points ...
      iloveclim_om_c(1,:) = iloveclim_om_c(iloveclim_nix-1,:)
      iloveclim_om_c(iloveclim_nix,:) = iloveclim_om_c(2,:)
      iloveclim_om_n(1,:) = iloveclim_om_c(iloveclim_nix-1,:)
      iloveclim_om_n(iloveclim_nix,:) = iloveclim_om_c(2,:)
      iloveclim_om_p(1,:) = iloveclim_om_c(iloveclim_nix-1,:)
      iloveclim_om_p(iloveclim_nix,:) = iloveclim_om_c(2,:)
      iloveclim_om_ro2(1,:) = iloveclim_om_c(iloveclim_nix-1,:)
      iloveclim_om_ro2(iloveclim_nix,:) = iloveclim_om_c(2,:)



      iloveclim_setup_init = .TRUE.



      RETURN

!-----------------------------------------------------------------------
      END SUBROUTINE setup_iloveclim_medusa_xchange
!-----------------------------------------------------------------------


!=======================================================================      
      END MODULE MOD_ILOVECLIM_SETUP
!=======================================================================      
