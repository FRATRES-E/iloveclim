!-----|--1--------2---------3---------4---------5---------6---------7-|
!      This module implements the Insolation Temperature Melt (ITM)
!       parametrisation, here taken from Robinson, 2011, PhD Dissert.
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche, Emiel Spanier
!      Date   : 04 mars 2014
!      Derniere modification : 04 mars 2014
!-----|--1--------2---------3---------4---------5---------6---------7-|

       MODULE ablation_PDD


       IMPLICIT NONE

! --- afq pdd parameters, from Fausto et al.
       REAL, PRIVATE, PARAMETER :: Cice_warm   = 7./1000.
       REAL, PRIVATE, PARAMETER :: Cice_cold   = 15./1000.
       REAL, PRIVATE, PARAMETER :: Csnow_warm  = 3./1000.
       REAL, PRIVATE, PARAMETER :: Csnow_cold  = 3./1000.
       REAL, PRIVATE, PARAMETER :: T_warm      = 10
       REAL, PRIVATE, PARAMETER :: T_cold      = -1
       REAL, PRIVATE, PARAMETER :: Sigma_low   = 1.574
       REAL, PRIVATE, PARAMETER :: Sigma_slope = 0.0012224

       
       CONTAINS

         
!-----|--1--------2---------3---------4---------5---------6---------7-|
       SUBROUTINE daily_melt_pdd(melt_rate,alti_virtlevels,surf_air_temp,precip_snow)

!$ USE OMP_LIB

       REAL,DIMENSION(:,:,:),INTENT(IN) ::  surf_air_temp
       REAL,DIMENSION(:,:,:),INTENT(IN) ::  precip_snow
       REAL,DIMENSION(:),INTENT(IN) ::  alti_virtlevels
! --- afq Melt rate (m.day^-1)
       REAL,DIMENSION(:,:,:),INTENT(OUT) :: melt_rate

! --- afq Local variables
       REAL,DIMENSION(UBOUND(surf_air_temp,dim=1),UBOUND(surf_air_temp,dim=2),UBOUND(alti_virtlevels,dim=1)) :: Cice_3d, Csnow_3d
       REAL,DIMENSION(UBOUND(alti_virtlevels,dim=1)) :: Sigma_1d
       REAL,DIMENSION(UBOUND(alti_virtlevels,dim=1)) :: S22
       REAL,DIMENSION(UBOUND(surf_air_temp,dim=1),UBOUND(surf_air_temp,dim=2),UBOUND(alti_virtlevels,dim=1)) :: pdd
       REAL,DIMENSION(UBOUND(surf_air_temp,dim=1),UBOUND(surf_air_temp,dim=2),UBOUND(alti_virtlevels,dim=1)) :: pds
       REAL,DIMENSION(UBOUND(surf_air_temp,dim=1),UBOUND(surf_air_temp,dim=2),UBOUND(alti_virtlevels,dim=1)) :: pdsi
       REAL,DIMENSION(UBOUND(alti_virtlevels,dim=1)) :: pddct
       REAL,DIMENSION(UBOUND(surf_air_temp,dim=1),UBOUND(surf_air_temp,dim=2),UBOUND(alti_virtlevels,dim=1)) :: SImax

       REAL, PARAMETER :: dtp = 2.   !increment of temperature for PDD integration
       REAL, PARAMETER :: pi_l = 3.1415926
       
       REAL :: temp, summ
       INTEGER :: i,j,nlev,k

       
! --- afq melt coefficients depend on local temperature
!$OMP PARALLEL WORKSHARE
       WHERE (surf_air_temp(:,:,:).GE.T_warm)
          Cice_3d(:,:,:)=Cice_warm
          Csnow_3d(:,:,:)=Csnow_warm
       ELSEWHERE ((T_cold.LE.surf_air_temp(:,:,:)).AND.(surf_air_temp(:,:,:).LT.T_warm))
          Cice_3d(:,:,:)=Cice_warm+((Cice_cold-Cice_warm)/((T_warm-T_cold)**3))*((T_warm-surf_air_temp(:,:,:))**3)
          Csnow_3d(:,:,:)=Csnow_warm
       ELSEWHERE
          Cice_3d(:,:,:)=Cice_cold
          Csnow_3d(:,:,:)=Csnow_cold
       ENDWHERE

! --- afq temperature variability depends on altitude (should be computed only once!)
       Sigma_1d(:) = Sigma_low + Sigma_slope * alti_virtlevels(:)
       S22(:) = 0.5 / Sigma_1d(:) / Sigma_1d(:)

! --- afq integration constant pddct, check the units (daily pdd)
       pddct(:) = dtp / sigma_1d(:) / SQRT(2.*pi_l)
!$OMP END PARALLEL WORKSHARE
       
! --- afq csi, percentage of fresh snow that can refreeze (should be computed only once!)
!       CSI_1d(:)=min(max((alti_virtlevels(:)-800.)*0.000833,0.),1.)
! note afq marion dufresne, refreezing still in GRISLI for now...       

!$OMP PARALLEL
!$OMP DO PRIVATE(summ,temp,k)
       DO nlev = LBOUND(alti_virtlevels,dim=1), UBOUND(alti_virtlevels,dim=1)
          DO i= LBOUND(surf_air_temp,dim=1), UBOUND(surf_air_temp,dim=1)
             DO j= LBOUND(surf_air_temp,dim=2), UBOUND(surf_air_temp,dim=2)
                summ=0.0
                temp=0.0          ! variable d'integration
                k=1
                average_day: DO WHILE ((temp.LE.surf_air_temp(i,j,nlev)+2.5*sigma_1D(nlev)).AND.k.LE.50) ! pdd d'un jour par mois
                   summ=summ+temp*EXP(-((temp-surf_air_temp(i,j,nlev))*(temp-surf_air_temp(i,j,nlev)))*s22(nlev))
                   temp=temp+dtp     ! integration step
                   k= k+1
                END DO average_day
                pdd(i,j,nlev)=summ*pddct(nlev)   ! pdd of only one day
             END DO
          END DO
       END DO
!$OMP END DO       

!$OMP WORKSHARE
       ! pds: number of pdds required to melt the fresh snow
       pds(:,:,:)=precip_snow(:,:,:)/Csnow_3D(:,:,:)
       ! SImax: maximum super imposed ice
       SImax(:,:,:)=precip_snow(:,:,:)*0.6 !CSI_1D(:) afq marion dufresne (csi in grisli is 0.6 in standard vers.)
       ! psdi: number of pdds required to melt the fresh snow and the super imposed ice
       pdsi(:,:,:)=SImax(:,:,:)/Cice_3D(:,:,:)

       WHERE (pdd(:,:,:).LE.pds(:,:,:))
          ! case1: we don't have enough pdd to melt the fresh snow (we try and part refreezes)
          !BM(:,:)=ACC(:,:)-PDD(:,:)*Csnow_2D*(1-CSI_2D(:,:))
          melt_rate(:,:,:) = pdd(:,:,:) * Csnow_3D(:,:,:) * 0.4
       ENDWHERE
       WHERE ((pds(:,:,:).LT.pdd(:,:,:)).AND.(pdd(:,:,:).LE.pds(:,:,:)+pdsi(:,:,:)))
          ! case2: we melt all the fresh snow and part of the superimposed ice
          !BM(:,:)=SIMAX(:,:)-(PDD(:,:)-PDS(:,:))*Cice_2D
          melt_rate(:,:,:) = pds(:,:,:) * Csnow_3d(:,:,:) * 0.4 + (pdd(:,:,:)-pds(:,:,:)) * Cice_3D(:,:,:)
       ENDWHERE
       WHERE (pds(:,:,:)+pdsi(:,:,:).LE.pdd(:,:,:))
          ! case3: we melt all the fresh snow + part of the old ice
          !BM(:,:)=(PDS(:,:)+PDSI(:,:)-PDD(:,:))*Cice_2D
          melt_rate(:,:,:) = pds(:,:,:)*Csnow_3d(:,:,:) + (pdd(:,:,:)-pds(:,:,:)-pdsi(:,:,:)) * Cice_3D(:,:,:)
       ENDWHERE
!$OMP END WORKSHARE
!$OMP END PARALLEL

!       write(*,*) "T profile"
!       write(*,*) surf_air_temp(2,24,:)
       
!       write(*,*) "snow profile"
!       write(*,*) precip_snow(2,24,:)
       
!       write(*,*) "PDD profile"
!       write(*,*) pdd(2,24,:)  
       
       RETURN

       END SUBROUTINE daily_melt_pdd

!-----|--1--------2---------3---------4---------5---------6---------7-|

       END MODULE ablation_PDD
