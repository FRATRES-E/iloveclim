!=======================================================================
!  ec_fourier_mod.f90
!
!  Transformees de Fourier multiples (advection spectrale ECBilt) sur la
!  base FFTW3. Remplace l'interface heritee NAG (C06FPF / C06FQF) par une
!  interface Fortran moderne, sans les arguments fantomes TRIG/WORK ni le
!  INIT/IFAIL de NAG.
!
!  Operations (M = nlat sequences de longueur N = nlon, grille agg(nlat,nlon)) :
!     ec_fourier_grid2spec(agg) : grille -> coefficients (ancien C06FPF)
!     ec_fourier_spec2grid(agg) : coefficients -> grille (ancien C06FQF)
!
!  Le stockage "half-complex" dans agg est INCHANGE (colonnes mr=m+1 pour
!  les parties reelles, mi=nlon+1-m pour les imaginaires) : il est utilise
!  tel quel par les transformees de Legendre de atmdyn_mod. Ce module ne
!  fait que la partie Fourier, exactement comme avant, mais proprement.
!
!  Equivalences (validees a ~1e-16, cf. banc de non-regression) :
!     grid2spec = FFTW_R2HC * (1/sqrt(N))
!     spec2grid = FFTW_HC2R * (1/sqrt(N)), parties imaginaires niees en entree
!  Stockage NAG entrelace X(j+k*M) <-> FFTW many: stride=M, dist=1, howmany=M.
!
!  Planification : FFTW_MEASURE, avec "wisdom" persiste par resolution dans
!  un repertoire (defaut run/inputdata/FFTW3). Si le wisdom n'existe pas, il
!  est mesure puis ecrit ; sinon il est recharge (planification quasi
!  instantanee). Cf. ec_fourier_init.
!
!  Compatible -fdefault-real-8 : tous les reels sont real(C_DOUBLE) (kind 8).
!=======================================================================
module ec_fourier_mod
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'
  private

  public :: ec_fourier_init, ec_fourier_final
  public :: ec_fourier_grid2spec, ec_fourier_spec2grid

  ! --- etat du module ---
  logical,        save :: initialised = .false.
  integer,        save :: mlat = -1, nlon = -1
  type(C_PTR),    save :: plan_fwd = C_NULL_PTR   ! R2HC
  type(C_PTR),    save :: plan_bwd = C_NULL_PTR   ! HC2R
  real(C_DOUBLE), save :: rnorm = 1.0_C_DOUBLE    ! 1/sqrt(N)

contains

  !---------------------------------------------------------------------
  !  Initialisation : a appeler une fois (par resolution) avant les
  !  transformees. m = nlat, n = nlon. wisdom_dir optionnel.
  !---------------------------------------------------------------------
  subroutine ec_fourier_init(m, n, wisdom_dir)
    integer,          intent(in)           :: m, n
    character(len=*), intent(in), optional :: wisdom_dir
    real(C_DOUBLE), allocatable :: scratch(:)
    integer(C_INT)           :: rank, nn(1), howmany, istride, idist, planflag
    integer(C_FFTW_R2R_KIND) :: kf(1), kb(1)
    character(len=512) :: wdir, wfile
    logical :: have_wisdom

    ! Re-init eventuelle (changement de resolution) : on detruit l'ancien etat.
    if (initialised) call ec_fourier_final()

    mlat = m ; nlon = n
    rnorm = 1.0_C_DOUBLE / sqrt(real(n, C_DOUBLE))

    wdir = 'run/inputdata/FFTW3'
    if (present(wisdom_dir)) wdir = trim(wisdom_dir)

    ! Nom de fichier wisdom specifique a la resolution (m x n) : un jeu de
    ! plans par couple (nlat,nlon), conformement au besoin multi-resolution.
    write(wfile,'(A,"/wisdom_",I0,"x",I0,".dat")') trim(wdir), m, n

    have_wisdom = load_wisdom(trim(wfile))
    if (have_wisdom) then
       planflag = FFTW_MEASURE + FFTW_WISDOM_ONLY  ! reutilise le plan mesure
    else
       planflag = FFTW_MEASURE
    end if

    allocate(scratch(m*n)) ; scratch = 0.0_C_DOUBLE
    rank = 1 ; nn(1) = n ; howmany = m
    istride = m ; idist = 1
    kf(1) = FFTW_R2HC ; kb(1) = FFTW_HC2R

    plan_fwd = fftw_plan_many_r2r(rank, nn, howmany, &
                                  scratch, nn, istride, idist, &
                                  scratch, nn, istride, idist, &
                                  kf, planflag)
    plan_bwd = fftw_plan_many_r2r(rank, nn, howmany, &
                                  scratch, nn, istride, idist, &
                                  scratch, nn, istride, idist, &
                                  kb, planflag)

    ! Si WISDOM_ONLY a echoue (plan nul faute de wisdom complet), on
    ! retombe sur une vraie mesure puis on sauvegarde.
    if (.not. c_associated(plan_fwd) .or. .not. c_associated(plan_bwd)) then
       if (c_associated(plan_fwd)) call fftw_destroy_plan(plan_fwd)
       if (c_associated(plan_bwd)) call fftw_destroy_plan(plan_bwd)
       plan_fwd = fftw_plan_many_r2r(rank, nn, howmany, &
                                     scratch, nn, istride, idist, &
                                     scratch, nn, istride, idist, &
                                     kf, FFTW_MEASURE)
       plan_bwd = fftw_plan_many_r2r(rank, nn, howmany, &
                                     scratch, nn, istride, idist, &
                                     scratch, nn, istride, idist, &
                                     kb, FFTW_MEASURE)
       have_wisdom = .false.
    end if
    deallocate(scratch)

    if (.not. c_associated(plan_fwd) .or. .not. c_associated(plan_bwd)) &
       stop 'ec_fourier_init: echec de creation des plans FFTW'

    ! Sauvegarder le wisdom nouvellement mesure (best-effort).
    if (.not. have_wisdom) call save_wisdom(trim(wdir), trim(wfile))

    initialised = .true.
  end subroutine ec_fourier_init

  !---------------------------------------------------------------------
  subroutine ec_fourier_final()
    if (c_associated(plan_fwd)) call fftw_destroy_plan(plan_fwd)
    if (c_associated(plan_bwd)) call fftw_destroy_plan(plan_bwd)
    plan_fwd = C_NULL_PTR ; plan_bwd = C_NULL_PTR
    initialised = .false.
  end subroutine ec_fourier_final

  !---------------------------------------------------------------------
  !  Grille -> coefficients de Fourier (ancien C06FPF). In-place sur agg.
  !---------------------------------------------------------------------
  subroutine ec_fourier_grid2spec(agg)
    real(C_DOUBLE), intent(inout) :: agg(mlat*nlon)
    call check_ready('ec_fourier_grid2spec')
    call fftw_execute_r2r(plan_fwd, agg, agg)
    agg = agg * rnorm
  end subroutine ec_fourier_grid2spec

  !---------------------------------------------------------------------
  !  Coefficients de Fourier -> grille (ancien C06FQF). In-place sur agg.
  !  La transformee Hermitienne NAG = HC2R avec parties imaginaires niees.
  !---------------------------------------------------------------------
  subroutine ec_fourier_spec2grid(agg)
    real(C_DOUBLE), intent(inout) :: agg(mlat*nlon)
    integer :: j, k
    call check_ready('ec_fourier_spec2grid')
    ! negation des positions imaginaires (frequences 0-based nlon/2+1..nlon-1)
    do k = nlon/2+1, nlon-1
       do j = 1, mlat
          agg(k*mlat + j) = -agg(k*mlat + j)
       end do
    end do
    call fftw_execute_r2r(plan_bwd, agg, agg)
    agg = agg * rnorm
  end subroutine ec_fourier_spec2grid

  !====================== utilitaires prives ==========================

  subroutine check_ready(who)
    character(len=*), intent(in) :: who
    if (.not. initialised) then
       write(*,'(A)') 'ec_fourier_mod: '//who//' appele avant ec_fourier_init'
       stop 1
    end if
  end subroutine check_ready

  ! Charge le wisdom depuis un fichier. .true. si succes.
  logical function load_wisdom(path) result(ok)
    character(len=*), intent(in) :: path
    integer(C_INT) :: ret
    logical :: exists
    ok = .false.
    inquire(file=path, exist=exists)
    if (.not. exists) return
    ret = fftw_import_wisdom_from_filename(path//C_NULL_CHAR)
    ok = (ret /= 0_C_INT)
  end function load_wisdom

  ! Sauve le wisdom courant. Cree le repertoire au besoin (best-effort).
  subroutine save_wisdom(dir, path)
    character(len=*), intent(in) :: dir, path
    integer(C_INT) :: ret
    integer        :: ios
    logical        :: exists
    ! creation du repertoire : on tente, sans faire echouer le run si KO.
    inquire(file=dir, exist=exists)
    if (.not. exists) call execute_command_line('mkdir -p '//dir, wait=.true., exitstat=ios)
    ret = fftw_export_wisdom_to_filename(path//C_NULL_CHAR)
    ! ret==0 -> echec d'ecriture : non fatal (on aura juste a re-mesurer).
  end subroutine save_wisdom

end module ec_fourier_mod
