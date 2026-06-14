!=======================================================================
!  c06_fftw.f90
!
!  Remplacement open-source (FFTW3) des routines NAG
!     C06FPF : M transformees reelles -> Hermitiennes  (forward)
!     C06FQF : M transformees Hermitiennes -> reelles  (backward/inverse)
!
!  Objectif : substitution "drop-in" dans ECBilt. Memes signatures,
!  meme format de stockage (half-complex), meme normalisation 1/sqrt(N).
!
!  Stockage NAG (par sequence j, longueur N, M sequences entrelacees) :
!     X(j + k*M)  est le k-ieme echantillon de la sequence j
!     -> stride entre echantillons d'une sequence = M
!     -> distance entre sequences                 = 1
!  Ce schema correspond exactement aux parametres FFTW "advanced/guru" :
!     n = N, howmany = M, istride = ostride = M, idist = odist = 1
!
!  Format spectral half-complex (identique NAG <-> FFTW R2HC/HC2R) :
!     positions 0..N/2   : parties reelles   r_0 r_1 ... r_{N/2}
!     positions N/2+1..N-1 : parties imag.   i_{(N-1)/2} ... i_1  (inversees)
!
!  Normalisation : NAG applique 1/sqrt(N) en forward ET en backward.
!  FFTW n'applique rien. On multiplie donc par 1/sqrt(N) dans les deux sens.
!=======================================================================
module c06_fftw_mod
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'

  ! Cache de plans : on evite de recreer un plan a chaque appel.
  ! Indexe sur (M,N). ECBilt n'utilise qu'un seul couple (nlat,nlon),
  ! mais on gere un petit cache par robustesse.
  integer, parameter :: MAXP = 8
  type plancache_t
     integer            :: m = -1, n = -1
     type(C_PTR)        :: pf = C_NULL_PTR   ! plan forward  r2hc
     type(C_PTR)        :: pb = C_NULL_PTR   ! plan backward hc2r
  end type plancache_t
  type(plancache_t), save :: cache(MAXP)
  integer, save :: ncache = 0

contains

  !---------------------------------------------------------------------
  subroutine c06_get_plans(m, n, pf, pb)
    integer, intent(in)      :: m, n
    type(C_PTR), intent(out) :: pf, pb
    integer :: i
    real(C_DOUBLE), allocatable :: scratch(:)
    integer(C_INT) :: rank, nn(1), howmany, istride, idist
    integer(C_FFTW_R2R_KIND) :: kf(1), kb(1)

    do i = 1, ncache
       if (cache(i)%m == m .and. cache(i)%n == n) then
          pf = cache(i)%pf ; pb = cache(i)%pb ; return
       end if
    end do

    ! Creation des plans (une seule fois par couple m,n)
    allocate(scratch(m*n))
    scratch = 0.0_C_DOUBLE
    rank    = 1
    nn(1)   = n
    howmany = m
    istride = m        ! echantillons d'une sequence espaces de M
    idist   = 1        ! sequences contigues (stride 1 entre elles)
    kf(1)   = FFTW_R2HC
    kb(1)   = FFTW_HC2R

    pf = fftw_plan_many_r2r(rank, nn, howmany, &
                            scratch, nn, istride, idist, &
                            scratch, nn, istride, idist, &
                            kf, FFTW_ESTIMATE)
    pb = fftw_plan_many_r2r(rank, nn, howmany, &
                            scratch, nn, istride, idist, &
                            scratch, nn, istride, idist, &
                            kb, FFTW_ESTIMATE)
    deallocate(scratch)

    ncache = ncache + 1
    if (ncache > MAXP) stop 'c06_fftw: plan cache overflow'
    cache(ncache)%m = m ; cache(ncache)%n = n
    cache(ncache)%pf = pf ; cache(ncache)%pb = pb
  end subroutine c06_get_plans

  !---------------------------------------------------------------------
  !  C06FPF : forward, reel -> half-complex.  INIT ignore (plans en cache).
  !---------------------------------------------------------------------
  subroutine C06FPF(M, N, X, INIT, TRIG, WORK, IFAIL)
    integer,          intent(in)    :: M, N
    real(C_DOUBLE),   intent(inout) :: X(M*N)
    character(len=1), intent(in)    :: INIT
    real(C_DOUBLE),   intent(inout) :: TRIG(2*N)
    real(C_DOUBLE),   intent(inout) :: WORK(M*N)
    integer,          intent(inout) :: IFAIL
    type(C_PTR) :: pf, pb
    real(C_DOUBLE) :: scale

    if (M < 1) then ; IFAIL = 1 ; return ; end if
    if (N < 1) then ; IFAIL = 2 ; return ; end if

    call c06_get_plans(M, N, pf, pb)

    ! INIT='i'/'I' : appel d'initialisation pur cote NAG. Les plans sont
    ! deja construits ; on ne touche pas X (chez NAG le resultat sur le
    ! tableau bidon est sans interet). On sort proprement.
    if (INIT == 'i' .or. INIT == 'I') then
       IFAIL = 0 ; return
    end if

    call fftw_execute_r2r(pf, X, X)
    scale = 1.0_C_DOUBLE / sqrt(real(N, C_DOUBLE))
    X = X * scale
    IFAIL = 0
  end subroutine C06FPF

  !---------------------------------------------------------------------
  !  C06FQF : backward, half-complex -> reel. INIT ignore (plans en cache).
  !---------------------------------------------------------------------
  subroutine C06FQF(M, N, X, INIT, TRIG, WORK, IFAIL)
    integer,          intent(in)    :: M, N
    real(C_DOUBLE),   intent(inout) :: X(M*N)
    character(len=1), intent(in)    :: INIT
    real(C_DOUBLE),   intent(inout) :: TRIG(2*N)
    real(C_DOUBLE),   intent(inout) :: WORK(M*N)
    integer,          intent(inout) :: IFAIL
    type(C_PTR) :: pf, pb
    real(C_DOUBLE) :: scale

    if (M < 1) then ; IFAIL = 1 ; return ; end if
    if (N < 1) then ; IFAIL = 2 ; return ; end if

    call c06_get_plans(M, N, pf, pb)

    if (INIT == 'i' .or. INIT == 'I') then
       IFAIL = 0 ; return
    end if

    ! C06FQF (transformee Hermitienne NAG) = HC2R de FFTW avec les
    ! parties imaginaires d'entree niees. Positions imaginaires =
    ! frequences 0-based N/2+1 .. N-1 -> X(k*M+j).
    call c06_negate_imag(M, N, X)
    call fftw_execute_r2r(pb, X, X)
    scale = 1.0_C_DOUBLE / sqrt(real(N, C_DOUBLE))
    X = X * scale
    IFAIL = 0
  end subroutine C06FQF

  !---------------------------------------------------------------------
  subroutine c06_negate_imag(M, N, X)
    integer,        intent(in)    :: M, N
    real(C_DOUBLE), intent(inout) :: X(M*N)
    integer :: j, k
    do k = N/2+1, N-1
       do j = 1, M
          X(k*M + j) = -X(k*M + j)
       end do
    end do
  end subroutine c06_negate_imag

end module c06_fftw_mod
