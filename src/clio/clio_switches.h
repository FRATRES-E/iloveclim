! started the coding of a new generation NetCDF output module
! EXPERIMENTAL, expect rough edges!
#define CLIO_OUT_NEWGEN 0

!   to obtain additional printing conditioned on the first iteration
!      (1 = activated)
#define CTRL_FIRST_ITER 0

!   replacement of ltest in the CLIO param file in the future
!   for now, check that actual ltest is the same value as this one
#define L_TEST 3

!   replacement for the nitrap parameter (restoring fluxes)
!   also to check that the value here is identical to run.param
#define NIT_RAP 0

! 0 CLIO standard vertical mixing
! 1 use energy-constrained 3D map of tidal diffusivity
!   /!\ This option requires file 'tidal_mixing.nc' in inputdata/clio.
!   /!\ This option also needs avkb and avnub to be reduced to zero in run.param at levels deeper than 300 m.
!
#define TIDEMIX 0

#define XSLOP 1
#define I_COUPL 1

#define forced_winds 0
