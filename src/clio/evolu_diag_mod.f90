!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!   Copyright 2026 iLOVECLIM / LUDUS coding group

!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at

!       http://www.apache.org/licenses/LICENSE-2.0

!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#include "choixcomposantes.h"

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: [evolu_diag_mod]
!
!>     @brief Chronological diagnostics for the CLIO ocean: builds the "evolu" output file.
!
!      DESCRIPTION:
!>     This module merges the former fixed-form routines informe.f (header/registration) and inforun.f (per-timestep
!!     computation) into a single procedural module. Merging removes a long-standing hazard: column NAMES (titvar) and
!!     column VALUES (vinfor) used to be produced in two separate files, each with its own hand-maintained running index
!!     "nv". Nothing guaranteed that titvar(k) described the quantity stored in vinfor(k); a single mismatched increment
!!     silently mislabelled the scientific output. Here a single in-module registry is the one source of truth for the
!!     column ordering, and inforun writes values through idx("NAME"), so a name and its value can no longer drift apart.
!!
!!       informe(nn99)  : registers every diagnostic (fixes the column order once), computes averaging volumes,
!!                        opens "evolu" and writes its header. Called once at init, exactly as before.
!!       inforun(...)   : per-timestep computation. [ Temps 2 -- not in this file yet ]
!!
!>     @note dmr&clo -- Temps 1 (this file) is a CONSTANT-OUTPUT refactor: the "evolu" file must remain bit-for-bit
!!           identical to the legacy code for the default configuration (nsmax=2, FRAZER_ARCTIC=0, OCYCC=0, ISM=0,
!!           ISOOCN=0). Physical-constant substitutions from global_constants_mod (tK_zero_C, vol_mass_dens_wat, ...)
!!           are only FLAGGED in comments here, never applied, because they could alter the last bit. They belong to a
!!           later, separately-validated "constant harmonisation" step.
!>     @note dmr&clo -- The #if(0) blocks and the FRAZER_ARCTIC branches of the legacy inforun.f are NOT part of Temps 1.
!!           They are reintegrated later, each as its own chantier, behind explicit CPP guards, with their bug fixes.
!>     @note dmr&clo -- Legacy hand-jump #1 "nv = nv0 + nvhsf + 1" after the three HSF columns is a no-op: geogra.f sets
!!           nvhsf = 3 unconditionally, so the jump lands exactly where three plain register() calls already leave us.
!!           Dropped; nvhsf stays USE-associated (owned by geogra) and is only cross-checked here.
!>     @note dmr&clo -- Legacy hand-jump #2 "nv = nv0 + 3*nsmax + 1" in the T/S/A/B block WROTE 12 titles then rewound
!!           the counter, so in the default config the last 6 (A,A1-o,|A-o|,B,B1-o,|B-o|) were written and then OVER-
!!           WRITTEN by |w|,|u|,|v|,K.E,... and never appeared in evolu. Verified against a real evolu header: the
!!           surviving block-C columns are exactly T-c,T1-o,|T-o|,S-30,S1-o,|S-o| (3 per tracer x nsmax=2). We generate
!!           only those, via a loop over ns, and the "+1" of the legacy formula is recognised as an indexing artefact
!!           (from a stray nv=nv+1 before T-c), NOT a real column. This removes the silent-overwrite bug at its root.
!>     @note dmr&clo -- OUT OF SCOPE (known, deliberately untouched): a given nsmax value does not map to a unique
!!           tracer combination in CLIO. Naming tracers 3+ correctly is a CLIO-core design question, not a diagnostics
!!           one. Until it is resolved, block C hard-stops for nsmax>2 rather than emit ambiguous/overwritten columns.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module evolu_diag_mod

       use global_constants_mod, only: dblp=>dp, ip
       use para_mod,             only: nchsep_p => nchsep

       implicit none

       private

! dmr&clo --- public contract: the two legacy entry points, unchanged from the callers' point of view.
       public :: informe, inforun

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Averaging policy per column, replacing the scattered "if (ntmoy.ge.1) ktsum=1" / "if (ntmoy.eq.2) ..." lines
!           of the legacy. register() resolves the policy against ntmoy at registration time, keeping the averaging
!           intent LOCAL to each column. Two legacy patterns remain genuinely block-wide (a retroactive do-loop over a
!           whole block when ntmoy>=1 or ==2); those are reproduced explicitly where they occur, not folded in here.
!             AVG_NEVER  : ktsum stays 0                                 [legacy: no ktsum line]
!             AVG_ALWAYS : ktsum = 1 unconditionally                     [legacy: ktsum(nv) = 1]
!             AVG_TMOY1  : ktsum = 1 iff ntmoy >= 1                      [legacy: if (ntmoy.ge.1) ktsum = 1]
!             AVG_TMOY2  : ktsum = 1 iff ntmoy == 2                      [legacy: if (ntmoy.eq.2) ktsum = 1]
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       integer(ip), parameter :: AVG_NEVER  = 0
       integer(ip), parameter :: AVG_ALWAYS = 1
       integer(ip), parameter :: AVG_TMOY1  = 2
       integer(ip), parameter :: AVG_TMOY2  = 3

! dmr&clo --- Direction flag for the generic strait/section scalar-flux kernel (Option B). Used in Temps 2.
       integer(ip), parameter :: DIR_MERID = 1
       integer(ip), parameter :: DIR_ZONAL = 2

! dmr&clo --- Flux mode selector for face_scalar_flux (see kernel header).
       integer(ip), parameter :: FLUX_TOTAL      = 1
       integer(ip), parameter :: FLUX_ADV        = 2
       integer(ip), parameter :: FLUX_DIF        = 3
       integer(ip), parameter :: FLUX_ADV_LEGACY = 4     ! dmr&clo ANOMALY, marked for removal (reserve R1)

! dmr&clo --- persistent state for the Atlantic salt-content tendency (dtsaltA), formerly a save local in inforun
!             (reserve R6: move to restart if inter-segment continuity is ever required).
       real(dblp), save :: saltatlantic_old_sav = 0._dblp
       logical,    save :: first_dtsalt         = .true.

! dmr&clo --- module-private running counter for the registry; the single writer of the column ordering.
       integer(ip), save :: nv_reg = 0

! dmr&clo --- IMMUTABLE canonical name table (design decision "alpha"). idx("NAME") resolves against THIS, never against
!             titvar. Rationale: legacy informe rewrites titvar after the header (prefixes '  ' and truncates to
!             nchinf=5 chars), so titvar cannot serve as a durable lookup table. reg_name keeps the full names, is set
!             once by register(), and is never mutated -- so idx() works at any point in the run and is insensitive to
!             nchinf. This preserves the by-name design for inforun while leaving titvar's legacy behaviour untouched
!             (header stays bit-identical; no side effect on other readers of titvar such as moyen/binout).
       character(len=nchsep_p), dimension(:), allocatable, save :: reg_name

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   REGISTRY HELPERS (private)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! dmr&clo --- register one diagnostic: append its name, set its averaging flag from the policy + current ntmoy, return
!             the attributed index. This is the ONLY place that advances the column counter.
       integer(ip) function register(name, policy) result(iv)

         use bloc_mod,         only: ntmoy
         use para_mod,         only: ninfmx
         use reper_mod,        only: titvar, ktsum
         use newunit_clio_mod, only: clio3_out_id

         character(len=*), intent(in) :: name
         integer(ip),      intent(in) :: policy

         nv_reg = nv_reg + 1
         iv     = nv_reg

         if (iv > ninfmx) then
           write(clio3_out_id,*) 'STOP evolu_diag: registry overflow, ninfmx = ', ninfmx, ' name = ', name
           error stop 'evolu_diag: registry overflow'
         endif

         titvar(iv)   = name
         reg_name(iv) = name          ! dmr&clo canonical, immutable copy for idx() (see module-level note "alpha")

         select case (policy)
           case (AVG_ALWAYS)
             ktsum(iv) = 1
           case (AVG_TMOY1)
             ktsum(iv) = merge(1, 0, ntmoy >= 1)
           case (AVG_TMOY2)
             ktsum(iv) = merge(1, 0, ntmoy == 2)
           case default            ! AVG_NEVER
             ktsum(iv) = 0
         end select

       end function register

! dmr&clo --- resolve a registered name to its column index; hard-stop on an unknown name (catches typos at once).
!             Searches the IMMUTABLE reg_name table (not titvar), so it works at any point in the run regardless of the
!             legacy post-header rewrite of titvar. Linear search is fine: inforun resolves each name into a local
!             integer once, outside the timestep loop.
       integer(ip) function idx(name) result(iv)

         use newunit_clio_mod, only: clio3_out_id

         character(len=*), intent(in) :: name
         integer(ip) :: k

         if (.not. allocated(reg_name)) then
           write(clio3_out_id,*) 'STOP evolu_diag: idx() called before informe() registered the columns. name = ', name
           error stop 'evolu_diag: idx before registration'
         endif

         do k = 1, nv_reg
           if (trim(reg_name(k)) == trim(name)) then
             iv = k
             return
           endif
         enddo

         write(clio3_out_id,*) 'STOP evolu_diag: unknown diagnostic name in idx() = ', name
         error stop 'evolu_diag: unknown diagnostic name'

       end function idx

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   ENTRY POINT 1 -- informe: register the column order, compute averaging volumes, open + header "evolu".
!           Behaviour identical to legacy informe.f for the default configuration. This is where the column ORDER is
!           decided, once and for all, through the register() sequence below.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subroutine informe(nn99)

         use const_mod,        only: epsil
         use para0_mod,        only: imax, jmax, kmax, nsmax
         use para_mod,         only: nchinf, nbsmax, nchsep, ninfmx
         use bloc0_mod,        only: dx, dy, js1, js2, ks1, ks2, nitrun, spvr    &
                                   , is1, iu1, dz, dzw, z, zw, scalr, is2, iu2    &
                                   , daeta, fqajc, fss
         use bloc_mod,         only: ninfo, nstart, ntmoy, refexp, zsurf, zvolo  &
                                   , zvols, zvolv, zvolw, ctmi, zsurfs, zsurfo, zsurfv
         use ice_mod,          only: vwx
         use reper_mod,        only: dlat, fmtw, jmnor, jmsud, knor, ksud, nferme &
                                   , nocean, ylat1, zmdeta, nvinfo, nvhsf, titvar &
                                   , ktsum, scalwr, vinfom
         use newunit_clio_mod, only: clio3_out_id, mouchard_id, evolu_id, testevolu_id

         implicit none

         integer(ip), intent(in) :: nn99

! dmr&clo --- local variables
         character(len=nchinf) :: titinf
         character(len=8)      :: fmtinf, cc8
         character(len=30)     :: fmtr, fmtitr

         integer(ip) :: i, j, k, kz, n, nfrinf, nn, nn0, nninfo, nnntot, ns, nv, nv0
         integer(ip) :: run_param_id

         real(dblp), dimension(kmax) :: surfs, surfo, surfv
         real(dblp) :: ddxy, volo, vols, volv, volw, xx1, xxx, yylatn, yylats, zlim1

 1000    format(A30,1x,E13.6,1x,I6,1x,I6,1x,I6,1x,I6)
 1111    format(3(F7.1,1X,F7.3,1X),I3,A)

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1 ) Initialisation, ouverture et ecriture de l'entete du fichier.
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

! dmr&clo --- ktsum is now set per-column by register(); we still clear the whole array so any unregistered tail
!             stays 0 (defensive, matches legacy pre-zeroing).
         do nv = 1, ninfmx
           ktsum(nv) = 0
         enddo

! dmr&clo --- reset the registry counter and (re)allocate the immutable canonical-name table (re-entry safe).
         nv_reg = 0
         if (allocated(reg_name)) deallocate(reg_name)
         allocate(reg_name(ninfmx))
         reg_name(:) = ''

!--facteur d'echelle pour "D.eta" (unite = 10e-6 m/s) :
         zmdeta = 1.D+6

!--lecture des 3 1eres lignes du fichier "run.param" :
         open(newunit=run_param_id, file='run.param', status='old')
         read(run_param_id,'(A6)')       refexp
         read(run_param_id,*)
         read(run_param_id,'(A8,1X,I3)') fmtinf, nfrinf
         close(run_param_id)

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! dmr&clo --- REGISTRATION OF THE COLUMN ORDER.
!             One register() call per column, in the exact legacy order. Returned indices are ignored (order is what
!             matters); inforun fetches what it needs via idx("NAME"). Anchors (nv0, nvhsf, nocean) are captured from
!             nv_reg as the sequence passes them.
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!-- block A : iteration / time / global scalars (legacy columns 1..6) -----------------------------------------------
         nv = register('NoIt' , AVG_NEVER )     ! 1
         nv = register('T yr' , AVG_TMOY1 )     ! 2
         nv = register('EgAjC', AVG_ALWAYS)     ! 3
         nv = register('V_AjC', AVG_ALWAYS)     ! 4
         nv = register('D.Eta', AVG_ALWAYS)     ! 5
         nv = register('M.Eta', AVG_TMOY2 )     ! 6   (legacy alt name 'EgEta' kept in source history)

!-- block B : HSF heat fluxes + basin budgets (legacy columns 7..32 in default config) ------------------------------
         nv0 = nv_reg                            ! legacy nv0 anchor (start of the ntmoy>=1 block closed below)

         nv = register('DrHSF', AVG_NEVER)      ! 7   Drake     heat flux
         nv = register('InHSF', AVG_NEVER)      ! 8   Indonesia heat flux
         nv = register('BeHSF', AVG_NEVER)      ! 9   Bering    heat flux
! dmr&clo --- legacy did "nv = nv0 + nvhsf + 1" here; with nvhsf==3 (geogra) it is a no-op and is dropped.
!             Cross-check the invariant so a future change to nvhsf is caught loudly instead of silently shifting cols.
         if (nv_reg /= nv0 + nvhsf) then
           write(clio3_out_id,*) 'STOP evolu_diag: HSF block size /= nvhsf. nv_reg,nv0,nvhsf =', nv_reg, nv0, nvhsf
           error stop 'evolu_diag: HSF/nvhsf mismatch'
         endif

         nv = register('DANS', AVG_NEVER)       ! 10
         nv = register('DANN', AVG_NEVER)       ! 11
         nv = register('ISCS', AVG_NEVER)       ! 12
         nv = register('ISCN', AVG_NEVER)       ! 13
         nv = register('FRAS', AVG_NEVER)       ! 14
         nv = register('FRAN', AVG_NEVER)       ! 15
         nv = register('PNWS', AVG_NEVER)
         nv = register('PNWN', AVG_NEVER)
#if ( FRAZER_ARCTIC == 1 )
! dmr&clo --- Temps 2 reintegration point: Frazer Davies basin columns. Guarded, off by default.
         nv = register('BASS', AVG_NEVER)
         nv = register('BASN', AVG_NEVER)
#endif
         nv = register('ADGIN', AVG_NEVER)
         nv = register('ADPro', AVG_NEVER)
         nv = register('ADOut', AVG_NEVER)      ! 20
         nv = register('AABpr', AVG_NEVER)
         nv = register('AABex', AVG_NEVER)
         nv = register('AABat', AVG_NEVER)
         nv = register('Fc30A', AVG_NEVER)
         nv = register('Fs30A', AVG_NEVER)
! dmr&clo --- names are stored in FULL (legacy form); the header format (A5, nchinf=5) truncates them on WRITE to
!             'Fsa30','Fsd30','dtsal','Mov30' as seen in real output. idx() matches the full name via reg_name, so
!             truncation on display never breaks resolution.
         nv = register('Fsa30A', AVG_NEVER)     ! 26  (header shows truncated 'Fsa30')
         nv = register('Fsd30A', AVG_NEVER)     !     (header 'Fsd30')
         nv = register('FsaNA', AVG_NEVER)
         nv = register('FsdNA', AVG_NEVER)
         nv = register('dtsaltA', AVG_NEVER)    !     (header 'dtsal')
         nv = register('Mov30A', AVG_NEVER)     !     (header 'Mov30')
         nv = register('Fsber', AVG_NEVER)      ! 32
#if ( FRAZER_ARCTIC == 1 )
! dmr&clo --- Temps 2 reintegration point: Frazer Davies strait-flux columns.
         nv = register('BSFF', AVG_NEVER)
         nv = register('BSFP', AVG_NEVER)
         nv = register('BSFN', AVG_NEVER)
         nv = register('FSFF', AVG_NEVER)
         nv = register('FSFP', AVG_NEVER)
         nv = register('FSFN', AVG_NEVER)
         nv = register('BAFF', AVG_NEVER)
         nv = register('BAFP', AVG_NEVER)
         nv = register('BAFN', AVG_NEVER)
         nv = register('ADVE', AVG_NEVER)
         nv = register('DIFF', AVG_NEVER)
         nv = register('AOW',  AVG_NEVER)
         nv = register('AOE',  AVG_NEVER)
         nv = register('AOFC', AVG_NEVER)
         nv = register('CAFF', AVG_NEVER)
         nv = register('CAFP', AVG_NEVER)
         nv = register('CAFN', AVG_NEVER)
         nv = register('KSFF', AVG_NEVER)
         nv = register('KSFP', AVG_NEVER)
         nv = register('KSFN', AVG_NEVER)
#endif

! dmr&clo --- legacy: "if (ntmoy.ge.1) then; do nn=1+nv0,nv; ktsum(nn)=1; enddo; endif"  (block-B wide, reproduced).
         if (ntmoy >= 1) then
           do nn = 1 + nv0, nv_reg
             ktsum(nn) = 1
           enddo
         endif

!-- block C : basin-mean scalars, 3 columns per tracer (T-c/T1-o/|T-o|, S-30/S1-o/|S-o|, ...) -----------------------
! dmr&clo --- See module header note on legacy hand-jump #2. We emit exactly 3*nsmax columns (the surviving ones),
!             NOT the 12 written-then-overwritten legacy names. Bit-verified against a real evolu header for nsmax=2.
! dmr&clo --- OUT OF SCOPE guard: naming of tracers 3+ is a CLIO-core question (nsmax does not map 1:1 to a tracer set).
         if (nsmax > 2) then
           write(clio3_out_id,*) 'STOP evolu_diag: block C only defined for nsmax<=2 (T,S). nsmax = ', nsmax
           write(clio3_out_id,*) '  Naming of tracers 3+ is a CLIO-core design question, deliberately out of scope.'
           error stop 'evolu_diag: block C undefined for nsmax>2'
         endif

         nv0 = nv_reg                            ! anchor for the ntmoy==2 block closed after block D

         do ns = 1, nsmax
           if (ns == 1) then
             nv = register('T-c' , AVG_NEVER)   ! basin-mean pot. temperature, offset by scalwr(1)
             scalwr(1) = -273.15                ! dmr&clo TODO Temps 3: -> -tK_zero_C. Not now (bit-exact).
             nv = register('T1-o', AVG_NEVER)   ! level-1 value
             nv = register('|T-o|', AVG_NEVER)  ! |anomaly|
           else                                  ! ns == 2
             nv = register('S-30', AVG_NEVER)   ! basin-mean salinity, offset by scalwr(2)
             scalwr(2) = -30.
             nv = register('S1-o', AVG_NEVER)
             nv = register('|S-o|', AVG_NEVER)
           endif
         enddo

!-- block D : |w| |u| |v| K.E and the per-level anomaly families ----------------------------------------------------
         nv = register('|w|', AVG_TMOY2)
         nv = register('|u|', AVG_TMOY2)
         nv = register('|v|', AVG_TMOY2)
         nv = register('K.E', AVG_TMOY2)

         if (ninfo >= 0) then
!- Moyenne des valeurs absolues des differences :
           do k = ks1, ks2
             write(cc8,'(A3,I2)') 'T-o', k
             nv = register(cc8, AVG_NEVER)
           enddo
           do k = ks1, ks2
             write(cc8,'(A3,I2)') 'S-o', k
             nv = register(cc8, AVG_NEVER)
           enddo
         else
           write(cc8,'(A3,I2)') 'T-o', 1
           nv = register(cc8, AVG_NEVER)
!- Moyenne des differences :
           do k = 1 + ks1, ks2
             write(cc8,'(A1,I2,A2)') 'T', k, '-o'
             nv = register(cc8, AVG_NEVER)
           enddo
           write(cc8,'(A3,I2)') 'S-o', 1
           nv = register(cc8, AVG_NEVER)
           do k = 1 + ks1, ks2
             write(cc8,'(A1,I2,A2)') 'S', k, '-o'
             nv = register(cc8, AVG_NEVER)
           enddo
         endif

! dmr&clo --- legacy: "if (ntmoy.eq.2) then; do nn=1+nv0,nv; ktsum(nn)=1; enddo; endif"  (blocks C+D wide).
         if (ntmoy == 2) then
           do nn = 1 + nv0, nv_reg
             ktsum(nn) = 1
           enddo
         endif

!-- block E : sea-ice / iceberg / thermal-expansion tail -----------------------------------------------------------
         nocean = nv_reg                         ! legacy nocean anchor

         nv = register('AIEFN', AVG_NEVER)
         nv = register('AIEFS', AVG_NEVER)
         nv = register('A15N', AVG_NEVER)
         nv = register('A15S', AVG_NEVER)
         nv = register('A85N', AVG_NEVER)
         nv = register('A85S', AVG_NEVER)
         nv = register('ALEN', AVG_NEVER)
         nv = register('ALES', AVG_NEVER)
         nv = register('VOLN', AVG_NEVER)
         nv = register('VOLS', AVG_NEVER)
         nv = register('VONN', AVG_NEVER)
         nv = register('VONS', AVG_NEVER)
         nv = register('ECGN', AVG_NEVER)
         nv = register('ECGS', AVG_NEVER)
         nv = register('FRAG', AVG_NEVER)
#if ( FRAZER_ARCTIC == 1 )
         nv = register('FSIF', AVG_NEVER)
#endif
         nv = register('SPNG', AVG_NEVER)
#if ( FRAZER_ARCTIC == 1 )
         nv = register('KSIF', AVG_NEVER)
#endif
         nv = register('BERG', AVG_NEVER)
#if ( FRAZER_ARCTIC == 1 )
         nv = register('BSIF', AVG_NEVER)
         nv = register('CAIF', AVG_NEVER)
         nv = register('BAIF', AVG_NEVER)
#endif
         nv = register('ThEx', AVG_NEVER)
         nv = register('ISMM', AVG_NEVER)
         nv = register('IcbN', AVG_NEVER)
         nv = register('IcbS', AVG_NEVER)

! dmr&clo --- legacy: "if (ntmoy.ge.1) then; do nn=1+nocean,nv; ktsum(nn)=1; enddo; endif"
         if (ntmoy >= 1) then
           do nn = 1 + nocean, nv_reg
             ktsum(nn) = 1
           enddo
         endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

         nvinfo = nv_reg
         if (nvinfo > ninfmx) then
           write(clio3_out_id,*) 'Arret ! Depassement Nombre Max de variables traitees par la routine "informe"'
           error stop 'evolu_diag: nvinfo > ninfmx'
         endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--Definition et Ecriture de l'entete :
!- nombre d'enregistrements :
         nninfo = abs(ninfo)
         nn0    = (nstart - 1) / nninfo
         if (nstart == 1) nn0 = -1
         nferme = (nstart - 1 + nitrun) / nninfo
         nnntot = nferme - nn0
         nn0    = nninfo * (1 + nn0)
         nferme = nninfo * nferme

!- definition des formats :
         write(fmtw,'(A,I3,A,I1,A)')   '(', nfrinf, '(A', nchsep, ','//fmtinf//'))'
         write(fmtr,'(A,I3,A,I1,A)')   '(', nfrinf, '(', nchsep, 'X,'//fmtinf//'))'
         write(fmtitr,'(A,I3,A,I1,A)') '(', ninfmx, 'A', nchinf, ')'

!--Ouverture du fichier "evolu" :
         open(newunit=evolu_id, file='outputdata/ocean/evolu', status='unknown')
#if ( ISM >= 2 )
         open(newunit=testevolu_id, file='outputdata/ocean/testevolu', form='formatted', status='unknown')
#endif
!- ecriture de 2 lignes d entete :
         write(evolu_id,1000) fmtr, spvr, nvinfo, nnntot, 0, nfrinf
#if ( ISM >= 2 )
         write(testevolu_id,1000) fmtr, spvr, nvinfo, nnntot, 0, nfrinf
#endif
         xxx = 0.001 * DFLOAT(nninfo)
         xx1 = 0.001 * DFLOAT(nn0)
         write(cc8,'(A,I1)') ' ntmoy=', ntmoy
         write(evolu_id,1111) DFLOAT(nchinf), 0., xx1, xxx, 0., 0., 0, cc8
#if ( ISM >= 2 )
         write(testevolu_id,1111) DFLOAT(nchinf), 0., xx1, xxx, 0., 0., 0, cc8
#endif

!- ecriture de 2 lignes de titre :
         write(evolu_id,'(A,I8,A,I8,A,I5)') 'Evolution chronologique - Experience '//refexp  &
                                          //'   de', nn0, ' a', nferme, ' pas', nninfo
         write(evolu_id,fmtitr) (titvar(nv), nv = 1, nvinfo)
#if ( ISM >= 2 )
         write(testevolu_id,'(A,I8,A,I8,A,I5)') 'Evolution chronologique - Experience '//refexp  &
                                              //'   de', nn0, ' a', nferme, ' pas', nninfo
         write(testevolu_id,fmtitr) (titvar(nv), nv = 1, nvinfo)
#endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--preparation de "titvar" pour l ecriture parmi les valeurs numeriques :
         do nv = 2, nvinfo
           titinf     = titvar(nv)(:nchinf)
           titvar(nv) = '  '//titinf
         enddo

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! initialisation of vwx
         do j = 1, 57
           do k = 1, kmax + 1
             do n = 0, nbsmax
               vwx(j,k,n) = 0.0
             enddo
           enddo
         enddo
         write(mouchard_id,*) 'vwx', vwx(43,1,3)

!--calcul des differents "volumes" intervenant dans les moyennes :
         ddxy = dx * dy
         vols = 0.
         volo = 0.
         volv = 0.
         volw = 0.
         do k = ks1, ks2
           surfs(k) = 0.
           surfo(k) = 0.
           surfv(k) = 0.
           do j = js1, js2
             do i = is1(j), is2(j)
               surfs(k) = surfs(k) + ctmi(i,j,k,0)
               surfo(k) = surfo(k) + min( ctmi(i,j,k,0), (scalr(i,j,k,1)-spvr) )
             enddo
             do i = iu1(j), iu2(j)
               surfv(k) = surfv(k) + ctmi(i,j,k,1)
             enddo
           enddo
           vols = vols + surfs(k) * dz(k)
           volo = volo + surfo(k) * dz(k)
           volv = volv + surfv(k) * dz(k)
           if (k < ks2) volw = volw + surfs(k) * dzw(k+1)
         enddo
         zvols = 0.
         zvolo = 0.
         zvolv = 0.
         zvolw = 0.
         if (vols > epsil) zvols = 1. / vols
         if (volo > epsil) zvolo = 1. / volo
         if (volv > epsil) zvolv = 1. / volv
         if (volw > epsil) zvolw = 1. / volw
         vols = vols * ddxy
         volo = volo * ddxy
         volv = volv * ddxy
         volw = volw * ddxy
         do k = ks1, ks2
           if (surfs(k) > epsil) then
             zsurfs(k) = 1. / surfs(k)
           else
             zsurfs(k) = 0.
           endif
           if (surfo(k) > epsil) then
             zsurfo(k) = 1. / surfo(k)
           else
             zsurfo(k) = 0.
           endif
           if (surfv(k) > epsil) then
             zsurfv(k) = 1. / surfv(k)
           else
             zsurfv(k) = 0.
           endif
           surfs(k) = surfs(k) * ddxy
           surfo(k) = surfo(k) * ddxy
           surfv(k) = surfv(k) * ddxy
         enddo
         zsurf  = zsurfs(ks2)
         zmdeta = zmdeta * zsurf

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--initialisation des variables utilisees uniquement pour "informe":
         do j = 1, jmax
           do i = 1, imax
             daeta(i,j) = 0.
           enddo
         enddo

!--Remise a zero de fss, Energ. Aj.Conv & Freq.Aj.Conv (avant 1ere it) :
         do k = 1, kmax
           do j = 1, jmax
             do i = 1, imax
               fqajc(i,j,k) = 0.
             enddo
           enddo
         enddo
         do ns = 0, nsmax
           do j = 1, jmax
             do i = 1, imax
               fss(i,j,ns) = 0.
             enddo
           enddo
         enddo

!--Initialisation of the arrays for the accumulation
         do nv = 1, nvinfo
           vinfom(nv) = 0.
         enddo

!--Definition of the limits for the downsloping budget
         ksud  = ks2
         zlim1 = -1000.0
         do kz = ks2, ks1, -1
           if (z(kz) >= zlim1) ksud = kz
         enddo
         knor  = ks2
         zlim1 = -450.0
         do kz = ks2, ks1, -1
           if (z(kz) >= zlim1) knor = kz
         enddo
         yylats = -55.
         jmsud  = nint( (yylats-ylat1)/dlat ) + 1
         yylatn = 70.
         jmnor  = nint( (yylatn-ylat1)/dlat ) + 1

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

         if (nn99 == 2) then
           write(mouchard_id,*) 'nvinfo, nocean, ntmoy :', nvinfo, nocean, ntmoy
           write(mouchard_id,*) 'zvols , vols :', zvols, vols
           write(mouchard_id,*) 'zvolo , volo :', zvolo, volo
           write(mouchard_id,*) 'zvolv , volv :', zvolv, volv
           write(mouchard_id,*) 'zvolw , volw :', zvolw, volw
           write(mouchard_id,*) 'zsurfs, surfs (ks2) :', zsurfs(ks2), surfs(ks2)
           write(mouchard_id,'(1P6E12.5)') (zsurfs(k), k=ks1, ks2)
           write(mouchard_id,'(1P6E12.5)') ( surfs(k), k=ks1, ks2)
           write(mouchard_id,*) 'zsurfo, surfo (ks2) :', zsurfo(ks2), surfo(ks2)
           write(mouchard_id,'(1P6E12.5)') (zsurfo(k), k=ks1, ks2)
           write(mouchard_id,'(1P6E12.5)') ( surfo(k), k=ks1, ks2)
           write(mouchard_id,*) 'zsurfv, surfv (ks2) :', zsurfv(ks2), surfv(ks2)
           write(mouchard_id,'(1P6E12.5)') (zsurfv(k), k=ks1, ks2)
           write(mouchard_id,'(1P6E12.5)') ( surfv(k), k=ks1, ks2)
           write(mouchard_id,'(A,2I4,2F9.2)') ' downsloping jmnor,knor,lat_N,zw :'  &
                          , jmnor, knor, ylat1+(jmnor-1)*dlat, zw(knor)
           write(mouchard_id,'(A,2I4,2F9.2)') ' downsloping jmsud,ksud,lat_S,zw :'  &
                          , jmsud, ksud, ylat1+(jmsud-1)*dlat, zw(ksud)
           write(mouchard_id,*)
           write(mouchard_id,'(A,I4,2A,I4,A)') 'File "evolu" :', nnntot, ' time ',  &
             'records of', nvinfo, ' (=nvinfo) var. ; nv,title,ktsum :'
           write(mouchard_id,'(4(A,I3,1X,A,I2,A1))') (' no=', n, titvar(n), ktsum(n), ',', n=1, nvinfo)
           write(mouchard_id,*)
         endif

         return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine informe -
       end subroutine informe
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   ENTRY POINT 2 -- inforun: per-timestep diagnostics. Assembled from D2 lots 1-4.
!           Approach 1: every vinfor write is addressed BY NAME through indices resolved once (idx_init) into local
!           save'd integers / small arrays. No running "nv" arithmetic for column addressing. The registry (informe)
!           is the single source of truth for the ordering, so a value can no longer land in the wrong column.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subroutine inforun(nn99)

         use const_mod,    only: cpo, gpes, one, rho0, svrdrp, yeaday, zero
         use para0_mod,    only: imax, jmax, kmax, nsmax
         use para_mod,     only: nbsmax
         use bloc0_mod,    only: dx, dy, iberp, jberp, jeq, js1, js2, ju1     &
              , ju2, ks1, ks2, ku2, nsav, spvr, tpstot, fqajc, daeta         &
              , is1, is2, dzw, eta, ub, vb, u, v, w, dz, z, scal, scalr       &
              , scal0, b, tmu, iu1, iu2, ahs, tms, unsdy, unsdx, dts, alphmi
         use bloc_mod,     only: aire, ctmi, cmx, cmy, smy, smx, cmxy, koutpu  &
              , ninfo, nstart, ntmoy, numit, zsurf, zvols, zvolo, zvolv, zvolw &
              , zsurfo
         use isoslope_mod, only: viso, uiso
         use ice_mod,      only: ficebergn, ficebergs, toticesm, xlg         &
              , albq, hgbq, hnbq, ddtb, vwx
         use dynami_mod,   only: bound, vg, ug, dxs1
         use reper_mod,    only: dlat, fmtw, ndhsf, nferme, nocean, nvhsf     &
              , nvinfo, ylat1, iehsf, ishsf, jshsf, jehsf, iszon, iezon       &
              , titvar, vinfor, vinfom, scalwr, ktsum, zmdeta
         use ipcc_output_mod,  only: moc, thex, tmc
         use newunit_clio_mod, only: clio3_out_id, mouchard_id, evolu_id, testevolu_id

         implicit none

         integer(ip), intent(in) :: nn99

! dmr&clo --- common block preserved verbatim from the legacy (iceberg/ISM heat coupling).
         real(dblp), dimension(imax,jmax) :: heaticbism
         real(dblp)                       :: vicbismn, vicbisms
         common /icbism/ heaticbism, vicbismn, vicbisms

! dmr&clo --- local scalars / small work arrays
         real(dblp),  dimension(3)   :: vvk
         real(dblp),  dimension(0:1) :: zitsum
         logical                     :: flgout

         integer(ip) :: i, ii, ii1, ii2, j, jj, jj1, jj2, k, kk, n, nhsf      &
              , nn, nninfo, nniter, nnp, ns, nv
         real(dblp)  :: ctmobs, factke, phiy, ss, ssc2, sszonm, sumk, therma  &
              , vber, vbord, vv, vv2, v2cd2, vvn, vvp, vvpsm, yy, zz          &
              , saltatlantic, volatlantic, saltatlantic_new, saltatlantic_dt  &
              , vvtfa, vvtfd, vvtba, vvtbd, vvtca, vvtcd, vva, vvd

! dmr&clo --- one-time resolved column indices (approach 1). All save'd; resolved on first call.
         logical,     save :: idx_init = .false.
         integer(ip), save :: i_noit, i_tyr, i_egajc, i_vajc, i_deta, i_meta
         integer(ip), save, dimension(:), allocatable :: idx_hsf
         integer(ip), save :: i_dans, i_dann, i_iscs, i_iscn, i_fras, i_fran, i_pnws, i_pnwn
         integer(ip), save :: i_adgin, i_adpro, i_adout, i_aabpr, i_aabex, i_aabat
         integer(ip), save :: i_fc30a, i_fs30a, i_fsa30, i_fsd30, i_fsana, i_fsdna
         integer(ip), save :: i_dtsalt, i_mov30, i_fsber
         integer(ip), save, dimension(:), allocatable :: i_scmean, i_sclvl1, i_scanom
         integer(ip), save, dimension(:), allocatable :: i_to_lvl, i_so_lvl
         integer(ip), save :: i_absw, i_absu, i_absv, i_ke
         integer(ip), save, dimension(7) :: i_ice_n, i_ice_s
         integer(ip), save :: i_frag, i_spng, i_berg, i_thex, i_ismm, i_icbn, i_icbs
         integer(ip), save :: nm1n, nm2n
         integer(ip)       :: itr, m
         character(len=8)  :: cc8

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! dmr&clo  0) One-time resolution of all column indices used by inforun, by NAME.
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
         if (.not. idx_init) then
           i_noit  = idx('NoIt') ; i_tyr   = idx('T yr')
           i_egajc = idx('EgAjC'); i_vajc  = idx('V_AjC')
           i_deta  = idx('D.Eta'); i_meta  = idx('M.Eta')

           allocate(idx_hsf(nvhsf))
           idx_hsf(1) = idx('DrHSF') ; idx_hsf(2) = idx('InHSF') ; idx_hsf(3) = idx('BeHSF')
           do itr = 2, nvhsf
             if (idx_hsf(itr) /= idx_hsf(itr-1) + 1) then
               write(clio3_out_id,*) 'STOP inforun: HSF indices not contiguous at ', itr
               error stop 'inforun: HSF idx not contiguous'
             endif
           enddo

           i_dans = idx('DANS') ; i_dann = idx('DANN')
           i_iscs = idx('ISCS') ; i_iscn = idx('ISCN')
           i_fras = idx('FRAS') ; i_fran = idx('FRAN')
           i_pnws = idx('PNWS') ; i_pnwn = idx('PNWN')

           i_adgin = idx('ADGIN'); i_adpro = idx('ADPro'); i_adout = idx('ADOut')
           i_aabpr = idx('AABpr'); i_aabex = idx('AABex'); i_aabat = idx('AABat')

           i_fc30a = idx('Fc30A'); i_fs30a = idx('Fs30A')
           i_fsa30 = idx('Fsa30A'); i_fsd30 = idx('Fsd30A')
           i_fsana = idx('FsaNA'); i_fsdna = idx('FsdNA')
           i_dtsalt = idx('dtsaltA'); i_mov30 = idx('Mov30A'); i_fsber = idx('Fsber')

           allocate(i_scmean(nsmax), i_sclvl1(nsmax), i_scanom(nsmax))
           i_scmean(1) = idx('T-c') ; i_sclvl1(1) = idx('T1-o') ; i_scanom(1) = idx('|T-o|')
           i_scmean(2) = idx('S-30'); i_sclvl1(2) = idx('S1-o') ; i_scanom(2) = idx('|S-o|')

           allocate(i_to_lvl(ks1:ks2), i_so_lvl(ks1:ks2))
           if (ninfo >= 0) then
             do k = ks1, ks2
               write(cc8,'(A3,I2)') 'T-o', k ; i_to_lvl(k) = idx(cc8)
               write(cc8,'(A3,I2)') 'S-o', k ; i_so_lvl(k) = idx(cc8)
             enddo
           else
             write(cc8,'(A3,I2)') 'T-o', 1 ; i_to_lvl(ks1) = idx(cc8)
             do k = 1+ks1, ks2
               write(cc8,'(A1,I2,A2)') 'T', k, '-o' ; i_to_lvl(k) = idx(cc8)
             enddo
             write(cc8,'(A3,I2)') 'S-o', 1 ; i_so_lvl(ks1) = idx(cc8)
             do k = 1+ks1, ks2
               write(cc8,'(A1,I2,A2)') 'S', k, '-o' ; i_so_lvl(k) = idx(cc8)
             enddo
           endif

           i_absw = idx('|w|') ; i_absu = idx('|u|') ; i_absv = idx('|v|') ; i_ke = idx('K.E')

           i_ice_n(1) = idx('AIEFN'); i_ice_n(2) = idx('A15N'); i_ice_n(3) = idx('A85N')
           i_ice_n(4) = idx('ALEN') ; i_ice_n(5) = idx('VOLN'); i_ice_n(6) = idx('VONN')
           i_ice_n(7) = idx('ECGN')
           i_ice_s(1) = idx('AIEFS'); i_ice_s(2) = idx('A15S'); i_ice_s(3) = idx('A85S')
           i_ice_s(4) = idx('ALES') ; i_ice_s(5) = idx('VOLS'); i_ice_s(6) = idx('VONS')
           i_ice_s(7) = idx('ECGS')

           i_frag = idx('FRAG'); i_spng = idx('SPNG'); i_berg = idx('BERG')
           i_thex = idx('ThEx'); i_ismm = idx('ISMM'); i_icbn = idx('IcbN'); i_icbs = idx('IcbS')

! dmr&clo --- averaging anchors (reserve R8: asymmetric, reproduced faithfully). nm1n=PNWN, nm2n=AABat.
           nm1n = idx('PNWN') ; nm2n = idx('AABat')

           idx_init = .true.
         endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1 ) Prepare & start filling "vinfor".
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
         flgout = mod(numit,ninfo) == 0 .or. numit == 1

         do nv = 1, nvinfo
           vinfor(nv) = 0.
         enddo

         vinfor(i_noit) = DFLOAT(numit)
         vinfor(i_tyr)  = tpstot / ( 86400. * yeaday )

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) Global variables
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
         if (flgout) then
!--2.1 Convective-adjustment energy + elevation change
           do j = js1, js2
             do i = is1(j), is2(j)
               vinfor(i_egajc) = vinfor(i_egajc) + ctmi(i,j,ks2,0) * fqajc(i,j,1)
               vinfor(i_deta)  = vinfor(i_deta)  + ctmi(i,j,ks2,0) * daeta(i,j)
               daeta(i,j) = 0.
             enddo
           enddo
           vinfor(i_egajc) = vinfor(i_egajc) * zsurf
           vinfor(i_deta)  = vinfor(i_deta)  * zmdeta

!--AjC frequency over the ocean volume -> V_AjC
           do k = ks1+1, ks2
             sumk = 0.
             kk = k - 1
             do j = js1, js2
               do i = is1(j), is2(j)
                 sumk = sumk + ctmi(i,j,kk,0) * fqajc(i,j,k)
               enddo
             enddo
             vinfor(i_vajc) = vinfor(i_vajc) + sumk * dzw(k)
           enddo
           vinfor(i_vajc) = vinfor(i_vajc) * zvolw * 1000.
         endif

         if (flgout .or. ntmoy == 2) then
!--2.2 Basin-mean elevation -> M.Eta
           do j = js1, js2
             do i = is1(j), is2(j)
               vinfor(i_meta) = vinfor(i_meta) + ctmi(i,j,ks2,0) * eta(i,j)
             enddo
           enddo
           vinfor(i_meta) = vinfor(i_meta) * zsurf
         endif

!--2.3 Sverdrup transport: Drake / Indonesia / Bering -> DrHSF/InHSF/BeHSF
         do nhsf = 1, nvhsf
           nn = idx_hsf(nhsf)
           if (iehsf(nhsf) == 0) then
             ii = ishsf(nhsf)
             do j = jshsf(nhsf), jehsf(nhsf)
               vinfor(nn) = vinfor(nn) + cmy(ii,j,3) * ub(ii,j)
             enddo
             vinfor(nn) = vinfor(nn) * dy * svrdrp
           else
             jj = jshsf(nhsf)
             do i = ishsf(nhsf), iehsf(nhsf)
               vinfor(nn) = vinfor(nn) + cmx(i,jj,3) * vb(i,jj)
             enddo
             vinfor(nn) = vinfor(nn) * dx * svrdrp
           endif
         enddo

!--2.4 N/S transports at straits: DANS/DANN, ISCS/ISCN, FRAS/FRAN, PNWS/PNWN
         do nhsf = 4, min(7, ndhsf)
           select case (nhsf)
             case (4) ; nn = i_dans
             case (5) ; nn = i_iscs
             case (6) ; nn = i_fras
             case (7) ; nn = i_pnws
           end select
           if (iehsf(nhsf) == 0) then
             i = ishsf(nhsf)
             do k = ks1, ks2
               vvn = 0. ; vvp = 0.
               do j = jshsf(nhsf), jehsf(nhsf)
                 vvn = vvn + tmu(i,j,k) * cmy(i,j,3) * min(zero, u(i,j,k))
                 vvp = vvp + tmu(i,j,k) * cmy(i,j,3) * max(zero, u(i,j,k))
               enddo
               vinfor(nn)   = vinfor(nn)   + vvn * dz(k)
               vinfor(nn+1) = vinfor(nn+1) + vvp * dz(k)
             enddo
             vinfor(nn)   = vinfor(nn)   * dy * svrdrp
             vinfor(nn+1) = vinfor(nn+1) * dy * svrdrp
           else
             j = jshsf(nhsf)
             do k = ks1, ks2
               vvn = 0. ; vvp = 0.
               do i = ishsf(nhsf), iehsf(nhsf)
                 vvn = vvn + tmu(i,j,k) * cmx(i,j,3) * min(zero, v(i,j,k))
                 vvp = vvp + tmu(i,j,k) * cmx(i,j,3) * max(zero, v(i,j,k))
               enddo
               vinfor(nn)   = vinfor(nn)   + vvn * dz(k)
               vinfor(nn+1) = vinfor(nn+1) + vvp * dz(k)
             enddo
             vinfor(nn)   = vinfor(nn)   * dx * svrdrp
             vinfor(nn+1) = vinfor(nn+1) * dx * svrdrp
           endif
         enddo

!--2.5 THC circulation via vwx (living version). Legacy 'zold' branch dropped (dead code).
!- ADGIN
         yy = 68. ; jj1 = 1 + nint( (yy-ylat1)/dlat + 0.5 )
         yy = 75. ; jj2 = 1 + nint( (yy-ylat1)/dlat + 0.5 )
         zz = -500.0
         vinfor(i_adgin) = 0.0
         do jj = jj1, jj2
           do k = ks1, ks2
             if (z(k) > zz) cycle
             vinfor(i_adgin) = max( vinfor(i_adgin), vwx(jj,k,3) )
           enddo
         enddo
         vinfor(i_adgin) = vinfor(i_adgin) * dx * svrdrp
!- ADPro (exported as moc)
         yy = 45. ; jj1 = 1 + nint( (yy-ylat1)/dlat + 0.5 )
         yy = 75. ; jj2 = 1 + nint( (yy-ylat1)/dlat + 0.5 )
         zz = -500.0
         vinfor(i_adpro) = 0.0
         do jj = jj1, jj2
           do k = ks1, ks2
             if (z(k) > zz) cycle
             vinfor(i_adpro) = max( vinfor(i_adpro), vwx(jj,k,3) )
           enddo
         enddo
         moc = vinfor(i_adpro)
!- ADOut / AABat at 20 S
         yy = -20.0 ; jj = 1 + nint( (yy-ylat1)/dlat + 0.5 )
         zz = -500.0
         vinfor(i_adout) = 0.0 ; vinfor(i_aabat) = 0.0
         do k = ks1, ks2
           if (z(k) > zz) cycle
           vinfor(i_adout) = max( vinfor(i_adout), vwx(jj,k,3) )
           vinfor(i_aabat) = min( vinfor(i_aabat), vwx(jj,k,3) )
         enddo
!- AABpr
         yy = -70.0 ; jj1 = 1 + nint( (yy-ylat1)/dlat + 0.5 )
         yy = -60.0 ; jj2 = 1 + nint( (yy-ylat1)/dlat + 0.5 )
         zz = -500.0
         vinfor(i_aabpr) = 0.0
         do jj = jj1, jj2
           do k = ks1, ks2
             if (z(k) > zz) cycle
             vinfor(i_aabpr) = min( vinfor(i_aabpr), vwx(jj,k,0) )
           enddo
         enddo
!- AABex at 30 S
         yy = -30.0 ; jj = 1 + nint( (yy-ylat1)/dlat + 0.5 )
         zz = -500.0
         vinfor(i_aabex) = 0.0
         do k = ks1, ks2
           if (z(k) > zz) cycle
           vinfor(i_aabex) = min( vinfor(i_aabex), vwx(jj,k,0) )
         enddo

!--2.6a Heat flux at 30 S (Fc30A)
         yy = -30.0 ; jj = 1 + nint( (yy-ylat1)/dlat + 0.5 )
         vvpsm = 0.0 ; ssc2 = 2.0 * scal0(ks1,1)
         do k = ks1, ks2
           vv = 0.0
           do i = iszon(jj,nbsmax), iezon(jj,nbsmax)
             vv = vv + face_scalar_flux(DIR_MERID, FLUX_TOTAL, i, jj, k, 1, ssc2)
           enddo
           vvpsm = vvpsm + vv*dz(k)
         enddo
         vinfor(i_fc30a) = vvpsm * dx * rho0*cpo * 1.D-15

!--Atlantic average salinity (reference for Mov30/dtsaltA)
         saltatlantic = 0.0 ; volatlantic = 0.0
         do j = 1, jmax     ! dmr&clo R9 fix: was 'do j=1,imax' (legacy) -- j indexes latitude (jmax), not imax. The
!                           !   out-of-bounds iterations j>jmax read empty iszon/iezon ranges, so this changes nothing
!                           !   numerically (verified: bit-identical output) while removing the OOB access.
           do i = iszon(j,nbsmax), iezon(j,nbsmax)
             do k = 1, kmax
               saltatlantic = saltatlantic + aire(i,j)*dz(k)*tms(i,j,k)*scal(i,j,k,2)
               volatlantic  = volatlantic  + aire(i,j)*dz(k)*tms(i,j,k)
             enddo
           enddo
         enddo
         saltatlantic = saltatlantic / volatlantic

!--Salt flux at 30 S (Fs30A)
         vvpsm = 0.0 ; ssc2 = 2.0 * scal0(ks1,2)
         do k = ks1, ks2
           vv = 0.0
           do i = iszon(jj,nbsmax), iezon(jj,nbsmax)
             vv = vv + face_scalar_flux(DIR_MERID, FLUX_TOTAL, i, jj, k, 2, ssc2)
           enddo
           vvpsm = vvpsm + vv*dz(k)
         enddo
         vinfor(i_fs30a) = vvpsm * dx

!--Advective salt flux at 30 S (Fsa30A) -- FLUX_ADV_LEGACY anomaly (reserve R1)
         vvpsm = 0.0
         do k = ks1, ks2
           vv = 0.0
           do i = iszon(jj,nbsmax), iezon(jj,nbsmax)
             vv = vv + face_scalar_flux(DIR_MERID, FLUX_ADV_LEGACY, i, jj, k, 2, ssc2)
           enddo
           vvpsm = vvpsm + vv*dz(k)
         enddo
         vinfor(i_fsa30) = vvpsm * dx

!--Diffusive salt flux at 30 S (Fsd30A)
         vvpsm = 0.0
         do k = ks1, ks2
           vv = 0.0
           do i = iszon(jj,nbsmax), iezon(jj,nbsmax)
             vv = vv + face_scalar_flux(DIR_MERID, FLUX_DIF, i, jj, k, 2, ssc2)
           enddo
           vvpsm = vvpsm + vv*dz(k)
         enddo
         vinfor(i_fsd30) = vvpsm * dx

!--Advective & diffusive salt fluxes at the northern boundary (Fram+Barents+CAA)
         vvtfa = 0.0 ; vvtfd = 0.0 ; vvtba = 0.0 ; vvtbd = 0.0 ; vvtca = 0.0 ; vvtcd = 0.0
         do k = ks1, ks2
!- Fram (meridional) ii=106,107 j=55
           vva = 0.0 ; vvd = 0.0 ; jj = 55
           do ii = 106, 107
             vva = vva + face_scalar_flux(DIR_MERID, FLUX_ADV_LEGACY, ii, jj, k, 2, 0._dblp)
             vvd = vvd + face_scalar_flux(DIR_MERID, FLUX_DIF,        ii, jj, k, 2, 0._dblp)
           enddo
           vvtfa = vvtfa + vva*dz(k) ; vvtfd = vvtfd + vvd*dz(k)*(-1)
!- Barents (zonal) ii=110 j=54,56
           vva = 0.0 ; vvd = 0.0 ; ii = 110
           do jj = 54, 56
             vva = vva + face_scalar_flux(DIR_ZONAL, FLUX_ADV_LEGACY, ii, jj, k, 2, 0._dblp)
             vvd = vvd + face_scalar_flux(DIR_ZONAL, FLUX_DIF,        ii, jj, k, 2, 0._dblp)
           enddo
           vvtba = vvtba + vva*dz(k) ; vvtbd = vvtbd + vvd*dz(k)*(-1)
!- CAA (meridional) ii=101,102 j=57
           vva = 0.0 ; vvd = 0.0 ; jj = 57
           do ii = 101, 102
             vva = vva + face_scalar_flux(DIR_MERID, FLUX_ADV_LEGACY, ii, jj, k, 2, 0._dblp)
             vvd = vvd + face_scalar_flux(DIR_MERID, FLUX_DIF,        ii, jj, k, 2, 0._dblp)
           enddo
           vvtca = vvtca + vva*dz(k) ; vvtcd = vvtcd + vvd*dz(k)*(-1)
         enddo
         vinfor(i_fsana) = (vvtfa*dx) + (vvtba*dy) + (vvtca*dx)
         vinfor(i_fsdna) = (vvtfd*dx) + (vvtbd*dy) + (vvtcd*dx)

!--Atlantic salt-content tendency (dtsaltA)
         saltatlantic_new = 0.0
         do j = 1, jmax     ! dmr&clo R9 fix: was 'do j=1,imax' (legacy). Same reasoning as the saltatlantic loop above.
           do i = iszon(j,nbsmax), iezon(j,nbsmax)
             do k = 1, kmax
               saltatlantic_new = saltatlantic_new + aire(i,j)*dz(k)*tms(i,j,k)*scal(i,j,k,2)
             enddo
           enddo
         enddo
         if (first_dtsalt) then
           saltatlantic_dt = 0.0 ; first_dtsalt = .false.
         else
           saltatlantic_dt = (saltatlantic_new - saltatlantic_old_sav) / ddtb
         endif
         saltatlantic_old_sav = saltatlantic_new
         vinfor(i_dtsalt) = saltatlantic_dt

!--Mov at 30 S (explicit)
         yy = -30.0 ; jj = 1 + nint( (yy - ylat1) / dlat + 0.5 )
         vvpsm = 0.0
         do k = ks1, ks2
           ss = 0.0
           do i = iszon(jj,nbsmax), iezon(jj,nbsmax)
             ss = ss + 0.5 * ( scal(i,jj-1,k,2) + scal(i,jj,k,2) )
           enddo
           sszonm = ss / ( iezon(jj,nbsmax) - iszon(jj,nbsmax) )
           vv = 0.0
           do i = iszon(jj,nbsmax), iezon(jj,nbsmax)
             vv2   = 0.5 * ( v(i,jj,k) + v(i+1,jj,k) ) + viso(i,jj,k)
             v2cd2 = 0.5 * cmx(i,jj,2) * vv2
             phiy  = tms(i,jj,k) * tms(i,jj-1,k) * ( v2cd2 * (sszonm - saltatlantic) )
             vv    = vv + phiy
           enddo
           vvpsm = vvpsm + vv*dz(k)
         enddo
         vinfor(i_mov30) = vvpsm * dx * (-1/saltatlantic) / 1e6

!--Salt flux at Bering (Fsber), surface level only
         ii = 102 ; jj = 65 ; ssc2 = 2.0 * scal0(ks1,2)
         vvpsm = face_scalar_flux(DIR_MERID, FLUX_TOTAL, ii, jj, ks2, 2, ssc2) * dz(ks2)
         vinfor(i_fsber) = vvpsm * dx

!--2.6 Block C : basin-mean scalars + per-level |anomaly| families (only when flgout .or. ntmoy==2)
         if (flgout .or. ntmoy == 2) then
           nnp = 3
           if (ninfo < 0) nnp = 2
           do ns = 1, nsmax
             k = ks2
             vvk(1) = 0. ; vvk(2) = 0. ; vvk(3) = 0.
             do j = js1, js2
               do i = is1(j), is2(j)
                 ctmobs = min( ctmi(i,j,k,0), (scalr(i,j,k,ns)-spvr) )
                 vvk(1) = vvk(1) + ctmi(i,j,k,0) * scal(i,j,k,ns)
                 vv2    = ctmobs * (scal(i,j,k,ns) - scalr(i,j,k,ns))
                 vvk(2) = vvk(2) + vv2
                 vvk(3) = vvk(3) + abs(vv2)
               enddo
             enddo
             vinfor(i_scmean(ns)) = vvk(1) * dz(k)
             vinfor(i_scanom(ns)) = vvk(3) * dz(k)
             vinfor(i_sclvl1(ns)) = vvk(2) * zsurfo(k)
             if (ns <= 2) then
               m = (ks2+1) - k
               if (ns == 1) vinfor(i_to_lvl(m)) = vvk(3) * zsurfo(k)
               if (ns == 2) vinfor(i_so_lvl(m)) = vvk(3) * zsurfo(k)
             endif
             do k = ks2-1, ks1, -1
               vvk(1) = 0. ; vvk(2) = 0. ; vvk(3) = 0.
               do j = js1, js2
                 do i = is1(j), is2(j)
                   ctmobs = min( ctmi(i,j,k,0), (scalr(i,j,k,ns)-spvr) )
                   vvk(1) = vvk(1) + ctmi(i,j,k,0) * scal(i,j,k,ns)
                   vv2    = ctmobs * (scal(i,j,k,ns) - scalr(i,j,k,ns))
                   vvk(2) = vvk(2) + vv2
                   vvk(3) = vvk(3) + abs(vv2)
                 enddo
               enddo
               vinfor(i_scmean(ns)) = vinfor(i_scmean(ns)) + vvk(1) * dz(k)
               vinfor(i_scanom(ns)) = vinfor(i_scanom(ns)) + vvk(3) * dz(k)
               if (ns <= 2) then
                 m = (ks2+1) - k
                 if (ns == 1) vinfor(i_to_lvl(m)) = vvk(nnp) * zsurfo(k)
                 if (ns == 2) vinfor(i_so_lvl(m)) = vvk(nnp) * zsurfo(k)
               endif
             enddo
             vinfor(i_scmean(ns)) = vinfor(i_scmean(ns)) * zvols + scalwr(ns)
             if (ns == 1) tmc = vinfor(i_scmean(ns)) * cpo * rho0
             vinfor(i_scanom(ns)) = vinfor(i_scanom(ns)) * zvolo
           enddo

!--Basin mean of |w|
           vinfor(i_absw) = 0.
           do k = ks1+1, ks2
             vvk(1) = 0.
             do j = js1, js2
               do i = is1(j), is2(j)
                 vvk(1) = vvk(1) + ctmi(i,j,k-1,0) * abs(w(i,j,k))
               enddo
             enddo
             vinfor(i_absw) = vinfor(i_absw) + vvk(1) * dzw(k)
           enddo
           vinfor(i_absw) = vinfor(i_absw) * zvolw

!--Basin mean of |u|, |v|, K.E
           factke = rho0 * 0.5
           vinfor(i_absu) = 0. ; vinfor(i_absv) = 0. ; vinfor(i_ke) = 0.
           do k = ks1, ks2
             vvk(1) = 0. ; vvk(2) = 0. ; vvk(3) = 0.
             do j = ju1, ju2
               do i = iu1(j), iu2(j)
                 vvk(1) = vvk(1) + cmxy(i,j,3) * abs(u(i,j,k))
                 vvk(2) = vvk(2) + cmxy(i,j,3) * abs(v(i,j,k))
                 vvk(3) = vvk(3) + cmxy(i,j,3) * ( u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) )
               enddo
             enddo
             vvk(3) = vvk(3) * factke
             vinfor(i_absu) = vinfor(i_absu) + vvk(1) * dz(k)
             vinfor(i_absv) = vinfor(i_absv) + vvk(2) * dz(k)
             vinfor(i_ke)   = vinfor(i_ke)   + vvk(3) * dz(k)
           enddo
           vinfor(i_absu) = vinfor(i_absu) * zvolv
           vinfor(i_absv) = vinfor(i_absv) * zvolv
           vinfor(i_ke)   = vinfor(i_ke)   * zvolv
         endif

!--2.7 Sea-ice / freshwater / icebergs / thermal expansion (every call)
         call sea_ice_hemisphere(jeq, js2,  i_ice_n)
         call sea_ice_hemisphere(js1, jeq-1, i_ice_s)

#if ( HRCLIO == 0 )
         j = 56 ; ii1 = 106 ; ii2 = 108
#elif ( HRCLIO == 1 )
         j = 109 ; ii1 = 209 ; ii2 = 215
#endif
         vinfor(i_frag) = ice_freshwater_flux(j, ii1, ii2)

#if ( HRCLIO == 0 )
         j = 55 ; ii1 = 109 ; ii2 = 114
#elif ( HRCLIO == 1 )
         j = 108 ; ii1 = 216 ; ii2 = 225
#endif
         vinfor(i_spng) = ice_freshwater_flux(j, ii1, ii2)

         vbord = 1.0 + (1.0 - bound)
         vber  = vg(iberp,jberp)*tmu(iberp,jberp,ku2) / vbord * 1e-6
         jj    = jberp - max(0, int(sign(one,vber)))
         vinfor(i_berg) = vber*dxs1(iberp,jberp-1)*hgbq(iberp,jj)*(1.0-albq(iberp,jj))          &
                        + vber*dxs1(iberp-1,jberp-1)*hgbq(iberp-1,jj)*(1.0-albq(iberp-1,jj))

         therma = 0.0
         do k = ks1, ks2
           do j = js1, js2
             do i = is1(j), is2(j)
               therma = therma + dz(k)*ctmi(i,j,k,0)*b(i,j,k)
             enddo
           enddo
         enddo
         vinfor(i_thex) = (-1)*therma*zsurf/gpes
         thex = vinfor(i_thex)

         vinfor(i_ismm) = toticesm
         vinfor(i_icbn) = ficebergn/xlg*360.
#if ( ISM == 1 )
         if (flgism) vinfor(i_icbn) = vinfor(i_icbn) + vicbismn
#endif
         vinfor(i_icbs) = ficebergs/xlg*360.
#if ( ISM == 1 )
         if (flgism) vinfor(i_icbs) = vinfor(i_icbs) + vicbisms
#endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  4) Accumulation
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
         do nv = 2, nvinfo
           vinfom(nv) = vinfom(nv) + vinfor(nv)
         enddo

         if (.not. flgout) return

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  5) Averaging + write
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
         nninfo = abs(ninfo)
         nniter = numit - nstart + 1
         zitsum(0) = 1.
         if (nniter == 0) then
           zitsum(1) = 0.
         elseif (nniter < nninfo) then
           zitsum(1) = 1. / DFLOAT(nniter)
         else
           zitsum(1) = 1. / DFLOAT(nninfo)
         endif

         do nv = 2, nm1n-1
           vinfor(nv) = vinfom(nv) * zitsum(ktsum(nv))
           vinfom(nv) = 0.0
         enddo
         do nv = nm2n, nvinfo
           vinfor(nv) = vinfom(nv) * zitsum(ktsum(nv))
           vinfom(nv) = 0.0
         enddo
         do nv = nm1n, nm2n-1
           vinfom(nv) = 0.0
         enddo

         write(evolu_id,fmtw) (titvar(nv), vinfor(nv), nv=1,nvinfo)
#if ( ISM >= 2 )
         write(testevolu_id,fmtw) (titvar(nv), vinfor(nv), nv=1,nvinfo)
#endif
         if (numit == nferme) then
           write(clio3_out_id,*) ' evolu (in inforun), numit,nferme "', numit, nferme
           write(evolu_id,*)
           close(evolu_id)
#if ( ISM >= 2 )
           close(testevolu_id)
#endif
         endif

!--6) Reset accumulators
         if (koutpu <= 1 .or. mod(numit,nsav) /= 0) then
           do k = 1, kmax
             do j = 1, jmax
               do i = 1, imax
                 fqajc(i,j,k) = 0.
               enddo
             enddo
           enddo
         endif

         return
       end subroutine inforun
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   PRIVATE HELPERS (flux kernel + factored motifs)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! dmr&clo --- face_scalar_flux: scalar flux through ONE cell face, for tracer ns, returned as phiy per `mode`.
!             `direction` selects the meridional (v/cmx/smy/unsdy, neighbour j-1) vs zonal (u/cmy/smx/unsdx, neighbour
!             i-1) stencil. Modes: FLUX_TOTAL/ADV/DIF and FLUX_ADV_LEGACY (anomaly, reserve R1). No accumulation, no
!             scaling, no sign -- those stay at the call site. Verified bit-identical to the legacy per-face formulas.
       real(dblp) function face_scalar_flux(direction, mode, i, j, k, ns, ssc2) result(phiy)

         use const_mod,    only: one
         use bloc0_mod,    only: v, u, unsdy, unsdx, dts, alphmi, scal, ahs, tms
         use bloc_mod,     only: cmx, cmy, smy, smx, cmxy
         use isoslope_mod, only: viso, uiso
         use newunit_clio_mod, only: clio3_out_id

         integer(ip), intent(in) :: direction, mode, i, j, k, ns
         real(dblp),  intent(in) :: ssc2
         real(dblp) :: vv2, tcd2, cn, aalpha, ccdif, diffc, tmsprod, sL, sR

         if (direction == DIR_MERID) then
           vv2     = 0.5 * ( v(i,j,k) + v(i+1,j,k) ) + viso(i,j,k)
           tcd2    = 0.5 * cmx(i,j,2) * vv2
           cn      = smy(i,j,2) * unsdy * dts(k) * vv2
           ccdif   = ahs(k) * unsdy
           tmsprod = tms(i,j,k) * tms(i,j-1,k)
           sL      = scal(i,j-1,k,ns)
           sR      = scal(i,j,k,ns)
         else
           vv2     = 0.5 * ( u(i,j,k) + u(i,j+1,k) ) + uiso(i,j,k)
           tcd2    = 0.5 * cmy(i,j,2) * vv2
           cn      = smx(i,j,2) * unsdx * dts(k) * vv2
           ccdif   = ahs(k) * unsdx
           tmsprod = tms(i,j,k) * tms(i-1,j,k)
           sL      = scal(i-1,j,k,ns)
           sR      = scal(i,j,k,ns)
         endif

         aalpha = min( one, abs(cn) + alphmi(k) )
         diffc  = aalpha * abs(tcd2) + cmxy(i,j,2) * ccdif

         select case (mode)
           case (FLUX_TOTAL)
             phiy = tmsprod * ( tcd2*(sL + sR - ssc2) + diffc*(sL - sR) )
           case (FLUX_ADV)
             phiy = tmsprod * ( tcd2*(sL + sR - ssc2) )
           case (FLUX_DIF)
             phiy = tmsprod * ( diffc*(sL - sR) )
           case (FLUX_ADV_LEGACY)
             phiy = tmsprod * ( tcd2 * 0.5*(sL + sR) )
           case default
             write(clio3_out_id,*) 'STOP face_scalar_flux: unknown mode = ', mode
             error stop 'face_scalar_flux: unknown mode'
         end select

       end function face_scalar_flux

! dmr&clo --- ice_freshwater_flux: freshwater carried in sea ice across a strait (Fram, Barents/Kara). Factored from
!             the copied vfram/vbord/dxs1/hnbq/hgbq blocks. Sums i=ii1-1..ii2 at row j.
       real(dblp) function ice_freshwater_flux(j, ii1, ii2) result(ff)
         use const_mod,  only: one
         use bloc0_mod,  only: tmu, ks2
         use dynami_mod, only: bound, vg, dxs1
         use ice_mod,    only: hnbq, hgbq, albq
         integer(ip), intent(in) :: j, ii1, ii2
         real(dblp)  :: vbord, vfram
         integer(ip) :: i, ip1, jj
         ff    = 0.0
         vbord = 1.0 + (1.0 - bound)
         do i = ii1-1, ii2
           ip1   = i + 1
           vfram = ( vg(i,j)*tmu(i,j,ks2) + vg(ip1,j)*tmu(ip1,j,ks2) )      &
                 / ( max( tmu(i,j,ks2) + tmu(ip1,j,ks2), vbord ) )
           jj    = j - max(0, int(sign(one, vfram)))
           ff    = ff + vfram*dxs1(i,j-1)*1e-6 * (hnbq(i,jj)*0.33 + hgbq(i,jj)*0.9) * (1.0-albq(i,jj))
         enddo
       end function ice_freshwater_flux

! dmr&clo --- sea_ice_hemisphere: fills the 7 sea-ice columns for rows j0..j1 into icol(1:7):
!               1 AIEF 2 A15 3 A85 4 ALE 5 VOL 6 VON 7 ECG(RMS drift speed, finalized by sqrt(drift/VOL)).
       subroutine sea_ice_hemisphere(j0, j1, icol)
         use const_mod,  only: one
         use bloc0_mod,  only: ks2, is1, is2, tms
         use bloc_mod,   only: aire
         use ice_mod,    only: albq, hgbq, hnbq
         use dynami_mod, only: ug, vg
         use reper_mod,  only: vinfor
         integer(ip), intent(in) :: j0, j1
         integer(ip), dimension(7), intent(in) :: icol
         integer(ip) :: i, j
         real(dblp)  :: fice
         do i = 1, 7
           vinfor(icol(i)) = 0.
         enddo
         do j = j0, j1
           do i = is1(j), is2(j)
             if (tms(i,j,ks2) == 1) then
               fice = 1.0 - albq(i,j)
               vinfor(icol(1)) = vinfor(icol(1)) + fice*aire(i,j)
               if (albq(i,j) < 0.15) vinfor(icol(2)) = vinfor(icol(2)) + aire(i,j)
               if (albq(i,j) < 0.85) vinfor(icol(3)) = vinfor(icol(3)) + aire(i,j)
               vinfor(icol(4)) = vinfor(icol(4)) + albq(i,j)*aire(i,j)*fice / max(fice, 1e-6*one)
               vinfor(icol(5)) = vinfor(icol(5)) + fice*aire(i,j)*hgbq(i,j)
               vinfor(icol(6)) = vinfor(icol(6)) + fice*aire(i,j)*hnbq(i,j)
               vinfor(icol(7)) = vinfor(icol(7)) + fice*aire(i,j)*hgbq(i,j)*( ug(i,j)*ug(i,j) + vg(i,j)*vg(i,j) )
             endif
           enddo
         enddo
         vinfor(icol(7)) = sqrt( vinfor(icol(7)) / max( vinfor(icol(5)), 1e-5*one ) )
       end subroutine sea_ice_hemisphere

      end module evolu_diag_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
