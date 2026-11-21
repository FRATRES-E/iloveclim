# iLOVECLIM — Progress Summary
### June–August 2026 (last quarter in a nutshell)

---

## `main` branch

> **63 commits, 4 contributors: a quarter of deep code modernization, plus new ocean tracers and land-carbon feedbacks.** Existing configurations should reproduce prior results — with two exceptions flagged below.

### ⚠ Action required
- **VECODE / carbon cycle needs a reinstall** before use (module rewrite — see below).
- **Verify ¹⁴C** if you rely on it: a restart-bug fix shifts results by a few percent. Carbon itself is binary-identical.

### What changed

**Code modernization — Fortran-77 → Fortran-90 modules** `[ready · non-regression validated]`
- **ECBilt** (atmosphere) modularized in a five-step June campaign: `atmdyn0`, `atmphys`, `initial0`, `ecbilt0` are now proper F90 modules; hardcoded declarations cleaned up and legacy `COMATM` common-block switches removed.
- **VECODE + carbon cycle** fully rewritten into modules (`veget.f90`, `eco2_mod.f90`, `c14_prod_mod.f90`, `carbone_co2_mod.f90`…). Carbon is binary-identical; a ¹⁴C restart bug was fixed. `[needs reinstall]`
- **CLIO `defgrid`** split into three F90 modules.

**New & refined biogeochemistry**
- **Argon** added as an ocean tracer `[opt-in: ARGON flag]` — initialization, ocean–atm fluxes, transport, extra advected tracer, NewGen output. New ocean–atmosphere isotopic fractionations alongside.
- **FROG soil-carbon feedback** onto atmospheric CO₂, plus **silica cycle** work and oxygen-flux debugging/homogenization. `[ready]`
- **Permafrost** bug fixed; several carbon-flux corrections.

**Output & infrastructure**
- **`evolu` diagnostics written directly as netCDF** — robustness-focused rewrite plus controlling namelist. `[experimental · off by default · namelists copied manually]` A `MOV30`/`FSA30` computation bug was fixed en route.
- **ERA5 offline wind forcing** now compiles in the standard version; absolute paths for wind/source files.
- Standard group-testing job files added; `i-operator` submodule bumped (latest installs GRISLI as a library under `ifort`).

### Contributor activity

| Contributor | Commits | Lines (+/−) | Files | Main focus |
|---|---|---|---|---|
| **dmr** | 36 | +15,428 / −16,049 | 77 | ECBilt & VECODE/carbon modernization, evolu netCDF, defgrid |
| **clermonte** | 17 | +457 / −62 | 14 | Argon tracer, isotopic fractionations, oxygen output |
| **nbouttes** | 9 | +232 / −60 | 15 | FROG soil-carbon feedback, silica cycle, permafrost & oxygen fixes |
| **MathisVI** | 1 | +3 / −0 | 1 | Missing LGM flag in `emic.f` |

Two complementary modes this quarter: one contributor clearing technical debt across the whole codebase (net lines *down* despite 77 files reworked — legacy `.f` replaced by leaner modules), two adding new science on top. *Line churn rewards big rewrites and undercounts subtle fixes — read it as activity, not value. Excludes 4 merge commits; code directories only (`src/`, `src_cpl/`, `tools/`).*

---

## `devel` branch — A Glimpse of What's Coming

*Most of `devel` this quarter is `main` merged forward, so results are consistent across both branches. Two developments, though, are `devel`-only — early foundations for the next release.*

**A new modular coupler** `[in progress]` — The monolithic `emic.f` main routine (~1000 lines) is being replaced by a proper component-based coupling framework: a component registry, a coupler core, per-component modules (FROG already ported), and a coupling timer. OCYCC/CYCC is partially ported and Climate Indices are wired into the new coupler. This is the structural backbone for cleaner, more maintainable model coupling going forward.

**NAG → FFTW3** `[landed on devel]` — The proprietary NAG dependency for Fourier transforms is gone (~3300 lines removed), replaced by the open FFTW3 library. This removes a licensing barrier and should make the model easier to build and share across machines.

**In short:** the near-term theme on `devel` is *architecture* — a modern coupling structure and one fewer proprietary dependency — laying groundwork rather than changing science. Expect these to firm up before merging to `main`.
