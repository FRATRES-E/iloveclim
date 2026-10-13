# Standalone tests of the FV transport modules (src/fvtrans), ROADMAP step 1.1.
#
#   make                      build the test programs
#   make check COEF=<path>    build and run them; COEF is ECBilt's inputdata/coef.dat
#   make suite                build the evaluation suite fvt_suite (needs netCDF-Fortran, nf-config)
#   ./run_suite.sh <outdir>   run it at T21, T42, T85 and plot (python3 with numpy, matplotlib, netCDF4)
#
# Flags follow the model's gfortran.options where they matter (double precision default, big-endian unformatted
# with 4-byte record markers), plus the style sheet's -std=f2008 -Wall -Wextra.

ifeq ($(origin FC),default)
FC      := gfortran
endif
TOP     := ../../..
SRC_FVT := $(TOP)/src/fvtrans
SRC_GC  := $(TOP)/src/lbm/global_constants_mod.f90
INC     := -I$(TOP)/defs/ecbilt_clio/global_flags -I$(TOP)/src/clio

FFLAGS  ?= -O2 -g -cpp -std=f2008 -Wall -Wextra -ffree-line-length-132 -fimplicit-none \
           -fconvert=big-endian -frecord-marker=4 -fcheck=all -ffpe-trap=invalid,zero,overflow
COEF    ?= $(TOP)/inputs/ecbilt/coef.dat

OBJS    := global_constants_mod.o fvt_legendre_mod.o fvt_grid_mod.o fvt_ffsl_mod.o fvt_specflux_mod.o
PROGS   := test_fvt_grid test_fvt_ffsl test_fvt_specflux test_fvt_specref

# evaluation suite: optimised, no run-time checks, objects kept apart from the test build
NCFFLAGS := $(shell nf-config --fflags 2>/dev/null)
NCLIBS   := $(shell nf-config --flibs 2>/dev/null)
SFLAGS   ?= -O3 -cpp -std=f2008 -Wall -Wextra -ffree-line-length-132 -fimplicit-none
SDIR     := suite_build
SSRC     := $(SRC_GC) $(SRC_FVT)/fvt_legendre_mod.f90 $(SRC_FVT)/fvt_grid_mod.f90 $(SRC_FVT)/fvt_ffsl_mod.f90 \
            fvt_testcases_mod.f90 fvt_ncout_mod.f90 fvt_specref_mod.f90

.PHONY: all check suite clean

all: $(PROGS)

global_constants_mod.o: $(SRC_GC)
	$(FC) $(FFLAGS) -std=gnu -c $< -o $@

%.o: $(SRC_FVT)/%.f90
	$(FC) $(FFLAGS) $(INC) -c $< -o $@

fvt_legendre_mod.o: global_constants_mod.o
fvt_grid_mod.o:     fvt_legendre_mod.o
fvt_ffsl_mod.o:     fvt_grid_mod.o
fvt_specflux_mod.o: fvt_grid_mod.o

# the tests compare reals bitwise on purpose (exact poles, symmetry)
test_fvt_grid: test_fvt_grid.f90 $(OBJS)
	$(FC) $(FFLAGS) -Wno-compare-reals $(INC) $< $(OBJS) -o $@

test_fvt_ffsl: test_fvt_ffsl.f90 $(OBJS)
	$(FC) $(FFLAGS) $(INC) $< $(OBJS) -o $@

# ECBilt's NAG FFT (legacy fixed form), needed to reproduce the model transforms
nag.o: $(TOP)/src/ecbilt/nag.f
	$(FC) -O2 -cpp -std=legacy -w $(INC) -c $< -o $@

test_fvt_specflux: test_fvt_specflux.f90 $(OBJS) nag.o
	$(FC) $(FFLAGS) $(INC) $< $(OBJS) nag.o -o $@

# spectral baseline (current ECBilt chain), a test-side module
fvt_specref_mod.o: fvt_specref_mod.f90 fvt_grid_mod.o
	$(FC) $(FFLAGS) $(INC) -c $< -o $@

test_fvt_specref: test_fvt_specref.f90 $(OBJS) fvt_specref_mod.o nag.o
	$(FC) $(FFLAGS) $(INC) $< $(OBJS) fvt_specref_mod.o nag.o -o $@

check: all
	./test_fvt_grid $(COEF)
	./test_fvt_ffsl
	./test_fvt_specflux $(COEF)
	./test_fvt_specref

suite: fvt_suite

fvt_suite: fvt_suite.f90 $(SSRC)
	mkdir -p $(SDIR)
	cd $(SDIR) && for f in $(SSRC); do $(FC) $(SFLAGS) $(NCFFLAGS) -I../$(TOP)/defs/ecbilt_clio/global_flags \
	  -I../$(TOP)/src/clio -c ../$$f || exit 1; done
	cd $(SDIR) && $(FC) -O2 -cpp -std=legacy -w -I../$(TOP)/defs/ecbilt_clio/global_flags -I../$(TOP)/src/clio \
	  -c ../$(TOP)/src/ecbilt/nag.f
	# compiled inside $(SDIR) so that only its modules are visible (not stale ones from the test build here)
	cd $(SDIR) && $(FC) $(SFLAGS) $(NCFFLAGS) -I../$(TOP)/defs/ecbilt_clio/global_flags -I../$(TOP)/src/clio \
	  ../$< *.o $(NCLIBS) -o ../$@

clean:
	rm -rf *.o *.mod $(PROGS) fvt_suite $(SDIR)
