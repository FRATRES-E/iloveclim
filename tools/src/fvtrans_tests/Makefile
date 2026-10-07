# Standalone tests of the FV transport modules (src/fvtrans), ROADMAP step 1.1.
#
#   make                      build the test programs
#   make check COEF=<path>    build and run them; COEF is ECBilt's inputdata/coef.dat
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

OBJS    := global_constants_mod.o fvt_legendre_mod.o fvt_grid_mod.o
PROGS   := test_fvt_grid

.PHONY: all check clean

all: $(PROGS)

global_constants_mod.o: $(SRC_GC)
	$(FC) $(FFLAGS) -std=gnu -c $< -o $@

%.o: $(SRC_FVT)/%.f90
	$(FC) $(FFLAGS) $(INC) -c $< -o $@

fvt_legendre_mod.o: global_constants_mod.o
fvt_grid_mod.o:     fvt_legendre_mod.o

# the tests compare reals bitwise on purpose (exact poles, symmetry)
test_fvt_grid: test_fvt_grid.f90 $(OBJS)
	$(FC) $(FFLAGS) -Wno-compare-reals $(INC) $< $(OBJS) -o $@

check: all
	./test_fvt_grid $(COEF)

clean:
	rm -f *.o *.mod $(PROGS)
