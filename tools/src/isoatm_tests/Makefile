# Standalone tests of the water isotope modules (src/isoatm), ROADMAP step 2b.
#
#   make          build the test programs
#   make check    build and run them (exit code 1 on failure)
#
# The isotope code is compiled with WISOATM = 1: a copy of choixcomposantes.h with that setting is made in inc/.
# Flags as in ../fvtrans_tests (no double-precision promotion: the modules use explicit kinds).

ifeq ($(origin FC),default)
FC      := gfortran
endif
TOP     := ../../..
SRC_GC  := $(TOP)/src/lbm/global_constants_mod.f90
SRC_AT  := $(TOP)/src/ecbilt/comatm.f90
SRC_ISO := $(TOP)/src/isoatm
FLAGDIR := $(TOP)/defs/ecbilt_clio/global_flags
INC     := -Iinc -I$(TOP)/src/clio

FFLAGS  ?= -O2 -g -cpp -std=f2008 -Wall -Wextra -ffree-line-length-132 -fimplicit-none \
           -fcheck=all -ffpe-trap=invalid,zero,overflow

OBJS    := global_constants_mod.o comatm.o iso_alphas.o
PROGS   := test_iso_alphas

.PHONY: all check clean

all: $(PROGS)

inc/choixcomposantes.h: $(FLAGDIR)/choixcomposantes.h
	mkdir -p inc
	cp $(FLAGDIR)/*.h inc/
	sed -i 's/^#define WISOATM 0/#define WISOATM 1/' $@

global_constants_mod.o: $(SRC_GC)
	$(FC) $(FFLAGS) -std=gnu -c $< -o $@

comatm.o: $(SRC_AT) inc/choixcomposantes.h global_constants_mod.o
	$(FC) $(FFLAGS) -std=gnu -Wno-unused $(INC) -c $< -o $@

iso_alphas.o: $(SRC_ISO)/iso_alphas.f90 inc/choixcomposantes.h comatm.o
	$(FC) $(FFLAGS) $(INC) -c $< -o $@

test_iso_alphas: test_iso_alphas.f90 $(OBJS)
	$(FC) $(FFLAGS) $(INC) $^ -o $@

check: $(PROGS)
	@for p in $(PROGS); do ./$$p || exit 1; done

clean:
	rm -rf *.o *.mod inc $(PROGS)
