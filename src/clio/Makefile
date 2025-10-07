# Generated automatically from Makefile.src.in by configure.
#
#       $Id: Makefile.src.in,v 1.5 2001/08/01 10:55:33 severijn Exp $
#

srcdir          = .

top_srcdir      = ..


include make.macros-gfortran


EXTRACFLAGS     =
EXTRACPPFLAGS   =
# EXTRAFFLAGS     = -I${FOPTSEP}. -I${FOPTSEP}${top_srcdir}/include ${NETCDFINCLUDE} ${UDUNITSINCLUDE}
EXTRALDFLAGS    =
EXTRALIBS       = ${NETCDFLIB} ${UDUNITSLIB}


FSRCS   = \
        comemic_mod.f90                 \
        para0_mod.f90                   \
        iceberg_coupl_mod.f90           \
        bloc0_mod.f90                   \
        para_mod.f90                    \
        bloc_mod.f90                    \
        const_mod.f90                   \
        datadc_mod.f90                  \
        densit_mod.f90                  \
        dynami_mod.f90                  \
        ice_mod.f90                     \
        iso_param_mod.f90               \
        isoslope_mod.f90                \
        reper_mod.f90                   \
        thermo_mod.f90                  \
        trace_mod.f90                   \
        type_mod.f90                    \
        moment_mod.f90                  \
        varno_mod.f90                   \
        vareq_mod.f90                   \
        netcdf_out_mod.f90              \
        brines_mod.f90                  \
        geogra_mod.f90                  \
        correct_mod.f90                 \
        raccord_mod.f90                 \
        staocb_mod.f90                  \
        informe_mod.f90                 \
        inforun_mod.f90                 \
        forcat_mod.f90                  \
        dynamifunc_mod.f90              \
        icdyna_mod.f90                  \
        icdadv_mod.f90                  \
        init_isoocn_mod.f90             \
        clio_utils_mod.f90              \
        redrunb_mod.f90                 \
        etat_mod.f90                    \
        seaice_mod.f90                  \
        thersf_mod.f90                  \
        diftur_mod.f90                  \
        engtur_mod.f90                  \
        conti3d_mod.f90                 \
        flucor_mod.f90                  \
        isopycnal_mod.f90               \
        slopez_mod.f90                  \
        rahmflx_mod.f90                 \
        scalars_mod.f90                 \
        uv_funcs_mod.f90                \
        barot_mod.f90                   \
        clio_routines_mod.f90           \
        clio_mod.f90


TFSRCS  =

CMD     =
OBJS    = ${FSRCS:.f90=.o}
MODS    = ${FSRCS:.f90=.mod}

LIBRARY = libclio_4_0_0.a
LIBOBJS = ${FSRCS:.f90=.o}


distfiles = \
        ${FSRCS}        \
        ${TFSRCS}

dist_garbage = \
        Makefile

garbage = \
        ${LIBOBJS:.o=.l}        \
        ${LIBOBJS:.o=.lst}      \
        ${LIBOBJS:.o=.T}        \
        ${OBJS:.o=.l}           \
        ${OBJS:.o=.lst}         \
        ${OBJS:.o=.T}

more_garbage =


include make.rules


run:
		@echo Nothing to do for run

depend:
		${SED} -e '/^# DO NOT DELETE THIS LINE -- MAKE DEPEND DEPENDS ON IT/,$$d' < Makefile > Makefile.new
		echo '# DO NOT DELETE THIS LINE -- MAKE DEPEND DEPENDS ON IT' >> Makefile.new
		${FMAKEDEPEND} ${FSRCS} ${TFSRCS} >> Makefile.new
		mv Makefile Makefile.old
		mv Makefile.new Makefile

install:
		${INSTALL_DATA} ${LIBRARY} ${libdir}
		@for header in *.${MODEXT}; do  \
                echo Installing $$header ...;           \
                ${INSTALL_DATA} $$header ${includedir}; \
        done

install-strip: install

main: library ${CMD}

uninstall:
		@echo Not yet implemented


#
# Dependencies
#

${CMD}: ${OBJS} ${LIBRARY}
		${LD} ${LDFLAGS} -o $@ ${OBJS} $(MODS) ${LIBRARY} ${LIBS}
#
# The End of All Things (Op. Cit.)
