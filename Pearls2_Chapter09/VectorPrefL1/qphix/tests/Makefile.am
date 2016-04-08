INCFLAGS=-I@top_srcdir@/include -I@top_builddir@/include

AM_CXXFLAGS=$(INCFLAGS) @QPHIX_AUX_CFLAGS@
AM_LDFLAGS=-L@top_builddir@/lib @LDFLAGS@ @QPHIX_AUX_LDFLAGS@
LDADD = -lqphix_solver @LIBS@ @QPHIX_AUX_LIBS@

check_PROGRAMS = time_dslash_noqdp
if QPHIX_BUILD_WITH_QDP
check_PROGRAMS += t_dslash
endif

if QPHIX_BUILD_CLOVER
check_PROGRAMS += time_clov_noqdp
if QPHIX_BUILD_WITH_QDP
check_PROGRAMS +=  t_clov_dslash t_clov_invert_from_file
endif
endif

# check_PROGRAMS += test_blas

bin_PROGRAMS = $(check_PROGRAMS)

FRAMEWORK_SRCS =  unittest.h \
	dslashm_w.h \
	dslashm_w.cc \
	reunit.h \
	invbicgstab_test.h \
	reunit.cc

t_dslash_SOURCES = $(FRAMEWORK_SRCS) \
	testDslashFull.h \
	testDslashFull.cc \
	t_dslash.cc

t_dslash_DEPENDENCIES = build_lib

time_dslash_noqdp_SOURCES= \
	time_dslash_noqdp.cc \
	timeDslashNoQDP.h \
	timeDslashNoQDP.cc

if QPHIX_BUILD_CLOVER
t_clov_dslash_SOURCES = $(FRAMEWORK_SRCS) \
	testClovDslashFull.h \
	testClovDslashFull.cc \
	t_clov_dslash.cc \
	mesfield.h \
	mesfield.cc \
	clover_term_qdp_w.h \
	clover_fermact_params_w.h

t_clov_invert_from_file_SOURCES = $(FRAMEWORK_SRCS) \
	testClovInvertFromFile.h \
	testClovInvertFromFile.cc \
	t_clov_invert_from_file.cc \
	mesfield.h \
	mesfield.cc \
	clover_term_qdp_w.h \
	clover_fermact_params_w.h

time_clov_noqdp_SOURCES= \
	time_clov_noqdp.cc \
	timeClovNoQDP.h \
	timeClovNoQDP.cc

t_clov_dslash_DEPENDENCIES = build_lib
t_clov_invert_from_file_DEPENDENCIES = build_lib
time_clov_noqdp_DEPENDENCIES = build_lib

endif

#test_blas_SOURCES=$(FRAMEWORK_SRCS) \
#	test_blas.cc \
#	testBlas.h \
#	testBlas.cc
#
#test_blas_DEPENDENCIES= build_lib

DEPENDENCIES = build_lib
${bin_PROGRAMS}: ${DEPENDENCIES}
${check_PROGRAMS}: ${DEPENDENCIES}
${EXTRA_PROGRAMS}: ${DEPENDENCIES}

.PHONY: build_lib
	cd @top_builddir@/lib; $(MAKE)
