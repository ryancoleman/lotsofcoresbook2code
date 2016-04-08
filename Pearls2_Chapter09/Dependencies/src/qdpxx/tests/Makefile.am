#
# Boiler Plate stuff: Autoconf Substitutions etc
#
## Should only need @QMP_CFLAGS@, @QMP_LDFLAGS@ and @QMP_LIBS@
## n.b. Only defined for some parallel architectures
AM_CXXFLAGS = -I@top_srcdir@/include -I@top_builddir@/include \
              -I@top_srcdir@/other_libs/qio/include \
	      -I@top_builddir@/other_libs/qio/include \
              -I@top_srcdir@/other_libs/qio/other_libs/c-lime/include \
              -I@top_builddir@/other_libs/qio/other_libs/c-lime/include \
              -I@top_srcdir@/other_libs/xpath_reader/include \
              @BAGEL_QDP_CXXFLAGS@ @HDF5_CXXFLAGS@ @LIBXML2_CXXFLAGS@ @QMP_CFLAGS@ @QMT_CXXFLAGS@

if BUILD_LIBINTRIN
AM_CXXFLAGS += -I@top_srcdir@/other_libs/@LIBINTRIN_DIR@/include
AM_CXXFLAGS += -I@top_builddir@/other_libs/@LIBINTRIN_DIR@/include
endif


AM_LDFLAGS = -L@top_builddir@/lib \
             -L@top_builddir@/other_libs/qio/lib \
             -L@top_builddir@/other_libs/qio/other_libs/c-lime/lib \
             -L@top_builddir@/other_libs/xpath_reader/lib \
             @LDFLAGS@ @BAGEL_QDP_LDFLAGS@ @HDF5_LDFLAGS@ @QMP_LDFLAGS@ @QMT_LDFLAGS@

if BUILD_LIBINTRIN
AM_LDFLAGS += -L@top_builddir@/other_libs/libintrin/lib
endif

LDADD = -lqdp -lXPathReader -lxmlWriter -lqio -llime @BAGEL_QDP_LIBS@ @HDF5_LIBS@ @LIBXML2_LIBS@ @QMP_LIBS@ @QMT_LIBS@ @LIBS@

if BUILD_LIBINTRIN
LDADD += -lintrin
endif

#
# Local Headers
#

#
# The programs to build
# 
check_PROGRAMS = test_vaxpy_double time_vaxpy_double test_matmat_double test_cmul time_matmat_double


# The program and its dependencies
test_HDRS=unittest.h \
	testvol.h

test_vaxpy_double_SOURCES = $(test_HDRS) \
	testVaypxDouble.h \
	testVaypxDouble.cc \
	testVaxpyDouble.h \
	testVaxpyDouble.cc \
	testVaxmyDouble.h \
	testVaxmyDouble.cc \
	testVaxpbyDouble.h \
	testVaxpbyDouble.cc \
	testVScalDouble.h \
	testVScalDouble.cc \
	testLocalSumSqDouble.h \
	testLocalSumSqDouble.cc \
	testLocalVcdotRealDouble.h \
	testLocalVcdotRealDouble.cc \
	testLocalVcdotDouble.h \
	testLocalVcdotDouble.cc \
	test_vaxpy_double.cc
test_vaxpy_double_DEPENDENCIES = build_libs


time_vaxpy_double_SOURCES = $(test_HDRS) \
	timeVaxpyDouble.h \
	timeVaxpyDouble.cc \
	timeVaxpby.cc \
	time_vaxpy_double.cc

time_vaxpy_double_DEPENDENCIES = build_libs

test_matmat_double_SOURCES = test_matmat_double.cc \
	testMatScalMultDouble.h \
	testMatScalMultDouble.cc \
	testMatPeqMatDouble.h \
	testMatPeqMatDouble.cc \
	testMatEqMatMatDouble.h \
	testMatEqMatMatDouble.cc \
	testMatEqMatHermDouble.h \
	testMatEqMatHermDouble.cc \
	testMatEqHermMatDouble.h \
	testMatEqHermMatDouble.cc \
	testMatEqHermHermDouble.h \
	testMatEqHermHermDouble.cc 

test_cmul_SOURCES = test_cmul.cc \
	testCMulDouble.h \
	testCMulDouble.cc

test_matmat_double_DEPENDENCIES = build_libs

time_matmat_double_SOURCES = time_matmat_double.cc \
	timeMatEqMatMatDouble.h \
	timeMatEqMatMatDouble.cc \
	timeMatEqMatHermDouble.h \
	timeMatEqMatHermDouble.cc \
	timeMatEqHermMatDouble.h \
	timeMatEqHermMatDouble.cc \
	timeMatEqHermHermDouble.h \
	timeMatEqHermHermDouble.cc

time_matmat_double_DEPENDENCIES = build_libs
# build lib is a target that goes tot he build dir of the library and 
# does a make to make sure all those dependencies are OK. In order
# for it to be done every time, we have to make it a 'phony' target
.PHONY: build_libs
build_libs:
	cd @top_builddir@/lib ; $(MAKE)
