TOPDIR=@top_srcdir@
BUILDDIR=@top_builddir@

AM_CXXFLAGS= -I$(TOPDIR)/include -I$(BUILDDIR)/include @QPHIX_AUX_CFLAGS@

lib_LIBRARIES = libqphix_solver.a

libqphix_solver_a_SOURCES = \
	print_utils.cc

if QPHIX_BUILD_QPX
libqphix_solver_a_SOURCES += bgq_threadbind.cc
else
libqphix_solver_a_SOURCES += generic_threadbind.cc
endif
