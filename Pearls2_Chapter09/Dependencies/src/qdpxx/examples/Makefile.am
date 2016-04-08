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
AM_CXXFLAGS += -I@top_srcdir@/other_libs/@LIBINTRIN_DIR@/include \
               -I@top_builddir@/other_libs/@LIBINTRIN_DIR@/include
endif

if BUILD_FILEDB
AM_CXXFLAGS += -I@top_srcdir@/other_libs/filedb/src \
               -I@top_srcdir@/other_libs/filedb/filehash
endif

AM_LDFLAGS = -L@top_builddir@/lib

if QDP_USE_LIBXML2
AM_LDFLAGS += -L@top_builddir@/other_libs/qio/lib \
             -L@top_builddir@/other_libs/qio/other_libs/c-lime/lib \
             -L@top_builddir@/other_libs/xpath_reader/lib
endif
AM_LDFLAGS += @LDFLAGS@ @BAGEL_QDP_LDFLAGS@ @HDF5_LDFLAGS@ @QMP_LDFLAGS@ @QMT_LDFLAGS@

if BUILD_LIBINTRIN
AM_LDFLAGS += -L@top_builddir@/other_libs/libintrin/lib
endif

if BUILD_FILEDB
AM_LDFLAGS += -L@top_builddir@/other_libs/filedb/src \
              -L@top_builddir@/other_libs/filedb/filehash
endif

LDADD = -lqdp

if QDP_USE_HDF5
LDADD += @HDF5_LIBS@
endif

if QDP_USE_LIBXML2
LDADD += -lXPathReader -lxmlWriter -lqio -llime
endif
LDADD += @BAGEL_QDP_LIBS@ @LIBXML2_LIBS@ @QMP_LIBS@ @QMT_LIBS@ @LIBS@

if BUILD_LIBINTRIN
LDADD += -lintrin
endif

if BUILD_FILEDB
LDADD += -lfiledb -lfilehash
endif

#
# Local Headers
#
HDRS=examples.h 

#
# The programs to build
# 
check_PROGRAMS = t_skeleton t_io t_mesplq t_db \
      t_xml t_entry t_nersc t_shift t_exotic t_basic t_qio \
      t_cugauge t_transpose_spin t_partfile t_su3 \
      t_map_obj_disk t_map_obj_memory t_clov_force

EXTRA_PROGRAMS  = t_qio_factory t_gsum t_iprod


if BUILD_WILSON_EXAMPLES
check_PROGRAMS += t_dslashm t_formfac t_spectrum t_qdp t_linalg 
EXTRA_PROGRAMS += t_subtype t_foo t_blas t_cblas t_blas_g5 t_blas_g5_2 t_blas_g5_3 t_spinproj t_spinproj2
endif

EXTRA_PROGRAMS += nersc2ildg test_ildglat lhpc2ildg

if BUILD_STAGGERED_EXAMPLES
check_PROGRAMS += 
endif

#
# The program and its dependencies
t_skeleton_SOURCES = t_skeleton.cc
t_foo_SOURCES = t_foo.cc $(HDRS)
t_subtype_SOURCES = t_subtype.cc $(HDRS)
t_xml_SOURCES = t_xml.cc $(HDRS)
t_qio_SOURCES = t_qio.cc $(HDRS)
t_qio_DEPENDENCIES = build_lib rebuild_other_libs
t_qio_factory_SOURCES = t_qio_factory.cc $(HDRS)
t_qio_factory_DEPENDENCIES = build_lib rebuild_other_libs
t_partfile_SOURCES = t_partfile.cc $(HDRS) \
	unittest.h testOpenPartFile.h
t_partfile_DEPENDENCIES = build_lib rebuild_other_libs

t_basic_SOURCES = t_basic.cc
t_shift_SOURCES = t_shift.cc $(HDRS)
t_io_SOURCES = t_io.cc $(HDRS)
t_entry_SOURCES = t_entry.cc $(HDRS)
t_nersc_SOURCES = t_nersc.cc mesplq.cc reunit.cc $(HDRS)
t_dslashm_SOURCES = t_dslashm.cc dslashm_w.cc $(HDRS)
t_mesplq_SOURCES = t_mesplq.cc mesplq.cc reunit.cc $(HDRS)
t_exotic_SOURCES = t_exotic.cc rgauge.cc taproj.cc expm12.cc reunit.cc $(HDRS)
t_formfac_SOURCES = t_formfac.cc formfac_w.cc $(HDRS)
t_spectrum_SOURCES = t_spectrum.cc baryon_w.cc mesons_w.cc mesplq.cc $(HDRS)
t_qdp_SOURCES =  t_qdp.cc formfac_w.cc dslashm_w.cc baryon_w.cc \
	 mesons_w.cc mesplq.cc reunit.cc $(HDRS)
t_linalg_SOURCES =  t_linalg.cc linalg1.cc linalg.h

t_blas_SOURCES = t_blas.cc blas1.cc blas1.h
t_blas_DEPENDENCIES = build_lib

t_gsum_SOURCES = t_gsum.cc
t_gsum_DEPENDENCIES = build_lib

t_iprod_SOURCES = t_iprod.cc
t_iprod_DEPENDENCIES = build_lib

t_cblas_SOURCES= t_cblas.cc cblas1.cc cblas1.h 

t_clov_force_SOURCES=t_clov_force.cc reunit.cc $(HDRS)
t_clov_force_DEPENDENCIES= build_lib

t_db_SOURCES = t_db.cc $(HDRS)
t_map_obj_disk_SOURCES = t_map_obj_disk.cc $(HDRS)
t_map_obj_memory_SOURCES = t_map_obj_memory.cc $(HDRS)

t_blas_g5_SOURCES = t_blas_g5.cc $(HDRS)
t_blas_g5_2_SOURCES = t_blas_g5_2.cc $(HDRS)
t_blas_g5_3_SOURCES = t_blas_g5_3.cc $(HDRS)
t_spinproj_SOURCES = t_spinproj.cc 
t_spinproj2_SOURCES = t_spinproj2.cc 

t_cugauge_SOURCES = t_cugauge.cc mesplq.cc $(HDRS) 
t_transpose_spin_SOURCES = t_transpose_spin.cc $(HDRS)

nersc2ildg_SOURCES = nersc2ildg.cc $(HDRS) mesplq.cc
test_ildglat_SOURCES = test_ildglat.cc $(HDRS) mesplq.cc
test_ildglat_DEPENDENCIES = build_lib

t_su3_SOURCES = t_su3.cc 
t_su3_DEPENDENCIES = build_lib

lhpc2ildg_SOURCES = lhpc2ildg.cc $(HDRS) mesplq.cc
lhpc2ildg_DEPENDENCIES = build_lib

# build lib is a target that goes tot he build dir of the library and 
# does a make to make sure all those dependencies are OK. In order
# for it to be done every time, we have to make it a 'phony' target
.PHONY: build_lib
build_lib:
	cd @top_builddir@/lib ; $(MAKE)

.PHONY: rebuild_other_libs
rebuild_other_libs:
	cd @top_builddir@/other_libs ; $(MAKE)
