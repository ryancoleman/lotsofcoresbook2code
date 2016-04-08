# Have to add top_builddir/include to pick up made qdp_config_internal.h

INCFLAGS = -I@top_srcdir@/include -I@top_builddir@/include \
           -I@top_srcdir@/other_libs/qio/include \
           -I@top_builddir@/other_libs/qio/include \
           -I@top_srcdir@/other_libs/qio/other_libs/c-lime/include \
           -I@top_builddir@/other_libs/qio/other_libs/c-lime/include

if QDP_USE_LIBXML2
INCFLAGS += -I@top_srcdir@/other_libs/xpath_reader/include
endif

if BUILD_LIBINTRIN
INCFLAGS += -I@top_srcdir@/other_libs/@LIBINTRIN_DIR@/include -I@top_builddir@/other_libs/@LIBINTRIN_DIR@/include
endif

AM_CXXFLAGS = $(INCFLAGS) @BAGEL_QDP_CXXFLAGS@ @HDF5_CXXFLAGS@ @LIBXML2_CXXFLAGS@ @QMP_CFLAGS@  @QMT_CXXFLAGS@

lib_LIBRARIES = libqdp.a

## Base set of sources and conditionally included sources
libqdp_a_SOURCES = qdp_map.cc qdp_subset.cc qdp_random.cc \
	qdp_layout.cc qdp_io.cc qdp_byteorder.cc qdp_util.cc \
	qdp_stdio.cc \
        qdp_profile.cc qdp_strnlen.cc qdp_crc32.cc \
        qdp_stopwatch.cc \
        qdp_rannyu.cc

if QDP_USE_LIBXML2
libqdp_a_SOURCES += qdp_xmlio.cc qdp_iogauge.cc qdp_qdpio.cc qdp_qio_strings.cc qdp_map_obj_disk.cc
endif

if QDP_USE_HDF5
libqdp_a_SOURCES += qdp_hdf5.cc
endif

# Scalar	
if ARCH_SCALAR
libqdp_a_SOURCES += qdp_scalar_init.cc qdp_scalar_layout.cc \
	qdp_scalar_specific.cc qdp_scalarsite_specific.cc
endif

# Parallel-scalar
if ARCH_PARSCALAR
libqdp_a_SOURCES += qdp_parscalar_init.cc qdp_parscalar_layout.cc \
	qdp_parscalar_specific.cc \
	qdp_scalarsite_specific.cc
endif

# Optimized code using sse extensions
if QDP_USE_SCALAR_SSE
libqdp_a_SOURCES += \
	scalarsite_sse/qdp_scalarsite_sse.cc \
	scalarsite_sse/qdp_scalarsite_linalg_double.cc \
	scalarsite_sse/sse_blas_local_sumsq_double.cc \
	scalarsite_sse/sse_blas_local_vcdot_double.cc \
	scalarsite_sse/sse_blas_local_vcdot_real_double.cc \
	scalarsite_sse/sse_blas_vaxmbyz4_double.cc \
	scalarsite_sse/sse_blas_vaxmyz4_double.cc \
	scalarsite_sse/sse_blas_vaxpbyz4_double.cc \
	scalarsite_sse/sse_blas_vaxpy4_double.cc \
	scalarsite_sse/sse_blas_vaypx4_double.cc \
	scalarsite_sse/sse_blas_vscal4_double.cc \
	scalarsite_sse/sse_linalg_m_eq_scal_m_double.cc \
	scalarsite_sse/sse_linalg_m_peq_m_double.cc \
	scalarsite_sse/sse_linalg_m_eq_mm_double.cc \
	scalarsite_sse/sse_linalg_m_eq_mh_double.cc \
	scalarsite_sse/sse_linalg_m_eq_hm_double.cc \
	scalarsite_sse/sse_linalg_m_eq_hh_double.cc
endif


# Scalar with Vector extensions
if ARCH_SCALARVEC
libqdp_a_SOURCES += qdp_scalarvec_init.cc qdp_scalarvec_layout.cc \
	qdp_scalarvec_specific.cc qdp_scalarvecsite_specific.cc
endif

# Parallel scalar with Vector extensions
if ARCH_PARSCALARVEC
libqdp_a_SOURCES += qdp_parscalarvec_init.cc qdp_parscalarvec_layout.cc \
	qdp_parscalarvec_specific.cc qdp_scalarvecsite_specific.cc
endif

# Optimized code using sse extensions
if QDP_USE_SCALARVEC_SSE
libqdp_a_SOURCES += qdp_scalarvecsite_sse.cc
endif
 
if BUILD_DEFAULT_ALLOCATOR
libqdp_a_SOURCES += qdp_default_allocator.cc
endif 

if BUILD_BGQ_THREADBIND
libqdp_a_SOURCES += qdp_bgq_threadbind.cc
else
libqdp_a_SOURCES += qdp_generic_threadbind.cc
endif
