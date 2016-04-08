## Process this file with automake to produce Makefile.in

AM_CFLAGS  = -I@top_srcdir@/include -I@top_srcdir@/other_libs/c-lime/include \
             -I@top_builddir@/include -I@top_builddir@/other_libs/c-lime/include @QMP_CFLAGS@
AM_LDFLAGS = -L../lib -L../other_libs/c-lime/lib @QMP_LDFLAGS@
LDADD      = -lqio -llime @QMP_LIBS@ @LIBS@ -lm

# programs for all architectures
check_PROGRAMS  = 
bin_PROGRAMS =

if USING_QMP
check_PROGRAMS += qio-test1         \
                  qio-test2         \
                  qio-test3         \
                  qio-test4         \
                  qio-test5         \
                  qio-test6         \
                  qio-test7         \
                  qio-test8
endif

# Only scalar arch
if ARCH_SCALAR
bin_PROGRAMS += qio-convert-mesh-singlefs
bin_PROGRAMS += qio-convert-mesh-pfs
bin_PROGRAMS += qio-convert-mesh-ppfs
bin_PROGRAMS += qio-convert-nersc
check_PROGRAMS += 
endif

# Only parallel arch
if ARCH_PARSCALAR
check_PROGRAMS +=
bin_PROGRAMS += qio-copy-mesh-ppfs
endif

EXTRA_DIST = qio-host-test.sh layout_test layout_test_ppfs binary_host_test

ADD_MPP_SOURCE = qio-test.c qio-test-util.c layout_hyper.c qio-test.h
ADD_MESH_SOURCE =  qio-convert-mesh.c qio-mesh-utilities.c layout_hyper_mesh.c qio-convert-mesh.h
ADD_COPY_SOURCE = qio-mesh-utilities.c qio-convert-mesh.h

qio_test1_SOURCES = qio-test1.c ${ADD_MPP_SOURCE}
qio_test2_SOURCES = qio-test2.c ${ADD_MPP_SOURCE}
qio_test3_SOURCES = qio-test3.c ${ADD_MPP_SOURCE}
qio_test4_SOURCES = qio-test4.c ${ADD_MPP_SOURCE}
qio_test5_SOURCES = qio-test5.c ${ADD_MPP_SOURCE}
qio_test6_SOURCES = qio-test6.c ${ADD_MPP_SOURCE}
qio_test7_SOURCES = qio-test7.c ${ADD_MPP_SOURCE}
qio_test8_SOURCES = qio-test8.c ${ADD_MPP_SOURCE}
qio_convert_mesh_singlefs_SOURCES = qio-convert-mesh-singlefs.c ${ADD_MESH_SOURCE}
qio_convert_mesh_pfs_SOURCES = qio-convert-mesh-pfs.c ${ADD_MESH_SOURCE}
qio_convert_mesh_ppfs_SOURCES = qio-convert-mesh-ppfs.c ${ADD_MESH_SOURCE}
qio_convert_nersc_SOURCES = qio-convert-nersc.c
qio_copy_mesh_ppfs_SOURCES = qio-copy-mesh-ppfs.c ${ADD_COPY_SOURCE}

DEPENDENCIES = ../lib/libqio.a ../other_libs/c-lime/lib/liblime.a
${check_PROGRAMS}: ${DEPENDENCIES}
${bin_PROGRAMS}: ${DEPENDENCIES}
