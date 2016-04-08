
AM_CFLAGS = -I@top_srcdir@/include -I@top_srcdir@/other_libs/c-lime/include \
            -I@top_builddir@/include -I@top_builddir@/other_libs/c-lime/include @QMP_CFLAGS@

lib_LIBRARIES = libqio.a

QIO_SRCS = \
   qio/QIO_close_read.c \
   qio/QIO_close_write.c \
   qio/QIO_info.c \
   qio/QIO_info_ildg_record.c \
   qio/QIO_info_private.c \
   qio/QIO_info_usqcd_gauge.c \
   qio/QIO_info_usqcd_ksprop.c \
   qio/QIO_info_usqcd_prop.c \
   qio/QIO_next_record.c \
   qio/QIO_open_read.c \
   qio/QIO_open_write.c \
   qio/QIO_read.c \
   qio/QIO_read_record_data.c \
   qio/QIO_read_record_info.c \
   qio/QIO_string.c \
   qio/QIO_utils.c \
   qio/QIO_write.c \
   qio/QIO_host_file_conversion.c \
   qio/QIO_host_utils.c

DML_GENERIC = \
   dml/DML_crc32.c \
   dml/DML_utils.c

DML_PARSCALAR = ${OBJECTS} dml/DML_parscalar.c dml/DML_route.c
DML_SCALAR = ${OBJECTS} dml/DML_scalar.c

LRL_SRCS = lrl/LRL_main.c

GENERIC_SRCS = $(QIO_SRCS) $(DML_GENERIC) $(LRL_SRCS)

if ARCH_SCALAR
  libqio_a_SOURCES = $(GENERIC_SRCS) $(DML_SCALAR)
  EXTRA_libqio_a_SOURCES = $(DML_PARSCALAR)
endif

if ARCH_PARSCALAR
  libqio_a_SOURCES = $(GENERIC_SRCS) $(DML_PARSCALAR)
  EXTRA_libqio_a_SOURCES = $(DML_SCALAR)
endif
