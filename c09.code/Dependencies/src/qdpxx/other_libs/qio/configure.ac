AC_PREREQ(2.53)
AC_INIT(qio, 2.4.2, detar@physics.utah.edu)
AC_CONFIG_AUX_DIR(config)

AC_CONFIG_SRCDIR([include/qio.h])
AC_CONFIG_HEADER([include/qio_config_internal.h])

AM_INIT_AUTOMAKE

dnl for saving the configure options
qio_conf_opts=

dnl
dnl --enable-parallel-io
dnl
AC_ARG_ENABLE(parallel-io,
  AC_HELP_STRING([--enable-parallel-io],
    [Only supported for parallel file systems with a parallel fseek]),
  [ac_parallel_io=1
   qio_conf_opts="$qio_conf_opts --enable-parallel-io"],
  [ac_parallel_io=0]
)

dnl
dnl --enable-parallel-arch argument
dnl
AC_ARG_ENABLE(parallel-arch,
  AC_HELP_STRING([--enable-parallel-arch=<arch>],
    [Build QIO for parallel architecture <arch>]),
  AC_SUBST(PARALLEL_ARCH, [$enableval]),
  AC_SUBST(PARALLEL_ARCH, [parscalar])
)

dnl
dnl --enable-qmp-route
dnl
AC_ARG_ENABLE(qmp-route,
  AC_HELP_STRING([--enable-qmp-route],
     [Enable usage of QMP_route if it exists, default is off]),
  [ac_qmp_route=$enableval
   qio_conf_opts="$qio_conf_opts --enable-qmp-route=$enableval"],
  [ac_qmp_route="no"] 
)

dnl
dnl --with-qmp=
dnl
AC_ARG_WITH(qmp,
  AC_HELP_STRING(
    [--with-qmp=DIR],
    [Build QIO on top of the QMP installed in DIR.]
  ),
  [QMP_HOME="${with_qmp}"]
)

dnl
dnl Enable buffering of partition output
dnl This can go away if the DML_partition_out is fixed to work
dnl on grid machines
dnl
AC_ARG_ENABLE(dml-output-buffering,
  AC_HELP_STRING(
     [--enable-dml-output-buffering],
     [Buffering enabled in partition output]),
   [AC_MSG_NOTICE([Enabling DML output buffering])
    AC_DEFINE([QIO_USE_DML_OUT_BUFFERING],[], [Enables DML output buffering of sites])],
   [AC_MSG_NOTICE([Disabled DML output buffering. Size only affects single site of data])]
)

dnl
dnl Enable changing DML buffersize
dnl
AC_ARG_ENABLE(dml-bufsize, 
   AC_HELP_STRING(
     [--enable-dml-bufsize=X],
     [Set the DML buffer size to X bytes]
   ),
  [AC_MSG_NOTICE([Setting DML_BUF_BYTES to $enableval])
   AC_DEFINE_UNQUOTED([QIO_DML_BUF_BYTES], [$enableval], [Maximum DML Buffer Size in bytes])
   qio_conf_opts="$qio_conf_opts --enable-dml-bufsize=$enableval"
  ],
  [AC_MSG_NOTICE([Setting DML_BUF_BYTES to 262144 bytes (256K) ]) 
   AC_DEFINE([QIO_DML_BUF_BYTES], [262144], [Maximum DML buffer size in bytes])
  ]
)

dnl
dnl Faster but not necessarily safter QMP route workaround
dnl
AC_ARG_ENABLE(fast-route,
   AC_HELP_STRING(
     [--enable-fast-route],
     [When not using QMP_route, use the slow version of DML route]
   ),
   [AC_MSG_NOTICE([Enabling Fast DML_route from J. Osborn])
    AC_DEFINE([QIO_USE_FAST_ROUTE],[], [Enable J. Osborns Fast DML route])
    qio_conf_opts="$qio_conf_opts --enable-fast-route"],
   [AC_MSG_NOTICE([Using Slow Route from Balint]) ]
)

dnl
dnl QIO only supports scalar and parscalar (basically, not-parallel and
dnl parallel) for other cases, map it to one of these.
dnl
case ${PARALLEL_ARCH} in
  scalarvec)
    AC_MSG_NOTICE([Mapping architecture ${PARALLEL_ARCH} to  scalar])
    AC_SUBST(LOCAL_PARALLEL_ARCH, [scalar])
    ;;
  parscalarvec)
    AC_MSG_NOTICE([Mapping architecture ${PARALLEL_ARCH} to  parscalar])
    AC_SUBST(LOCAL_PARALLEL_ARCH, [parscalar])
    ;;
  scalar)
    AC_SUBST(LOCAL_PARALLEL_ARCH, [scalar])
    ;;
  parscalar)
    AC_SUBST(LOCAL_PARALLEL_ARCH, [parscalar])
    ;;
  *)
    AC_MSG_ERROR([ Unknown value for --enable-parallel-arch ])
    ;;
esac

dnl
dnl Now check appropriate cases
dnl
case ${LOCAL_PARALLEL_ARCH} in
  scalar)
    AC_SUBST(ARCHDEF_SUBSTITUTION, ARCH_SCALAR)
    ;;
  parscalar)
    AC_SUBST(ARCHDEF_SUBSTITUTION, ARCH_PARSCALAR)
    ;;
  *)
    AC_MSG_ERROR([ Unknown value for --enable-parallel-arch ])
    ;;
esac

dnl The test code that can be run in scalar mode expects a scalar
dnl version of QMP.  So we always check for QMP.

AC_MSG_NOTICE([Checking for QMP ])

if test "X${QMP_HOME}X" = "XX" ; then
  AC_PATH_PROG([QMP_CONFIG], [qmp-config], [])
else
  AC_PATH_PROG([QMP_CONFIG], [qmp-config], [], [${QMP_HOME}/bin ${PATH}])
fi

if test "X${QMP_CONFIG}X" != "XX" ; then
  AC_MSG_NOTICE([Found QMP configuration program ${QMP_CONFIG}])
  # use the QMP CC unless it is explicitly set by the user
  if ! set | grep -q '^CC=' ; then
    export CC="`${QMP_CONFIG} --cc`"
    AC_SUBST(CC,  "$CC")
    AC_MSG_NOTICE([QMP compiler: ${CC}])
  fi
  AC_SUBST(QMP_CFLAGS,  "`${QMP_CONFIG} --cflags`")
  AC_MSG_NOTICE([QMP compile flags: ${QMP_CFLAGS}])
  AC_SUBST(QMP_LDFLAGS, "`${QMP_CONFIG} --ldflags`")
  AC_MSG_NOTICE([QMP linking flags: ${QMP_LDFLAGS}])
  AC_SUBST(QMP_LIBS,    "`${QMP_CONFIG} --libs`")
  AC_MSG_NOTICE([QMP libraries flags: ${QMP_LIBS}])
else 
  AC_MSG_WARN([QMP configuration program qmp-config not found.])
  AC_MSG_WARN([Set environment variables QMP_CFLAGS QMP_LDFAGS QMP_LIBS
    before configure])
fi

# Set up large file support
dnl
dnl --enable-largefile
dnl
AC_ARG_ENABLE(largefile,
  AC_HELP_STRING(
    [--enable-largefile],
    [ For files larger than 2 GB ]
  ),
  [AC_MSG_NOTICE([Enabling large file support.])
   AC_SYS_LARGEFILE
   qio_conf_opts="$qio_conf_opts --enable-largefile"]
)

dnl save the configure options
AC_SUBST(QIO_CONFIGURE_OPTIONS, [$qio_conf_opts])

# Checks for programs.
AC_PROG_CC
AC_PROG_RANLIB
AC_CHECK_TOOL(AR, ar, [ar])

dnl
dnl verify we can compile and link against QMP, if needed
dnl
PAC_QMP_LINK_CC_FUNC(
  [${QMP_CFLAGS}],
  [${QMP_LDFLAGS}],
  [${QMP_LIBS}],
  ,
  ,
  [qmp_link_ok=yes],
  [qmp_link_ok=no]
)

AC_MSG_CHECKING([if we can compile/link of a simple QMP C program])
case ${LOCAL_PARALLEL_ARCH} in
  parscalar)
    if test "X${qmp_link_ok}X" = "XyesX" ; then
      AC_MSG_RESULT(yes)
    else
      AC_MSG_RESULT(no)
      AC_MSG_ERROR([Cannot compile/link a basic QMP C program!
        Check QMP_CFLAGS, QMP_LDFLAGS, QMP_LIBS.])
    fi
    ;;
  *)
    AC_MSG_RESULT(ignored)
    ;;
esac

AC_SUBST(QDP_CFLAGS)
AC_SUBST(QDP_LDFLAGS)
AC_SUBST(QLA_CFLAGS)
AC_SUBST(QLA_LDFLAGS)

# Check for the fseeko functions
AC_FUNC_FSEEKO

# Checks for header files.
## AC_HEADER_STDC
## AC_CHECK_HEADERS([stdlib.h string.h strings.h])
AC_CHECK_HEADERS([malloc.h])

# Checks for type configuration
AC_CHECK_HEADERS([stdint.h inttypes.h],
[ dnl for inttypes.h and stdint.h for uint_xxx types
  dnl if successful check for the actual types too
  AC_CHECK_TYPES([uint16_t, uint32_t, uint64_t],
                 [],
                 [AC_MSG_ERROR([stdint.h found but either uint16_t, uint32_t or uint64_t not found]) ]
                )
],
[ 
  dnl no inttypes.h or stdint.h found check common unsigned types
  dnl for sizes and make appropriate decisions in the lime_fixed_types.h file
  AC_CHECK_SIZEOF(unsigned char)
  AC_CHECK_SIZEOF(unsigned short)
  AC_CHECK_SIZEOF(unsigned int)
  AC_CHECK_SIZEOF(unsigned long)
  AC_CHECK_SIZEOF(unsigned long long)
]
)

# Check for QMP route on parscalar
if test "X${LOCAL_PARALLEL_ARCH}X" = "XparscalarX" ; then 
   QMP_BKUP_CFLAGS="${CFLAGS}"
   QMP_BKUP_LDFLAGS="${LDFLAGS}"
   QMP_BKUP_LIBS="${LIBS}"
   AC_LANG_SAVE
   CFLAGS="${CFLAGS} ${QMP_CFLAGS}"
   LDFLAGS="${LDFLAGS} ${QMP_LDFLAGS}"
   LIBS="${LIBS} ${QMP_LIBS}"
   AC_CHECK_FUNCS(QMP_route)
   CFLAGS="${QMP_BKUP_CFLAGS}"
   LDFLAGS="${QMP_BKUP_LDFLAGS}"
   LIBS="${LIBS}"
   AC_LANG_RESTORE
fi

# Turn on parallel-io flag
if test ${ac_parallel_io} -eq 1; then
  AC_DEFINE_UNQUOTED(QIO_USE_PARALLEL_READ, ${ac_parallel_io}, [Enable parallel file-system read])
  AC_DEFINE_UNQUOTED(QIO_USE_PARALLEL_WRITE, ${ac_parallel_io}, [Enable parallel file-system write])
  AC_MSG_NOTICE([Parallel file-system extensions enabled]);
fi

if test "X${ac_qmp_route}X" == "XyesX"; then
  AC_MSG_NOTICE([Enabling use of QMP route if it exists])
  AC_DEFINE_UNQUOTED(QIO_USE_QMP_ROUTE, [1], [Enable QMP_route if it exists])
else
  AC_MSG_NOTICE([Disabling use of QMP route - even if it exists])
fi

# Now have all the options... Do some configuration 
AM_CONDITIONAL(ARCH_SCALAR,    [test "X${LOCAL_PARALLEL_ARCH}X" = "XscalarX"])
AM_CONDITIONAL(ARCH_PARSCALAR, [test "X${LOCAL_PARALLEL_ARCH}X" = "XparscalarX"])
AM_CONDITIONAL(USING_QMP, [test "X${QMP_CONFIG}X" != "XX"])

# pass configuration to distcheck target
DISTCHECK_CONFIGURE_FLAGS="--enable-parallel-arch=$PARALLEL_ARCH CC=\"$CC\" CFLAGS=\"$CFLAGS\""
if test -n "$QMP_HOME"; then
  DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-qmp=$QMP_HOME"
fi
AC_SUBST(DISTCHECK_CONFIGURE_FLAGS)

# export CC and CFLAGS so that the other_libs will use the same
# export CC CFLAGS
# Configure subdir packages
AC_CONFIG_SUBDIRS(other_libs/c-lime)

# Configure makefiles
AC_CONFIG_FILES(Makefile)
AC_CONFIG_FILES(other_libs/Makefile)
AC_CONFIG_FILES(include/Makefile)
AC_CONFIG_FILES(lib/Makefile)
AC_CONFIG_FILES(bin/Makefile)
AC_CONFIG_FILES(bin/qio-config)
AC_CONFIG_FILES(examples/Makefile)
AC_CONFIG_FILES(doc/Makefile)
AC_OUTPUT
