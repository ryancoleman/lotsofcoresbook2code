AC_INIT(xpath_reader, 1.4.0, bj@ph.ed.ac.uk)
AC_CONFIG_AUX_DIR(config)
AC_PROG_CXX(g++)
AC_PROG_RANLIB()
AC_CHECK_TOOL(AR,ar,[ar])
AM_INIT_AUTOMAKE

AC_CONFIG_HEADERS([include/xpath_reader_config_internal.h])
dnl now libxml2
AC_ARG_WITH(libxml2,
  AC_HELP_STRING(
    [--with-libxml2=DIR],
    [Build on top of libxml2 where libxml2 is installed in DIR]
  ),
  [LIBXML2_HOME="$with_libxml2"]
)

echo ${LIBXML2_HOME}

if test "X${LIBXML2_HOME}X" = "XX" ; then
  AC_PATH_PROG(LIBXML2_CONFIG, [xml2-config], [])
else
  AC_PATH_PROG(LIBXML2_CONFIG, [xml2-config], [], [${LIBXML2_HOME}/bin:${PATH}])
fi

if test "X${LIBXML2_CONFIG}X" = "XX" ; then
  AC_MSG_ERROR([libxml2 configuration program xml2-config not found.])
fi

AC_MSG_NOTICE([Found libxml2 configuration program ${LIBXML2XX_CONFIG}])
AC_SUBST(LIBXML2_CXXFLAGS, "`${LIBXML2_CONFIG} --cflags`")
AC_MSG_NOTICE([libxml2 compile flags: ${LIBXML2_CXXFLAGS}])
AC_SUBST(LIBXML2_LIBS,     "`${LIBXML2_CONFIG} --libs`")
AC_MSG_NOTICE([libxml2 libraries flags: ${LIBXML2_LIBS}])

dnl This function is defined in acinclude.m4
dnl it tries to try to compile the program 
AC_MSG_CHECKING([if we can compile/link a simple libxml2 program])
PAC_LIBXML2_LINK_CXX_FUNC(
  ${LIBXML2_CXXFLAGS},
  ${LIBXML2_LIBS},
  ,
  ,
  [libxml2_link_ok=yes],
  [libxml2_link_ok=no]
)

if test "X${libxml2_link_ok}X" = "XyesX" ; then 
  AC_MSG_RESULT(yes)
else 
  AC_MSG_RESULT(no)
  AC_MSG_ERROR([Cannot compile/link a program with libxml2.
    Use --with-libxml2=<dir> to select a working version.])
fi 

AC_CONFIG_FILES(Makefile)
AC_CONFIG_FILES(lib/Makefile)
AC_CONFIG_FILES(include/Makefile)
AC_CONFIG_FILES(examples/Makefile)
AC_OUTPUT
