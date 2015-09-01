# configure.ac file for FileDB package
AC_INIT([filedb],[1.2],[filedb authors <chen@jlab.org>],[filedb])
AC_COPYRIGHT([Jefferson Science Assoicates, LLC. Jefferson Lab])
AC_MSG_NOTICE([Configuring filedb])

AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE([-Wall -Werror])

AC_CONFIG_HEADERS([src/filedb_config_internal.h])

## Get compiling host information
BUILD_OS=linux
AC_CANONICAL_BUILD
# Check whether this is a linux machine
case $build in
     *linux*)
     AC_MSG_NOTICE([This is a linux machine.])
     ;;
     *)
     AC_MSG_WARN([This package has not been tried on your platform $build !])
     BUILD_OS=$build_os
     ;;
esac

AC_PROG_CC

AC_PROG_CXX

AC_LANG(C,C++)

AM_PROG_AS
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

AC_PROG_INSTALL

AC_PROG_RANLIB

AC_CHECK_TOOL(AR, ar, [ar])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_C_VOLATILE

if test "X{$LDFLAGS}X" != "XX"
then
	SAVELDFLAGS=$LDFLAGS
fi

if test "X{$CPPFLAGS}X" != "XX"
then
	SAVECPPFLAGS=$CPPFLAGS
fi

# Check big/litte endian
FILEDB_BYTE_ORDER=UNKNOWN
AC_C_BIGENDIAN([FILEDB_BYTE_ORDER=BIG_ENDIAN],[FILEDB_BYTE_ORDER=LITTLE_ENDIAN],[FILEDB_BYTE_ORDER=UNKNOWN],[FILEDB_BYTE_ORDER=UNIV])

if test "$FILEDB_BYTE_ORDER" = BIG_ENDIAN
then
	AC_MSG_NOTICE([Compile FileDB for big endian machines.])
	AC_SUBST(FILEDB_BYTE_ORDER_DEF,[-D_FILEDB_BIG_ENDIAN])
elif test "$FILEDB_BYTE_ORDER" = LITTLE_ENDIAN
then
	AC_MSG_NOTICE([Compile FileDB for little endian machines.])
	AC_SUBST(FILEDB_BYTE_ORDER_DEF,[-D_FILEDB_LITTLE_ENDIAN])
else
	AC_MSG_ERROR([Cannot determine endian type of this machine.])
fi

# Turn on massive debug information
FILEDB_DEBUG_ON=NO
AC_ARG_ENABLE([debug-hashdb],
AC_HELP_STRING([--enable-debug-hashdb],[turn on low level debug verbose messages.]),
[if test "$enableval" = yes; then FILEDB_DEBUG_ON=YES;fi])

if test "$FILEDB_DEBUG_ON" = YES
then
	AC_MSG_NOTICE([Low level hash database debug message is on.])
	AC_SUBST(FILEHASH_DEBUG,[-D_FFDB_DEBUG])
fi

FILEDB_TEMPLATE_COMPILE_TEST=NO
AC_ARG_ENABLE([template-compile-test],
AC_HELP_STRING([--enable-template-compile-test],[Test compile template headers through dumb instantiations.]),
[if test "$enableval" = yes; then FILEDB_TEMPLATE_COMPILE_TEST=YES;fi])

if test "$FILEDB_TEMPLATE_COMPILE_TEST" = YES
then
	AC_MSG_NOTICE([Testing compile template headers is on.])
	AC_SUBST(FILEDB_TEMP_COMP_TEST,[-D_FILEDB_COMPILE_TEST])
fi

AC_CONFIG_FILES(Makefile)
AC_CONFIG_FILES(filehash/Makefile)
AC_CONFIG_FILES(src/Makefile)

AC_OUTPUT
