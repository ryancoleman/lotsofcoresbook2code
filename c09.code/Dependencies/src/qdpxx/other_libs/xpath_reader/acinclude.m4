dnl PAC_LIBXML2_LINK_CXX_FUNC(
dnl   LIBXML2_CXXFLAGS,
dnl   LIBXML2_LIBS,
dnl   LIBXML2_VARS,
dnl   LIBXML2_FUNC,
dnl   [action if working],
dnl   [action if not working]
dnl )
dnl
dnl  LIBXML2_CXXFLAGS for the necessary includes paths (-I)
dnl  LIBXML2_LIBS     for the libraries (-l<lib> etc)
dnl  LIBXML2_VARS     for the declaration of variables needed
dnl                   to call LIBXML2_FUNC code fragment
dnl  LIBXML2_FUNC     for the body of a C++ function call or even general code
dnl                 fragment on which to run a compile/link test.
dnl                 If LIBXML2_VARS and LIBXML2_FUNC are empty, a basic test
dnl                 of compiling and linking a LIBXML2 program is run.
dnl
AC_DEFUN(
  PAC_LIBXML2_LINK_CXX_FUNC,
  [
dnl - set local parallel compiler environments
dnl - so input variables can be CXXFLAGS, LDFLAGS or LIBS
    pac_LIBXML2_CXXFLAGS="$1"
    pac_LIBXML2_LIBS="$2"
    AC_LANG_SAVE
    AC_LANG_CPLUSPLUS
dnl - save the original environment
    pac_saved_CXXFLAGS="$CXXFLAGS"
    pac_saved_LDFLAGS="$LDFLAGS"
    pac_saved_LIBS="$LIBS"
dnl - set the parallel compiler environment
    CXXFLAGS="$CXXFLAGS $pac_LIBXML2_CXXFLAGS"
    LDFLAGS="$LDFLAGS $pac_LIBXML2_LDFLAGS"
    LIBS="$LIBS $pac_LIBXML2_LIBS"
    AC_TRY_LINK(
      [
        #include <libxml/xmlmemory.h>
	#include <libxml/parser.h>
      ], [
        int argc ; char **argv ;
        xmlDocPtr doc;
	char *docname="foo";	
	doc = xmlParseFile(docname);
        $3 ;
        $4 ;
      ],
      [pac_libxml2_working=yes],
      [pac_libxml2_working=no]
    )
    CXXFLAGS="$pac_saved_CXXFLAGS"
    LDFLAGS="$pac_saved_LDFLAGS"
    LIBS="$pac_saved_LIBS"
    AC_LANG_RESTORE
    if test "X${pac_libxml2_working}X" = "XyesX" ; then
       ifelse([$5],,:,[$5])
    else
       ifelse([$6],,:,[$6])
    fi
  ]
)
