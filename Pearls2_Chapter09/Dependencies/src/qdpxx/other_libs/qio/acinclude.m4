dnl Balint Joo, 13/12/2002
dnl George T. Fleming, 03/03/2003
dnl
dnl Stole this from mpich-1.2.4/mpe
dnl
dnl PAC_QMP_LINK_CC_FUNC(
dnl   QMP_CFLAGS,
dnl   QMP_LDFLAGS,
dnl   QMP_LIBS,
dnl   QMP_VARS,
dnl   QMP_FUNC,
dnl   [action if working],
dnl   [action if not working]
dnl )
dnl
dnl  QMP_CFLAGS   is the include option (-I) for QMP includes
dnl  QMP_LDFLAGS  is the link path (-L) option for QMP libraries
dnl  QMP_LIBS     is the library (-l) option for QMP libaries
dnl  QMP_VARS     is the the declaration of variables needed to call QMP_FUNC
dnl  QMP_FUNC     is the body of QMP function call to be checked for existence
dnl               e.g.  QMP_VARS="QMP_u32_t foo;"
dnl                     QMP_FUNC="foo = QMP_get_SMP_count();"
dnl               if QMP_FUNC is empty, assume linking with basic QMP program.
dnl               i.e. check if QMP definitions are valid
dnl
AC_DEFUN(
  PAC_QMP_LINK_CC_FUNC,
  [
dnl - set local parallel compiler environments
dnl   so input variables can be CFLAGS, LDFLAGS or LIBS
    pac_QMP_CFLAGS="$1"
    pac_QMP_LDFLAGS="$2"
    pac_QMP_LIBS="$3"
    AC_LANG_SAVE
    AC_LANG_C
dnl - save the original environment
    pac_saved_CFLAGS="$CFLAGS"
    pac_saved_LDFLAGS="$LDFLAGS"
    pac_saved_LIBS="$LIBS"
dnl - set the parallel compiler environment
    CFLAGS="$CFLAGS $pac_QMP_CFLAGS"
    LDFLAGS="$LDFLAGS $pac_QMP_LDFLAGS"
    LIBS="$LIBS $pac_QMP_LIBS"
    AC_TRY_LINK(
      [#include "qmp.h"],
      [
        int argc ; char **argv ; QMP_thread_level_t provided;
        $4 ;
        QMP_init_msg_passing(&argc, &argv, QMP_THREAD_SINGLE, &provided ) ;
        $5 ;
        QMP_finalize_msg_passing() ;
      ],
      [pac_qmp_working=yes],
      [pac_qmp_working=no]
    )
    CFLAGS="$pac_saved_CFLAGS"
    LDFLAGS="$pac_saved_LDFLAGS"
    LIBS="$pac_saved_LIBS"
    AC_LANG_RESTORE
    if test "X${pac_qmp_working}X" = "XyesX" ; then
       ifelse([$6],,:,[$6])
    else
       ifelse([$7],,:,[$7])
    fi
  ]
)
