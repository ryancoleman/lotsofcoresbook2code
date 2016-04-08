dnl
m4_pattern_allow(QDP_AC_ALIGNMENT_SIZE)dnl
dnl
dnl Balint Joo, 13/12/2002
dnl George T. Fleming, 03/03/2003
dnl
dnl Stole this from mpich-1.2.4/mpe
dnl
dnl PAC_QMP_LINK_CXX_FUNC(
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
dnl               if QMP_FUNC is empty, assume linking with basic MPI program.
dnl               i.e. check if QMP definitions are valid
dnl
AC_DEFUN(
  PAC_QMP_LINK_CXX_FUNC,
  [
dnl - set local parallel compiler environments
dnl   so input variables can be CFLAGS, LDFLAGS or LIBS
    pac_QMP_CFLAGS="$1"
    pac_QMP_LDFLAGS="$2"
    pac_QMP_LIBS="$3"
	AC_LANG_PUSH([C++])
dnl - save the original environment
    pac_saved_CXXFLAGS="$CXXFLAGS"
    pac_saved_LDFLAGS="$LDFLAGS"
    pac_saved_LIBS="$LIBS"
dnl - set the parallel compiler environment
    CXXFLAGS="$CXXFLAGS $pac_QMP_CFLAGS"
    LDFLAGS="$LDFLAGS $pac_QMP_LDFLAGS"
    LIBS="$LIBS $pac_QMP_LIBS"
    AC_TRY_LINK(
      [#include "qmp.h"],
      [
        int argc ; char **argv ;
        QMP_thread_level_t prv;
        $4 ;
        QMP_init_msg_passing(&argc, &argv, QMP_THREAD_SINGLE, &prv) ;
        $5 ;
        QMP_finalize_msg_passing() ;
      ],
      [pac_qmp_working=yes],
      [pac_qmp_working=no]
    )
    CXXFLAGS="$pac_saved_CXXFLAGS"
    LDFLAGS="$pac_saved_LDFLAGS"
    LIBS="$pac_saved_LIBS"
    AC_LANG_POP([C++])
    if test "X${pac_qmp_working}X" = "XyesX" ; then
       ifelse([$6],,:,[$6])
    else
       ifelse([$7],,:,[$7])
    fi
  ]
)


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
dnl  LIBXML2_FUNC     for the body of a QDP++ function call or even general code
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
    AC_LANG_PUSH([C++])
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
    AC_LANG_POP([C++])
    if test "X${pac_libxml2_working}X" = "XyesX" ; then
       ifelse([$5],,:,[$5])
    else
       ifelse([$6],,:,[$6])
    fi
  ]
)

dnl Balint Joo, 13/12/2002
dnl George T. Fleming, 03/03/2003
dnl
dnl Stole this from mpich-1.2.4/mpe
dnl
dnl PAC_BAGEL_QDP_LINK_CXX_FUNC(
dnl   BAGEL_QDP_CFLAGS,
dnl   BAGEL_QDP_LDFLAGS,
dnl   BAGEL_QDP_LIBS,
dnl   BAGEL_QDP_VARS,
dnl   BAGEL_QDP_FUNC,
dnl   [action if working],
dnl   [action if not working]
dnl )
dnl
dnl  BAGEL_QDP_CFLAGS   is the include option (-I) for BAGEL_QDP includes
dnl  BAGEL_QDP_LDFLAGS  is the link path (-L) option for BAGEL_QDP libraries
dnl  BAGEL_QDP_LIBS     is the library (-l) option for BAGEL_QDP libaries
dnl  BAGEL_QDP_VARS     is the the declaration of variables needed to call BAGEL_QDP_FUNC
dnl  BAGEL_QDP_FUNC     is the body of BAGEL_QDP function call to be checked for existence
dnl               e.g.  BAGEL_QDP_VARS="BAGEL_QDP_u32_t foo;"
dnl                     BAGEL_QDP_FUNC="foo = BAGEL_QDP_get_SMP_count();"
dnl               if BAGEL_QDP_FUNC is empty, assume linking with basic MPI program.
dnl               i.e. check if BAGEL_QDP definitions are valid
dnl
AC_DEFUN(
  PAC_BAGEL_QDP_LINK_CXX_FUNC,
  [
dnl - set local parallel compiler environments
dnl   so input variables can be CFLAGS, LDFLAGS or LIBS
    pac_BAGEL_QDP_CFLAGS="$1"
    pac_BAGEL_QDP_LDFLAGS="$2"
    pac_BAGEL_QDP_LIBS="$3"
    AC_LANG_PUSH([C++])
dnl - save the original environment
    pac_saved_CXXFLAGS="$CXXFLAGS"
    pac_saved_LDFLAGS="$LDFLAGS"
    pac_saved_LIBS="$LIBS"
dnl - set the parallel compiler environment
    CXXFLAGS="$CXXFLAGS $pac_BAGEL_QDP_CFLAGS"
    LDFLAGS="$LDFLAGS $pac_BAGEL_QDP_LDFLAGS"
    LIBS="$LIBS $pac_BAGEL_QDP_LIBS"
    AC_TRY_LINK(
      [
        #include "bagel_qdp.h"
      ],
      [
        int argc ; char **argv ;
	BAGELQDPFloat *xptr;
	BAGELQDPFloat *aptr;
	BAGELQDPFloat *zptr;
	BAGELQDPFloat *yptr;
	unsigned int n3vec;
	qdp_vaxpy3(zptr,aptr,xptr,yptr,n3vec);
      ],
      [pac_bagel_qdp_working=yes],
      [pac_bagel_qdp_working=no]
    )
    CXXFLAGS="$pac_saved_CXXFLAGS"
    LDFLAGS="$pac_saved_LDFLAGS"
    LIBS="$pac_saved_LIBS"
    AC_LANG_POP([C++])
    if test "X${pac_bagel_qdp_working}X" = "XyesX" ; then
       ifelse([$6],,:,[$6])
    else
       ifelse([$7],,:,[$7])
    fi
  ]
)

dnl Abhinav Sarje, 10/09/2013
dnl Andre Walker-Loud, 10/09/2013
dnl
dnl PAC_HDF5_LINK_C_FUNC(
dnl   HDF5_CFLAGS,
dnl   HDF5_LIBS,
dnl   HDF5_VARS,
dnl   HDF5_FUNC,
dnl   [action if working],
dnl   [action if not working]
dnl )
dnl
dnl  HDF5_CFLAGS for the necessary includes paths (-I)
dnl  HDF5_LIBS     for the libraries (-l<lib> etc)
dnl  HDF5_VARS     for the declaration of variables needed
dnl                   to call HDF5_FUNC code fragment
dnl  HDF5_FUNC     for the general code
dnl                 fragment on which to run a compile/link test.
dnl                 If HDF5_VARS and HDF5_FUNC are empty, a basic test
dnl                 of compiling and linking a HDF5 program is run.
dnl
AC_DEFUN(
  PAC_HDF5_LINK_C_FUNC,
  [
dnl - set local parallel compiler environments
dnl - so input variables can be CFLAGS, LDFLAGS or LIBS
    pac_HDF5_CFLAGS="$1"
    pac_HDF5_CXXFLAGS="$2"
    pac_HDF5_LIBS="$3"
    AC_LANG_PUSH([C])
dnl - save the original environment
    pac_saved_CFLAGS="$CFLAGS"
    pac_saved_CXXFLAGS="$CXXFLAGS"
    pac_saved_LDFLAGS="$LDFLAGS"
    pac_saved_LIBS="$LIBS"
dnl - set the parallel compiler environment
    CFLAGS="$CFLAGS $pac_HDF5_CFLAGS"
    CXXFLAGS="$CXXFLAGS $pac_HDF5_CXXFLAGS"
    LDFLAGS="$LDFLAGS $pac_HDF5_LDFLAGS"
    LIBS="$LIBS $pac_HDF5_LIBS"
    AC_TRY_LINK(
      [
        #include <hdf5.h>
      ], [
        int argc ; char **argv ;
                hid_t prod_id, file_id;
                prod_id = H5Pcreate(H5P_FILE_ACCESS);
                file_id = H5Fcreate("test.hd5", H5F_ACC_TRUNC, H5P_DEFAULT, prod_id);
                H5Pclose(prod_id);
                H5Fclose(file_id);
      ],
      [pac_hdf5_working=yes],
      [pac_hdf5_working=no]
    )
    CFLAGS="$pac_saved_CFLAGS"
    CXXFLAGS="$pac_saved_CXXFLAGS"
    LDFLAGS="$pac_saved_LDFLAGS"
    LIBS="$pac_saved_LIBS"
    AC_LANG_POP([C])
    if test "X${pac_hdf5_working}X" = "XyesX" ; then
       ifelse([$6],,:,[$6])
    else
       ifelse([$7],,:,[$7])
    fi
  ]
)
