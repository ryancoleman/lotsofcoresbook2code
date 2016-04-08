AC_INIT(cpp_wilson_dslash, 1.0.0, bjoo@jlab.org)
AC_CONFIG_AUX_DIR(config)


AC_CANONICAL_HOST

AC_PROG_CC
AC_PROG_CXX
AC_PROG_RANLIB

dnl initialise automake
AM_INIT_AUTOMAKE

dnl put defines in include sse/config.h
AC_CONFIG_HEADERS([include/qphix/qphix_config_internal.h])
dnl dnl **************************************************************
dnl dnl ***                 Get User Options                       ***
dnl dnl **************************************************************

dnl learn the parallel architecture
AC_ARG_ENABLE(parallel-arch,
  AC_HELP_STRING([--enable-parallel-arch=<arch>],
    [Build Assembler for parallel architecture <arch>]),
  [ PARALLEL_ARCH="${enableval}"]
)

AC_ARG_ENABLE(testing,
   AC_HELP_STRING([--enable-testing],
   [Enable make check to automatically run tests, rather than just build them]),   [ do_tests="${enableval}" ]
)

AC_ARG_ENABLE(soalen,
   AC_HELP_STRING([--enable-soalen=N],
     [Build using SOALEN=N where N=4,8,16 for MIC N=4,8 for AVX]),
     [SOALEN=${enableval}],
     [SOALEN=8]
)

AC_ARG_ENABLE(cean,
    AC_HELP_STRING([--enable-cean],
      [ Enable Cilk Array Notation extentions ]),
      [ USE_CEAN=${enableval}],
      [ USE_CEAN="no" ]
)

AC_ARG_ENABLE(mm_malloc,
	AC_HELP_STRING([--enable-mm-malloc],
	  [ Use _mm_malloc to allocate aligned memory. Otherwise use posix memalign. (default: enabled) ]),
	  [ USE_MM_MALLOC=${enableval} ],
	  [ USE_MM_MALLOC="yes" ]
)

AC_ARG_ENABLE(clover,
   AC_HELP_STRING([--enable-clover],
     [ Enable building of clover term. default is no ]),
     [ build_clover=${enableval}],
     [ build_clover="no" ]
)

dnl check for QMP -- needed for parallel configures
AC_ARG_WITH(qmp,
  AC_HELP_STRING(
    [--with-qmp=DIR],
    [Build Assembler on top of QMP, where QMP is installed in DIR.]
  ),
  [QMP_HOME="${with_qmp}"]
)

AC_ARG_ENABLE(fake-comms,
	AC_HELP_STRING(
	[--enable-fake-comms],
	[Allocate comms buffers and point send buffers to recv. Face processing but no actual comms]
	),
	[ USE_FAKE_COMMS=${enableval} ],
	[ USE_FAKE_COMMS="no" ]
)


dnl Use Jie Chens threads library
AC_ARG_WITH(qmt,
   AC_HELP_STRING(
    [--with-qmt=<DIR>],
    [Enable use of QMT Threads library installed in DIR]
   ),
   [ qmt_enabled="yes"
     QMT_HOME=${with_qmt} ],
   [ qmt_enabled="no" ]
)

dnl Use OpenMP threading
AC_ARG_ENABLE(openmp,
   AC_HELP_STRING(
    [--enable-openmp],
    [Enable building of OpenMP dispatcher. Pass flags necessary for OpenMP on CFLAGS]
   ),
   [ omp_enabled="yes" ],
   [ omp_enabled="no" ]
)


dnl with QDP -- QDP will override everything else
dnl since a parallel QDP will know its architecture
dnl precition and have all the flags for QMP set already
dnl in its CXXFLAGS etc etc
dnl 
dnl This flag and all the other flags are complimentary
dnl to each other therefore. Also --with-qdp is also 
dnl passed down the command line when compiling Chroma,
dnl so a sub-configure will also work
dnl we can also put the QDP++ packers in here too...
AC_ARG_WITH(qdp,
  AC_HELP_STRING(
     [--with-qdp=DIR],
     [ Build Assembler with QDP++, where QDP++ is installed in DIR. ]
  ),
  [ QDPXX_HOME="${with_qdp}"
    QDPXX_GIVEN="yes"
  ],
  [ QDPXX_GIVEN="no" ]
)

AC_ARG_ENABLE(proc,
   AC_HELP_STRING(
	[--enable-proc=proc],
	[ Build code specializations for proc. proc=NONE,SSE,AVX,MIC ]
   ),
   [ proc="${enable_proc}" ],
   [ proc="NONE"]
)

dnl dnl ******************** DONE WITH USER INPUT *********************
dnl dnl ***        Now Configuration                               ****
dnl dnl ***************************************************************
dnl 
dnl     ***************************************************************
dnl     ***       There are two distinct branches depending         ***
dnl     ***      on whether we are compiling with QDP++ or not      ***
dnl     ***************************************************************

dnl 


if test "X${QDPXX_GIVEN}X" = "XyesX"; then 
     dnl Find the QDP++ config program 
     if test "X${QDPXX_HOME}X" = "XX" ; 
     then
	AC_PATH_PROG(QDPXX_CONFIG, [qdp++-config], [])
     else
	AC_PATH_PROG(QDPXX_CONFIG, [qdp++-config], [], [${QDPXX_HOME}/bin:${PATH}])
     fi

     if test "X${QDPXX_CONFIG}X" = "XX" ; then
	AC_MSG_ERROR([QDP++ configuration program qdp++-config not found.])
     fi

     AC_MSG_NOTICE([Found QDP++ configuration program ${QDPXX_CONFIG}])
     AC_MSG_CHECKING([if QDP++ uses checkerboarded layout])
     layout=`${QDPXX_CONFIG} --layout`
     if test "X${layout}X" == "Xcb2X" -o "X${layout}X" == "Xcb3dX" ; then 
	AC_MSG_RESULT([yes])
     else
        AC_MSG_RESULT([no])
	AC_MSG_ERROR([This Dslash Operator only works for checkerboarded QDP++. Sorry])
     fi


     AC_SUBST(QPHIX_AUX_CFLAGS, "`${QDPXX_CONFIG} --cxxflags`")
     AC_MSG_NOTICE([QDP++ compile flags: ${QPHIX_AUX_CFLAGS}])

     AC_SUBST(QPHIX_AUX_LDFLAGS,  "`${QDPXX_CONFIG} --ldflags`")
     AC_MSG_NOTICE([QDP++ linking flags: ${QPHIX_AUX_LDFLAGS}])

     AC_SUBST(QPHIX_AUX_LIBS,     "`${QDPXX_CONFIG} --libs`")
     AC_MSG_NOTICE([QDP++ libraries flags: ${QPHIX_AUX_LIBS}])

     QDPXX_ND="`${QDPXX_CONFIG} --Nd`"
     AC_MSG_NOTICE([QDP++ spacetime dimension: ${QDPXX_ND}])
     
     QDPXX_NC="`${QDPXX_CONFIG} --Nc`"
     AC_MSG_NOTICE([QDP++ number of colors: ${QDPXX_NC}])

     QDPXX_NS="`${QDPXX_CONFIG} --Ns`"
     AC_MSG_NOTICE([QDP++ number of spin components: ${QDPXX_NS}])

     PARALLEL_ARCH="`${QDPXX_CONFIG} --parallel-arch`"

     AC_MSG_NOTICE([QDP++ has parallel arch: ${PARALLEL_ARCH} ] )

     dnl This function is defined in acinclude.m4
     dnl it tries to try to compile the program 
     PAC_QDPXX_LINK_CXX_FUNC(
	${QPHIX_AUX_CFLAGS},
	${QPHIX_AUX_LDFLAGS},
	${QPHIX_AUX_LIBS},
	,
	,
	[qdpxx_link_ok=yes],
	[qdpxx_link_ok=no]
     )

     AC_MSG_CHECKING([if we can compile/link a simple QDP++ program])
     if test "X${qdpxx_link_ok}X" = "XyesX" ; then 
	AC_MSG_RESULT(yes)
     else 
	AC_MSG_RESULT(no)
	AC_MSG_ERROR([Cannot compile/link a program with QDP++.
    Use --with-qdp++=<dir> to select a working version.])
    fi 



    dnl *********************************************************
    dnl ***              end of section with QDP++            ***
    dnl *********************************************************
else 
     dnl ***************************************************************
     dnl *** Section without QDP++ -- Need to check for QMP possibly ***
     dnl ***************************************************************
     case "${PARALLEL_ARCH}" in
	  scalar)
		AC_MSG_NOTICE([ Scalar Build! Not checking for QMP ])
		AC_SUBST(QPHIX_AUX_CFLAGS)
		AC_SUBST(QPHIX_AUX_LDFLAGS)
		AC_SUBST(QPHIX_AUX_LIBS)
		;;

	  parscalar)
		AC_MSG_NOTICE([ Parscalar build! Checking for QMP ])

		dnl look for qmp-config program
		if test "X${QMP_HOME}X" = "XX" ; then
		   AC_PATH_PROG([QMP_CONFIG], [qmp-config], [])
		else
		   AC_PATH_PROG([QMP_CONFIG], [qmp-config], [], [${QMP_HOME}/bin:${PATH}])
		fi

		dnl look for flags needed for QMP
		if test "X${QMP_CONFIG}X" != "XX" ; then
		   AC_MSG_NOTICE([Found QMP configuration program ${QMP_CONFIG}])
		   AC_SUBST(QPHIX_AUX_CFLAGS,  "`${QMP_CONFIG} --cflags`")
		   AC_MSG_NOTICE([QMP compile flags: ${QPHIX_AUX_CFLAGS}])

		   AC_SUBST(QPHIX_AUX_LDFLAGS, "`${QMP_CONFIG} --ldflags`")
		   AC_MSG_NOTICE([QMP linking flags: ${QPHIX_AUX_LDFLAGS}])

		   AC_SUBST(QPHIX_AUX_LIBS,    "`${QMP_CONFIG} --libs`")
		   AC_MSG_NOTICE([QMP libraries flags: ${QPHIX_AUX_LIBS}])
		else 
		   AC_MSG_WARN([QMP configuration program qmp-config not found.])
		   AC_MSG_WARN([Set environment variables QMP_CFLAGS QMP_LDFAGS QMP_LIBS before configure])
		fi

		dnl verify we can compile and link against QMP, if needed
		AC_PROG_CXX
		PAC_QMP_LINK_CC_FUNC(
			[${QPHIX_AUX_CFLAGS}],
			[${QPHIX_AUX_LDFLAGS}],
			[${QPHIX_AUX_LIBS}],
			,
			,
			[qmp_link_ok=yes],
			[qmp_link_ok=no]
		)
		
		AC_MSG_CHECKING([if we can compile/link of a simple QMP C program])
		if test "X${qmp_link_ok}X" = "XyesX" ; then
		   AC_MSG_RESULT(yes)
		else
		   AC_MSG_RESULT(no)
		   AC_MSG_ERROR([Cannot compile/link a basic QMP C program! Check QMP_CFLAGS, QMP_LDFLAGS, QMP_LIBS.])
		fi
		;;
	*)
		AC_MSG_ERROR([ Unknown value for --enable-parallel-arch ])
		;;
	esac
fi


dnl ************************************************************************
dnl **** Compile Defines to select right stuff. Always SSE now             *
dnl ************************************************************************


case ${PARALLEL_ARCH} in
scalar)
	AC_DEFINE_UNQUOTED(QPHIX_SCALAR, [1], [Scalar Arch])
	if test "X${USE_FAKE_COMMS}X" = "XyesX"; 
	then
	   AC_MSG_NOTICE([ Going to use Fake comms ])	
	   AC_DEFINE_UNQUOTED(QPHIX_DO_COMMS,[1], [Do the actual comms])
   	   AC_DEFINE_UNQUOTED(QPHIX_FAKE_COMMS, [1], [Fake Comms])
	else
	   AC_MSG_NOTICE([ This version does no comms])
        fi
	;;
parscalar)
	AC_DEFINE_UNQUOTED(QPHIX_PARSCALAR, [1], [Parscalar Arch])
	if test "X${USE_FAKE_COMMS}X" = "XyesX"; 
	then	
      	   AC_MSG_ERROR([Cannot use parscalar arch with fake comms])
	else
	   AC_MSG_NOTICE([ Regular QMP COMMS ])
    	   AC_DEFINE_UNQUOTED(QPHIX_DO_COMMS,[1], [Do the actual comms])
	   AC_DEFINE_UNQUOTED(QPHIX_QMP_COMMS,[1], [QMP Comms])
	fi					   
	;;
*)
	AC_MSG_ERROR([Unsupported parallel arch])
	;;
esac

case ${SOALEN} in
1)
	AC_MSG_NOTICE([Building for SOALEN=1])
	AC_DEFINE_UNQUOTED(QPHIX_SOALEN,[1], [SOALEN is 1])
	;;
4)
	AC_MSG_NOTICE([Building for SOALEN=4])
	AC_DEFINE_UNQUOTED(QPHIX_SOALEN,[4], [SOALEN is 4])
	;;
8)
	AC_MSG_NOTICE([Building for SOALEN=8])
	AC_DEFINE_UNQUOTED(QPHIX_SOALEN,[8], [SOALEN is 8])
	;;
16)
	AC_MSG_NOTICE([Building for SOALEN=16])
	AC_DEFINE_UNQUOTED(QPHIX_SOALEN,[16], [SOALEN is 16])
	;;
*)
	AC_MSG_ERROR([Unsupported SOALEN])
	;;
esac

case ${USE_CEAN} in
yes) 
	AC_MSG_NOTICE([Enabling Cilk Notation Use])
	AC_DEFINE_UNQUOTED(QPHIX_USE_CEAN,[1], [CEAN is on])
	;;
*)
	;;
esac

case ${USE_MM_MALLOC} in
yes) 
	AC_MSG_NOTICE([Enabling _mm_malloc use])
	AC_DEFINE_UNQUOTED(QPHIX_USE_MM_MALLOC,[1], [MM_MALLOC is on])
	;;
*)
	AC_MSG_NOTICE([Not using _mm_malloc. Using posix memalign])
	;;
esac

AM_CONDITIONAL(QPHIX_SCALAR, [ test "X${PARALLEL_ARCH}X" == "XscalarX" ] )
AM_CONDITIONAL(QPHIX_PARSCALAR, [ test "X${PARALLEL_ARCH}X" == "XparscalarX" ] )
AM_CONDITIONAL(QPHIX_BUILD_WITH_QDP, [ test "X${QDPXX_GIVEN}X" == "XyesX" ])

if test "x${build_clover}x" == "xyesx";
then
        AC_MSG_NOTICE([Enabling Clover Term])
        AC_DEFINE([QPHIX_BUILD_CLOVER],[1], [ Build Clover Term ])
else
        AC_MSG_NOTICE([Not Building Clover Term])
fi

build_qpx="no"

case ${proc} in
AVX|avx) 
	AC_DEFINE([QPHIX_AVX_SOURCE],[1], [ Generate AVX Specializations ])
	AC_DEFINE([QPHIX_LLC_CACHE_ALIGN], [64], [ L2 Cache alignment ] )

	AC_MSG_NOTICE([Configuring AVX Specializations ])
	if test $SOALEN -gt 8;
	then 
	  AC_MSG_ERROR([Cannot have SOALEN more than 8 for AVX]);
	fi		       	    
;;
MIC|mic)
	AC_DEFINE([QPHIX_MIC_SOURCE],[1], [ Generate MIC Specializations ])
	AC_DEFINE([QPHIX_LLC_CACHE_ALIGN], [64], [ L2 Cache alignment ] )

	AC_MSG_NOTICE([Configuring MIC Specializations ])
;;
SCALAR|scalar)
	AC_DEFINE([QPHIX_SCALAR_SOURCE],[1], [ Generate Scalar Specializations ])
	AC_DEFINE([QPHIX_LLC_CACHE_ALIGN], [64], [ L2 Cache alignment ] )
	AC_MSG_NOTICE([Configuring Scalar Specializations ])
;;        
QPX|qpx)
	AC_DEFINE([QPHIX_QPX_SOURCE],[1], [ Generate QPX Specializations ])
	AC_DEFINE([QPHIX_LLC_CACHE_ALIGN], [128], [ L2 Cache alignment ] )
	AC_MSG_NOTICE([Configuring QPX Specializations ])
	if test $SOALEN -ne 4;
	then
	  AC_MSG_ERROR([Cannot have SOALEN not equal to 4 for QPX]);
	fi
	build_qpx="yes"
;;
*)
	AC_MSG_NOTICE([Not generating any code specializations])
	;;
esac

AM_CONDITIONAL(QPHIX_BUILD_CLOVER, [ test "x${build_clover}x" = "xyesx" ])
AM_CONDITIONAL(QPHIX_BUILD_QPX, [ test "x${build_qpx}x" = "xyesx"])

dnl dnl produce output
AC_CONFIG_FILES(Makefile)
AC_CONFIG_FILES(include/Makefile)
AC_CONFIG_FILES(lib/Makefile)
AC_CONFIG_FILES(tests/Makefile)
AC_OUTPUT
