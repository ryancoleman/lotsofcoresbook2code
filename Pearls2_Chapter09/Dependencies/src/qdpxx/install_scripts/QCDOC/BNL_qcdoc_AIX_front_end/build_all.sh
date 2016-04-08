#!/usr/bin/bash

# Get build functions
. ./build_functions.sh

HERE=/home/edwardsr/qcd
INSTALL_ROOT=/home/edwardsr/arch/v2.6.0-CJ
QOS=/qcdoc/sfw/qos/devel/v2.6.0-CJ/aix5.2f
PRECISION=double
HOST_SYS=powerpc-gnu-elf
BUILD_SYS=none
BAGEL_QDP=${HERE}/bagel_qdp
BAGEL_INSTALL_DIR=${INSTALL_ROOT}/bagel
BAGEL_CPU=ppc440
BAGEL_ALLOC=qalloc
BAGEL_COMM=qmp
BAGEL_WILSON_DIR=${INSTALL_ROOT}/bagel_${PRECISION}_${BAGEL_CPU}_${BAGEL_ALLOC}_${BAGEL_COMM}
BAGELQDP_INSTALLDIR=${INSTALL_ROOT}/bagel_qdp_${PRECISION}_${BAGEL_CPU}

LIBXML_SRCDIR=${HERE}/libxml2-2.6.6
LIBXML=${INSTALL_ROOT}/libxml2-2.6.6
QDP_PARALLEL_ARCH=parscalar
QDP_DO_EDRAM=yes
QDP_DO_BLAS=yes
CHROMA_DO_PAB_DSLASH=yes
CHROMA_DO_GMP=yes
CHROMA_GMPDIR=/qcdoc/sfw/packages/gmp/qos-2.5.9
# Munge directory names
QDP_INSTALLDIR=${INSTALL_ROOT}/qdp_${PRECISION}
if test "X${QDP_DO_EDRAM}X" == "XyesX";
then
	QDP_INSTALLDIR=${QDP_INSTALLDIR}_edram;
else
	QDP_INSTALLDIR=${QDP_INSTALLDIR}_ddr;
fi

if test "X${QDP_DO_BAGEL}X" == "XyesX";
then 
	QDP_INSTALLDIR=${QDP_INSTALLDIR}_blas;
fi

CHROMA_INSTALLDIR=${INSTALL_ROOT}/chroma_${PRECISION}
if test "X${QDP_DO_EDRAM}X" == "XyesX";
then
        CHROMA_INSTALLDIR=${CHROMA_INSTALLDIR}_edram;
else
        CHROMA_INSTALLDIR=${CHROMA_INSTALLDIR}_ddr;
fi
                                                                             
if test "X${CHROMA_DO_PAB_DSLASH}X" == "XyesX";
then
        CHROMA_INSTALLDIR=${CHROMA_INSTALLDIR}_pab_dslash;
fi

##
## The actual building is done here
##
##

## Build BAGEL
build_bagel ${HERE}/bagel ${BAGEL_INSTALL_DIR}

source ${QOS}/scripts/setup.sh
## Build Wilson Dslash

build_bagel_wilson_dslash ${HERE}/bagel_wilson_dslash \
                          ${BAGEL_WILSON_DIR} \
		          ${BAGEL_INSTALL_DIR} \
	                  ${PRECISION} ${BAGEL_COMM} ${BAGEL_ALLOC} \
	                  ${BAGEL_CPU} ${HOST_SYS} ${BUILD_SYS} ${QOS}

build_bagel_qdp ${HERE}/bagel_qdp \
                ${BAGEL_INSTALL_DIR} \
		${BAGELQDP_INSTALLDIR} \
                ${PRECISION} \
                ${BAGEL_CPU} \
                ${HOST_SYS} \
		${BUILD_SYS}

build_libxml ${LIBXML_SRCDIR} ${LIBXML} ${HOST_SYS} ${BUILD_SYS}

## Build QDP++
build_qdp  "${HERE}/qdp++" ${QDP_INSTALLDIR} ${QOS} ${LIBXML} ${PRECISION} ${QDP_DO_EDRAM} ${QDP_DO_BLAS} ${HOST_SYS} ${BUILD_SYS} ${BAGELQDP_INSTALLDIR}

## Build Chroma
build_chroma ${HERE}/chroma ${CHROMA_INSTALLDIR} ${QDP_INSTALLDIR} ${HOST_SYS} ${BUILD_SYS} ${CHROMA_DO_PAB_DSLASH} ${BAGEL_WILSON_DIR} ${CHROMA_DO_GMP} ${CHROMA_GMPDIR}

pushd ${INSTALL_ROOT}
find . -name "*" -type d -exec chmod ugo+rx {} \; -print
popd
