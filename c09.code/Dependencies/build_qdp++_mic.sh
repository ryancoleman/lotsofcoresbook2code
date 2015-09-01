#!/bin/bash
#
#################
# BUILD QMP
#################
source ../env/env_mic_comp_pref.sh

SRCDIR=${TOPDIR}/Dependencies/src
BUILDDIR=${TOPDIR}/Dependencies/build
INSTALLDIR=${TOPDIR}/Dependencies/install

pushd ${SRCDIR}/qdpxx
autoreconf
popd

pushd ${BUILDDIR}

if [ -d ./build_qdp++_mic ]; 
then 
  rm -rf ./build_qdp++_mic
fi

mkdir  ./build_qdp++_mic
cd ./build_qdp++_mic


${SRCDIR}/qdpxx/configure \
	--prefix=${INSTALLDIR}/qdp++_mic \
	--enable-precision=single \
        --enable-parallel-arch=scalar \
        --disable-filedb \
        --enable-largefile \
        --enable-parallel-io \
        --enable-alignment=64 \
	--with-libxml2=${INSTALLDIR}/libxml2_mic \
	CXXFLAGS="${PK_CXXFLAGS}" \
	CFLAGS="${PK_CFLAGS}" \
	CXX="${PK_CXX}" \
	CC="${PK_CC}" \
	--host=x86_64-linux-gnu --build=none \
	${OMPENABLE}

${MAKE}
${MAKE} install

popd
