#!/bin/bash
#
#################
# BUILD QMP
#################
source ../env/env_avx.sh

SRCDIR=${TOPDIR}/Dependencies/src
BUILDDIR=${TOPDIR}/Dependencies/build
INSTALLDIR=${TOPDIR}/Dependencies/install

pushd ${SRCDIR}/qdpxx
autoreconf
popd

pushd ${BUILDDIR}

if [ -d ./build_qdp++ ]; 
then 
  rm -rf ./build_qdp++
fi

mkdir  ./build_qdp++
cd ./build_qdp++


${SRCDIR}/qdpxx/configure \
	--prefix=${INSTALLDIR}/qdp++ \
	--enable-precision=single \
        --enable-parallel-arch=scalar \
        --disable-filedb \
        --enable-largefile \
        --enable-parallel-io \
        --enable-alignment=64 \
	--with-libxml2=${INSTALLDIR}/libxml2 \
	CXXFLAGS="${PK_CXXFLAGS}" \
	CFLAGS="${PK_CFLAGS}" \
	CXX="${PK_CXX}" \
	CC="${PK_CC}" \
	--host=x86_64-linux-gnu --build=none \
	${OMPENABLE}

${MAKE}
${MAKE} install

popd
