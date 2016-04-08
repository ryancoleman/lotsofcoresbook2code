#!/bin/bash
#
#################
# BUILD QMP
#################
source ../env/env_avx.sh
DEPDIR=${TOPDIR}/Dependencies/install
INSTALL=`pwd`/install
SRCDIR=`pwd`

if [ -d ./build_avx ]; 
then 
  rm -rf ./build_avx
fi

mkdir  ./build_avx
cd ./build_avx


${SRCDIR}/qphix/configure \
	--prefix=${INSTALLDIR}/dslash-avx-s8 \
	--with-qdp=${DEPDIR}/qdp++ \
	--enable-proc=AVX \
	--enable-soalen=8 \
	--enable-clover \
	--enable-cean \
	--enable-mm-malloc \
	CXXFLAGS="${PK_CXXFLAGS}" \
	CFLAGS="${PK_CFLAGS}" \
	CXX="${PK_CXX}" \
	CC="${PK_CC}" \
	--host=x86_64-linux-gnu --build=none \
	${OMPENABLE}

# LDFLAGS="-L/opt/intel/vtune_amplifier_xe/lib64" LIBS="-littnotify" \

${MAKE} clean
${MAKE}

