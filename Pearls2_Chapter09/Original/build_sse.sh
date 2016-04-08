#!/bin/bash
#
#################
# BUILD Chroma
#################
source ../env/env_avx.sh

DEPDIR=${TOPDIR}/Dependencies/install
INSTALL=`pwd`/install
SRCDIR=`pwd`

if [ -d ./build_avx_sse ]; 
then 
  rm -rf ./build_avx_sse
fi

mkdir  ./build_avx_sse
cd ./build_avx_sse

####
#  DISABLE C++ Dslash because of include file conflicts
###
${SRCDIR}/cpp_wilson_dslash/configure --prefix=${INSTALL} \
	--with-qdp=${DEPDIR}/qdp++ \
	--enable-sse2 --enable-sse3 \
        ${OMPENABLE} \
        CC="${PK_CC}"  CXX="${PK_CXX}" \
	CXXFLAGS="" CFLAGS="" \
        --host=x86_64-linux-gnu --build=none

${MAKE}
${MAKE} check
