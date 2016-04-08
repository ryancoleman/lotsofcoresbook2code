#!/bin/bash
#
#################
# BUILD Chroma
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

####
#  DISABLE C++ Dslash because of include file conflicts
###
${SRCDIR}/cpp_wilson_dslash/configure --prefix=${INSTALL} \
	--with-qdp=${DEPDIR}/qdp++ \
        ${OMPENABLE} \
        CC="${PK_CC}"  CXX="${PK_CXX}" \
	CXXFLAGS="" CFLAGS="" \
        --host=x86_64-linux-gnu --build=none

${MAKE}
${MAKE} check
