#!/bin/bash
#
#################
# BUILD Chroma
#################
source ../env/env_mic.sh

DEPDIR=${TOPDIR}/Dependencies/install
INSTALL=`pwd`/install
SRCDIR=`pwd`

if [ -d ./build_mic ]; 
then 
  rm -rf ./build_mic
fi

mkdir  ./build_mic
cd ./build_mic

####
#  DISABLE C++ Dslash because of include file conflicts
###
${SRCDIR}/cpp_wilson_dslash/configure --prefix=${INSTALL} \
	--with-qdp=${DEPDIR}/qdp++_mic \
        ${OMPENABLE} \
        CC="${PK_CC}"  CXX="${PK_CXX}" \
	CXXFLAGS="${PK_CXXFLAGS}" CFLAGS="${PK_CFLAGS}" \
        --host=x86_64-linux-gnu --build=none

# If you are compiling with the VTune ITT API add this line on to the end of the 
# last one with a continuation
#  LDFLAGS="-L/opt/intel/vtune_amplifier_xe/bin64/k1om" LIBS="-littnotify"

${MAKE}
${MAKE} check
