#!/bin/bash
#
#################
# BUILD QMP
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


${SRCDIR}/qphix/configure \
	--prefix=${INSTALLDIR}/dslash-mic-s8 \
	--with-qdp=${DEPDIR}/qdp++_mic \
	--enable-proc=scalar \
	--enable-soalen=1 \
	--enable-clover \
	--enable-cean \
	--enable-mm-malloc \
	CXXFLAGS="${PK_CXXFLAGS}" \
	CFLAGS="${PK_CFLAGS}" \
	CXX="${PK_CXX}" \
	CC="${PK_CC}" \
	--host=x86_64-linux-gnu --build=none \
	${OMPENABLE}

# If you want to use the VTune ITT API uncomment this line and add it on as 
# A continuation to the line above.
#
# LDFLAGS="-L/opt/intel/vtune_amplifier_xe/bin64/k1om" LIBS="-littnotify"

${MAKE} clean
${MAKE}

