#!/bin/bash

# Copyright (c) 2015 Matthias Noack (ma.noack.pr@gmail.com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

CC=icpc
HAM_PATH=thirdparty/ham
VC_PATH=thirdparty/vc
VECCLASS_PATH=thirdparty/vcl
VECCLASS_PATH_MIC=thirdparty/vcl_mic

VECLIB="VEC_INTEL"
#VECLIB="VEC_VC"
#VECLIB="VEC_VCL"
NUM_ITERATIONS=26 # including warmup below
NUM_WARMUP=1

OPTIONS="-std=c++11 -g -O3 -restrict -openmp -qopt-report=5" #-Wall
#OPTIONS="-std=c++11 -g -O3 -restrict -openmp -qopt-report=5 -DUSE_INITZERO" #-Wall
OPTIONS_MIC="$OPTIONS -mmic"
#OPTIONS_MIC="$OPTIONS_MIC -opt-prefetch-distance=6,1"
OPTIONS_HOST="$OPTIONS -xHost"

INCLUDE="-Iinclude -I${HAM_PATH}/include"
INCLUDE_HOST="$INCLUDE -I${VC_PATH}/include -I${VECCLASS_PATH}"
INCLUDE_MIC="$INCLUDE -I${VC_PATH}/include -I${VECCLASS_PATH_MIC}"

LIB="-lrt"
LIB_HOST="$LIB ${VC_PATH}/lib/libVc.a"
LIB_MIC="$LIB ${VC_PATH}/lib/libVc_MIC.a"

BUILD_DIR_HOST="bin"
BUILD_DIR_MIC="bin.mic"

FILES=( \
common.cpp \
kernel/commutator_reference.cpp \
kernel/commutator_omp_aosoa.cpp \
kernel/commutator_omp_aosoa_constants.cpp \
kernel/commutator_omp_aosoa_direct.cpp \
kernel/commutator_omp_aosoa_constants_direct.cpp \
kernel/commutator_omp_aosoa_constants_direct_perm.cpp \
kernel/commutator_omp_aosoa_constants_direct_perm2to3.cpp \
kernel/commutator_omp_aosoa_constants_direct_perm2to5.cpp \
kernel/commutator_omp_manual_aosoa.cpp \
kernel/commutator_omp_manual_aosoa_constants.cpp \
kernel/commutator_omp_manual_aosoa_constants_perm.cpp \
kernel/commutator_omp_manual_aosoa_direct.cpp \
kernel/commutator_omp_manual_aosoa_constants_direct.cpp \
kernel/commutator_omp_manual_aosoa_constants_direct_perm.cpp \
kernel/commutator_omp_manual_aosoa_constants_direct_unrollhints.cpp \
kernel/commutator_omp_manual_aosoa_constants_direct_perm_unrollhints.cpp \
)

# compile
build()
{
	local BUILD_DIR=$1
	local OPTIONS=$2
	local INCLUDE=$3
	local LIB=$4

	mkdir -p $BUILD_DIR
	
		
	for file in "${FILES[@]}"
	do
		local tmp=${file##*/}
		local NAME=${tmp%.*}
		$CC -c $OPTIONS $INCLUDE -o ${BUILD_DIR}/${NAME}.o src/${file} 
		OBJS=" $OBJS ${BUILD_DIR}/${NAME}.o"
	done

	echo $OBJS

#	$CC -c $OPTIONS $INCLUDE -o ${BUILD_DIR}/common.o src/common.cpp 
#	$CC -c $OPTIONS $INCLUDE -o ${BUILD_DIR}/commutator_reference.o src/kernel/commutator_reference.cpp
#	$CC -c $OPTIONS $INCLUDE -o ${BUILD_DIR}/commutator_omp_auto_final.o src/kernel/commutator_omp_auto_final.cpp
#	$CC -c $OPTIONS $INCLUDE -o ${BUILD_DIR}/commutator_omp_manual_aosoa_constants_direct.o src/kernel/commutator_omp_manual_aosoa_constants_direct.cpp
#	$CC -c $OPTIONS $INCLUDE -o ${BUILD_DIR}/commutator_omp_manual_final.o src/kernel/commutator_omp_manual_final.cpp
#	$CC $OPTIONS $INCLUDE -o ${BUILD_DIR}/benchmark_omp ${BUILD_DIR}/commutator_omp_manual_final.o ${BUILD_DIR}/commutator_omp_manual_aosoa_constants_direct.o  ${BUILD_DIR}/commutator_omp_auto_final.o ${BUILD_DIR}/commutator_reference.o ${BUILD_DIR}/common.o src/benchmark_omp.cpp $LIB

	$CC $OPTIONS $INCLUDE -o ${BUILD_DIR}/benchmark_omp $OBJS src/benchmark_omp.cpp $LIB
}

usage ()
{ 
	echo "Usage and defaults:";
	echo -e "\t-c\t Build CPU variant.";
	echo -e "\t-a\t Build MIC (Accelerator) variant."; 
	echo -e "\t-i ${NUM_ITERATIONS}\t Number of iterations (including warmups).";
	echo -e "\t-w ${NUM_WARMUP}\t Number of warmup iterations.";
	echo -e "\t-v ${VECLIB}\t Vector library: VEC_INTEL | VEC_VC | VEC_VCL";        
}

BUILT_SOMETHING=false

# evaluate command line
while getopts ":i:w:v:cah" opt; do
	case $opt in
	i) # iterations
		echo "Setting NUM_ITERATIONS to $OPTARG" >&2
		NUM_ITERATIONS=$OPTARG
		;;
	w) # warmup iterations
		echo "Setting NUM_WARMUP to $OPTARG" >&2
		NUM_WARMUP=$OPTARG
		;;
	v) # vec lib
		echo "Setting VECLIB to $OPTARG" >&2
		VECLIB=$OPTARG
		;;
	c) # CPU
		echo "Building for CPU" >&2
		BUILT_CPU=true
		;;
	a) # Accelerator
		echo "Building for Accelerator" >&2
		BUILT_ACC=true
		;;
	h) # usage
		usage
		exit 0
		;;
	:) # missing value
		echo "Option -$OPTARG requires an argument."
		usage
		exit 1
		;;
	\?)
		echo "Invalid option: -$OPTARG" >&2
		usage
		exit 1
		;;
  esac
done

if [ "$BUILT_CPU" = "true" ]
then
	build "$BUILD_DIR_HOST" "$OPTIONS_HOST -DNUM_ITERATIONS=${NUM_ITERATIONS} -DNUM_WARMUP=${NUM_WARMUP} -D${VECLIB}" "$INCLUDE_HOST" "$LIB_HOST"
	BUILT_SOMETHING=true
fi

if [ "$BUILT_ACC" = "true" ]
then
	build "$BUILD_DIR_MIC" "$OPTIONS_MIC -DNUM_ITERATIONS=${NUM_ITERATIONS} -DNUM_WARMUP=${NUM_WARMUP} -D${VECLIB}" "$INCLUDE_MIC" "$LIB_MIC"
	BUILT_SOMETHING=true
fi

if [ "$BUILT_SOMETHING" = "false" ]
then
	echo "Please use at least one of -c -a";
	usage
fi

# TODO: generate asm 

