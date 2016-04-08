#!/bin/bash

# Copyright (c) 2015 Matthias Noack (ma.noack.pr@gmail.com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

# NOTE:
# Variables to change: CC, OPENCL_PATH, OPTIONS

CC=icpc # also change OPTIONS
OPTIONS="-std=c++11 -Wall -g -O3 -xHost -restrict -openmp -DDIM=7"
#CC=g++ # NOTE: also change OPTIONS below
#OPTIONS="-std=c++11 -Wall -g -O3 -march=native -Drestrict=__restrict__ -fopenmp -DDIM=7"
OPENCL_PATH=/opt/intel/opencl/ # NOTE: change this to your installation directory
#OPENCL_PATH=/usr/local/cuda

NUM_ITERATIONS=26 # including warmup below
NUM_WARMUP=1
INTEL_PREFETCH_LEVEL=1 # sets -auto-prefetch-level= for OpenCL compilation

HAM_PATH=thirdparty/ham
CLU_PATH_INCLUDE=thirdparty/CLU.git
CLU_PATH_LIB=thirdparty/CLU/clu_runtime

INCLUDE="-Iinclude -I${HAM_PATH}/include/ -I${OPENCL_PATH}/include -I${CLU_PATH_INCLUDE}"
LIB="-lrt -L${OPENCL_PATH}/lib64 -lOpenCL ${CLU_PATH_LIB}/libclu_runtime.a"

BUILD_DIR=bin
mkdir -p $BUILD_DIR

build()
{
	local CONFIG=$1
	local SUFFIX=$2
	$CC -c $OPTIONS $CONFIG $INCLUDE -o ${BUILD_DIR}/common.o src/common.cpp 
	$CC -c $OPTIONS $CONFIG $INCLUDE -o ${BUILD_DIR}/commutator_reference.o src/kernel/commutator_reference.cpp

	$CC $OPTIONS $CONFIG $INCLUDE -o ${BUILD_DIR}/benchmark_ocl${SUFFIX} ${BUILD_DIR}/commutator_reference.o ${BUILD_DIR}/common.o src/benchmark_ocl.cpp $LIB
}

usage ()
{ 
	echo "Usage and defaults:";
	echo -e "\t-c\t Build CPU variant.";
	echo -e "\t-a\t Build MIC (Accelerator) variant."; 
	echo -e "\t-g\t Build GPU variant.";  
	echo -e "\t-i ${NUM_ITERATIONS}\t Number of iterations (including warmups).";
	echo -e "\t-w ${NUM_WARMUP}\t Number of warmup iterations.";
	echo -e "\t-p ${INTEL_PREFETCH_LEVEL}\t Value used for the Intel-specific OpenCL compiler option '-auto-prefetch-level='.";        
}

BUILT_SOMETHING=false

# evaluate command line
while getopts ":i:w:p:cagh" opt; do
	case $opt in
	i) # iterations
		echo "Setting NUM_ITERATIONS to $OPTARG" >&2
		NUM_ITERATIONS=$OPTARG
		;;
	w) # warmup iterations
		echo "Setting NUM_WARMUP to $OPTARG" >&2
		NUM_WARMUP=$OPTARG
		;;
	p) # Prefetch level (Intel only)
		echo "Setting INTEL_PREFETCH_LEVEL to $OPTARG" >&2
		INTEL_PREFETCH_LEVEL=$OPTARG
		;;
	c) # CPU
		echo "Building for CPU" >&2
		BUILT_CPU=true
		;;
	a) # Accelerator
		echo "Building for Accelerator" >&2
		BUILT_ACC=true
		;;
	g) # GPU
		echo "Building for GPU" >&2
		BUILT_GPU=true
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
	build "-DDEVICE_TYPE=CL_DEVICE_TYPE_CPU -DNUM_ITERATIONS=${NUM_ITERATIONS} -DNUM_WARMUP=${NUM_WARMUP} -DINTEL_PREFETCH_LEVEL=${INTEL_PREFETCH_LEVEL} -DVEC_LENGTH=4 -DVEC_LENGTH_AUTO=16" "_cpu"
	BUILT_SOMETHING=true
fi

if [ "$BUILT_ACC" = "true" ]
then
	build "-DDEVICE_TYPE=CL_DEVICE_TYPE_ACCELERATOR -DNUM_ITERATIONS=${NUM_ITERATIONS} -DNUM_WARMUP=${NUM_WARMUP} -DINTEL_PREFETCH_LEVEL=${INTEL_PREFETCH_LEVEL} -DVEC_LENGTH=8 -DVEC_LENGTH_AUTO=16" "_mic"
	BUILT_SOMETHING=true
fi

if [ "$BUILT_GPU" = "true" ]
then
	build "-DDEVICE_TYPE=CL_DEVICE_TYPE_GPU -DNUM_ITERATIONS=${NUM_ITERATIONS} -DNUM_WARMUP=${NUM_WARMUP} -DVEC_LENGTH=8 -DVEC_LENGTH_AUTO=16" "_gpu"
	BUILT_SOMETHING=true
fi

if [ "$BUILT_SOMETHING" = "false" ]
then
	echo "Please use at least one of -c -a -g";
	usage
fi

