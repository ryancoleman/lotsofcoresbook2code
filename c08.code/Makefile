# Copyright (c) 2015 Intel Corporation
# All rights reserved.
# 
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are met:
# 
#  * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#  * Neither the name of Intel Corporation
#    nor the names of its contributors may be used to endorse or promote
#    products derived from this software without specific prior written
#    permission.
# 
#    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#    ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#    DISCLAIMED. IN NO EVENT SHALL INTEL CORPORATION BE LIABLE FOR ANY
#    DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
#    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
#    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
#    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# Source files
#  original source file
#  	am_call_original.cpp
#
#  Optimized source file using #pragma simd
#       am_call.cpp
#
#  Optimized source file using C/C++ vector extension
#  	am_call_vector.cpp
#  
# Executive files
#  Original Version - Scalar code runs in one thread
#       am_call_original
#
#  Scalar Serial - Scalar, serial  Optimized code
#       am_call_ss
#
#  Vector Serial - Scalar Serial Version with #pragma SIMD vectorization
#       am_call_vs
#
#  Scalar Parallel - Scalar Serial code with OpenMP parallelization
#       am_call_sp
#
#  Vector Parallel - Vector Serial code with OpenMP parallelization
#       am_call_vp
#
#  Vector Extension Vectorization
#  	am_call_vector
#
#  KNC Vector Parallel -- vector parallel version for KNC target
#       am_call_knc
#
CXX=icpc
COMPILER_VERSION:="$(CXX)-$(shell $(CXX) --version | head -n1 | cut -d' ' -f3)"
CXXFLAGS=-O3 -ipo -restrict -fimf-precision=low -fimf-domain-exclusion=31  \
         -fno-alias -qopt-report=5 

LDFLAGS= -ltbbmalloc

all: am_call_ss am_call_vs am_call_sp am_call_vp am_call_original am_call_vector am_call_knc

am_call_ss: am_call.cpp 
	$(CXX) $(LDFLAGS) -xHOST -openmp-stubs -no-vec -DSCALAR -unroll0 $(CXXFLAGS) -DCOMPILER_VERSION=\"$(COMPILER_VERSION)\" -o $@ $^

am_call_vs: am_call.cpp 
	$(CXX) $(LDFLAGS) -xHOST -openmp-stubs -unroll0 $(CXXFLAGS) -DCOMPILER_VERSION=\"$(COMPILER_VERSION)\" -o $@ $^

am_call_sp: am_call.cpp 
	$(CXX) $(LDFLAGS) -xHOST -openmp -DSCALAR -no-vec -unroll0 $(CXXFLAGS) -DCOMPILER_VERSION=\"$(COMPILER_VERSION)\" -o $@ $^

am_call_vp: am_call.cpp 
	$(CXX) $(LDFLAGS) -xHOST -openmp $(CXXFLAGS) -DCOMPILER_VERSION=\"$(COMPILER_VERSION)\" -o $@ $^

am_call_original: am_call_original.cpp 
	$(CXX) $(LDFLAGS) -xHOST $(CXXFLAGS) -DCOMPILER_VERSION=\"$(COMPILER_VERSION)\" -o $@ $^

am_call_vector: am_call_vector.cpp 
	$(CXX) $(LDFLAGS) -xCORE_AVX2 $(CXXFLAGS) -DCOMPILER_VERSION=\"$(COMPILER_VERSION)\" -o $@ $^

am_call_knc: am_call.cpp 
	$(CXX) $(LDFLAGS) -mmic -openmp $(CXXFLAGS) -DCOMPILER_VERSION=\"$(COMPILER_VERSION)\" -o $@ $^
clean:
	rm -r -f am_call_ss am_call_vs am_call_sp am_call_vp am_call_original am_call_vector am_call_knc
