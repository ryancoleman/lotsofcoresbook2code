#    The MIT License (MIT)
#    
#    Copyright (c) 2015 OpenVec
#    
#    Permission is hereby granted, free of charge, to any person obtaining a copy
#    of this software and associated documentation files (the "Software"), to deal
#    in the Software without restriction, including without limitation the rights
#    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#    copies of the Software, and to permit persons to whom the Software is
#    furnished to do so, subject to the following conditions:
#    
#    The above copyright notice and this permission notice shall be included in all
#    copies or substantial portions of the Software.
#    
#    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#    SOFTWARE.
#    
#    
#    Authors:
#    Paulo Souza
#    Leonardo Borges
#    Cedric Andreolli
#    Philippe Thierry
#

#Default KNC
knc_target?=mic0

ifeq ($(arch),intel) # Intel compiler with best flags for the host invoking the compiler
CC=icc -std=c99
CPP=icpc
CFLAGS=-Wall -O3 -xHOST
OMP_FLAGS=-openmp

else # Intel KNC
ifeq ($(arch),knc)
CC=icc -std=c99
CPP=icpc
CFLAGS=-Wall -O3 -mmic
OMP_FLAGS=-openmp

else # ARM neon
ifeq ($(arch),neon)
CC=gcc -std=c99
CPP=g++
CFLAGS=-Wall -O3 -mfpu=neon -mtune=cortex-a9 -mcpu=cortex-a9
OMP_FLAGS=

else # GCC AVX
ifeq ($(arch),gccavx)
CC=gcc -std=c99
CPP=g++
CFLAGS=-Wall -O3 -mavx
OMP_FLAGS=

else # Default gcc with sse
CC=gcc -std=c99
CPP=g++
CFLAGS=-Wall -O3 -msse
OMP_FLAGS=
endif
endif
endif
endif

all:  selftest selftest.double cmp_to_zero benchmark_saxpy easy_conditional memcpy \
  reduction stencil_unroll benchmark_saxpy_unroll stencil_naive \
  load_store vector_tail test_floor_ceil floor
	echo "Done."

cleanall:	clean

ifeq ($(arch),knc)
run_selftest:
	scp cmp_to_zero.exec $(knc_target):/tmp;ssh $(knc_target) '/tmp/cmp_to_zero.exec | grep Passed'
	scp easy_conditional_cpp.exec $(knc_target):/tmp;ssh $(knc_target) '/tmp/easy_conditional_cpp.exec | grep Passed'
	scp easy_conditional.exec $(knc_target):/tmp;ssh $(knc_target) '/tmp/easy_conditional.exec | grep Passed'
	scp load_store_cpp.exec $(knc_target):/tmp;ssh $(knc_target) '/tmp/load_store_cpp.exec | grep Passed'
	scp memcpy.exec $(knc_target):/tmp;ssh $(knc_target) '/tmp/memcpy.exec | grep Passed'
	scp reduction.exec $(knc_target):/tmp;ssh $(knc_target) '/tmp/reduction.exec | grep Passed'
	scp selftest_cpp.exec $(knc_target):/tmp;ssh $(knc_target) '/tmp/selftest_cpp.exec | grep Passed'
	scp selftest.exec $(knc_target):/tmp;ssh $(knc_target) '/tmp/selftest.exec | grep Passed'
	scp selftest.double_cpp.exec $(knc_target):/tmp;ssh $(knc_target) '/tmp/selftest.double_cpp.exec | grep Passed'
	scp selftest.double.exec $(knc_target):/tmp;ssh $(knc_target) '/tmp/selftest.double.exec | grep Passed'
	scp vector_tail.exec $(knc_target):/tmp;ssh $(knc_target) '/tmp/vector_tail.exec | grep Passed'
	scp test_floor_ceil_cpp.exec $(knc_target):/tmp;ssh $(knc_target) '/tmp/test_floor_ceil_cpp.exec | grep Passed'

else
run_selftest:
	./cmp_to_zero.exec | grep Passed
	./easy_conditional_cpp.exec | grep Passed
	./easy_conditional.exec | grep Passed
	./load_store_cpp.exec | grep Passed
	./memcpy.exec | grep Passed
	./reduction.exec | grep Passed
	./selftest_cpp.exec | grep Passed
	./selftest.exec | grep Passed
	./selftest.double_cpp.exec | grep Passed
	./selftest.double.exec | grep Passed
	./vector_tail.exec | grep Passed
	./test_floor_ceil_cpp.exec | grep Passed
endif

clean:
	\rm -f *.exec *.o

selftest:
	$(CC) $(CFLAGS) selftest.c -lm -o selftest.exec
	$(CPP) $(CFLAGS) selftest.cpp -o selftest_cpp.exec

selftest.double:
	$(CC) $(CFLAGS) selftest.double.c -lm -o selftest.double.exec
	$(CPP) $(CFLAGS) selftest.double.cpp -o selftest.double_cpp.exec

cmp_to_zero:
	$(CC) $(CFLAGS) cmp_to_zero.c -o cmp_to_zero.exec

benchmark_saxpy:
	$(CC) $(CFLAGS) benchmark_saxpy.c -o benchmark_saxpy.exec
	$(CPP) $(CFLAGS) benchmark_saxpy.cpp -o benchmark_saxpy_cpp.exec

benchmark_saxpy_unroll:
	$(CC) $(CFLAGS) benchmark_saxpy_unroll.c -o benchmark_saxpy_unroll.exec
	$(CPP) $(CFLAGS) benchmark_saxpy_unroll.cpp -o benchmark_saxpy_unroll_cpp.exec

easy_conditional:
	$(CC) $(CFLAGS) easy_conditional.c -o easy_conditional.exec
	$(CPP) $(CFLAGS) easy_conditional.cpp -o easy_conditional_cpp.exec

memcpy:
	$(CC) $(CFLAGS) memcpy.c -o memcpy.exec

reduction:
	$(CC) $(CFLAGS) reduction.c -o reduction.exec

stencil_naive:
	$(CC) $(CFLAGS) $(OMP_FLAGS) stencil_naive.c -o stencil_naive.exec
	$(CPP) $(CFLAGS) $(OMP_FLAGS) stencil_naive.cpp -o stencil_naive_cpp.exec

stencil_unroll:
	$(CC) $(CFLAGS) $(OMP_FLAGS) stencil_unroll.c -o stencil_unroll.exec
	$(CPP) $(CFLAGS) $(OMP_FLAGS) stencil_unroll.cpp -o stencil_unroll_cpp.exec

floor:
	$(CPP) $(CFLAGS) floor.cpp -o floor_cpp.exec

load_store:
	$(CPP) $(CFLAGS) load_store.cpp -o load_store_cpp.exec

test_floor_ceil:
	$(CPP) $(CFLAGS) test_floor_ceil.cpp -o test_floor_ceil_cpp.exec

vector_tail:
	$(CC) $(CFLAGS) vector_tail.c -o vector_tail.exec
