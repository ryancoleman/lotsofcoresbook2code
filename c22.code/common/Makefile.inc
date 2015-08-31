################################################
## This file contains default settings for 

ROOT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

CINC:=-I$(ROOT_DIR)
COPT:=-g -O3 -fno-inline

CXXOPT:=
CLIB:=-g -lm -lrt
CXXLIB:=


ifeq ($(COMP),gnu)
        CXX:=g++
        CC:=gcc
		FC:=gfortran
else
	CXX:=icpc
	CC:=icc
	FC:=ifort
	COPT+=-qopt-report:5 -qopt-report-phase=vec
endif

CXXOPT+=-I$(TBBROOT)/include

ifeq ($(PLAT),knc)
	COPT+=-mmic -D TIMER_ITERATIONS=5000
	CLIB+=-mmic
	GETTIME_O:=../common/get_time_MIC.o
	CXXLIB+=-L$(TBBROOT)/lib/mic
	EXECMD:=micnativeloadex ./main
else
	COPT+=-m64
	GETTIME_O:=../common/get_time.o
	CXXLIB+=-L$(TBBROOT)/lib/intel64/gcc4.4
	EXECMD:=./main
endif

CXXLIB+=-ltbbmalloc -ltbb

ifneq ($(OMP),no)
	ifeq ($(CC),icc)
		COPT+=-qopenmp-simd -vecabi=legacy
	else
		COPT+=-fopenmp-simd
	endif
	COPT+=-D_OPENMP
endif
