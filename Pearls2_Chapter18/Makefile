# Copyright Intel Corporation 2013

SHELL = /bin/sh

#root directory of Intel compilers
ICC_ROOT=/opt/intel/composer_xe_2015.2.164
#root directory of Intel MPI
IMPI_ROOT=/opt/intel/impi/5.0.3.048/

#16 beta compilers
#ICC_ROOT=/opt/intel/compilers_and_libraries_2016.0.056/linux
#IMPI_ROOT=/opt/intel/compilers_and_libraries_2016.0.056/linux/mpi/

#compilers
MPICXX     = mpiicpc
MPICC      = mpiicc
CXX=icpc
CC=icc

# common compiler flags
ICFLAGS    =  -restrict -unroll  -O3 -qopenmp -g  -DINLINE=inline -DMUST_INLINE=inline -DHAVE_MKL
DEBUGFLAGS = -debug inline-debug-info
# host target option
HOSTFLAGS   = -xAVX
# mic target option 
MICFLAGS    = -mmic -fp-model source
# C++ flags
CXXFLAGS   = $(ICFLAGS) -std=c++11 $(DEBUGFLAGS)
# CC flags
CFLAGS     = $(ICFLAGS) -std=c99 $(DEBUGFLAGS)
# Flags to enable CREW on mic
CREWFLAGS   = -D__INTEL_CREW -mP2OPT_hpo_par_crew_codegen=T
OFFLOADFLAGS= -offload-attribute-target=mic -DMIC_OFFLOAD

#directory where binaries and libraries are located on mic
MICDIR     = /tmp

srcdir = .
IFLAGS = -I$(srcdir) -I$(MKLROOT)/include
LDFLAGS=-lpthread -lm
MKL_LIBS=-Wl,--start-group -L$(MKLROOT)/lib/intel64 -lmkl_cdft_core -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group

#static link to MKL libraries
MKL_MIC_LIBS=-Wl,--start-group $(MKLROOT)/lib/mic/libmkl_intel_lp64.a $(MKLROOT)/lib/mic/libmkl_intel_thread.a $(MKLROOT)/lib/mic/libmkl_core.a -Wl,--end-group
#MKL_MIC_LIBS=-Wl,--start-group -L$(MKLROOT)/lib/mic -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group

all: dgemm fft3d diffusion

dgemm: $(OBJS) $(srcdir)/dgemm.cpp
	$(MPICXX) $(IFLAGS) $(CXXFLAGS) $(HOSTFLAGS) -o host/dgemm $(srcdir)/dgemm.cpp $(MKL_LIBS) $(LDFLAGS)
	$(MPICXX) $(IFLAGS) $(CXXFLAGS) $(MICFLAGS) -o mic/dgemm $(srcdir)/dgemm.cpp $(MKL_MIC_LIBS) $(LDFLAGS)

fft3d: $(OBJS) $(srcdir)/fft3d.cpp
	$(MPICXX) $(IFLAGS) $(CXXFLAGS) $(HOSTFLAGS) -o host/fft3d $(srcdir)/fft3d.cpp $(MKL_LIBS) $(LDFLAGS)
	$(MPICXX) $(IFLAGS) $(CXXFLAGS) $(MICFLAGS) -o mic/fft3d $(srcdir)/fft3d.cpp $(MKL_MIC_LIBS) $(LDFLAGS)

diffusion: $(OBJS) $(srcdir)/diffusion.c
	$(CC) $(IFLAGS) $(CFLAGS) $(HOSTFLAGS) -o host/diffusion  $(srcdir)/diffusion.c 
	$(CC) $(IFLAGS) $(CFLAGS) $(CREWFLAGS) $(MICFLAGS) -o mic/diffusion  $(srcdir)/diffusion.c 

run_host: host/dgemm
	OMP_NUM_THREADS=4,1 host/dgemm
	OMP_NUM_THREADS=4,1 host/fft3d
	OMP_NUM_THREADS=16,1 host/diffusion mic nx=256 count=10 nf=10 out=

mic_dgemm: mic/dgemm mic/fft3d mic/diffusion
	@echo ""
	@echo "Testing dgemm using MPI=1, OpenMP=60 and MKL=4"
	@echo ""
	export I_MPI_MIC=enable
	ssh mic0 LD_LIBRARY_PATH=$(MICDIR) \
	OMP_PLACES=threads OMP_PROC_BIND=spread,close \
	OMP_NESTED=true OMP_NUM_THREADS=60,4 \
	KMP_HOT_TEAMS_MAX_LEVEL=2 KMP_HOT_TEAMS_MODE=1 /tmp/dgemm

mic_fft3d: mic/fft3d
	@echo ""
	@echo "Testing fft3d using MPI=1, OpenMP=60 and MKL=4"
	@echo ""
	ssh mic0 LD_LIBRARY_PATH=$(MICDIR) \
	OMP_PLACES=threads OMP_PROC_BIND=spread,close \
	OMP_NESTED=true OMP_NUM_THREADS=60,4 \
	KMP_HOT_TEAMS_MAX_LEVEL=2 KMP_HOT_TEAMS_MODE=1 /tmp/fft3d -g 64 -i 10

mic_diffusion: mic/diffusion
	@echo ""
	@echo "Testing diffusion mic using MPI=1, OpenMP=240"
	@echo ""
	ssh mic0 LD_LIBRARY_PATH=$(MICDIR) \
	OMP_PLACES=threads OMP_PROC_BIND=close \
	OMP_NUM_THREADS=240 /tmp/diffusion mic nx=240 count=10 nf=10 out=

miccopy:
	scp mic/* mic0:$(MICDIR)

miclibcopy:
	scp $(ICC_ROOT)/compiler/lib/mic/libiomp5.so mic0:$(MICDIR)
	scp $(IMPI_ROOT)/mic/bin/pmi_proxy mic0:/tmp
	scp $(IMPI_ROOT)/mic/lib/libmp*.s* mic0:/tmp/

	@echo ""
	@echo "REMEMBER TO export LD_LIBRARY_PATH=$(MICDIR) ON THE COPROCESSOR (if needed)"
	@echo ""

clean:
	@rm -rf mic/* host/*

# * Other names and brands may be claimed as the property of others.



