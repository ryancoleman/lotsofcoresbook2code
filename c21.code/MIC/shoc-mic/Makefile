# Basics
CC	= icc
CPP     = icc
LD	= icc
CXX     = icc

SHOC_COMMON	=./common
MICROOT		=./bin
BINDIR		=./bin
OBJDIR		=./obj

# Workload folders
VPATH=$(SHOC_COMMON) level0 md reduction scan triad spmv sort fft gemm s3d 

# Common objects
COMMON_OBJS 	= main.o Option.o OptionParser.o Timer.o ResultDatabase.o ProgressBar.o
COMMON_OBJFILES = $(addprefix $(OBJDIR)/, $(COMMON_OBJS))

# Workload objects
BENCH_OBJS = BusSpeedReadback.o BusSpeedDownload.o DeviceMemory.o MaxFlops.o \
	     MD.o S3D.o			\
	     Reduction.o Scan.o Triad.o \
	     Spmv.o GEMM.o FFT.o
BENCHMARKPROG = $(patsubst %.o,$(BINDIR)/%,$(BENCH_OBJS))

# Flags to include MIC MKL library
MIC_MKL_LIBS	 = -lmkl_core

# Flags to enable compiler reporting - Modify according detail level needs
REPORTING	 = # -vec-report0 -opt-report 1 -vec-report0 -opt-report-phase:offload -openmp-report

# Compiler flags
CFLAGS           = -O3 -openmp -parallel -intel-extensions -I$(SHOC_COMMON) -mP2OPT_hlo_pref_issue_second_level_prefetch=F -mP2OPT_hlo_pref_issue_first_level_prefetch=F	\
		   -offload-option,mic,compiler,"-mP2OPT_hpo_vec_check_dp_trip=F -fimf-precision=low -fimf-domain-exclusion=15" $(REPORTING)

# Workload specific compiler flags
STENCIL_CPPFLAGS =

# Linker flags
LDFLAGS 	= -mkl

$(OBJDIR)/%.o: %.cpp
	$(CC) -c $< $(CFLAGS) $(CXXFLAGS) -o $@

$(BINDIR)/%: $(OBJDIR)/%.o 
	$(CC) -o $@ $(CFLAGS) $(CXXFLAGS) $(LDFLAGS) $(COMMON_OBJFILES) $< $(LIBS) 
 
all : $(BENCHMARKPROG)
	for d in $(SUBDIRS); do (cd $$d; $(MAKE) ); done

$(BENCHMARKPROG) : $(COMMON_OBJFILES)

clean :
	rm -f $(OBJDIR)/*
	rm -f $(BINDIR)/*
	make -C ./stencil2d clean
	
# Individual Targets

# FFT
fft : $(BINDIR)/FFT

$(BINDIR)/FFT : $(COMMON_OBJFILES)

# GEMM
gemm : $(BINDIR)/GEMM

$(BINDIR)/GEMM : $(COMMON_OBJFILES)

# MD
md : $(BINDIR)/MD

$(BINDIR)/MD : $(COMMON_OBJFILES)

# Reduction
reduction : $(BINDIR)/Reduction

$(BINDIR)/Reduction : $(COMMON_OBJFILES)

# SCAN
scan: $(BINDIR)/Scan

$(BINDIR)/Scan : $(COMMON_OBJFILES)

# Sort
$(BINDIR)/Sort : $(COMMON_OBJFILES)

sort: $(BINDIR)/Sort

# S3D
$(BINDIR)/S3D : $(COMMON_OBJFILES)

s3d : $(BINDIR)/S3D

# Stencil 2D
stencil2d:
	make -C ./stencil2d

.PHONY: stencil2d clean all
