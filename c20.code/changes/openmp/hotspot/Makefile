# C compiler
CC = icc
CC_FLAGS = -g -fopenmp -O2 -offload-option,mic,compiler,"-no-opt-prefetch -qopt-report-phase=vec -qopt-report=5"
#CC_FLAGS = -g -fopenmp -O2 -offload-option,mic,compiler,"-qopt-report-phase=vec -qopt-report=5"
#CC_FLAGS = -g -fopenmp -O2

all: hotspot hotspot_offload 


hotspot: hotspot_openmp.cpp Makefile 
	$(CC) $(CC_FLAGS) hotspot_openmp.cpp -o hotspot 

hotspot_offload: hotspot_openmp.cpp Makefile
	$(CC) $(CC_FLAGS) -DOMP_OFFLOAD hotspot_openmp.cpp -o hotspot_offload

clean:
	rm -f hotspot hotspot_offload
