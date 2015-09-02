# C compiler
CC = icc
CC_FLAGS = -g -O3 -fopenmp -offload-option,mic,compiler,"-no-opt-prefetch -qopt-report-phase=vec -qopt-report=5" 

all: needle needle_offload

needle: 
	$(CC) $(CC_FLAGS) needle.cpp -o needle 

needle_offload:
	$(CC) $(CC_FLAGS) -DOMP_OFFLOAD needle.cpp -o needle_offload

clean:
	rm -f needle needle_offload
