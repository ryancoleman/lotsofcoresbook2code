CC	= icc #-g
CPP	= icc
LD	= icc
CXX	= icc

#flags for prefetch distance scaling
CFLAGS	= -O3 -intel-extensions  -finline -watch=mic-cmd -qopt-prefetch-distance=$(L2_DIST),$(L1_DIST) -DUSE_PARALLEL_VERSION
#flags for default prefetch distance
#CFLAGS	= -O3 -intel-extensions  -finline -watch=mic-cmd -DUSE_PARALLEL_VERSION

#CFLAGS	= -O3 -intel-extensions  -finline -watch=mic-cmd -opt-prefetch-distance=0,4 -DUSE_PARALLEL_VERSION
#LDFLAGS = -lrt
LDFLAGS = -lrt -qoffload-option,mic,ld,"-lpthread -rpath /opt/intel/mic/myo/lib/ -rpath /opt/intel/mic/coi/device-linux-release/lib/"
#-lscif -offload-build
#-offload-ldopts="-L$(MKLROOT)lib/mic;/opt/intel/composer_xe_2013.1.117/compiler/lib/mic/;$(MKLROOT)lib/intel64 -lscif"
ODIR	= ./bin
OBJ	= SmithWaterman.o DataLoader.o

CFLAGS += -qopt-report-phase=loop -qopt-report=5 -qopt-report-stdout
#CFLAGS += -qopt-report-phase=loop -qopt-report=5 -qopt-report-stdout

%.o: %.cpp
	 $(CC) -c -o $@ $< $(CFLAGS)	 

swmic: $(OBJ)
	 $(CC) -o $@ $^ $(LDFLAGS)

.PHONY: clean

clean:
		rm -rf *.o
		rm -rf *.out
		rm -rf swmic
