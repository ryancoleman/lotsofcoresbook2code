# FILE LICENSE TAG: SAMPLE 

# First, declare some useful macros

# Self-explanatory

COI_DIR ?= /usr/include/intel-coi

# We use the Composer XE compiler:

CC = icc -Wall -Werror-all

# Self-explanatory

COI_INC ?= -I$(COI_DIR)

# Self-explanatory

# Disabling some vectorization-related diagnostics as ComposerXE 2015 would abort compilation upon failing to
# # vectorize loops which include malloc() if -O3 or -O2 is used
COMMON_FLAGS ?= $(COI_INC) -I../../common -I../../../include -I../../../src/include -mkl -diag-disable 13368 -diag-disable 15527

# Two configurations can be built: DEBUG or RELEASE.  By default, RELEASE config is built.
# To create DEBUG config, make CFG=DEBUG

CFG = RELEASE
DEBUG_FLAGS = -g -O0 -D_DEBUG
RELEASE_FLAGS = -DNDEBUG -O3 

ifeq (RELEASE,$(CFG))
OPT = $(RELEASE_FLAGS)
else
OPT = $(DEBUG_FLAGS)
endif

# CC_FLAGS refers to both host and device compiler flags:

CC_FLAGS = $(COMMON_FLAGS) $(OPT) -Wl,--enable-new-dtags 

# Self-explanatory:

HOST_FLAGS = $(CC_FLAGS) -lcoi_host
DEV_FLAGS = $(CC_FLAGS) -lcoi_device -rdynamic

# The binaries built from the Makefile will go into either bin/host or bin/dev
# representing the host bits, and device bits respectively:

BIN_HOST = ../../../bin/host
BIN_DEV  = ../../../bin/dev

# HOST_CC means CC compiler for host bits
# DEV_CC  means CC compiler for device bits

HOST_CC ?= $(CC)       -openmp
DEV_CC  ?= $(CC) -mmic -openmp

# The all target - builds everything

all: $(BIN_HOST) $(BIN_DEV) $(BIN_HOST)/lu_tiled_hstreams $(BIN_DEV)/lu_sink_1.so

# Create the bin dirs:

$(BIN_HOST) $(BIN_DEV):
	mkdir -p $@

# Clean the build:

clean:
	rm -fr $(BIN_HOST)/lu_tiled_hstreams $(BIN_DEV)/lu_sink_1.so  *.o ../../common/*.o

######################################################################################

# End of preliminaries.

# Next define some macros for building libhstreams and the test app:

# Define CPP source file lists:

LU_HOST_CPP_SRCS     = lu_tile_hstreams.cpp ../../common/dtime.cpp ../../common/hStreams_custom.cpp ../../common/matrices_generator.cpp
LU_DEV_CPP_SRCS      = ../../common/hStreams_custom_kernels_sink.cpp

# Derive object file lists from CPP lists:

LU_HOST_CPP_OBJS     = $(LU_HOST_CPP_SRCS:.cpp=.host.o)
LU_DEV_CPP_OBJS      = $(LU_DEV_CPP_SRCS:.cpp=.dev.o)

# Building lu host object files uses this rule:

$(LU_HOST_CPP_OBJS): %.host.o: %.cpp
	$(HOST_CC) $(HOST_FLAGS) -c $^ -o $@

# Building lu device object files use this rule:

$(LU_DEV_CPP_OBJS): %.dev.o: %.cpp
	$(DEV_CC) $(DEV_FLAGS) -fPIC -shared -c $^ -o $@

# Building lu, for the host uses this rule:

$(BIN_HOST)/lu_tiled_hstreams: $(LU_HOST_CPP_OBJS)
	$(HOST_CC) $(HOST_FLAGS) -o $@ $^ -L$(BIN_HOST)  -lhstreams_source

# Building lu, for the device uses this rule:
$(BIN_DEV)/lu_sink_1.so: $(LU_DEV_CPP_OBJS)
	$(DEV_CC) $(DEV_FLAGS) -fPIC -shared -o $@ $^ -L$(BIN_DEV) -lhstreams_sink -Wl,-soname,lu_sink_1.so -static-intel

