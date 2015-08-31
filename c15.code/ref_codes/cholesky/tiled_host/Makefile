# FILE LICENSE TAG: SAMPLE 

# First, declare some useful macros

# We use the Composer XE compiler:

CC = icc -Wall -Werror-all

# Self-explanatory

# Disabling some vectorization-related diagnostics as ComposerXE 2015 would abort compilation upon failing to
# vectorize loops which include malloc() if -O3 or -O2 is used
COMMON_FLAGS ?= -I../../common  -mkl -diag-disable 13368 -diag-disable 15527

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

CC_FLAGS = $(COMMON_FLAGS) $(OPT) -Wl,--enable-new-dtags 

# Self-explanatory:

HOST_FLAGS = $(CC_FLAGS)

# The binaries built from the Makefile will go into bin/host

BIN_HOST = ../../../bin/host

# HOST_CC means CC compiler for host bits

HOST_CC ?= $(CC)       -openmp

# The all target

all: $(BIN_HOST) $(BIN_HOST)/cholesky_tiled_host

# Create the bin dir:

$(BIN_HOST):
	mkdir -p $@

# Clean the build:

clean:
	rm -fr $(BIN_HOST)/cholesky_tiled_host  *.o ../../common/*.o

######################################################################################

# End of preliminaries.

# Define CPP source file lists:

CHOLESKY_HOST_CPP_SRCS     = cholesky_tile.cpp ../../common/dtime.cpp ../../common/matrices_generator.cpp

# Derive object file lists from CPP lists:

CHOLESKY_HOST_CPP_OBJS     = $(CHOLESKY_HOST_CPP_SRCS:.cpp=.host.o)

# Building cholesky host object files uses this rule:

$(CHOLESKY_HOST_CPP_OBJS): %.host.o: %.cpp
	$(HOST_CC) $(HOST_FLAGS) -c $^ -o $@

# Building cholesky, for the host uses this rule:

$(BIN_HOST)/cholesky_tiled_host: $(CHOLESKY_HOST_CPP_OBJS)
	$(HOST_CC) $(HOST_FLAGS) -o $@ $^ 
