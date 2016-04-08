#### Edit the paths according to your build directories####
#FLANN_PATH= /home/guest1/visualsearch_orig/flann-1.8.4-src/
FLANN_PATH=../../flann-1.8.4-src
INTEL_PATH=/opt/intel/composer_xe_2015.3.187/compiler/lib/intel64
###########################################################

DEBUGFLAGS = -g -ggdb -O0 -pg
OPTFLAGS= -O3  
CFLAGS	= -Wno-write-strings $(OPTFLAGS) 
#CFLAGS	= -Wno-write-strings $(DEBUGFLAGS) 
LIBS	+= -lflann -lpthread
LIB_PATH	+= -L $(FLANN_PATH)/build/lib -Wl,-rpath,$(INTEL_PATH)
INC_PATH	+= -I  $(FLANN_PATH)/src/cpp

CC = icpc

SOURCE := server.cpp util.cpp
PROGS := server

OBJS = $(SOURCE:.cpp=.o)

ifneq ($(findstring $(MAKEFLAGS),s),s)
ifndef V
	QUIET_CC	= @echo '   ' CC $@;
	QUIET_LINK	= @echo '   ' LINK $@;
	QUIET_DEP	= @echo '   ' DEP $@;
endif
endif

all: $(PROGS) $(SCRIPTS) FORCE

.PHONY: all install clean
.PHONY: FORCE cscope

$(OBJS): $(SOURCE)
	$(CC) $(CFLAGS) $(INC_PATH) -c $(SOURCE)
	
server: $(OBJS)
	$(CC) $(LDFLAGS) $(CFLAGS) $(LIB_PATH) -o $@ $(OBJS) $(LIBS)

clean: FORCE
	-rm -f .depend $(OBJS) $(PROGS) core.* core 

