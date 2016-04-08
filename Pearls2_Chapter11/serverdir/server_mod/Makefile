
DEBUGFLAGS = -g -ggdb -O0 -pg
OPTFLAGS= -O3  
#CFLAGS	= -Wno-write-strings $(OPTFLAGS) 
CFLAGS	= -Wno-write-strings $(DEBUGFLAGS) 
LIBS	+= -lflann -lpthread
LIB_PATH	+= -L ~/intel-img/flann-1.8.4-cpu/build/lib
INC_PATH	+= -I  /home1/02857/prash/intel-img/flann-1.8.4-cpu/src/cpp
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

