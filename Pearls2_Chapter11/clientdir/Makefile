
DEBUGFLAGS = -g
#OPTFLAGS= -O3  
CFLAGS	= -Wno-write-strings $(DEBUGFLAGS)  
CC = icpc

SOURCE := client.cpp util.cpp
PROGS := client

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
	$(CC) $(CFLAGS) -c $(SOURCE)
	
client: $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS)

clean: FORCE
	-rm -f .depend $(OBJS) $(PROGS) core.* core 

