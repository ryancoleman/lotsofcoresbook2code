##
## Create binaries for section 2 which describes vectorization in C
##

include ../common/Makefile.inc

SRCFILES:=$(notdir $(wildcard ./*.c ))
OBJS:=#../common/get_time.o
OBJS+=$(subst .c,.o,$(SRCFILES))

ifneq ($(PLAT), knc)
	ifeq ($(CC), icc)
        	COPT+=-xSSE4.2 -axCORE-AVX2
	else
        	COPT+=-march=core-avx2
	endif
endif

build: main

all: main run

run: main
	$(EXECMD)

main:$(OBJS)
	$(CC) -o main $(OBJS) $(GETTIME_O) $(CLIB)
	
%.o: %.c $(INCFILES)
	$(CC) -c $(COPT) $(CINC) -o $@ $<
	
clean:
	rm *.o main -f
