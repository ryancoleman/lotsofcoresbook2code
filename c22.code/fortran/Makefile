##
## Create binaries for section 4 which describes vectorization in C
##

include ../common/Makefile.inc

SRCFILES:=$(notdir $(wildcard ./*.f90 ))
OBJS:= #../common/get_time.o
OBJS+=$(subst .f90,.o,$(SRCFILES))

ifneq ($(PLAT), knc)
	ifeq ($(CXX), icpc)
        	COPT+=-xHost
	else
        	COPT+=-march=core-avx2
	endif
endif

build: main

all: main run

run: main
	./main

main:$(OBJS) $(GETTIME_O)
	$(FC) -o main $(OBJS) $(GETTIME_O) $(CLIB)
	
%.o: %.f90
	$(FC) -c $(COPT) $(CXXOPT) -o $@ $<
	
clean:
	rm *.o main -rf
