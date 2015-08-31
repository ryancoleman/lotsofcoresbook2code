##
## Create binaries for section 4 which describes vectorization in C
##

include ../common/Makefile.inc

SRCFILES:=$(notdir $(wildcard ./*.cpp ))
OBJS:= #../common/get_time.o
OBJS+=$(subst .cpp,.o,$(SRCFILES))

ifneq ($(PLAT), knc)
	ifeq ($(CXX), icpc)
        	COPT+=-xSSE4.2 -axCORE-AVX2
	else
        	COPT+=-march=core-avx2
	endif
endif

build: main

all: main run

run: main
	./main

main:$(OBJS) $(GETTIME_O)
	$(CXX) -o main $(OBJS) $(GETTIME_O) $(CLIB) $(CXXLIB)
	
%.o: %.cpp $(INCFILES)
	$(CXX) -c $(COPT) $(CXXOPT) $(CINC) -o $@ $<
	
clean:
	rm *.o main -rf
