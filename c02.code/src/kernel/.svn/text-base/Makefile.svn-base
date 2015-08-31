include ../macros.make

SHELL = /bin/sh
FFLAGS = $(PHYSFLAGS) $(R8FLAG) $(OMPFLAG) $(MICFLAGS)
LIBS = $(GPTLLIB)

BASENAME = wsm6_kernel
SRC = $(BASENAME).F90
OBJ = $(BASENAME).o

.SUFFIXES : .o .f90 .F90

$(WSM6EXE) : $(SRC)
	$(CPP) $(CPP_FLAGS) $(SRC) > $(BASENAME).f90
	$(CHUNKCMD) $(BASENAME).f90 > $(BASENAME)_chunk.f90
	$(LITERALKCMD) $(BASENAME)_chunk.f90 > $(BASENAME)_kcmd.f90
	${FCX} -o $(OBJ) ${FFLAGS} -c $(BASENAME)_kcmd.f90
	$(FCX) -o $@ $(OMPFLAG) $(MICFLAGS) $(OBJ) $(GPTLOBJS) $(LIBS)

clean:
	$(RM) *.f90 *.o *.mod
	- find ./ -type l -exec rm  -f  {} \;

