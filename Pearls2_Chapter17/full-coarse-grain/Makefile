SHELL= /bin/bash
FC= mpiifort 
FCFLAGS= -openmp  -g -traceback -O3 -mt_mpi
FCFLAGS_MIC= -mmic -openmp -g -traceback -O3 -mt_mpi

SRC= stest_hybrid_coarse.f90 shado_xchange_EW_coarse.f90\
     shado_xchange_NS_coarse.f90 shado_xchange_TB_coarse.f90 

default: all
all:  exe exe_mic

exe: $(SRC)
	$(FC) $(FCFLAGS) $(SRC) -o stest_hybrid_coarse.x

exe_mic: $(SRC)
	$(FC) $(FCFLAGS_MIC) $(SRC) -o stest_hybrid_coarse.x.MIC

clean:
	/bin/rm -f *.o *.x *.x.MIC

