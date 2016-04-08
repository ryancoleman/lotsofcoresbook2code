SHELL= /bin/bash
FC= mpiifort 
FCFLAGS= -openmp -g -traceback -O3
FCFLAGS_MIC= -mmic -openmp -g -traceback -O3

SRC= stest_hybrid_part_coarse.f90 shado_xchange_EW_part_coarse.f90\
     shado_xchange_NS_part_coarse.f90 shado_xchange_TB_part_coarse.f90 

default: all
all:  exe exe_mic

exe: $(SRC)
	$(FC) $(FCFLAGS) $(SRC) -o stest_hybrid_part_coarse.x

exe_mic: $(SRC)
	$(FC) $(FCFLAGS_MIC) $(SRC) -o stest_hybrid_part_coarse.x.MIC

clean:
	/bin/rm -f *.o *.x *.x.MIC

