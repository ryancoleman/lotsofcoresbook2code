SHELL= /bin/bash
FC= mpiifort 
FCFLAGS= -openmp -g -traceback -O3
FCFLAGS_MIC= -mmic -openmp -g -traceback -O3

SRC= stest_hybrid_fine.f shado_xchange_EW_fine.f \
     shado_xchange_NS_fine.f shado_xchange_TB_fine.f 

default: all
all:  exe exe_mic

exe: $(SRC)
	$(FC) $(FCFLAGS) $(SRC) -o stest_hybrid_fine.x

exe_mic: $(SRC)
	$(FC) $(FCFLAGS_MIC) $(SRC) -o stest_hybrid_fine.x.MIC

clean:
	/bin/rm -f *.o *.x *.x.MIC

