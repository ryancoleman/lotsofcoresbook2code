################
# CPU

FC=ifort
FFLAGS = -O3 -openmp
FFLAGS += -w -ftz -fno-alias 
FFLAGS += -g -debug inline-debug-info
FFLAGS += -xHost
FFLAGS += -vec-report6 #-opt-report

FFLAGS += -fp-model precise -FR # WRF 3.4 default CPU options for Intel compiler
#FFLAGS += -auto -fp-model fast=1 -no-prec-div -no-prec-sqrt -auto -align array64byte # WRF 3.6 default options for CPU

#
################


################
# Xeon Phi

FC_MIC=ifort
FFLAGS_MIC = -O3 -openmp
FFLAGS_MIC += -w -ftz -fno-alias 
FFLAGS_MIC += -g -debug inline-debug-info
FFLAGS_MIC += -vec-report6 #-opt-report

FFLAGS_MIC += -mmic -fp-model precise -FR -auto -align array64byte -opt-assume-safe-padding # precise math
#FFLAGS_MIC += -mmic -fimf-precision=low -fimf-domain-exclusion=15 -opt-assume-safe-padding -opt-streaming-stores always -opt-streaming-cache-evict=0 -mP2OPT_hlo_pref_use_outer_strategy=F -fp-model fast=1 -no-prec-div -no-prec-sqrt -auto -align array64byte # WRF default options for MIC

#
################


all: gsfcgce_fortran gsfcgce_fortran_opt gsfcgce_fortran_MIC gsfcgce_fortran_opt_MIC

gsfcgce_fortran: module_model_constants.o module_mp_gsfcgce.o gsfcgce_driver.o
	$(FC) $(FFLAGS) $^ -o gsfcgce_fortran

gsfcgce_fortran_opt: module_model_constants.o module_mp_gsfcgce_opt.o gsfcgce_driver_opt.o
	$(FC) $(FFLAGS) $^ -o gsfcgce_fortran_opt

gsfcgce_driver.o: gsfcgce_driver.f90
	$(FC) -c $^ $(FFLAGS) -o $@

module_mp_gsfcgce.o: module_mp_gsfcgce.f90
	$(FC) -c $^ $(FFLAGS)

gsfcgce_driver_opt.o: gsfcgce_driver_opt.f90
	$(FC) -c $^ $(FFLAGS) -o $@

module_mp_gsfcgce_opt.o: module_mp_gsfcgce_opt.f90
	$(FC) -c $^ $(FFLAGS)

module_model_constants.o: module_model_constants.f90
	$(FC) -c $^ $(FFLAGS)


gsfcgce_fortran_MIC: module_model_constants_MIC.o module_mp_gsfcgce_MIC.o gsfcgce_driver_MIC.o
	$(FC_MIC) $(FFLAGS_MIC) $^ -o gsfcgce_fortran_MIC

gsfcgce_fortran_opt_MIC: module_model_constants_MIC.o module_mp_gsfcgce_opt_MIC.o gsfcgce_driver_opt_MIC.o
	$(FC_MIC) $(FFLAGS_MIC) $^ -o gsfcgce_fortran_opt_MIC

gsfcgce_driver_MIC.o: gsfcgce_driver.f90
	$(FC_MIC) -c $^ $(FFLAGS_MIC) -o $@

module_mp_gsfcgce_MIC.o: module_mp_gsfcgce.f90
	$(FC_MIC) -c $^ $(FFLAGS_MIC) -o $@

gsfcgce_driver_opt_MIC.o: gsfcgce_driver_opt.f90
	$(FC_MIC) -c $^ $(FFLAGS_MIC) -o $@

module_mp_gsfcgce_opt_MIC.o: module_mp_gsfcgce_opt.f90
	$(FC_MIC) -c $^ $(FFLAGS_MIC) -o $@

module_model_constants_MIC.o: module_model_constants.f90
	$(FC_MIC) -c $^ $(FFLAGS_MIC) -o $@

clean:
	rm -fr *.mod *.o *genmod* *.optrpt gsfcgce_fortran gsfcgce_fortran_opt gsfcgce_fortran_MIC gsfcgce_fortran_opt_MIC

