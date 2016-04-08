include make.inc

.PHONY: part-1
.PHONY: part-1-no-offload
.PHONY: part-2
.PHONY: clean

all:
	@$(MAKE) part-1
	@$(MAKE) part-2

part-1:
	@$(MAKE) -C part-1 veryclean
	@$(MAKE) -C part-1 version=1 offload=1
	@$(MAKE) -C part-1 clean
	@$(MAKE) -C part-1 version=2 offload=1
	@$(MAKE) -C part-1 clean
	@$(MAKE) -C part-1 version=3 offload=1
	@$(MAKE) -C part-1 clean
	@$(MAKE) -C part-1 version=4 offload=1
	@$(MAKE) -C part-1 clean
	@$(MAKE) -C part-1 version=5 offload=1
	@$(MAKE) -C part-1 clean
	@$(MAKE) -C part-1 version=6 offload=1
	@$(MAKE) -C part-1 clean
	@$(MAKE) -C part-1 version=7 offload=1
	@$(MAKE) -C part-1 clean
	@$(MAKE) -C part-1 version=8 offload=1
	@$(MAKE) -C part-1 clean
	@$(MAKE) -C part-1 version=9 offload=1
	@$(MAKE) -C part-1 clean
	@$(MAKE) -C part-1 version=10 offload=1
	@$(MAKE) -C part-1 clean

part-1-no-offload:
	@$(MAKE) -C part-1 veryclean
	@$(MAKE) -C part-1 version=1 offload=0
	@$(MAKE) -C part-1 clean
	@$(MAKE) -C part-1 version=2 offload=0
	@$(MAKE) -C part-1 clean
	@$(MAKE) -C part-1 version=3 offload=0
	@$(MAKE) -C part-1 clean
	@$(MAKE) -C part-1 version=4 offload=0
	@$(MAKE) -C part-1 clean
	@$(MAKE) -C part-1 version=5 offload=0
	@$(MAKE) -C part-1 clean
	@$(MAKE) -C part-1 version=6 offload=0
	@$(MAKE) -C part-1 clean
	@$(MAKE) -C part-1 version=7 offload=0
	@$(MAKE) -C part-1 clean
	@$(MAKE) -C part-1 version=8 offload=0
	@$(MAKE) -C part-1 clean
	@$(MAKE) -C part-1 version=9 offload=0
	@$(MAKE) -C part-1 clean
	@$(MAKE) -C part-1 version=10 offload=0
	@$(MAKE) -C part-1 clean

part-2:
	@$(MAKE) -C part-2 veryclean
	@$(MAKE) -C part-2 nesting=manual
	@$(MAKE) -C part-2 clean
	@$(MAKE) -C part-2 nesting=omp
	@$(MAKE) -C part-2 clean
	@$(MAKE) -C part-2 nesting=teams
	@$(MAKE) -C part-2 clean

clean:
	@$(MAKE) -C part-1 clean
	@$(MAKE) -C part-2 clean

veryclean:
	@$(MAKE) -C part-1 veryclean
	@$(MAKE) -C part-2 veryclean
#part-2:
#	make -C part-2
