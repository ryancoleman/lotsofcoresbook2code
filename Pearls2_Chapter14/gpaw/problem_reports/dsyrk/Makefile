.PHONY: all clean realclean

all: offload_performance_dsyrk offload_performance_dsyrk.mic

offload_performance_dsyrk: offload_performance_dsyrk.f90
	ifort -align array64byte -O3 -mkl -o offload_performance_dsyrk offload_performance_dsyrk.f90
offload_performance_dsyrk.mic: offload_performance_dsyrk.f90
	ifort -align array64byte -O3 -mkl -qoffload=none -mmic -o offload_performance_dsyrk.mic offload_performance_dsyrk.f90

clean:
	rm -f offload_performance_dsyrk offload_performance_dsyrk.mic

realclean: clean
	rm -f *~
