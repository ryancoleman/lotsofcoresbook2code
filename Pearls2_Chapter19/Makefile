all: thirdparty/ham run

thirdparty/ham: 
	echo "Making thirdparty"
	$(MAKE) -C thirdparty all

bin/benchmark_omp:
	sh make_omp.sh -c
bin/benchmark_ocl_cpu:
	sh make_ocl.sh -c
bin/benchmark_ocl_mic:
	sh make_ocl.sh -a
bin.mic/benchmark_omp:
	sh make_omp.sh -a

run: bin/benchmark_omp bin/benchmark_ocl_cpu bin/benchmark_ocl_mic bin.mic/benchmark_omp
	bin/benchmark_ocl_cpu
	bin/benchmark_ocl_mic
	./run_mic.sh bin.mic/benchmark_omp
	./bin/benchmark_omp

