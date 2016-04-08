all:
	make -C matMult
	make -C cholesky/tiled_host
	make -C cholesky/tiled_hstreams
	make -C lu/tiled_host
	make -C lu/tiled_hstreams
	cd matMult; echo "Running matMult"; ./runit.sh;
	cd cholesky/tiled_hstreams; echo "Running cholesky tiled_hstreams"; ./runit.sh;
	cd lu/tiled_hstreams; echo "Running lu tiled_hstreams"; ./runit.sh;

clean:
	make clean -C matMult
	make clean -C cholesky/tiled_host
	make clean -C cholesky/tiled_hstreams
	make clean -C lu/tiled_host
	make clean -C lu/tiled_hstreams
