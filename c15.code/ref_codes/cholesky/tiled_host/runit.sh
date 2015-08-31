# FILE LICENSE TAG: SAMPLE

cd ../../../bin/host

#scatter affinity set for host
export KMP_AFFINITY=scatter

./cholesky_tiled_host -m 4800 -t 6 -l col -i 5
./cholesky_tiled_host -m 4800 -t 6 -l row -i 5
