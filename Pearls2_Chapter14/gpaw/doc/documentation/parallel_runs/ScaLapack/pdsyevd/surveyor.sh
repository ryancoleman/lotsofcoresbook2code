qsub -A Gpaw -n 1 -t 10 --mode vn --env LD_LIBRARY_PATH="/bgsys/drivers/ppcfloor/gnu-linux/powerpc-bgp-linux/lib:$LD_LIBRARY_PATH" test.exe
