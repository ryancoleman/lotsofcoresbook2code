export NUM_TILES_X=4; export NUM_TILES_Y=60; export OMP_STACKSIZE=8m; export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME; ulimit -s unlimited; ./gsfcgce_fortran_MIC; ./gsfcgce_fortran_opt_MIC
