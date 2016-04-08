ulimit -s unlimited
export OMP_STACKSIZE=8m

export NUM_TILES_X=2
export NUM_TILES_Y=16

./gsfcgce_fortran
./gsfcgce_fortran_opt
