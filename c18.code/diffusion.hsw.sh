#!/bin/env bash

unset OMP_PROC_BIND
unset OMP_PLACES
unset KMP_AFFINITY

echo $LD_LIBRARY_PATH

appdir=host
app=host/diffusion

nx_list=(240 256 480 512)
count=(600 600 200 200)

#number of cores on HSW
ncores=28

i=0
for nx in "${nx_list[@]}"
do 

  #options passed to the executable
  load="nx=$nx count=${count[$i]} nf=10 out="

  export OMP_NESTED=true

  export KMP_HOT_TEAMS_MAX_LEVEL=2
  export KMP_HOT_TEAMS_MODE=1

  #Use nCores*nHT for mic, twoyz
  export OMP_NUM_THREADS=$ncores
  export OMP_PLACES=cores
  export OMP_PROC_BIND=close

  echo '=============================================='
  echo 'Tiled/Blocked/Collapse'
  $app mic $load
  rm $appdir/temp.tmp

  echo '=============================================='
  echo 'Explicity 2D partition YZ'
  $app twoyz $load
  rm $appdir/temp.tmp

  echo '=============================================='
  echo 'Explicity 2D partition YZ'
  $app twoyz $load np=4
  rm $appdir/temp.tmp

  #HT1 using two threads in the nested loop
  let nthreads=$ncores/2
  export OMP_PLACES=cores
  export OMP_NUM_THREADS=$nthreads,2
  export OMP_PROC_BIND=spread,close

  echo '=============================================='
  echo 'Nested OpenMP with a load-balance scheme'
  echo $OMP_NUM_THREADS,$OMP_PROC_BIND
  $app nested_lb $load
  rm $appdir/temp.tmp

  #HT2 using two threads in the nested loop
  export OMP_PLACES=threads
  export OMP_NUM_THREADS=$ncores,2
  export OMP_PROC_BIND=spread,close

  echo '=============================================='
  echo 'Nested OpenMP with a load-balance scheme'
  echo $OMP_NUM_THREADS,$OMP_PROC_BIND
  $app nested_lb $load
  rm $appdir/temp.tmp


  i=$((i+1))
done

