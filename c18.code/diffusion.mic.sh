#!/bin/env bash

unset OMP_PROC_BIND
unset OMP_PLACES
unset KMP_AFFINITY

echo $LD_LIBRARY_PATH

appdir=.
#app=../apps_xe15.2.164/diffusion-pearl
app=./diffusion

nx_list=(240 256 480 512)
count=(1200 1200 400 400)

ncores=60

export OMP_NESTED=true
export KMP_HOT_TEAMS_MAX_LEVEL=2
export KMP_HOT_TEAMS_MODE=1

i=0
for nx in "${nx_list[@]}"
do 

  #options passed to the executable
  load="nx=$nx count=${count[$i]} nf=10 out="

  #Use nCores*nHT for mic, twoyz
  let nthreads=4*$ncores
  export OMP_NUM_THREADS=$nthreads
  export OMP_PLACES=threads
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
  $app twoyz $load np=8
  rm $appdir/temp.tmp

  #Use NESTED OpenMP 4 for crew, nested_lib
  export OMP_PLACES=threads
  export OMP_NUM_THREADS=$ncores,4
  export OMP_PROC_BIND=spread,close

  echo '=============================================='
  echo 'Using Crew '
  $app crew $load
  rm $appdir/temp.tmp

  echo '=============================================='
  echo 'Nested OpenMP with a load-balance scheme'
  echo $OMP_NUM_THREADS,$OMP_PROC_BIND
  $app nested_lb $load
  rm $appdir/temp.tmp

  let nthreads=$ncores/2
  export OMP_NUM_THREADS=$nthreads,8
  export OMP_PROC_BIND=spread,close
  echo '=============================================='
  echo $OMP_NUM_THREADS,$OMP_PROC_BIND
  echo
  $app nested_lb np=8 $load
  rm $appdir/temp.tmp

  let nthreads=$ncores/2
  export OMP_NUM_THREADS=$nthreads,16
  export OMP_PROC_BIND=spread,close
  echo '=============================================='
  echo $OMP_NUM_THREADS,$OMP_PROC_BIND
  echo
  $app nested_lb np=16 $load
  rm $appdir/temp.tmp

  i=$((i+1))
done
