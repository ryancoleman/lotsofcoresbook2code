#!/bin/env bash
export OMP_NESTED=true
export MKL_DYNAMIC=false
unset OMP_PROC_BIND
unset OMP_PLACES
export OMP_STACKSIZE=16M

export KMP_HOT_TEAMS_MAX_LEVEL=2
#unset KMP_HOT_TEAMS_MAX_LEVEL
#unset KMP_HOT_TEAMS_MODE
export KMP_HOT_TEAMS_MODE=1

unset KMP_AFFINITY

app=./dgemm

title=gemm.p7120
np_max=240
for M in 128 256 512 1024;
do
  mkl=4
  export OMP_PROC_BIND=spread,close
  export MKL_NUM_THREADS=$mkl
  mtag=th$MKL_NUM_THREADS
  output=$title.$M.4ht.$mtag.out

  env  > $output
  for nth in 1 2 4 12 20 30 60;
  do
    let mpi=$np_max/$nth/$mkl
    export OMP_NUM_THREADS=$nth,$MKL_NUM_THREADS
    echo $M,$mpi,$nth,$MKL_NUM_THREADS
    mpirun -np $mpi $app $M &>> $output
    sleep 1
  done

 mkl=3                                                                 
  export OMP_PLACES=threads
  export OMP_PROC_BIND=spread,spread
  export MKL_NUM_THREADS=$mkl                                           
  mtag=th$MKL_NUM_THREADS
  output=$title.$M.3ht.$mtag.out                                        
  env  > $output                                                        
  for nth in 1 2 4 12 20 30 60;                                         
  do                                                                    
    let mpi=180/$nth/$mkl                                               
    export OMP_NUM_THREADS=$nth,$MKL_NUM_THREADS                        
    echo $M,$mpi,$nth,$MKL_NUM_THREADS                                  
    mpirun -np $mpi $app $M &>> $output                                 
    sleep 1                                                             
  done          

  mkl=2                               
  export OMP_PROC_BIND=spread,spread
  export MKL_NUM_THREADS=$mkl      
  mtag=th$MKL_NUM_THREADS
  output=$title.$M.2ht.$mtag.out     
  env  > $output                  
  for nth in 1 2 4 12 20 30 60; 
  do                            
    let mpi=120/$nth/$mkl   
    export OMP_NUM_THREADS=$nth,$MKL_NUM_THREADS
    echo $M,$mpi,$nth,$MKL_NUM_THREADS          
    mpirun -np $mpi $app $M &>> $output
    sleep 1                                                             
  done                      

  mkl=1                               
  export OMP_PROC_BIND=spread,spread
  export MKL_NUM_THREADS=$mkl      
  mtag=th$MKL_NUM_THREADS
  output=$title.$M.1ht.$mtag.out     
  env  > $output                  
  for nth in 1 2 4 12 20 30 60; 
  do                            
    let mpi=60/$nth/$mkl   
    export OMP_NUM_THREADS=$nth,$MKL_NUM_THREADS
    echo $M,$mpi,$nth,$MKL_NUM_THREADS          
    mpirun -np $mpi $app $M &>> $output
    sleep 1                                                             
  done                      
done

#for M in 128 256 512 1024 2048;
#do
#  mkl=8
#  export OMP_PROC_BIND=spread,close
#  export MKL_NUM_THREADS=$mkl
#  mtag=th$MKL_NUM_THREADS.15.1.133
#  output=$title.4ht.$M.$mtag.out
#
#  env  > $output
#  #for nth in 1 2 6 10 15 30;
#  for nth in 1 2 6 10 15 30;
#  do
#    let mpi=$np_max/$nth/$mkl
#    export OMP_NUM_THREADS=$nth,$MKL_NUM_THREADS
#    echo $M,$mpi,$nth,$MKL_NUM_THREADS
#    mpirun -np $mpi $app $M &>> $output
#    sleep 1
#  done
#
#  mkl=16
#  export OMP_PROC_BIND=spread,close
#  export MKL_NUM_THREADS=$mkl
#  mtag=th$MKL_NUM_THREADS.15.1.133
#  output=$title.4ht.$M.$mtag.out
#
#  env  > $output
#  for nth in 1 3 5 15;
#  do
#    let mpi=$np_max/$nth/$mkl
#    export OMP_NUM_THREADS=$nth,$MKL_NUM_THREADS
#    echo $M,$mpi,$nth,$MKL_NUM_THREADS
#    mpirun -np $mpi $app $M &>> $output
#    sleep 1
#  done
#done

