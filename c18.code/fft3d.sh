app=./fft3d

export OMP_NESTED=true
export MKL_DYNAMIC=false
export MKL_NUM_THREADS=16
unset OMP_PROC_BIND
unset OMP_PLACES

export KMP_HOT_TEAMS_MAX_LEVEL=2
export KMP_HOT_TEAMS_MODE=1

export OMP_PROC_BIND=spread,close
unset KMP_AFFINITY

nfft=1
np_max=240

ng=64
nsteps=32

omp=1
for nmkl in 4 8 16 20 40;
do
  let np=$np_max/$nmkl/$omp
  let howmany=$nfft
  export OMP_NUM_THREADS=$omp,$nmkl
  export MKL_NUM_THREADS=$nmkl
  let niters=$nsteps
  mpirun -np $np $app -i $niters -g $ng -s $howmany
done

omp=2
for nmkl in 4 8 12 30;
do
  let np=$np_max/$nmkl/$omp
  let howmany=$nfft
  export OMP_NUM_THREADS=$omp,$nmkl
  export MKL_NUM_THREADS=$nmkl
  let niters=$nsteps
  mpirun -np $np $app -i $niters -g $ng -s $howmany
done

omp=3
for nmkl in 4 8 16 20 40;
do
  let np=$np_max/$nmkl/$omp
  let howmany=$nfft
  export OMP_NUM_THREADS=$omp,$nmkl
  export MKL_NUM_THREADS=$nmkl
  let niters=$nsteps
  mpirun -np $np $app -i $niters -g $ng -s $howmany
done

ng=128
nsteps=10

omp=1
for nmkl in 4 8 16 20 40;
do
  let np=$np_max/$nmkl/$omp
  let howmany=$nfft
  export OMP_NUM_THREADS=$omp,$nmkl
  export MKL_NUM_THREADS=$nmkl
  let niters=$nsteps
  mpirun -np $np $app -i $niters -g $ng -s $howmany
done

omp=2
for nmkl in 4 8 12 30;
do
  let np=$np_max/$nmkl/$omp
  let howmany=$nfft
  export OMP_NUM_THREADS=$omp,$nmkl
  export MKL_NUM_THREADS=$nmkl
  let niters=$nsteps
  mpirun -np $np $app -i $niters -g $ng -s $howmany
done

omp=3
for nmkl in 4 8 16 20 40;
do
  let np=$np_max/$nmkl/$omp
  let howmany=$nfft
  export OMP_NUM_THREADS=$omp,$nmkl
  export MKL_NUM_THREADS=$nmkl
  let niters=$nsteps
  mpirun -np $np $app -i $niters -g $ng -s $howmany
done

ng=100
nsteps=20

omp=1
for nmkl in 4 8 16 20 40;
do
  let np=$np_max/$nmkl/$omp
  let howmany=$nfft
  export OMP_NUM_THREADS=$omp,$nmkl
  export MKL_NUM_THREADS=$nmkl
  let niters=$nsteps
  mpirun -np $np $app -i $niters -g $ng -s $howmany
done

omp=2
for nmkl in 4 8 12 30;
do
  let np=$np_max/$nmkl/$omp
  let howmany=$nfft
  export OMP_NUM_THREADS=$omp,$nmkl
  export MKL_NUM_THREADS=$nmkl
  let niters=$nsteps
  mpirun -np $np $app -i $niters -g $ng -s $howmany
done

omp=3
for nmkl in 4 8 16 20 40;
do
  let np=$np_max/$nmkl/$omp
  let howmany=$nfft
  export OMP_NUM_THREADS=$omp,$nmkl
  export MKL_NUM_THREADS=$nmkl
  let niters=$nsteps
  mpirun -np $np $app -i $niters -g $ng -s $howmany
done
