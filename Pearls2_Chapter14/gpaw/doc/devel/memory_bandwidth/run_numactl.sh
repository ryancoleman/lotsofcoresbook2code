#!/bin/sh

# Execute the GPAW code as one or more independent tasks

if test -z $NCORES;
then
   export NCORES=8
fi

if test -z $STARTCORES;
then
   export STARTCORES=0
fi

if test -z $CORES_PER_SOCKET;
then
   export CORES_PER_SOCKET=4
fi

if test -z $MACHINE;
then
   export MACHINE=TEST
fi

# export PYTHONPATH=~/gpaw.mkl:${PYTHONPATH}
export script=../../H2Al110.py

## CONFIGURE one of the following:

if [ -f /home/camp/modulefiles.sh ]; then
    . /home/camp/modulefiles.sh
    module load openmpi
fi

# Using the GCC compiler and the AMD ACML library
#. /usr/local/openmpi-1.2.5-gfortran/bin/mpivars-1.2.5.sh
#export LD_LIBRARY_PATH=/opt/acml-4.0.1/gfortran64/lib:${LD_LIBRARY_PATH}

# Using the Intel compiler and the Intel MKL library
#. /usr/local/openmpi-1.2.7.intel/bin/mpivars-1.2.7.sh
#. /opt/intel/cce/10.1.018/bin/iccvars.sh
#. /opt/intel/mkl/10.0.4.023/tools/environment/mklvarsem64t.sh
#. /opt/intel/fce/10.1.018/bin/ifortvars.sh
export OMP_NUM_THREADS=1

index=${STARTCORES}
while [ "$index" -le "$NCORES" ];
do
  if [ "$index" -eq 0 ];
      then
      p=1
  else
      p=$index
  fi
  #
  echo Benchmark for ${p} tasks started at `date`
  ( cd ${MACHINE}_py${p}_01; time mpiexec -np $p ../../taskit.BINDING.one.node 0 ${CORES_PER_SOCKET} python $script --runs=5 > out.txt )
  echo
  index=`expr $index + 2`
done
echo Benchmark runs ended at `date`
