#!/bin/sh
###These lines are for Moab
#MSUB -l nodes=1:ppn=8
#MSUB -l walltime=02:00:00
#MSUB -l partition=nehalem|harpertown|westmere

. /opt/gridware/modules-3.2.7/modules.sh
export APPS=`echo ~/CAMd`
export MODULEFILES="${APPS}/modulefiles"
module use --append ${MODULEFILES}
module load python
module load nose
module load numpy
module load python-ase
module load gpaw-setups
module load gpaw/0.9.0.8965-1

export OMP_NUM_THREADS=1

##### Running commands
cat $PBS_NODEFILE
NP=`cat $PBS_NODEFILE | wc -l`
mpirun -x PYTHONPATH -x GPAW_SETUP_PATH -x PATH -np $NP -machinefile $PBS_NODEFILE gpaw-python `which gpaw-test`
