#!/bin/bash --login

#BSUB -R "@@nmpi@@*{select[@@arch@@] span[ptile=@@ppn@@]}"
#BSUB -q @@queue@@
#BSUB -x
#BSUB -J @@shortname@@
#BSUB -o /home/mklemm/projects/gpaw/results/@@fullname@@.log

# get some more information about the system
echo ==========================================================================
which icpc
echo ==========================================================================
which mpirun
echo ==========================================================================
cpuinfo
echo ==========================================================================
cat /proc/cpuinfo
echo ==========================================================================
env
echo ==========================================================================
if [ "@@offload@@x" == "yesx" ]; then
    micinfo
else
    echo "No offload, Xeon only"
fi
echo ==========================================================================

echo Sourcing pymic and GPAW environment
export PYMIC_BASE=/home/mklemm/projects/gpaw/pymic
export MIC_ENV_PREFIX=PYMIC
export PYMIC_LD_LIBRARY_PATH=$MIC_LD_LIBRARY_PATH
export PYMIC_LIBRARY_PATH=$PYMIC_BASE/examples/double_it/:$PYMIC_BASE/examples/dgemm/:$PYMIC_BASE/examples/svd:$PYMIC_BASE/pymic:$PYMIC_BASE/benchmarks:$PYMIC_BASE/tests
export MIC_LD_LIBRARY_PATH=$MIC_LD_LIBRARY_PATH:$PYMIC_BASE/pymic
export PYTHONPATH=$PYTHONPATH:$PYMIC_BASE
export GPAW_BASE=/home/mklemm/projects/gpaw/gpaw
export LIBXCDIR=$GPAW_BASE/libxc-2.0.2/install
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GPAW_BASE/libxc-2.0.2/install/lib
export PYTHONPATH=$GPAW_BASE/ase-3.9.0-trunk:$PYTHONPATH
export PATH=$GPAW_BASE/ase-3.9.0-trunk/tools:$PATH
export GPAW_SETUP_PATH=$GPAW_BASE/gpaw-setups-0.9.9672
export PYTHONPATH=$GPAW_BASE:$PYTHONPATH
export PATH=$GPAW_BASE/build/bin.linux-x86_64-intel64-2.6:$GPAW_BASE/tools/:$PATH
export PYMIC_LD_LIBRARY_PATH=$PYMIC_LD_LIBRARY_PATH
export PYMIC_LIBRARY_PATH=$PYMIC_LIBRARY_PATH:$GPAW_BASE/gpaw/mic
echo ==========================================================================

# first test the pymic installation on the master node
echo "Testing pymic to make sure GPAW run correctly even when only using host"
cd $PYMIC_BASE
make tests
echo ==========================================================================

# setup the run directory
rundir=/lfs/lfs08/mklemm/gpaw/@@fullname@@-$LSB_JOBID
mkdir -p $rundir
time cp -av /home/mklemm/projects/gpaw/gpaw/benchmarks/* $rundir
cd $rundir
echo ==========================================================================

# launch the job
echo "Launching GPAW with @@benchmark@@.py"
export I_MPI_HYDRA_BOOTSTRAP=ssh
export I_MPI_DEBUG=5
export I_MPI_PIN=1
export I_MPI_DAPL_PROVIDER=ofa-v2-mlx4_0-1
export I_MPI_FABRICS=shm:dapl
if [ "@@debug@@x" == "yesx" ]; then
    export PYMIC_DEBUG=5
fi
gpaw-python @@benchmark@@.py --dry-run
cat $rundir/*@@benchmark@@*.txt
echo ==========================================================================
if [ "@@offload@@x" == "yesx" ]; then
    export GPAW_OFFLOAD=1
    time mpirun -prepend-rank -np @@nmpi@@ ./affinity_wrapper.sh @@ppn@@ gpaw-python @@benchmark@@.py
else
    export GPAW_OFFLOAD=0
    time mpirun -prepend-rank-np @@nmpi@@ gpaw-python @@benchmark@@.py
fi
echo ==========================================================================

# dump the output file of the job
cat $rundir/*@@benchmark@@*.txt
echo ==========================================================================
