#/bin/sh
#MSUB -n 8 # number of tasks
#MSUB -T 4600 # time
#MSUB -q hybrid # use hybrid for GPU
#MSUB -A paxxxx

set -x
cd ${BRIDGE_MSUB_PWD}
module use ${HOME}/CAMd/modulefiles
module load nose
# blas/lapack/atlas missing
export LD_LIBRARY_PATH=${HOME}/CAMd/atlas:${LD_LIBRARY_PATH}
module load numpy
module load scipy
module load python-ase
module load gpaw-setups/0.8.7929-1
#module load gpaw-setups/0.9.9672-1
module load gpaw
#ccc_mprun gpaw-python gpaw-0.9.0.8965/gpaw/test/2Al.py
ccc_mprun gpaw-python `which gpaw-test`
