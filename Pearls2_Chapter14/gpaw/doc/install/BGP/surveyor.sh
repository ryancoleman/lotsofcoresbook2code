type=Au_bulk3x3x3
cwd=`pwd`
acct=Gpaw
queue=default
time=60
nodes=64
mode=vn
mapping=ZYXT
# uncomment below for mapfile
# mapfile=BGMAP_128_4x4x4x8
# mapping=$mapfile
job=${type}_${nodes}_${mode}_${mapping}
input=${type}.py
tau=wrapper.py
scratch=/pvfs-surveyor/${USER}
install=/gpfs/home/${USER}
setups=/soft/apps/gpaw-setups-0.6.6300
aseversion=
gpawversion=
sharelibs=/lib:/bgsys/drivers/ppcfloor/comm/default/lib:/bgsys/drivers/ppcfloor/comm/sys/lib:/bgsys/drivers/ppcfloor/runtime/SPI
taulibs=/soft/apps/tau/tau-2.19.2/bgp/lib/bindings-bgptimers-gnu-mpi-python-pdt
bin=gpaw-python

rm -rf $scratch/$job
mkdir $scratch/$job
cp $input $scratch/$job
# copy mapfile
# cp $mapfile $scratch/$job 
# copy tau wrapper script
# cp $tau $scratch/$job
cd $scratch/$job

# pristine
qsub -A $acct -n $nodes -t $time -q $queue --mode $mode --env DCMF_EAGER=8388608:BG_MAPPING=$mapping:MPIRUN_ENABLE_TTY_REPORTING=0:OMP_NUM_THREADS=1:GPAW_SETUP_PATH=$setups:PYTHONPATH=${install}/gpaw${gpawversion}:${install}/ase${aseversion}:${taulibs}:$PYTHONPATH:LD_LIBRARY_PATH=${taulibs}:${sharelibs} ${install}/gpaw${gpawversion}/build/bin.linux-ppc64-2.6/${bin} ${input} --domain-decomposition=4,4,4 --state-parallelization=2 --sl_default=2,2,64

# TAU automatic profiling
# qsub -A $acct -n $nodes -t $time -q $queue --mode $mode --env TAU_VERBOSE=1:TAU_CALLPATH=0:TAU_CALLPATH_DEPTH=10:TAU_COMM_MATRIX=0:TAU_TRACK_MESSAGE=0:TAU_THROTTLE=0:TAU_COMPENSATE=1:TAU_METRICS=BGPTIMERS:DCMF_EAGER=8388608:BG_MAPPING=$mapping:MPIRUN_ENABLE_TTY_REPORTING=0:OMP_NUM_THREADS=1:GPAW_SETUP_PATH=$setups:PYTHONPATH=${install}/gpaw${gpawversion}:${install}/ase${aseversion}:${taulibs}:$PYTHONPATH:LD_LIBRARY_PATH=${taulibs}:${sharelibs} ${install}/gpaw${gpawversion}/build/bin.linux-ppc64-2.6/${bin} ${tau} --domain-decomposition=4,4,4 --state-parallelization=2 --sl_default=2,2,64