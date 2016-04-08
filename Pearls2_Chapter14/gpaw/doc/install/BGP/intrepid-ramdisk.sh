type=Au_bulk6x6x6
cwd=`pwd`
acct=Gpaw
queue=prod
time=45
nodes=8192
mode=vn
# uncomment below for mapfile
# mapfile=BGMAP_8x8x4x16 # 1024 nodes
# mapping=$mapfile
mapping=ZYXT 
input=${type}.py
scratch=/intrepid-fs0/users/${USER}/persistent
install=/gpfs/home/${USER}
setups=/soft/apps/gpaw-setups-0.6.6300
gpawversion=7800
sharelibs=/lib:/gpaw/V1R4M2
job=${type}_${nodes}_${mode}_${mapping}_r${gpawversion}
taulibs=/gpaw/tau-2.19.2
bin=gpaw-python
kernel=gpaw 

rm -rf $scratch/$job
mkdir $scratch/$job
cp $input $scratch/$job
cp $mapfile $scratch/$job 
# cp $pos $scratch/$job
cd $scratch/$job

# Many of the values for the DMCF variables below are defaults
# pristine
qsub --kernel $kernel -A $acct -n $nodes -t $time -q $queue --mode $mode --env DCMF_EAGER=8388608:BG_MAPPING=$mapping:MPIRUN_ENABLE_TTY_REPORTING=0:OMP_NUM_THREADS=1:GPAW_SETUP_PATH=$setups:PYTHONHOME=/gpaw:LD_LIBRARY_PATH=$sharelibs ${install}/gpaw-r${gpawversion}/build/bin.linux-ppc64-2.6/${bin} ${type}.py --domain-decomposition=8,8,4 --state-parallelization=16

# TAU manual instrumentation
# qsub --kernel $kernel -A $acct -n $nodes -t $time -q $queue --mode $mode --env TAU_VERBOSE=1:TAU_CALLPATH=0:TAU_CALLPATH_DEPTH=10:TAU_COMM_MATRIX=0:TAU_TRACK_MESSAGE=0:TAU_THROTTLE=0:TAU_COMPENSATE=1:TAU_METRICS=BGPTIMERS:DCMF_EAGER=8388608:BG_MAPPING=$mapping:MPIRUN_ENABLE_TTY_REPORTING=0:OMP_NUM_THREADS=1:GPAW_SETUP_PATH=$setups:PYTHONHOME=/gpaw:PYTHONPATH=${taulibs}:LD_LIBRARY_PATH=$sharelibs:${taulibs} ${install}/gpaw-r${gpawversion}/build/bin.linux-ppc64-2.6/${bin} ${type}.py --domain-decomposition=8,8,16 --state-parallelization=32
