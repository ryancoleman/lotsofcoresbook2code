type=b256H2O
cwd=`pwd`
acct=Gpaw
queue=default
time=60
nodes=64
mode=smp
mapping=ZYXT
# mapfile=BGMAP_128_4x4x4x8
# mapping=$mapfile
job=${type}_${nodes}_${mode}_${mapping}
input=${type}.py
# pos=Au102_revised.xyz
scratch=/pvfs-surveyor/${USER}
install=/soft/apps

rm -rf $scratch/$job
mkdir $scratch/$job
cp $input $scratch/$job
# cp $mapfile $scratch/$job
# cp $pos $scratch/$job
cd $scratch/$job

qsub -A $acct -n $nodes -t $time -q $queue --mode $mode --env BG_MAPPING=$mapping:MPIRUN_ENABLE_TTY_REPORTING=0:OMP_NUM_THREADS=1:GPAW_SETUP_PATH=$GPAW_SETUP_PATH:PYTHONPATH=${install}/gpaw-r6000:${install}/ase-r1428:$PYTHONPATH:LD_LIBRARY_PATH=$CN_LD_LIBRARY_PATH ${install}/gpaw-r6000/build/bin.linux-ppc64-2.6/gpaw-python ${type}.py --domain-decomposition=4,4,4 --state-parallelization=1 --sl_diagonalize=4,4,64
