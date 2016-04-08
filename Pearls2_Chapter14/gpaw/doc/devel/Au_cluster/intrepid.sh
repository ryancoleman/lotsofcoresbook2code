type=Au_cluster
cwd=`pwd`
acct=Gpaw
queue=prod
time=90
nodes=512
mode=vn
mapping=ZYXT
# mapfile=BGMAP_128_4x4x4x8
# mapping=$mapfile
job=${type}_${nodes}_${mode}_${mapping}
input=${type}.py
pos=Au102_revised.xyz
scratch=/intrepid-fs0/users/${USER}/persistent
install=/soft/apps

rm -rf $scratch/$job
mkdir $scratch/$job
cp $input $scratch/$job
# cp $mapfile $scratch/$job
cp $pos $scratch/$job
cd $scratch/$job

qsub -A $acct -n $nodes -t $time -q $queue --mode $mode --env BG_MAPPING=$mapping:MPIRUN_ENABLE_TTY_REPORTING=0:OMP_NUM_THREADS=1:GPAW_SETUP_PATH=$GPAW_SETUP_PATH:PYTHONPATH=${install}/gpaw-r6000:${install}/ase-r1438:$PYTHONPATH:LD_LIBRARY_PATH=$CN_LD_LIBRARY_PATH ${install}/gpaw-r6000/build/bin.linux-ppc64-2.6/gpaw-python ${type}.py --domain-decomposition=8,8,8 --state-parallelization=4 --sl_diagonalize=5,5,64
