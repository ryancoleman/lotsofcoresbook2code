#!/bin/bash

# get some information about the job
ppn=$1
shift
rank=$PMI_RANK
nmpi=$PMI_SIZE

# number of devices in the system
ndev=2

# number of cores per device
nphcores=61
nphcores=$((nphcores - 1))

# number of threads per physical core
tpc=4

# ranks per device
rpd=$((ppn / ndev))
if [ "$rpd" == "0" ]; then
    rpd=1
fi

# physical cores per device
ncores=$((nphcores / rpd))

# partition number of the current rank on its device
partition=$((rank % rpd))

# offset for the current rank
offset=$((ncores * partition))

# build core selection string
select="${ncores}c,${tpc}t,${offset}o"

# fire up the actual run
log="affinity-`printf %03d $rank`.log"
rm -f $log
echo "host `hostname` rank `printf %03d $rank` - $select " |& tee -a $log
env | grep PYMIC |& tee -a $log
PYMIC_KMP_AFFINITY=compact,verbose PYMIC_KMP_PLACE_THREADS=$select $@ |& tee -a $log
