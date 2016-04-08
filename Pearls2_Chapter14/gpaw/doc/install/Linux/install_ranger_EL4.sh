setenv APPS "/share/home/01067/tg803307/CAMd"
setenv MODULEFILES "${APPS}/modulefiles"

# build packages

set nose_version=0.11.3
tar zxf nose-${nose_version}.tar.gz
cd nose-${nose_version}
python setup.py install --root=${APPS}/nose-${nose_version}-1
cd ..

set numpy_version=1.5.0
tar zxf numpy-${numpy_version}.tar.gz
cd numpy-${numpy_version}
python setup.py install --root=${APPS}/numpy-${numpy_version}-1
cd ..

set ase_version=3.4.1.1765
tar zxf python-ase-${ase_version}.tar.gz

set gpaw_version=0.7.2.6974
tar zxf gpaw-${gpaw_version}.tar.gz

set gpaw_setups_version=0.6.6300
tar zxf gpaw-setups-${gpaw_setups_version}.tar.gz

mkdir -p ${MODULEFILES}/nose
cat <<EOF > ${MODULEFILES}/nose/${nose_version}-1
#%Module1.0
set apps_path ${APPS}
prereq python/2.5.2
prepend-path    PATH                \$apps_path/nose-${nose_version}-1/opt/apps/python/python-2.5.2/bin/
prepend-path    PYTHONPATH          \$apps_path/nose-${nose_version}-1/opt/apps/python/python-2.5.2/lib/python2.5/site-packages/
unset apps_path
EOF

mkdir -p ${MODULEFILES}/numpy
cat <<EOF > ${MODULEFILES}/numpy/${numpy_version}-1
#%Module1.0
set apps_path ${APPS}
prereq python/2.5.2
prereq nose/0.11.3-1
prepend-path    PATH                \$apps_path/numpy-${numpy_version}-1/opt/apps/python/python-2.5.2/bin
prepend-path    PYTHONPATH          \$apps_path/numpy-${numpy_version}-1/opt/apps/python/python-2.5.2/lib/python2.5/site-packages/
unset apps_path
EOF

mkdir -p ${MODULEFILES}/campos-ase3
cat <<EOF > ${MODULEFILES}/campos-ase3/${ase_version}-1
#%Module1.0
set apps_path ${APPS}
prereq python/2.5.2
prereq nose/0.11.3-1
prereq numpy/1.5.0-1
prepend-path    PATH                \$apps_path/python-ase-${ase_version}/tools
prepend-path    PYTHONPATH          \$apps_path/python-ase-${ase_version}/
unset apps_path
EOF

mkdir -p ${MODULEFILES}/campos-gpaw-setups
cat <<EOF > ${MODULEFILES}/campos-gpaw-setups/${gpaw_setups_version}-1
#%Module1.0
set apps_path ${APPS}
prepend-path    GPAW_SETUP_PATH     \$apps_path/gpaw-setups-${gpaw_setups_version}
unset apps_path
EOF

mkdir -p ${MODULEFILES}/campos-gpaw
cat <<EOF > ${MODULEFILES}/campos-gpaw/${gpaw_version}-1
#%Module1.0
set apps_path ${APPS}
prereq python/2.5.2
prereq nose/0.11.3-1
prereq numpy/1.5.0-1
prereq campos-ase3
prereq campos-gpaw-setups
prereq gcc/4.4.5
prereq openmpi/1.3b
prereq mkl/10.0
prepend-path    PATH                \$apps_path/gpaw-${gpaw_version}/tools
prepend-path    PATH                \$apps_path/gpaw-${gpaw_version}/build/bin.linux-x86_64-2.5/
prepend-path    PYTHONPATH          \$apps_path/gpaw-${gpaw_version}/
prepend-path    PYTHONPATH          \$apps_path/gpaw-${gpaw_version}/build/lib.linux-x86_64-2.5/
setenv OMP_NUM_THREADS 1
unset apps_path
EOF

module avail

module use --append ${MODULEFILES}
module load python/2.5.2
module load nose/0.11.3-1
module load numpy/1.5.0-1
module load campos-ase3
# test numpy
python -c "import numpy; numpy.test()"
# test ase
mkdir -p testase
cd testase
testase.py --no-display >& testase.log
cd ..
# build gpaw
cd gpaw-${gpaw_version}
python setup.py build_ext --customize=../customize_ranger_EL4.py --remove-default-flags >& build_ext.log
cd ..
module load campos-gpaw-setups
module load gcc/4.4.5
module load openmpi/1.3b
module load mkl/10.0
module load campos-gpaw
mkdir -p testgpaw
cd testgpaw
mpiexec -np 4 gpaw-python `which gpaw-test` >& tee testgpaw.log
