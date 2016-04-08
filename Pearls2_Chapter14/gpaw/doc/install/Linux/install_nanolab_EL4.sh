export APPS="/home/karsten/CAMd"
export MODULEFILES="${APPS}/modulefiles"

# build packages

numpy_version=1.1.1
tar zxf numpy-${numpy_version}.tar.gz
cd  numpy-${numpy_version}
# disable compiling with atlas
sed -i "s/_lib_atlas =.*/_lib_atlas = ['ignore_atlas']/g" numpy/distutils/system_info.py
python setup.py install --root=${APPS}/numpy-${numpy_version}-1
cd ..

ase_version=3.4.1.1765
tar zxf python-ase-${ase_version}.tar.gz
# patch cPickle
sed -i "s/cPickle as //g" python-ase-${ase_version}/ase/io/trajectory.py

gpaw_version=0.7.2.6974
tar zxf gpaw-${gpaw_version}.tar.gz

gpaw_setups_version=0.6.6300
tar zxf gpaw-setups-${gpaw_setups_version}.tar.gz

mkdir -p ${MODULEFILES}/numpy
cat <<EOF > ${MODULEFILES}/numpy/${numpy_version}-1
#%Module1.0
set apps_path ${APPS}
prepend-path    PATH                \$apps_path/numpy-${numpy_version}-1/usr/bin
prepend-path    PYTHONPATH          \$apps_path/numpy-${numpy_version}-1/usr/lib/python2.3/site-packages/
unset apps_path
EOF

mkdir -p ${MODULEFILES}/campos-ase3
cat <<EOF > ${MODULEFILES}/campos-ase3/${ase_version}-1
#%Module1.0
set apps_path ${APPS}
prereq numpy
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
prereq numpy
prereq campos-ase3
prereq campos-gpaw-setups
prereq intel_compilers/11.1
prereq openmpi/1.3.3
prepend-path    PATH                \$apps_path/gpaw-${gpaw_version}/tools
prepend-path    PATH                \$apps_path/gpaw-${gpaw_version}/build/bin.linux-i686-2.3
prepend-path    PYTHONPATH          \$apps_path/gpaw-${gpaw_version}/
prepend-path    PYTHONPATH          \$apps_path/gpaw-${gpaw_version}/build/lib.linux-i686-2.3
setenv OMP_NUM_THREADS 1
unset apps_path
EOF

module avail

module use --append ${MODULEFILES}
module load numpy
module load campos-ase3
# test numpy
python -c "import numpy; numpy.test()"
# test ase
mkdir -p testase
cd testase
testase.py --no-display 2>&1 | tee testase.log
cd ..
# build gpaw
cd gpaw-${gpaw_version}
python setup.py build_ext --customize=../customize_nanolab_EL4.py --remove-default-flags
cd ..
module load campos-gpaw-setups
module load intel_compilers/11.1
module load openmpi/1.3.3
module load campos-gpaw
mkdir -p testgpaw
cd testgpaw
mpiexec -np 4 gpaw-python `which gpaw-test` 2>&1 | tee testgpaw.log
