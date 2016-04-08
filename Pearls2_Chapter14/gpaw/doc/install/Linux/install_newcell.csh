setenv APPS "/afs/crc.nd.edu/user/j/jbray2/CAMd"
setenv MODULEFILES "${APPS}/modulefiles"

set nose_version=0.11.3
wget http://python-nose.googlecode.com/files/nose-${nose_version}.tar.gz
tar zxf nose-${nose_version}.tar.gz
cd nose-${nose_version}
python setup.py install --root=${APPS}/nose-${nose_version}-1
cd ..

set numpy_version=1.4.1
wget http://downloads.sourceforge.net/numpy/numpy-${numpy_version}.tar.gz
tar zxf numpy-${numpy_version}.tar.gz
cd  numpy-${numpy_version}
python setup.py install --root=${APPS}/numpy-${numpy_version}-1
cd ..

set ase_version=3.4.1.1765
wget https://wiki.fysik.dtu.dk/ase-files/python-ase-${ase_version}.tar.gz
tar zxf python-ase-${ase_version}.tar.gz

set gpaw_version=0.7.2.6974
wget https://wiki.fysik.dtu.dk/gpaw-files/gpaw-${gpaw_version}.tar.gz
tar zxf gpaw-${gpaw_version}.tar.gz

set gpaw_setups_version=0.6.6300
wget http://wiki.fysik.dtu.dk/gpaw-files/gpaw-setups-${gpaw_setups_version}.tar.gz
tar zxf gpaw-setups-${gpaw_setups_version}.tar.gz

mkdir -p ${MODULEFILES}/nose
cat <<EOF > ${MODULEFILES}/nose/${nose_version}-1
#%Module1.0
set apps_path ${APPS}
prepend-path    PATH                \$apps_path/nose-${nose_version}-1/usr/bin
prepend-path    PYTHONPATH          \$apps_path/nose-${nose_version}-1/usr/lib/python2.4/site-packages/
unset apps_path
EOF

mkdir -p ${MODULEFILES}/numpy
cat <<EOF > ${MODULEFILES}/numpy/${numpy_version}-1
#%Module1.0
set apps_path ${APPS}
prereq nose # nose
prepend-path    PATH                \$apps_path/numpy-${numpy_version}-1/usr/bin
prepend-path    PYTHONPATH          \$apps_path/numpy-${numpy_version}-1/usr/lib64/python2.4/site-packages/
unset apps_path
EOF

mkdir -p ${MODULEFILES}/campos-ase3
cat <<EOF > ${MODULEFILES}/campos-ase3/${ase_version}-1
#%Module1.0
set apps_path ${APPS}
prereq nose
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
prereq ompi/1.3.2-gnu
prereq numpy
prereq campos-ase3
prereq campos-gpaw-setups
prepend-path    PATH                \$apps_path/gpaw-${gpaw_version}/tools
prepend-path    PATH                \$apps_path/gpaw-${gpaw_version}/build/bin.linux-x86_64-2.4
prepend-path    PYTHONPATH          \$apps_path/gpaw-${gpaw_version}/
prepend-path    PYTHONPATH          \$apps_path/gpaw-${gpaw_version}/build/lib.linux-x86_64-2.4
setenv OMP_NUM_THREADS 1
unset apps_path
EOF

module avail

module use --append ${MODULEFILES}
module load nose
module load numpy
module load campos-ase3
# test numpy
python -c "import numpy; numpy.test()"
# test ase
mkdir -p tmp
cd tmp
testase.py --no-display
cd ..
# build gpaw
cd gpaw-${gpaw_version}
module load ompi/1.3.2-gnu
python setup.py build_ext --customize=../customize_newcell.py --remove-default-flags
cd ..
module load campos-gpaw-setups
module load campos-gpaw
