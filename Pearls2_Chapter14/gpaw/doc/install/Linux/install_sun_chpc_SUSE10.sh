export APPS=`echo ~/CAMd`
export MODULEFILES="${APPS}/modulefiles"

# build packages

python_version=2.7.3
wget http://www.python.org/ftp/python/${python_version}/Python-${python_version}.tgz
tar zxf Python-${python_version}.tgz
cd Python-${python_version}
./configure --prefix=${APPS}/Python-${python_version}-1
make 2>&1 | tee make.log
make install 2>&1 | tee make_install.log
cd ..

mkdir -p ${MODULEFILES}/python
cat <<EOF > ${MODULEFILES}/python/${python_version}-1
#%Module1.0
set apps_path ${APPS}
prepend-path    PATH                \$apps_path/Python-${python_version}-1/bin
unset apps_path
EOF

module use --append ${MODULEFILES}
module load python

export GPAW_PLATFORM=`python -c "from distutils import util, sysconfig; print util.get_platform()+'-'+sysconfig.get_python_version()"`
export PYTHONVERSION=`python -c "from distutils import sysconfig; print sysconfig.get_python_version()"`

nose_version=1.1.2
tar zxf nose-${nose_version}.tar.gz
cd nose-${nose_version}
python setup.py install --root=${APPS}/nose-${nose_version}-1
cd ..

numpy_version=1.6.1
tar zxf numpy-${numpy_version}.tar.gz
cd  numpy-${numpy_version}
# Use /usr/bin/gfortran and numpy's internal blas/lapack
sed -i "s/_lib_names = \['blas'\]/_lib_names = ['']/g"  numpy/distutils/system_info.py
sed -i "s/_lib_names = \['lapack'\]/_lib_names = ['']/g"  numpy/distutils/system_info.py
# avoid "Both g77 and gfortran runtimes linked in lapack_lite !" setting --fcompiler=gnu95
# note that this forces /usr/bin/gfortran to be used
python setup.py build --fcompiler=gnu95 2>&1 | tee build.log
python setup.py install --root=${APPS}/numpy-${numpy_version}-1 2>&1 | tee install.log
cd ..

ase_version=3.6.0.2515
tar zxf python-ase-${ase_version}.tar.gz

gpaw_version=0.9.0.8965
tar zxf gpaw-${gpaw_version}.tar.gz

gpaw_setups_version=0.8.7929
tar zxf gpaw-setups-${gpaw_setups_version}.tar.gz

mkdir -p ${MODULEFILES}/nose
cat <<EOF > ${MODULEFILES}/nose/${nose_version}-1
#%Module1.0
set apps_path ${APPS}
prereq python
prepend-path    PATH                \$apps_path/nose-${nose_version}-1${APPS}/Python-${python_version}-1/bin
prepend-path    PYTHONPATH          \$apps_path/nose-${nose_version}-1${APPS}/Python-${python_version}-1/lib/python${PYTHONVERSION}/site-packages/
unset apps_path
EOF

mkdir -p ${MODULEFILES}/numpy
cat <<EOF > ${MODULEFILES}/numpy/${numpy_version}-1
#%Module1.0
set apps_path ${APPS}
prereq nose
prepend-path    PATH                \$apps_path/numpy-${numpy_version}-1${APPS}/Python-${python_version}-1/bin
prepend-path    PYTHONPATH          \$apps_path/numpy-${numpy_version}-1${APPS}/Python-${python_version}-1/lib/python${PYTHONVERSION}/site-packages/
unset apps_path
EOF

mkdir -p ${MODULEFILES}/python-ase
cat <<EOF > ${MODULEFILES}/python-ase/${ase_version}-1
#%Module1.0
set apps_path ${APPS}
prereq numpy
prepend-path    PATH                \$apps_path/python-ase-${ase_version}/tools
prepend-path    PYTHONPATH          \$apps_path/python-ase-${ase_version}/
unset apps_path
EOF

mkdir -p ${MODULEFILES}/gpaw-setups
cat <<EOF > ${MODULEFILES}/gpaw-setups/${gpaw_setups_version}-1
#%Module1.0
set apps_path ${APPS}
prepend-path    GPAW_SETUP_PATH     \$apps_path/gpaw-setups-${gpaw_setups_version}
unset apps_path
EOF

mkdir -p ${MODULEFILES}/gpaw
cat <<EOF > ${MODULEFILES}/gpaw/${gpaw_version}-1
#%Module1.0
set apps_path ${APPS}
prereq python-ase
prereq gpaw-setups
prepend-path    PATH                \$apps_path/gpaw-${gpaw_version}/tools
prepend-path    PATH                \$apps_path/gpaw-${gpaw_version}/build/bin.${GPAW_PLATFORM}
prepend-path    PYTHONPATH          \$apps_path/gpaw-${gpaw_version}/
prepend-path    PYTHONPATH          \$apps_path/gpaw-${gpaw_version}/build/lib.${GPAW_PLATFORM}
setenv OMP_NUM_THREADS 1
unset apps_path
EOF

module load nose
module load numpy
# test numpy
python -c "import numpy; numpy.test()"

module load python-ase
# test ase
mkdir -p testase
cd testase
testase.py --no-display 2>&1 | tee testase.log
cd ..
# build gpaw
cd gpaw-${gpaw_version}
# disable fftw (something is wrong with /usr/local/lib/libfftw3.a)
sed -i 's/libfftw3\.so/libfftw3_not.so/' gpaw/fftw.py
python setup.py build_ext --customize=../customize_sun_chpc_SUSE10.py --remove-default-flags 2>&1 | tee build_ext.log
cd ..
module load gpaw-setups
module load gpaw
mkdir -p testgpaw
cd testgpaw
mpiexec -np 4 gpaw-python `which gpaw-test` 2>&1 | tee testgpaw.log
