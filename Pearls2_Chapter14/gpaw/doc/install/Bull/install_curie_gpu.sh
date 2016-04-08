export APPS=`readlink -f ~/CAMd`
export CAMD_MODULEFILES="${APPS}/modulefiles"

export GPAW_PLATFORM=`python -c "from distutils import util, sysconfig; print util.get_platform()+'-'+sysconfig.get_python_version()"`
export PYTHONVERSION=`python -c "from distutils import sysconfig; print sysconfig.get_python_version()"`

# build packages

nose_version=1.1.2
tar zxf nose-${nose_version}.tar.gz
cd nose-${nose_version}
python setup.py install --root=${APPS}/nose-${nose_version}-1
cd ..

mkdir -p ${CAMD_MODULEFILES}/nose
cat <<EOF > ${CAMD_MODULEFILES}/nose/${nose_version}-1
#%Module1.0
set apps_path ${APPS}
prepend-path    PATH                \$apps_path/nose-${nose_version}-1/usr/bin
prepend-path    PYTHONPATH          \$apps_path/nose-${nose_version}-1/usr/lib/python${PYTHONVERSION}/site-packages/
unset apps_path
EOF

numpy_version=1.5.1
tar zxf numpy-${numpy_version}.tar.gz
cd  numpy-${numpy_version}
# atlas on curie is built without -fPIC
# /usr/bin/ld: /usr/local/atlas-3.9.72/lib/libcblas.a(cblas_dgemm.o): relocation R_X86_64_32 against `.rodata.str1.8' can not be used when making a shared object; recompile with -fPIC
# and atlas-devel is not installed!
# so hack!
ln -s /usr/lib64/atlas/libatlas.so.3.0 libatlas.so
ln -s /usr/lib64/atlas/libcblas.so.3.0 libcblas.so
ln -s /usr/lib64/atlas/libclapack.so.3 libclapack.so
ln -s /usr/lib64/atlas/libf77blas.so.3 libf77blas.so
ln -s /usr/lib64/atlas/liblapack.so.3 liblapack.so
echo "[DEFAULT]" > site.cfg
echo "library_dirs = $PWD" >> site.cfg
echo "include_dirs = /usr/local/atlas-3.9.72/include" >> site.cfg
# avoid "Both g77 and gfortran runtimes linked in lapack_lite !" setting --fcompiler=gnu95
# note that this forces /usr/bin/gfortran to be used
python setup.py build --fcompiler=gnu95 2>&1 | tee build.log
python setup.py install --root=${APPS}/numpy-${numpy_version}-1 2>&1 | tee install.log
cd ..

mkdir -p ${CAMD_MODULEFILES}/numpy
cat <<EOF > ${CAMD_MODULEFILES}/numpy/${numpy_version}-1
#%Module1.0
set apps_path ${APPS}
prereq nose
prepend-path    PATH                \$apps_path/numpy-${numpy_version}-1/usr/bin
prepend-path    PYTHONPATH          \$apps_path/numpy-${numpy_version}-1/usr/lib64/python${PYTHONVERSION}/site-packages/
unset apps_path
EOF

# the atlas is missing on the hybrid (only?) nodes so hack again!
mkdir atlas && cd atlas
cp -p /usr/lib64/atlas/liblapack.so.3 .
cp -p /usr/lib64/atlas/libf77blas.so.3 .
cp -p /usr/lib64/atlas/libcblas.so.3 .
cp -p /usr/lib64/atlas/libatlas.so.3 .
cd ..

module use $CAMD_MODULEFILES
module load nose
module load numpy  # scipy build needs numpy!

scipy_version=0.9.0
tar zxf scipy-${scipy_version}.tar.gz
cd  scipy-${scipy_version}
# avoid g77 - leads to Segmentation faults
# note that this forces /usr/bin/gfortran to be used
python setup.py build --fcompiler=gnu95 2>&1 | tee build.log
python setup.py install --root=${APPS}/scipy-${scipy_version}-1 2>&1 | tee install.log
cd ..

mkdir -p ${CAMD_MODULEFILES}/scipy
cat <<EOF > ${CAMD_MODULEFILES}/scipy/${scipy_version}-1
#%Module1.0
set apps_path ${APPS}
prereq nose
prepend-path    PATH                \$apps_path/scipy-${scipy_version}-1/usr/bin
prepend-path    PYTHONPATH          \$apps_path/scipy-${scipy_version}-1/usr/lib64/python${PYTHONVERSION}/site-packages/
unset apps_path
EOF

ase_version=3.7.0.3168
tar zxf python-ase-${ase_version}.tar.gz

mkdir -p ${CAMD_MODULEFILES}/python-ase
cat <<EOF > ${CAMD_MODULEFILES}/python-ase/${ase_version}-1
#%Module1.0
set apps_path ${APPS}
prereq numpy
prepend-path    PATH                \$apps_path/python-ase-${ase_version}/tools
prepend-path    PYTHONPATH          \$apps_path/python-ase-${ase_version}/
unset apps_path
EOF

gpaw_setups_version=0.9.9672
tar zxf gpaw-setups-${gpaw_setups_version}.tar.gz

mkdir -p ${CAMD_MODULEFILES}/gpaw-setups
cat <<EOF > ${CAMD_MODULEFILES}/gpaw-setups/${gpaw_setups_version}-1
#%Module1.0
set apps_path ${APPS}
prepend-path    GPAW_SETUP_PATH     \$apps_path/gpaw-setups-${gpaw_setups_version}
unset apps_path
EOF

gpaw_version=0.9.0.8965
tar zxf gpaw-${gpaw_version}.tar.gz

mkdir -p ${CAMD_MODULEFILES}/gpaw
cat <<EOF > ${CAMD_MODULEFILES}/gpaw/${gpaw_version}-1
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
module load scipy
# test numpy and scipy
python -c "import numpy; numpy.test()"
python -c "import scipy; scipy.test()"

# test ase
module load python-ase
mkdir -p testase
cd testase
testase.py --no-display 2>&1 | tee testase.log
cd ..

# build gpaw
cd gpaw-${gpaw_version}
module load cuda
# if on rpa-gpu-expt branch
rm -f c/cukernels.o
cd c
nvcc -arch sm_20 -c cukernels.cu -Xcompiler -fPIC
cd ..
# wget https://svn.fysik.dtu.dk/projects/gpaw/trunk/config.py  # fixed in trunk
python setup.py build_ext --customize=../customize_curie_gpu.py --remove-default-flags 2>&1 | tee build_ext.log
cd ..
module load gpaw-setups
module load gpaw
