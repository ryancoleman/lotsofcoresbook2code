#!/bin/sh

export APPS="/home/firegam/CAMd"
export MODULEFILES="${APPS}/modulefiles"

INSTALL_DACAPO=False

# the following packages are needed
#sudo apt-get install python-all-dev
#sudo apt-get install python-numpy
#sudo apt-get install python-scientific
#sudo apt-get install python-matplotlib

#sudo apt-get install netcdf-dev # for dacapo
#sudo apt-get install netcdf-bin # for dacapo
#sudo apt-get install fftw-dev # for dacapo
#sudo apt-get install gfortran # for dacapo

# build packages

openmpi_version=1.4.3
tar jxf openmpi-${openmpi_version}.tar.bz2
cd openmpi-${openmpi_version}
./configure --prefix=${APPS}/openmpi-${openmpi_version}-1
make && make install
cd ..

nose_version=0.11.3
tar zxf nose-${nose_version}.tar.gz
cd nose-${nose_version}
python setup.py install --root=${APPS}/nose-${nose_version}-1
cd ..

numpy_version=1.5.0
tar zxf numpy-${numpy_version}.tar.gz
cd  numpy-${numpy_version}
# disable compiling with atlas
sed -i "s/_lib_atlas =.*/_lib_atlas = ['ignore_atlas']/g" numpy/distutils/system_info.py
python setup.py install --root=${APPS}/numpy-${numpy_version}-1
cd ..

ase_version=3.4.1.1765
tar zxf python-ase-${ase_version}.tar.gz

gpaw_version=0.7.2.6974
tar zxf gpaw-${gpaw_version}.tar.gz

gpaw_setups_version=0.6.6300
tar zxf gpaw-setups-${gpaw_setups_version}.tar.gz

if [ "${INSTALL_DACAPO}" = "True" ];
    then
    dacapo_version=2.7.16
    tar zxf campos-dacapo-${dacapo_version}.tar.gz
    cd campos-dacapo-${dacapo_version}
    export FFTW=/usr/lib64
    export NETCDF=/usr/lib64
    export BLASLAPACK='${NETCDF}/libnetcdff.a ${NETCDF}/libnetcdf.a -L/opt/intel/Compiler/11.1/072/mkl/lib/em64t/ -lmkl_lapack -lmkl_core -lmkl_sequential -lmkl_gf_lp64'
    export MPIDIR=/home/firegam/CAMd/openmpi-1.4.3-1
    export MPI_LIBDIR=${MPIDIR}/lib
    export MPI_BINDIR=${MPIDIR}/bin
    export MPI_INCLUDEDIR=${MPIDIR}/include
    cd src
    # patch -f90: f951: error: unrecognized command line option "-f90=gfortran"
    sed -i 's#-f90=$(GFORTRAN_FNOSECOND_UNDERSCORE_FC90)##g' Makefile
    make gfortran_fnosecond_underscore
    make gfortran_fnosecond_underscore MP=mpi
    cd gfortran_fnosecond_underscore_serial
    ln -s dacapo.run dacapo_${dacapo_version}-1_serial.run
    cd ..
    cd gfortran_fnosecond_underscore_mpi
    ln -s dacapo.run dacapo_${dacapo_version}-1_mpi.run
    cd ..
    cd ..
    wget https://svn.fysik.dtu.dk/projects/ase/trunk/ase/calculators/jacapo/tools/dacapo.run
    chmod ugo+rx dacapo.run
    cd ..

    dacapo_psp_version=1
    tar zxf campos-dacapo-pseudopotentials-${dacapo_psp_version}.tar.gz
    mkdir psp-1
    cp -pr campos-dacapo-pseudopotentials-1/psp/*/*/*.pseudo psp-1
fi

. ./set_env_armageddon.sh

# test numpy
python -c "import numpy; numpy.test()"
# test ase
mkdir -p testase
cd testase
testase.py --no-display 2>&1 | tee testase.log
cd ..
# build gpaw
cd gpaw-${gpaw_version}
python setup.py build_ext --customize=../customize_armageddon.py --remove-default-flags
cd ..
mkdir -p testgpaw
cd testgpaw
mpiexec -np 4 gpaw-python `which gpaw-test` 2>&1 | tee testgpaw.log

if [ "${INSTALL_DACAPO}" == "True" ];
    then
    mkdir -p testjacapo
    cd testjacapo
    # from https://wiki.fysik.dtu.dk/ase/ase/calculators/jacapo.html
cat <<EOF > CO.py
#!/usr/bin/env python
from ase import *
from ase.structure import molecule
from ase.calculators.jacapo import *

CO = molecule('CO')
CO.set_cell([6,6,6])
CO.center()

calc = Jacapo(nc='CO.nc',
          atoms=CO,
          pw=300,
          nbands=8)

print CO.get_potential_energy()

EOF

# test on 4 cores
python -c "for i in range(4): print 'localhost'" > pbs_nodefile && PBS_NODEFILE=pbs_nodefile python CO.py
fi
