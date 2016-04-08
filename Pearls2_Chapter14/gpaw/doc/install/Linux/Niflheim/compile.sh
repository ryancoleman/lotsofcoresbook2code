#!/bin/sh

if test -z $GPAW_HOME;
    then
    echo "Error: \$GPAW_HOME variable not set"
    exit 1
fi
if test -e $GPAW_HOME/build;
    then
    echo "Error: \$GPAW_HOME/build directory exists - please remove it manually in order to proceed"
    exit 1
fi

rm -rf $GPAW_HOME/build/
#echo "source /home/camp/modulefiles.sh&& module purge&& module load open64/4.2.3-0&& module load NUMPY&& cd $GPAW_HOME&& python setup.py --remove-default-flags --customize=./doc/install/Linux/Niflheim/el5-opteron-open64-acml-4.4.0-acml-4.4.0-hdf-SL-2.0.1.py build_ext 2>&1 | tee compile-el5-opteron-open64-acml-4.4.0-acml-4.4.0-hdf-SL-2.0.1.log" | ssh fjorm bash
echo "source /home/opt/modulefiles/modulefiles_el6.sh&& module purge&& module load NUMPY/1.7.1-1&& cd $GPAW_HOME&& python setup.py --remove-default-flags --customize=./doc/install/Linux/Niflheim/el6-x3455-tm-gfortran-openmpi-1.6.3-openblaso-0.2.8.1-sl-hdf5-1.8.10.py build_ext 2>&1 | tee compile-el6-x3455-tm-gfortran-openmpi-1.6.3-openblaso-0.2.8.1-sl-hdf5-1.8.10.log" | ssh slid bash
echo "source /home/opt/modulefiles/modulefiles_el6.sh&& module purge&& module load NUMPY/1.7.1-1&& cd $GPAW_HOME&& python setup.py --remove-default-flags --customize=./doc/install/Linux/Niflheim/el6-dl160g6-tm-gfortran-openmpi-1.6.3-openblaso-0.2.8.1-sl-hdf5-1.8.10.py build_ext 2>&1 | tee compile-el6-dl160g6-tm-gfortran-openmpi-1.6.3-openblaso-0.2.8.1-sl-hdf5-1.8.10.log" | ssh muspel bash
echo "source /home/opt/modulefiles/modulefiles_el6.sh&& module purge&& module load NUMPY/1.7.1-1&& cd $GPAW_HOME&& python setup.py --remove-default-flags --customize=./doc/install/Linux/Niflheim/el6-sl230s-tm-gfortran-openmpi-1.6.3-openblaso-0.2.8.1-sl-hdf5-1.8.10.py build_ext 2>&1 | tee compile-el6-sl230s-tm-gfortran-openmpi-1.6.3-openblaso-0.2.8.1-sl-hdf5-1.8.10.log" | ssh surt bash
# TAU
#echo "source /home/camp/modulefiles.sh&& module purge&& module load NUMPY&& module load TAU&&cd $GPAW_HOME&& python setup.py --remove-default-flags --customize=./doc/install/Linux/Niflheim/el5-opteron-gcc43-goto2-1.13-acml-4.4.0-TAU.py build_ext 2>&1 | tee compile-el5-opteron-gcc43-goto2-1.13-acml-4.4.0-TAU.log" | ssh fjorm bash

