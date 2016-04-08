#!/bin/bash
export ibmcmp_base=/soft/compilers/ibmcmp-nov2012
python=/soft/apps/python/scalable-python-2.6.7-cnk-gcc/bin/python
#python=/bgsys/tools/Python-2.6/bin/python
echo $python
${python} setup.py clean 
# manually remove _hdf5.so due to bug in build system
rm -f build/lib.linux-ppc64-2.6/_hdf5.so
# rm -f build/temp.linux-ppc64-2.6/*.o
# compile _gpaw.so
# ${python} setup.py build_ext --customize=customize_mira_xlc_serial.py >& build_xlc_serial.log
# compile gpaw-python
${python} setup.py build_ext --ignore-numpy --customize=customize_mira_xlc_mpi.py 2>&1 | tee build_xlc_mpi.log
