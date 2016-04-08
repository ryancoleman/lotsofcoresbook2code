#!/bin/sh
export CC=/bgsys/drivers/ppcfloor/gnu-linux/bin/powerpc64-bgq-linux-gcc
export BASECFLAGS="-fno-strict-aliasing"
export LD_LIBRARY_PATH=/bgsys/drivers/ppcfloor/gnu-linux/lib64
export PYTHONHOME=/soft/apps/python/scalable-python-2.6.7-cnk-gcc
export PYTHON=${PYTHONHOME}/bin/python
# root=/soft/apps/python/scalable-python-2.6.7-cnk-gcc
buildir=build

rm -rf ${builddir}
# ${PYTHON} setup.py install --root="$root" 2>&1 | tee numpy-1.3.0.log.mira
${PYTHON} setup.py install 2>&1 | tee numpy-1.3.0.log.mira