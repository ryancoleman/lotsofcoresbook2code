.. _vsc.univie:

==========
vsc.univie
==========

The vsc.univie machine is a cluster of dual socket, hexa-core Intel Xeon X5650
2.67 GHz processors with 2 GB of memory per core.

Instructions assume **bash**, installation under ``${HOME}/opt``.

Setup the root directory::

  mkdir -p ${HOME}/opt
  cd ${HOME}/opt

Set the versions::

  export nose=0.11.3
  # Warning: version 1.6.0 seems inconsistent about C-, Fortran-contiguous
  # http://mail.scipy.org/pipermail/numpy-discussion/2011-July/057557.html
  export numpy=1.5.1
  export scipy=0.9.0

  export acml=4.0.1

  export ase=3.5.1.2175
  export gpaw=0.8.0.8092
  export setups=0.8.7929
  
and create sh startup script::

  cat <<EOF > ${HOME}/opt/campos.sh
  #!/bin/sh
  #
  export GPAW_PLATFORM=`python -c "from distutils import util, sysconfig; print util.get_platform()+'-'+sysconfig.get_python_version()"`
  #
  export LD_LIBRARY_PATH=\${HOME}/opt/acml-${acml}/gfortran64/lib:\${LD_LIBRARY_PATH}
  export LD_LIBRARY_PATH=\${HOME}/opt/CBLAS.acml-${acml}/lib:\${LD_LIBRARY_PATH} # if cblas used
  #
  export PYTHONPATH=\${HOME}/opt/nose-${nose}-1/usr/lib/python2.4/site-packages:\${PYTHONPATH}
  export PATH=\${HOME}/opt/nose-${nose}-1/usr/bin:\${PATH}
  #
  export PYTHONPATH=\${HOME}/opt/numpy-${numpy}-1/usr/lib64/python2.4/site-packages:\${PYTHONPATH}
  export PATH=\${HOME}/opt/numpy-${numpy}-1/usr/bin:\${PATH}
  #
  export PYTHONPATH=\${HOME}/opt/scipy-${scipy}-1/usr/lib64/python2.4/site-packages:\${PYTHONPATH}
  export PATH=\${HOME}/opt/scipy-${scipy}-1/usr/bin:\${PATH}
  #
  export PYTHONPATH=\${HOME}/opt/python-ase-${ase}:\${PYTHONPATH}
  export PATH=\${HOME}/opt/python-ase-${ase}/tools:\${PATH}
  #
  export GPAW_SETUP_PATH=\${HOME}/opt/gpaw-setups-${setups}
  #
  export GPAW_HOME=\${HOME}/opt/gpaw-${gpaw}
  export PYTHONPATH=\${GPAW_HOME}:\${PYTHONPATH}
  export PYTHONPATH=\${GPAW_HOME}/build/lib.${GPAW_PLATFORM}:\${PYTHONPATH}
  export PATH=\${GPAW_HOME}/build/bin.${GPAW_PLATFORM}:\${PATH}
  export PATH=\${GPAW_HOME}/tools:\${PATH}
  EOF

Download and install acml::

  acml-${acml} # download
  cd acml-${acml}
  tar zxf acml-*.tgz && tar zxf contents-acml-*.tgz

**Note**: numpy with acml dotblas Segmentation Faults (well, for some
versions on numpy, etc?)  for
:file:`gpaw/test/numpy_core_multiarray_dot.py` or
:file:`gpaw/test/gemm.py`.  Still there is no performance improvement
for :file:`gpaw/test/gemm.py` (if it works), even if case of dynamic
linking of cblas/acml - check with ldd that _dotblas.so is linked to
both acml and cblas.

This is how you can download and install cblas::

  wget http://www.netlib.org/blas/blast-forum/cblas.tgz
  tar zxf cblas.tar.gz && mv -f CBLAS CBLAS.acml-${acml} && cd CBLAS.acml-${acml}

  cp -p Makefile.LINUX Makefile.in
  # fix Makefile.in
  export PLAT=LINUX
  export BLLIB=${HOME}/opt/acml-${acml}/gfortran64/lib/libacml.a
  export CC=gcc
  export FC=gfortran
  export CFLAGS='-O3 -funroll-all-loops -DADD_ -fPIC'
  export FFLAGS='-O3 -funroll-all-loops -DADD_ -fPIC'
  #
  sed -i "s<^PLAT =.*<PLAT = ${PLAT}<" Makefile.in
  sed -i "s<^BLLIB =.*<BLLIB = ${BLLIB}<" Makefile.in
  sed -i "s<^CC =.*<CC = ${CC}<" Makefile.in
  sed -i "s<^FC =.*<FC = ${FC}<" Makefile.in
  sed -i "s<^LOADER =.*<LOADER = \$(FC) -lpthread<" Makefile.in
  sed -i "s<^CFLAGS =.*<CFLAGS = ${CFLAGS}<" Makefile.in
  sed -i "s<^FFLAGS =.*<FFLAGS = ${FFLAGS}<" Makefile.in

  # for dynamic library add the following to the Makefile (before cleanall:)
  # Remember TAB in the Makefile!
  shared: alllib
        ( mkdir tmp && cd tmp && cp $(CBLIB) . && ar x $(CBLIB) && $(CC) -shared -o libcblas.so.1.0.0 *.o -Wl,-soname=libcblas.so.1 && cp -p libcblas.so.1.0.0 ../lib && cd ../lib && ln -s libcblas.so.1.0.0 libcblas.so.1 && ln -s libcblas.so.1.0.0 libcblas.so )

  make all 2>&1 | tee make_all.log
  make shared 2>&1 | tee make_shared.log

  # create link: numpy needs all the libraries in one directory
  # separate directories in site.cfg do not work
  cd lib && ln -s cblas_${PLAT}.a libcblas.a
  ln -s ${HOME}/opt/acml-${acml}/gfortran64/lib/libacml.a .
  ln -s ${HOME}/opt/acml-${acml}/gfortran64/lib/libacml_mv.a .
  # if dynamic library needed
  ln -s ${HOME}/opt/acml-${acml}/gfortran64/lib/libacml.so .
  ln -s ${HOME}/opt/acml-${acml}/gfortran64/lib/libacml_mv.so .
  cd ../..

Build nose/numpy/scipy::

  wget --no-check-certificate https://downloads.sourceforge.net/project/numpy/NumPy/${numpy}/numpy-${numpy}.tar.gz
  wget --no-check-certificate https://downloads.sourceforge.net/project/scipy/scipy/${scipy}/scipy-${scipy}.tar.gz
  wget http://python-nose.googlecode.com/files/nose-${nose}.tar.gz
  tar zxf nose-${nose}.tar.gz
  tar zxf numpy-${numpy}.tar.gz
  tar zxf scipy-${scipy}.tar.gz
  cd nose-${nose}
  python setup.py install --root=${HOME}/opt/nose-${nose}-1 2>&1 | tee install.log

use the following ``site.cfg`` to build numpy without cblas (that's safer)::

  cat <<EOF > ${HOME}/opt/numpy-${numpy}/site.cfg
  [DEFAULT]
  library_dirs = /usr/lib64:${HOME}/opt/acml-${acml}/gfortran64/lib
  include_dirs = ${HOME}/opt/acml-${acml}/gfortran64/lib/../include:/usr/include/suitesparse
  [blas]
  libraries = acml
  library_dirs = ${HOME}/opt/acml-${acml}/gfortran64/lib
  [lapack]
  libraries = acml, gfortran
  library_dirs = ${HOME}/opt/acml-${acml}/gfortran64/lib
  EOF

and this one with cblas based on acml::

  cat <<EOF > ${HOME}/opt/numpy-${numpy}/site.cfg
  [DEFAULT]
  library_dirs = /usr/lib64:${HOME}/opt/CBLAS.acml-${acml}/lib
  include_dirs = ${HOME}/opt/acml-${acml}/gfortran64/lib/../include:/usr/include/suitesparse:${HOME}/opt/CBLAS.acml-${acml}/include
  [blas]
  libraries = acml
  library_dirs = ${HOME}/opt/CBLAS.acml-${acml}/lib
  [lapack]
  libraries = acml, gfortran
  library_dirs = ${HOME}/opt/CBLAS.acml-${acml}/lib
  EOF

**Note**: the ``site.cfg`` file is used only to specify directories, libraries
from ``site.cfg`` are ignored by numpy. Moreover numpy needs all the libraries in one directory,
separate directories in site.cfg do not work.

continue with::

  cd ../numpy-${numpy}
  # force numpy to use internal blas for dotblas + acml, note the double quotes!
  sed -i "s/_lib_names = \['blas'\]/_lib_names = ['acml']/g"  numpy/distutils/system_info.py
  sed -i "s/_lib_names = \['lapack'\]/_lib_names = ['acml']/g"  numpy/distutils/system_info.py
  # or with dotblas acml (through cblas) - seems not working or Segmentation Faults
  sed -i "s<_lib_mkl = .*<_lib_mkl = ['acml','cblas']<" numpy/distutils/system_info.py
  sed -i "s<\['mkl_lapack32','mkl_lapack64'\]<['acml','gfortran']<" numpy/distutils/system_info.py
  sed -i "s<l = 'mkl'<l = 'acml'<" numpy/distutils/system_info.py

  # avoid "Both g77 and gfortran runtimes linked in lapack_lite !" setting --fcompiler=gnu95
  python setup.py build --fcompiler=gnu95 2>&1 | tee build.log
  python setup.py install --root=${HOME}/opt/numpy-${numpy}-1 2>&1 | tee install.log
  cd ..
  source ${HOME}/opt/campos.sh
  python -c "import numpy; numpy.test()"

  cd scipy-${scipy}
  python setup.py config_fc --fcompiler=gfortran install --root=${HOME}/opt/scipy-${scipy}-1 2>&1 | tee install.log
  cd ..
  python -c "import scipy; scipy.test()"

Make sure that you have the right mpicc::

  which mpicc
  /usr/mpi/qlogic/bin/mpicc

Install ASE/GPAW::

  wget https://wiki.fysik.dtu.dk/ase-files/python-ase-${ase}.tar.gz
  wget https://wiki.fysik.dtu.dk/gpaw-files/gpaw-${gpaw}.tar.gz
  wget http://wiki.fysik.dtu.dk/gpaw-files/gpaw-setups-${setups}.tar.gz
  tar zxf python-ase-${ase}.tar.gz
  tar zxf gpaw-${gpaw}.tar.gz
  tar zxf gpaw-setups-${setups}.tar.gz
  mkdir testase && cd testase && testase.py 2>&1 | tee ../testase.log
  wget https://svn.fysik.dtu.dk/projects/gpaw/trunk/doc/install/Linux/customize_vsc_univie.py
  cd ../gpaw-${gpaw}
  python setup.py --remove-default-flags --customize=../customize_vsc_univie.py build_ext 2>&1 | tee build_ext.log

The :file:`customize_vsc_univie.py` looks like:

.. literalinclude:: customize_vsc_univie.py

GPAW tests :file:`gpaw-test` can be submitted like this::

  qsub run.sh

where :file:`run.sh` looks like this::

  #!/bin/sh

  #$ -pe mpich 8
  #$ -V
  #$ -M my.name@example.at
  #$ -m be
  #$ -l h_rt=00:50:00

  if [ -z "${PYTHONPATH}" ]
  then
      export PYTHONPATH=""
  fi

  source ${HOME}/opt/campos.sh

  export OMP_NUM_THREADS=1

  mpirun -m $TMPDIR/machines -np $NSLOTS gpaw-python `which gpaw-test`

Please make sure that your jobs do not run multi-threaded, e.g. for a
job running on ``node02`` do from a login node::

  ssh node02 ps -fL

you should see **1** in the **NLWP** column. Numbers higher then **1**
mean multi-threaded job.

It's convenient to customize as described on the :ref:`parallel_runs` page.
