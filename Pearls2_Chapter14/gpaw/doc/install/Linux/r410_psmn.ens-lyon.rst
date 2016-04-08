.. _r410_psmn.ens_lyon:

==================
r410 psmn.ens-lyon
==================

Information about the machine https://www.psmn.ens-lyon.fr.

Instructions assume **cshrc**, installation under *${HOME}/softs*.

Enable the basic settings::

  source /usr/local/modeles/Cshrc
  # The next line causes gcc 4.2.1 to be loaded from /usr/local/bin: crazy!
  # source /usr/local/modeles/openmpi-1.4.1-intel-11.1.069-mkl
  # You must NOT do that:
  # when building GPAW, and other apllications using /usr/bin/gcc
  # remove it temporarily from ~/.tcshrc and logout/login again!

Setup the root directory::

  setenv HOMESOFTS ${HOME}/softs

  mkdir -p ${HOMESOFTS}
  cd ${HOMESOFTS}

  setenv GPAW_PLATFORM `python -c "from distutils import util, sysconfig; print util.get_platform()+'-'+sysconfig.get_python_version()"`

Set the versions::

  set nose=0.11.3
  # Warning: version 1.6.0 seems inconsistent about C-, Fortran-contiguous
  # http://mail.scipy.org/pipermail/numpy-discussion/2011-July/057557.html
  set numpy=1.5.1
  set scipy=0.9.0

  set acml=4.0.1

  set ase=3.5.1.2175
  set gpaw=0.8.0.8092
  set setups=0.8.7929
  
and create startup scripts::

  cat <<EOF > ${HOMESOFTS}/acml-${acml}-1.cshrc
  setenv LD_LIBRARY_PATH \${HOMESOFTS}/acml-${acml}/gfortran64/lib:\${LD_LIBRARY_PATH}
  EOF

  #cat <<EOF > ${HOMESOFTS}/acml-${acml}-1.cshrc
  #setenv LD_LIBRARY_PATH /softs/acml/acml${acml}/gfortran/gfortran64/lib:\${LD_LIBRARY_PATH}
  #EOF

  cat <<EOF > ${HOMESOFTS}/nose-${nose}-1.cshrc
  setenv PYTHONPATH \${HOMESOFTS}/nose-${nose}-1/usr/lib/python2.4/site-packages:\${PYTHONPATH}
  setenv PATH \${HOMESOFTS}/nose-${nose}-1/usr/bin:\${PATH}
  EOF

  cat <<EOF > ${HOMESOFTS}/numpy-${numpy}-1.cshrc
  setenv PYTHONPATH \${HOMESOFTS}/numpy-${numpy}-1/usr/lib64/python2.4/site-packages:\${PYTHONPATH}
  setenv PATH \${HOMESOFTS}/numpy-${numpy}-1/usr/bin:\${PATH}
  EOF

  cat <<EOF > ${HOMESOFTS}/scipy-${scipy}-1.cshrc
  setenv PYTHONPATH \${HOMESOFTS}/scipy-${scipy}-1/usr/lib64/python2.4/site-packages:\${PYTHONPATH}
  setenv PATH \${HOMESOFTS}/scipy-${scipy}-1/usr/bin:\${PATH}
  EOF

  cat <<EOF > ${HOMESOFTS}/python-ase-${ase}.cshrc
  setenv PYTHONPATH \${HOMESOFTS}/python-ase-${ase}:\${PYTHONPATH}
  setenv PATH \${HOMESOFTS}/python-ase-${ase}/tools:\${PATH}
  EOF

  cat <<EOF > ${HOMESOFTS}/gpaw-setups-${setups}.cshrc
  setenv GPAW_SETUP_PATH \${HOMESOFTS}/gpaw-setups-${setups}
  EOF

  cat <<EOF > ${HOMESOFTS}/gpaw-${gpaw}.cshrc
  setenv GPAW_PLATFORM `python -c "from distutils import util, sysconfig; print util.get_platform()+'-'+sysconfig.get_python_version()"`
  setenv GPAW_HOME \${HOMESOFTS}/gpaw-${gpaw}
  setenv PYTHONPATH \${GPAW_HOME}:\${PYTHONPATH}
  setenv PYTHONPATH \${GPAW_HOME}/build/lib.${GPAW_PLATFORM}:\${PYTHONPATH}
  setenv PATH \${GPAW_HOME}/build/bin.${GPAW_PLATFORM}:\${PATH}
  setenv PATH \${GPAW_HOME}/tools:\${PATH}
  EOF

  cat <<EOF > ${HOMESOFTS}/campos.cshrc
  #!/bin/tcsh
  #
  source ${HOMESOFTS}/acml-${acml}-1.cshrc
  #
  if ! ($?PYTHONPATH) then
      setenv PYTHONPATH ""
  endif
  #
  setenv HOMESOFTS \${HOME}/softs
  #
  source \${HOMESOFTS}/nose-${nose}-1.cshrc
  source \${HOMESOFTS}/numpy-${numpy}-1.cshrc
  source \${HOMESOFTS}/scipy-${scipy}-1.cshrc
  #
  source \${HOMESOFTS}/python-ase-${ase}.cshrc
  #
  source \${HOMESOFTS}/gpaw-setups-${setups}.cshrc
  source \${HOMESOFTS}/gpaw-${gpaw}.cshrc
  EOF

Make sure that you have the right compiler::

  % which gfortran  # gfortran 4.1.2
  /usr/bin/gfortran

Build nose/numpy/scipy::

  wget --no-check-certificate https://downloads.sourceforge.net/project/numpy/NumPy/${numpy}/numpy-${numpy}.tar.gz
  wget --no-check-certificate https://downloads.sourceforge.net/project/scipy/scipy/${scipy}/scipy-${scipy}.tar.gz
  wget http://python-nose.googlecode.com/files/nose-${nose}.tar.gz
  tar zxf nose-${nose}.tar.gz
  tar zxf numpy-${numpy}.tar.gz
  tar zxf scipy-${scipy}.tar.gz
  cd nose-${nose}
  python setup.py install --root=${HOMESOFTS}/nose-${nose}-1 >& install.log
  cd ..

continue with::

  cd numpy-${numpy}
  # Warning: numpy with gfortran44 fails under gpaw-python
  # ImportError: numpy.core.multiarray failed to import
  # Use /usr/bin/gfortran and numpy's internal blas/lapack
  sed -i "s/_lib_names = \['blas'\]/_lib_names = ['']/g"  numpy/distutils/system_info.py
  sed -i "s/_lib_names = \['lapack'\]/_lib_names = ['']/g"  numpy/distutils/system_info.py
  # avoid "Both g77 and gfortran runtimes linked in lapack_lite !" setting --fcompiler=gnu95
  # note that this forces /usr/bin/gfortran to be used
  python setup.py build --fcompiler=gnu95 >& build.log
  python setup.py install --root=${HOMESOFTS}/numpy-${numpy}-1 >& install.log
  cd ..
  source ${HOMESOFTS}/campos.cshrc
  python -c "import numpy; numpy.test()"

  # scipy cannot be installed due to missing suitesparse-devel and suitesparse
  #cd scipy-${scipy}
  #python setup.py config_fc --fcompiler=gfortran install --root=${HOMESOFTS}/scipy-${scipy}-1 2>&1 | tee install.log
  #cd ..
  #python -c "import scipy; scipy.test()"

Install ASE/GPAW::

  wget https://wiki.fysik.dtu.dk/ase-files/python-ase-${ase}.tar.gz
  wget https://wiki.fysik.dtu.dk/gpaw-files/gpaw-${gpaw}.tar.gz
  wget http://wiki.fysik.dtu.dk/gpaw-files/gpaw-setups-${setups}.tar.gz
  tar zxf python-ase-${ase}.tar.gz
  tar zxf gpaw-${gpaw}.tar.gz
  tar zxf gpaw-setups-${setups}.tar.gz
  mkdir testase && cd testase && testase.py >& ../testase.log
  wget https://svn.fysik.dtu.dk/projects/gpaw/trunk/doc/install/Linux/customize_r410_psmn.py
  cd ../gpaw-${gpaw}
  python setup.py --remove-default-flags --customize=../customize_r410_psmn.py build_ext >& build_ext.log

The :file:`customize_r410_psmn.py` looks like:

.. literalinclude:: customize_r410_psmn.py

GPAW tests :file:`gpaw-test` can be submitted like this::

  qsub -cwd -pe mpi8 8 run.tcsh

where :file:`run.tcsh` looks like this::

  #!/bin/tcsh

  setenv HOSTFILE $TMPDIR/machines

  setenv MPIPREFIX /softs/openmpi-gnu

  source ${HOME}/softs/campos.cshrc

  setenv OMP_NUM_THREADS 1

  ${MPIPREFIX}/bin/mpirun -v -prefix ${MPIPREFIX} -mca pls_rsh_agent "ssh" -mca btl openib,tcp,self -mca btl_tcp_if_include eth1,eth0 -np $NSLOTS `which gpaw-python` `which gpaw-test`

Only some nodes have infiniband, on those without ignore the errors::

  OpenIB on host r410lin???.ens-lyon.fr was unable to find any HCAs

Please make sure that your jobs do not run multi-threaded, e.g. for a
job running on ``dl165lin7`` do from a login node::

  ssh dl165lin7 ps -fL

you should see **1** in the **NLWP** column. Numbers higher then **1**
mean multi-threaded job.

It's convenient to customize as described on the :ref:`parallel_runs` page.
