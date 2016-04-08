.. _Niflheim:

========
Niflheim
========

Information about the Niflheim cluster can be found at
`<https://wiki.fysik.dtu.dk/niflheim>`_.

Preferably use the system default installation of GPAW and setups
- in this case simply enable the modules as described at https://wiki.fysik.dtu.dk/niflheim/Installed_software#enabling-modules, load the default GPAW environment with **module load GPAW**, and submit GPAW jobs with::

   gpaw-qsub -l nodes=... script.py

See https://wiki.fysik.dtu.dk/niflheim/Batch_jobs for details.

.. note::

   the name of the python script comes last, and after it
   any options to GPAW (e.g. --debug).

If you decide to install a development version of GPAW, this is what you do:

0. Install a separate GPAW version for each project. It is important to keep
   track of the GPAW/setups versions used during the project
   and include these details in the publication.
   The instructions below are given for a project called `devel`,
   i.e. the installation path will be :file:`/home/niflheim/$USER/devel/gpaw`.

   **Note** that this path appears in few places in the tools described on
   this page - make sure you change those occurences if installing
   in under a different directory. The instructions assume **bash**.

1. On the ``servcamd`` filesystem (login on your workstation)
   go to a directory on the Niflheim filesystem.

2. Checkout the GPAW source. Make **sure** that
   :file:`/home/niflheim/$USER/devel/gpaw` does **not** exist,
   before running the checkout!::

     svn checkout https://svn.fysik.dtu.dk/projects/gpaw/trunk /home/niflheim/$USER/devel/gpaw

3. Set the :envvar:`GPAW_HOME` environment variable::

     export GPAW_HOME=${HOME}/devel/gpaw

4. To compile the GPAW C-code, run the shell script
   :svn:`~doc/install/Linux/Niflheim/compile.sh` 
   from your gpaw directory, i.e.::

     cd ${GPAW_HOME}
     sh ./doc/install/Linux/Niflheim/compile.sh

   If you have login passwords active,
   this will force you to type your password four times. It is
   possible to remove the need for typing passwords on internal CAMd systems,
   using the procedure described at
   https://wiki.fysik.dtu.dk/it/SshWithoutPassword.

5. Set the relevant environment variables poiting to the compiled GPAW version.
   There are at least two ways to achive that:

    5a. **Prepend** :envvar:`PYTHONPATH` and :envvar:`PATH` environment variables in your login scripts.
    This method provides you with a default GPAW version (`devel`).
    Add to /home/niflheim/$USER/.bashrc::

	if [ "`echo $FYS_PLATFORM`" == "AMD-Opteron-el5" ]; then # fjorm
	    module load GPAW
	    export GPAW_PLATFORM="linux-x86_64-opteron-2.4"
	fi
	if test -n "`echo $FYS_PLATFORM | grep el6`"; then
	    module load GPAW
	    export GPAW_PLATFORM="linux-x86_64-`echo $FYS_PLATFORM | sed 's/-el6//'`-2.6"
	fi
	# GPAW_HOME must be set after loading the GPAW module!
	export GPAW_HOME=${HOME}/devel/gpaw
	export PATH=${GPAW_HOME}/build/bin.${GPAW_PLATFORM}:${PATH}
	export PATH=${GPAW_HOME}/tools:${PATH}
	export PYTHONPATH=${GPAW_HOME}:${PYTHONPATH}
	export PYTHONPATH=${GPAW_HOME}/build/lib.${GPAW_PLATFORM}:${PYTHONPATH}

	if test -n "`echo $FYS_PLATFORM | grep el6`"; then
	# http://docs.python.org/2/using/cmdline.html#envvar-PYTHONDONTWRITEBYTECODE
	    export PYTHONDONTWRITEBYTECODE=1  # disable creation of pyc files
	    module load NUMPY/1.7.1-1
	    module load SCIPY/0.12.0-1
	fi

    **Warning**: from the moment you save settings in
    :file:`/home/niflheim/$USER/.bashrc`, your new jobs
    (also those already waiting in the queue)
    will start using the new version.
    The jobs already running are not affected.
    Consider making such changes with no jobs in the queue.

    5b. Define a batch system submission function in bash.
    This is useful for installing different GPAW versions for different projects:

    - create the following bash script :file:`/home/niflheim/$USER/devel/gpaw/qsub.sh`::

	#!/bin/sh

	if [ -r "/home/camp/modulefiles.sh" ]; then
	    source /home/camp/modulefiles.sh
	fi
	if [ -r "/home/opt/modulefiles/modulefiles_el6.sh" ]; then
	    source /home/opt/modulefiles/modulefiles_el6.sh
	fi

	if [ "`echo $FYS_PLATFORM`" == "AMD-Opteron-el5" ]; then # fjorm
	    module load GPAW
	    export GPAW_PLATFORM="linux-x86_64-opteron-2.4"
	fi
	if test -n "`echo $FYS_PLATFORM | grep el6`"; then
	    module load GPAW
	    export GPAW_PLATFORM="linux-x86_64-`echo $FYS_PLATFORM | sed 's/-el6//'`-2.6"
	fi
	# GPAW_HOME must be set after loading the GPAW module!
	export GPAW_HOME=${HOME}/devel/gpaw
	export PATH=${GPAW_HOME}/build/bin.${GPAW_PLATFORM}:${PATH}
	export PATH=${GPAW_HOME}/tools:${PATH}
	export PYTHONPATH=${GPAW_HOME}:${PYTHONPATH}
	export PYTHONPATH=${GPAW_HOME}/build/lib.${GPAW_PLATFORM}:${PYTHONPATH}

	if test -n "`echo $FYS_PLATFORM | grep el6`"; then
	# http://docs.python.org/2/using/cmdline.html#envvar-PYTHONDONTWRITEBYTECODE
	    export PYTHONDONTWRITEBYTECODE=1  # disable creation of pyc files
	    module load NUMPY/1.7.1-1
	    module load SCIPY/0.12.0-1
	fi

	mpiexec gpaw-python "$name"

      Modify this file if needed (if you need different ASE/setups, etc)!

    - define the corresponding function in :file:`/home/niflheim/$USER/.bashrc`::

	 gpaw-qsub-devel ()
	 {
	 name="$1"
	 shift
	 qsub $@ -v name=$name ${HOME}/devel/gpaw/qsub.sh
	 }

    When submitting jobs specify the python script first!::

	gpaw-qsub-devel script.py -l nodes=...

6. If you prefer to use a personal setup's directory follow
   :ref:`installationguide_setup_files`.

When updating the gpaw code in the future:

- Go to the :envvar:`GPAW_HOME` directory and run::

    svn up

- If any of the c-code changed during the update repeat step 4.

.. note::

   Please ask the Niflheim's support staff to verify that gpaw-python runs single-threaded, e.g. for a job running on ``p024`` do from ``audhumbla``::

     ssh p024 ps -fL

   Numbers higher then **1** in the **NLWP** column mean multi-threaded job.

   In case of openmpi it is necessary to set the :envvar:`OMP_NUM_THREADS` variable::

     setenv OMP_NUM_THREADS 1 # [t]csh
     export OMP_NUM_THREADS=1 # [ba]sh

.. note::

   When setting any environment variables please do **not**
   overwrite the system default :envvar:`PATH`, :envvar:`PYTHONPATH`,
   nor :envvar:`GPAW_SETUP_PATH` environment variables.
   When setting the environment variables **prepend** them, i.e.:

   - using csh/tcsh::

       setenv PATH ${HOME}/bin:${PATH}

   - using bash::

       export PATH=${HOME}/bin:${PATH}
