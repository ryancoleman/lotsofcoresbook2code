.. _developer_installation:

======================
Developer installation
======================

The :ref:`installationguide_standard` will copy all the
Python files to the standard place for Python modules (something like
:file:`/usr/lib/python2.5/site-packages/gpaw`) or to
:file:`{<my-directory>}/lib/python/gpaw` if you used the
:file:`--home={<my-directory>}` option.  As a developer, you will want
Python to use the files from the svn checkout that you are hacking on.

Do the following:

  * Checkout the :ref:`latest_development_release` (if you
    need access to svn updates)
    or download a tarball of :ref:`latest_stable_release`.

    .. note::

       The instructions assume :envvar:`GPAW_HOME` is set to *${HOME}/gpaw*.

  * Go to the :file:`gpaw` directory::

     [~]$ cd gpaw

  * build :ref:`c_extension`::

     [gpaw]$ python setup.py build_ext 2>&1 | tee build_ext.log

    This will build two things:

    * :file:`_gpaw.so`:  A shared library for serial calculations containing
      GPAW's C-extension module.  The module will be in
      :file:`${GPAW_HOME}/build/lib.{<platform>}-{<python_version>}/`.

      .. note::

         Here is how to find the *<platform>* and *<python_version>* variables::
       
           python -c "from distutils import util; print util.get_platform()"
           python -c "from distutils import sysconfig; print sysconfig.get_python_version()"

      For example *<platform>* could be *linux-x86_64*, and
      *<python_version>* could be *2.5*.

    * :file:`gpaw-python`: A special Python interpreter for parallel
      calculations.  The interpreter has GPAW's C-code built in.  The
      :file:`gpaw-python` executable is located
      in :file:`${GPAW_HOME}/build/bin.{<platform>}-{<python_version>}/`.

      .. note::

         The :file:`gpaw-python` interpreter will be built only if
         :file:`setup.py` finds an ``mpicc`` compiler.

      See :ref:`parallel_installation` for more details about parallel runs.

  * Prepend :file:`${GPAW_HOME}` and :file:`${GPAW_HOME}/build/lib.{<platform>}-{<python_version>}/`
    onto your :envvar:`PYTHONPATH` and
    :file:`${GPAW_HOME}/build/bin.{<platform>}-{<python_version>}:${GPAW_HOME}/tools` onto
    :envvar:`PATH`, e.g. put into :file:`~/.tcshrc`::

     setenv GPAW_HOME ${HOME}/gpaw
     setenv GPAW_PLATFORM `python -c "from distutils import util, sysconfig; print util.get_platform()+'-'+sysconfig.get_python_version()"`
     setenv PYTHONPATH ${GPAW_HOME}:${PYTHONPATH}
     setenv PYTHONPATH ${GPAW_HOME}/build/lib.${GPAW_PLATFORM}:${PYTHONPATH}
     setenv PATH ${GPAW_HOME}/build/bin.${GPAW_PLATFORM}:${GPAW_HOME}/tools:${PATH}

    or if you use bash, put these lines into :file:`~/.bashrc`::

     export GPAW_HOME=${HOME}/gpaw
     export GPAW_PLATFORM=`python -c "from distutils import util, sysconfig; print util.get_platform()+'-'+sysconfig.get_python_version()"`
     export PYTHONPATH=${GPAW_HOME}:${PYTHONPATH}
     export PYTHONPATH=${GPAW_HOME}/build/lib.${GPAW_PLATFORM}:${PYTHONPATH}
     export PATH=${GPAW_HOME}/build/bin.${GPAW_PLATFORM}:${GPAW_HOME}/tools:${PATH}

  * continue on :ref:`installationguide_developer`.
