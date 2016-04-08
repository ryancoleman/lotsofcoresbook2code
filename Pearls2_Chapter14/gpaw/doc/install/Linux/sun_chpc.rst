.. _sun_chpc:

==============
sun.chpc.ac.za
==============

Here you find information about the the system http://www.chpc.ac.za/sun/.

The installation of user's packages on SUSE 10 **login01**,
64-bit described below uses
`modules <http://modules.sourceforge.net/>`_, and assumes ``sh`` shell:

- packages are installed under ``~/CAMd``::

   mkdir ~/CAMd
   cd ~/CAMd

- module files are located under ``~/CAMd/modulefiles``::

   mkdir ~/CAMd/modulefiles

- download the :svn:`~doc/install/Linux/customize_sun_chpc_SUSE10.py` file:

  .. literalinclude:: customize_sun_chpc_SUSE10.py

- download packages with :svn:`~doc/install/Linux/download_sun_chpc.sh`:

  .. literalinclude:: download_sun_chpc.sh

- install packages, deploy modules and test with :svn:`~doc/install/Linux/install_sun_chpc_SUSE10.sh`:

  .. literalinclude:: install_sun_chpc_SUSE10.sh

  **Note** that every time you wish to install a new version of a package,
  and deploy new module file, better keep the old module file.

- submit the test job (jobs must be submitted from under *~/scratch*)::

   mqsub msub_sun_chpc.sh

  using the following :file:`msub_sun_chpc.sh`:

  .. literalinclude:: msub_sun_chpc.sh
