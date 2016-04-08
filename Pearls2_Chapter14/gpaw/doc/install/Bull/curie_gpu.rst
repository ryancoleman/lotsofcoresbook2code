.. _curie_gpu:

====================
curie.ccc.cea.fr GPU
====================

Here you find information about the the system
http://www-hpc.cea.fr/en/complexe/tgcc-curie.htm.

**Warning**: May 14 2013: rpa-gpu-expt branch fails to compile
due to cublasZdgmm missing in cuda-4.2.

The system runs Bull Linux.
The installation assumes *bash* shell:

- packages are installed under ``~/CAMd``::

   mkdir ~/CAMd
   cd ~/CAMd

- module files are located under ``~/CAMd/modulefiles``::

   mkdir ~/CAMd/modulefiles

- download the :svn:`~doc/install/Bull/customize_curie_gpu.py` file:

  .. literalinclude:: customize_curie_gpu.py

- download packages with :svn:`~doc/install/Bull/download_curie_gpu.sh`:

  .. literalinclude:: download_curie_gpu.sh

- install packages, deploy modules and test with
  :svn:`~doc/install/Bull/install_curie_gpu.sh`:

  .. literalinclude:: install_curie_gpu.sh

  **Note** that every time you wish to install a new version of a package,
  and deploy new module file, better keep the old module file.

- submit a test job::

   ccc_msub msub_curie_gpu.sh

  using the following :file:`msub_curie_gpu.sh`:

  .. literalinclude:: msub_curie_gpu.sh
