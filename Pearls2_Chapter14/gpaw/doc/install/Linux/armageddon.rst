.. _armageddon:

===========================
armageddon.chimfar.unimo.it
===========================

The installation of user's packages described below assumes *bash* shell:

- packages are installed under ``~/CAMd``::

   mkdir ~/CAMd
   cd ~/CAMd

- download the :svn:`~doc/install/Linux/customize_armageddon.py`
  and :svn:`~doc/install/Linux/set_env_armageddon.sh` files::

   wget https://svn.fysik.dtu.dk/projects/gpaw/trunk/doc/install/Linux/customize_armageddon.py
   wget https://svn.fysik.dtu.dk/projects/gpaw/trunk/doc/install/Linux/set_env_armageddon.sh

  .. literalinclude:: customize_armageddon.py

- download packages with :svn:`~doc/install/Linux/download_armageddon.sh`,
  buy running ``sh download_armageddon.sh``:

  .. literalinclude:: download_armageddon.sh

- install packages and test with :svn:`~doc/install/Linux/install_armageddon.sh`,
  buy running ``sh install_armageddon.sh``.

- enable packages with :svn:`~doc/install/Linux/set_env_armageddon.sh`,
  buy running ``. set_env_armageddon.sh`` (put this in *~/.bashrc*).
