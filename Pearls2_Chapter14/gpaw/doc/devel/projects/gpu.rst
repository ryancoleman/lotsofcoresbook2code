===
GPU
===

Aalto Effort
============

(Samuli Hakala, Ville Havu, Jussi Enkovaara (CSC) )
We have been implementing the most performance critical C-kernels
in the finite-difference mode to utilize GPUs. Implementation is done
using CUDA and cuBLAS libraries when possible, Python interface to CUDA,
PyCUDA is also utilized. Code supports using multiple GPUs with MPI.
First tests indicate promising speedups compared
to CPU-only code, and we are hoping to test larger systems (where
the benefits are expected to be larger) soon. Currently, we are extending the
CUDA implementation to real-time TDDFT.

Code is not in full production level yet.

Stanford/SUNCAT Effort
======================

(Lin Li, Jun Yan, Christopher O'Grady)

We have created a GPU version of the GPAW RPA code.  We observe a
speedup of 40 (GPU vs. CPU+GPU) making GPUs about a factor of 5 more
cost effective for this calculation.  It is available in the svn branch
rpa-gpu-expt.

Installation of cuda on Fedora 18 x86_64
========================================

First, cuda_5.0.35 does not support gcc version 4.7 and up,
either kernel 3.7 and up https://bugzilla.redhat.com/show_bug.cgi?format=multiple&id=902045, however there exist workarounds.

Proceed as follows:

0. boot a 3.6 kernel

1. yum -y install wget make gcc-c++ freeglut-devel libXi-devel libXmu-devel mesa-libGLU-devel

2. configure rpmfusion-nonfree http://rpmfusion.org/Configuration

3. yum -y install xorg-x11-drv-nvidia-libs

4. disable X::

     rm /etc/systemd/system/default.target&& ln -s /lib/systemd/system/multi-user.target /etc/systemd/system/default.target

5. select a 3.6 version (see https://bugzilla.redhat.com/show_bug.cgi?format=multiple&id=902045 ) of kernel in ``/etc/grub2.cfg`` and make sure the *linux* line contains::

     nouveau.modeset=0 rd.driver.blacklist=nouveau

   See http://www.pimp-my-rig.com/2012/01/unload-nouveau-install-nvidia-driver.html

6. disable kernel and nvidia updates in ``/etc/yum.conf``::

     exclude=kernel* *nvidia*

7. reboot

8. download and install cuda::

     wget http://developer.download.nvidia.com/compute/cuda/5_0/rel-update-1/installers/cuda_5.0.35_linux_64_fedora16-1.run
     sh cuda_5.0.35_linux_64_fedora16-1.run -override compiler

   Keep the recommended installation paths (``/usr/local/cuda-5.0``, ...),
   and after the installation create a link::

     ln -s /usr/local/cuda /opt/cuda

9. add to ``~/bashrc``::

     export PATH=/opt/cuda/bin:${PATH}
     export LD_LIBRARY_PATH=/opt/cuda/lib64:$LD_LIBRARY_PATH

   and source it.

10. convert error into a warning in ``/usr/local/cuda-5.0/include/host_config.h``::

      // #error — unsupported GNU version! gcc 4.7 and up are not supported!
      #warning — unsupported GNU version! gcc 4.7 and up are not supported!


    See https://www.udacity.com/wiki/cs344/troubleshoot_gcc47

11. test::

      cd NVIDIA_CUDA-5.0_Samples/1_Utilities/deviceQuery
      make&& ./deviceQuery

12. if X does not run the ``/dev/nvidia*`` will not be created.
    The solution seems to be creating them by hand (``NvidiaDevCreator.sh``)::

      #!/bin/bash

      /sbin/modprobe nvidia

      if [ "$?" -eq 0 ]; then
        # Count the number of NVIDIA controllers found.
        NVDEVS=`lspci | grep -i NVIDIA`
        N3D=`echo "$NVDEVS" | grep "3D controller" | wc -l`
        NVGA=`echo "$NVDEVS" | grep "VGA compatible controller" | wc -l`

        N=`expr $N3D + $NVGA - 1`
        for i in `seq 0 $N`; do
          mknod -m 666 /dev/nvidia$i c 195 $i
        done

        mknod -m 666 /dev/nvidiactl c 195 255

      else
        exit 1
      fi

    Or switch back to X::

      cd /etc/systemd/system
      ln -s /lib/systemd/system/graphical.target default.target

13. install http://mathema.tician.de/software/pycuda (needed only for *Aalto Effort*)::

      cd&& git clone https://github.com/inducer/pycuda.git
      cd ~/pycuda
      PATH=$PATH:/opt/cuda/bin LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/cuda/lib64 python configure.py --update-user --boost-compiler=gcc
      /bin/mv -f ~/.aksetup-defaults.py siteconf.py
      sed -i "s/boost_python-py27/boost_python/" siteconf.py
      sed -i 's/boost_thread/boost_thread-mt/' siteconf.py
      sed -i "s#'\${CUDA_ROOT}/lib', ##" siteconf.py
      PATH=$PATH:/opt/cuda/bin LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/cuda/lib64 python setup.py install --root=~/pycuda-fc18-1

    and add to ``~/.bashrc``::

      export PYTHONPATH=~/pycuda-fc18-1/usr/lib64/python2.7/site-packages:${PYTHONPATH}

    May 7 2013: note that ``compyte`` has been removed from ``pucuda``,
    but the source of ``pucuda`` does not reflect that.
    Therefore ``git clone https://github.com/inducer/compyte.git``
    and create a link under ``pucuda`` install tree.
    In addition https://pypi.python.org/pypi/pytools,
    https://pypi.python.org/pypi/py,
    https://pypi.python.org/pypi/pytest are required by ``pucuda``.
