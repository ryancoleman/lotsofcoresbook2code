#!/bin/sh

export APPS="/home/firegam/CAMd"
export MODULEFILES="${APPS}/modulefiles"

INSTALL_DACAPO=False

# download packages
openmpi_version=1.4.3
wget http://www.open-mpi.org/software/ompi/v1.4/downloads/openmpi-${openmpi_version}.tar.bz2
nose_version=0.11.3
wget http://python-nose.googlecode.com/files/nose-${nose_version}.tar.gz
numpy_version=1.1.1
wget http://downloads.sourceforge.net/numpy/numpy-${numpy_version}.tar.gz
numpy_version=1.5.0
wget http://downloads.sourceforge.net/numpy/numpy-${numpy_version}.tar.gz
ase_version=3.4.1.1765
wget https://wiki.fysik.dtu.dk/ase-files/python-ase-${ase_version}.tar.gz
gpaw_version=0.7.2.6974
wget https://wiki.fysik.dtu.dk/gpaw-files/gpaw-${gpaw_version}.tar.gz
gpaw_setups_version=0.6.6300
wget http://wiki.fysik.dtu.dk/gpaw-files/gpaw-setups-${gpaw_setups_version}.tar.gz
if [ "${INSTALL_DACAPO}" = "True" ];
    then
    # dacapo
    dacapo_version=2.7.16
    wget "https://wiki.fysik.dtu.dk/dacapo/Installation?action=AttachFile&do=get&target=campos-dacapo-${dacapo_version}.tar.gz" -O campos-dacapo-${dacapo_version}.tar.gz
    dacapo_psp_version=1
    wget "https://wiki.fysik.dtu.dk/dacapo/Installation?action=AttachFile&do=get&target=campos-dacapo-pseudopotentials-${dacapo_psp_version}.tar.gz" -O campos-dacapo-pseudopotentials-${dacapo_psp_version}.tar.gz
fi
