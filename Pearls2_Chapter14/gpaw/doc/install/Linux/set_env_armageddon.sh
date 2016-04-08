#!/bin/sh

export APPS="/home/firegam/CAMd"

INSTALL_DACAPO=False

# openmpi
openmpi_version=1.4.3

PA=${APPS}/openmpi-${openmpi_version}-1/bin

case $PATH in
  *${PA}*)      ;;
  *?*)  PATH=${PA}:${PATH} ;;
  *)    PATH=${PA} ;;
esac
export PATH

PY=${APPS}/openmpi-${openmpi_version}-1/lib

case $LD_LIBRARY_PATH in
  *${PY}*)      ;;
  *?*)  LD_LIBRARY_PATH=${PY}:${LD_LIBRARY_PATH} ;;
  *)    LD_LIBRARY_PATH=${PY} ;;
esac
export LD_LIBRARY_PATH

# nose
nose_version=0.11.3

PA=${APPS}/nose-${nose_version}-1/usr/local/bin

case $PATH in
  *${PA}*)      ;;
  *?*)  PATH=${PA}:${PATH} ;;
  *)    PATH=${PA} ;;
esac
export PATH

PY=${APPS}/nose-${nose_version}-1/usr/local/lib/python2.6/dist-packages

case $PYTHONPATH in
  *${PY}*)      ;;
  *?*)  PYTHONPATH=${PY}:${PYTHONPATH} ;;
  *)    PYTHONPATH=${PY} ;;
esac
export PYTHONPATH

# numpy
numpy_version=1.5.0

PA=${APPS}/numpy-${numpy_version}-1/usr/local/bin

case $PATH in
  *${PA}*)      ;;
  *?*)  PATH=${PA}:${PATH} ;;
  *)    PATH=${PA} ;;
esac
export PATH

PY=${APPS}/numpy-${numpy_version}-1/usr/local/lib/python2.6/dist-packages

case $PYTHONPATH in
  *${PY}*)      ;;
  *?*)  PYTHONPATH=${PY}:${PYTHONPATH} ;;
  *)    PYTHONPATH=${PY} ;;
esac
export PYTHONPATH

# campos-ase3
ase_version=3.4.1.1765

PA=${APPS}/python-ase-${ase_version}/tools

case $PATH in
  *${PA}*)      ;;
  *?*)  PATH=${PA}:${PATH} ;;
  *)    PATH=${PA} ;;
esac
export PATH

PY=${APPS}/python-ase-${ase_version}/

case $PYTHONPATH in
  *${PY}*)      ;;
  *?*)  PYTHONPATH=${PY}:${PYTHONPATH} ;;
  *)    PYTHONPATH=${PY} ;;
esac
export PYTHONPATH

# campos-gpaw-setups
gpaw_setups_version=0.6.6300

PA=${APPS}/gpaw-setups-${gpaw_setups_version}

case $GPAW_SETUP_PATH in
  *${PA}*)      ;;
  *?*)  GPAW_SETUP_PATH=${PA}:${GPAW_SETUP_PATH} ;;
  *)    GPAW_SETUP_PATH=${PA} ;;
esac
export GPAW_SETUP_PATH

# campos-gpaw
gpaw_version=0.7.2.6974

PA=${APPS}/gpaw-${gpaw_version}/tools

case $PATH in
  *${PA}*)      ;;
  *?*)  PATH=${PA}:${PATH} ;;
  *)    PATH=${PA} ;;
esac
export PATH

PA=${APPS}/gpaw-${gpaw_version}/build/bin.linux-x86_64-2.6/

case $PATH in
  *${PA}*)      ;;
  *?*)  PATH=${PA}:${PATH} ;;
  *)    PATH=${PA} ;;
esac
export PATH

PY=${APPS}/gpaw-${gpaw_version}/

case $PYTHONPATH in
  *${PY}*)      ;;
  *?*)  PYTHONPATH=${PY}:${PYTHONPATH} ;;
  *)    PYTHONPATH=${PY} ;;
esac
export PYTHONPATH

PY=${APPS}/gpaw-${gpaw_version}/build/lib.linux-x86_64-2.6/

case $PYTHONPATH in
  *${PY}*)      ;;
  *?*)  PYTHONPATH=${PY}:${PYTHONPATH} ;;
  *)    PYTHONPATH=${PY} ;;
esac
export PYTHONPATH

if [ "${INSTALL_DACAPO}" = "True" ];
then
    # dacapo-pseudopotentials
    dacapo_psp_version=1

    PA=${APPS}/psp-${dacapo_psp_version}

    case $DACAPOPATH in
        *${PA}*)      ;;
        *?*)  DACAPOPATH=${PA}:${DACAPOPATH} ;;
        *)    DACAPOPATH=${PA} ;;
    esac
    export DACAPOPATH

    # dacapo
    dacapo_version=2.7.16

    PA=${APPS}/campos-dacapo-${dacapo_version}/src/gfortran_fnosecond_underscore_serial/dacapo_${dacapo_version}-1_serial.run

    case $DACAPOEXE_SERIAL in
        *${PA}*)      ;;
        *?*)  DACAPOEXE_SERIAL=${PA}:${DACAPOEXE_SERIAL} ;;
        *)    DACAPOEXE_SERIAL=${PA} ;;
    esac
    export DACAPOEXE_SERIAL

    PA=${APPS}/campos-dacapo-${dacapo_version}/src/gfortran_fnosecond_underscore_mpi/dacapo_${dacapo_version}-1_mpi.run

    case $DACAPOEXE_PARALLEL in
        *${PA}*)      ;;
        *?*)  DACAPOEXE_PARALLEL=${PA}:${DACAPOEXE_PARALLEL} ;;
        *)    DACAPOEXE_PARALLEL=${PA} ;;
    esac
    export DACAPOEXE_PARALLEL

    PA=${APPS}/campos-dacapo-${dacapo_version}/

    case $PATH in
        *${PA}*)      ;;
        *?*)  PATH=${PA}:${PATH} ;;
        *)    PATH=${PA} ;;
    esac
    export PATH

fi
