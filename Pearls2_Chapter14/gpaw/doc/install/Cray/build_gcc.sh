#!/bin/bash -x

export CC=cc
export CXX=g++
export MPICC=cc
export LINKFORSHARED='-Wl,-export-dynamic -dynamic'
export MPI_LINKFORSHARED='-Wl,-export-dynamic -dynamic'

install_dir='/some_path/scalable-python-gcc'
./configure --prefix=$install_dir --enable-mpi --disable-ipv6 2>&1 | tee loki-conf

module swap craype-interlagos craype-istanbul
module list
make 2>&1 | tee log-make
make install 2>&1 | tee log-inst

make clean

module swap craype-istanbul craype-interlagos
make mpi 2>&1 | tee log-make-mpi
cp ${install_dir}/bin/python .
make install-mpi 2>&1 | tee log-inst-mpi
