#!/bin/bash

export CC=cc
export CXX=g++
export MPICC=cc
# export XTPE_LINK_TYPE=dynamic
export LINKFORSHARED='-Wl,-export-dynamic -dynamic'
export MPI_LINKFORSHARED='-Wl,-export-dynamic -dynamic'

install_prefix=/appl/opt/python/scalable-gnu

# Make zlib built-in to the interpreter
sed -i -e 's/^#zlib.*/zlib zlibmodule.c -I\/usr\/include -L\/usr\/lib64 -lz/' Modules/Setup.dist


./configure --prefix=$install_prefix --enable-mpi --disable-ipv6 2>&1 | tee log-conf

make 2>&1 | tee log-make
make install 2>&1 | tee log-inst

make mpi 2>&1 | tee log-make-mpi
make install-mpi 2>&1 | tee log-inst-mpi
