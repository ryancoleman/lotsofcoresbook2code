#!/bin/sh

../configure --prefix=/usr/local/share/qdp++/parscalar-mpi --with-qmp=/usr/local/share/qmp/mpi --enable-parallel-arch=parscalar  --enable-sse2 CXXFLAGS="-O2 -finline-limit=50000 -march=pentium4" 
