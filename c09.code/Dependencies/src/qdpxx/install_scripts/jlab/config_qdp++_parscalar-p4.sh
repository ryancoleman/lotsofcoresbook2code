#!/bin/sh


../configure --prefix=/usr/local/qdp++/parscalar-p4 --with-qmp=/usr/local/qmp/mpich-p4-1.2.6 --enable-parallel-arch=parscalar --enable-sse2 CFLAGS="-std=c99 -pedantic" CXXFLAGS="-O2 -finline-limit=50000 -march=pentium4" 

# LDFLAGS="-static"

