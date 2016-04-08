#!/bin/sh


../configure  --prefix=/home/edwards/arch/qdp++/parscalar-gm-intel --with-qmp=/usr/local/qmp/mpich-gm-1.2.5 --enable-parallel-arch=parscalar CC=gcc CXX=/opt/intel_cc_80/bin/icc CXXFLAGS="-O" LDFLAGS="-static" 

# --enable-profiling
