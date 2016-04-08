#!/bin/sh

#../configure --prefix=/usr/local/qdp++/parscalar-gm-double --with-qmp=/usr/local/qmp/mpich-gm-1.2.5 --enable-parallel-arch=parscalar --enable-sse2 --enable-precision=double CXXFLAGS="-O2 -finline-limit=50000 -march=pentium4 -DUSE_REMOTE_QIO" LDFLAGS="-static -L/home/edwards/arch/lib" LIBS="-lremote -lunp_intel"

../configure --prefix=/usr/local/qdp++/parscalar-gm-double --with-qmp=/usr/local/qmp/mpich-gm-1.2.5 --enable-parallel-arch=parscalar --enable-sse2 --enable-precision=double CXXFLAGS="-O2 -finline-limit=50000 -march=pentium4" LDFLAGS="-static"
