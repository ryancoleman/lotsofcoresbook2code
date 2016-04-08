#!/bin/sh

#../configure --prefix=/usr/local/qdp++/parscalar-gigE --with-qmp=/usr/local/qmp/mvia --enable-parallel-arch=parscalar --enable-sse2 CFLAGS="-std=c99 -pedantic" CXXFLAGS="-O2 -finline-limit=50000 -march=pentium4 -DUSE_REMOTE_QIO" LDFLAGS="-static -L/home/edwards/arch/lib" LIBS="-lremote -lunp_intel"

../configure --prefix=/usr/local/qdp++/parscalar-gigE --with-qmp=/usr/local/qmp/mvia --enable-parallel-arch=parscalar --enable-sse2 CFLAGS="-std=c99 -pedantic" CXXFLAGS="-O2 -finline-limit=50000 -march=pentium4" LDFLAGS="-static"
