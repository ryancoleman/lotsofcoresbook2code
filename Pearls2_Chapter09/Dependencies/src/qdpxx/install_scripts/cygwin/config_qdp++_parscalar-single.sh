#!/bin/sh

../configure --prefix=/usr/local/share/qdp++/parscalar-single --with-qmp=/usr/local/share/qmp/single --enable-parallel-arch=parscalar  --enable-sse2 CXXFLAGS="-O2 -finline-limit=50000 -march=pentium4" 
