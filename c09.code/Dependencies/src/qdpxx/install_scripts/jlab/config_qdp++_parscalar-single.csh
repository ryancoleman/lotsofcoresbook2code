#!/bin/sh

../configure --prefix=$HOME/arch/qdp++/parscalar-single --with-qmp=/usr/local/qmp/qmp2-1-3/single --enable-parallel-arch=parscalar  --enable-sse2 CXXFLAGS='-O2 -finline-limit=50000 -march=pentium4' 
