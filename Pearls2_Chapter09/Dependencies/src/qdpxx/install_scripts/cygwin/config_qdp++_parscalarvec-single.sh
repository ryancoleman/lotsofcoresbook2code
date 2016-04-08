#!/bin/tcsh

../configure --prefix=/usr/local/share/qdp++/parscalarvec-single --with-qmp=/usr/local/share/qmp/single --enable-parallel-arch=parscalarvec  --enable-sse2 CXXFLAGS="-O2 -finline-limit=50000 -march=pentium4" 
