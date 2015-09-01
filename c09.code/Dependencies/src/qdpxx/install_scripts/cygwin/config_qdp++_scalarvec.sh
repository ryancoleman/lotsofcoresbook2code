#!/bin/sh

../configure --prefix=/usr/local/share/qdp++/scalarvec --enable-parallel-arch=scalarvec --enable-sse2 CXXFLAGS="-O2 -finline-limit=50000 -march=pentium4" 
