#!/bin/sh

../configure --prefix=/usr/local/share/qdp++/scalar --enable-parallel-arch=scalar --enable-sse2 CXXFLAGS="-O2 -finline-limit=50000 -march=pentium4" 
