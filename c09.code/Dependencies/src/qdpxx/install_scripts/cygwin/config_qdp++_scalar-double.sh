#!/bin/sh

../configure --prefix=/usr/local/share/qdp++/scalar-double --enable-parallel-arch=scalar --enable-sse2 --enable-precision=double CFLAGS="-O2 -msse -msse2 -march=pentium4" CXXFLAGS="-O2 -finline-limit=50000 -msse -march=pentium4" 
