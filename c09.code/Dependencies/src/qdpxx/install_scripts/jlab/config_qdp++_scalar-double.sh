#!/bin/sh

../configure --prefix=/usr/local/qdp++/scalar-double --enable-parallel-arch=scalar --enable-sse2 --enable-precision=double CXXFLAGS="-O2 -finline-limit=50000 -msse -march=pentium4" LDFLAGS=-static
