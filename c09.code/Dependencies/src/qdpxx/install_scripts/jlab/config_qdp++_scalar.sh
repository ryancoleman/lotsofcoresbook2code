#!/bin/sh

../configure --prefix=/usr/local/qdp++/scalar --enable-parallel-arch=scalar --enable-sse2 CXXFLAGS="-O2 -finline-limit=50000 -msse -march=pentium4" LDFLAGS=-static
