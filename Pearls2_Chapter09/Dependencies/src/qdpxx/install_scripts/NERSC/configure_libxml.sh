#!/usr/bin/bash

# Change here to set your own installation path
INSTALL="$HOME/Devel/SciDAC/install/libxml2-2.6.6-generic"

# No further changes below here
export CC=xlc_r
export CFLAGS="-O2"
./configure --prefix=${INSTALL} --disable-shared \
        --without-python \
        --without-readline \
        --without-threads \
        --without-history \
        --without-ftp \
        --without-http \
        --without-catalog \
        --without-docbook \
        --without-schemas \
        --without-zlib    \
        --without-iconv

