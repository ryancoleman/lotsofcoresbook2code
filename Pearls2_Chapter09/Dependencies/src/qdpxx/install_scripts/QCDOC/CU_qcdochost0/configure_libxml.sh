#!/usr/local/bin/bash

QOS=/qcdoc/sfw/qos/devel/QMP_v2/qcdochost0
XPREFIX=/qcdoc/local/gcc-3.2.2/bin/powerpc-gnu-elf-
SPECDIR=$QOS/quser/gcc-lib-user///
LIBDIR=$QOS/quser/gcc-lib-user///

. ${QOS}/scripts/setup.sh

./configure --prefix=/home/bj/install/QOS_QMP_v2/libxml-2.6.6 \
	CFLAGS="-O2" \
	--host=powerpc-gnu-elf \
	--build=sparc-sun-solaris2.9 \
	--without-zlib \
	--without-python \
	--without-readline \
	--without-threads \
	--without-history \
	--with-output \
	--without-writer \
	--without-reader \
	--without-ftp \
	--without-http \
	--without-pattern \
	--without-catalog \
	--without-docbook \
	--without-iconv \
	--without-schemas
