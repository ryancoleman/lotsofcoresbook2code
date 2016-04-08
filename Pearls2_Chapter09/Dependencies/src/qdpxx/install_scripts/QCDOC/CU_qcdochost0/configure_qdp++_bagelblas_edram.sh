#! /usr/local/bin/bash

../../qdp++/configure --enable-parallel-arch=parscalar \
	QMP_CFLAGS="-I$QOS/quser/include" \
        QMP_LDFLAGS="-L$QOS/quser" \
	QMP_LIBS="-lqcd_api" \
	CXXFLAGS="-I/home/bj/install/QOS_QMP_v2/bagel-1.3.2-single -O2 -finline-limit=50000" CFLAGS="-O2" \
        LDFLAGS="-L/home/bj/install/QOS_QMP_v2/bagel-1.3.2-single" \
	LIBS="-lwfmppc440s" \
	--enable-precision=single \
	--enable-qcdoc \
	--enable-qcdoc-edram \
        --disable-qmp-route \
	--with-libxml2=/home/bj/install/QOS_QMP_v2/libxml-2.6.6 \
	--host=powerpc-gnu-elf --build=sparc-sun-solaris2.9 \
	--prefix=/home/bj/install/QOS_QMP_v2/qdp_parscalar_single_blas
