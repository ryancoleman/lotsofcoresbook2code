#./configure --prefix=$HOME/scidac/build/qio --with-qmp=$QOS/quser \
# --host=qcdoc --build=none QMP_CFLAGS="-I$QOS/quser/include" \
# QMP_LDFAGS="-L$QOS/quser -L$LIBDIR -Xlinker" 
# QMP_LIBS="-lqcd_api" CC="env GCC_EXEC_PREFIX=$SPECDIR ${XPREFIX}gcc"

QOS=/qcdoc/sfw/qos/devel/QMP_v2/qcdochost0

../../qdp++/configure --enable-parallel-arch=parscalar \
	QMP_CFLAGS="-I$QOS/quser/include" \
        QMP_LDFLAGS="-L$QOS/quser" \
	QMP_LIBS="-lqcd_api" \
	CXXFLAGS="-O2 -finline-limit=50000" CFLAGS="-O2" \
	--enable-precision=single \
        --enable-qcdoc-edram \
        --disable-qmp-route \
	--with-libxml2=/home/bj/install/QOS_QMP_v2/libxml-2.6.6 \
	--host=powerpc-gnu-elf --build=sparc-sun-solaris2.9 \
	--prefix=/home/bj/install/qdp_parscalar_single
