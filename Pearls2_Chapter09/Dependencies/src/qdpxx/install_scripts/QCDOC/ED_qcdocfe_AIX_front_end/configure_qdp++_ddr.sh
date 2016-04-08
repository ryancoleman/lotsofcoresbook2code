# Do the configuration
#
# Source ./setup.sh first
#

../qdp++/configure --enable-parallel-arch=parscalar \
	QMP_CFLAGS="-I$QOS/quser/include" \
        QMP_LDFLAGS="-L$QOS/quser" \
	QMP_LIBS="-lqcd_api" \
	CXXFLAGS="-O2 -finline-limit=50000" CFLAGS="-O2" \
	--enable-precision=single \
        --disable-qmp-route \
	--with-libxml2=$LIBXML \
	--host=powerpc-gnu-elf --build=none \
	--prefix=/home/ed/bj/Devel/SciDAC/install/qdp_parscalar_single_ddr
