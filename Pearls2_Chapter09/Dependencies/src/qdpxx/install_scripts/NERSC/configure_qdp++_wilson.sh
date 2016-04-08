#! /usr/bin/bash

# The Root Install Directory 
export ROOT=/u1/bjoo/Devel/SciDAC/install

# The QDP target installation directory 
export QDP_HOME=$ROOT/qdp_parscalar_wilson

# The directory where libxml is installed
export LIBXML_HOME=$ROOT/libxml2-2.6.6-generic

# The directory where qmp is installed
export QMP_HOME=$ROOT/qmp_generic

# Make the template includes directroy -- Stupid XLC isms
if [ ! -d $ROOT/tempinc ];
then 
	echo $ROOT/tempinc does not exist -- creating
	mkdir -p $ROOT/tempinc
else
	echo $ROOT/tempinc exists - cleaning
	pushd $ROOT/tempinc
	rm -rf *
	popd
fi

./configure CC="mpcc_r" CXX="mpCC_r" CFLAGS="-O2 -qstrict -qfuncsect" \
	CXXFLAGS="-qarch=pwr3 -qtune=pwr3 -O2 -qtempinc=$ROOT/tempinc" \
	--prefix=$QDP_HOME \
	--enable-parallel-arch=parscalar \
	--enable-precision=single \
	--with-libxml2=$LIBXML_HOME \
	--with-qmp=$QMP_HOME \
	--host=powerpc-ibm-aix --build=none
