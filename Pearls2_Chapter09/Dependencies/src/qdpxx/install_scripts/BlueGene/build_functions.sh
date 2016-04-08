#!/usr/bin/bash

# Build QMP
function build_qmp { 
	qmpdir=$1
	qmp_installdir=$2

	echo Removing QMP Installdir: ${qmp_installdir}
	rm -rf ${qmp_installdir}

	qmp_builddir="./build_qmp"
	rm -rf ${qmp_builddir}
	mkdir -p ${qmp_builddir}
	pushd ${qmp_builddir}
	export PATH=/bgl/BlueLight/DRV521_2004-050113/ppc/blrts-gnu/bin:$PATH
        export PATH=$HOME/bin:$PATH
	command="${qmpdir}/configure "
	command=${command}" CC=\"powerpc-bgl-blrts-gnu-gcc \" "
	command=${command}" RANLIB=\"powerpc-bgl-blrts-gnu-ranlib\" "
	command=${command}" --with-qmp-comms-type=MPI "
	command=${command}" --with-qmp-comms-cflags=\"-I/bgl/BlueLight/DRV521_2004-050113/ppc/bglsys/include\" "
	command=${command}" --with-qmp-comms-ldflags=\"-L/bgl/BlueLight/DRV521_2004-050113/ppc/bglsys/lib  -L/bgl/BlueLight/DRV521_2004-050113/ppc/bglsys/lib\" "
	command=${command}" --with-qmp-comms-libs=\"-lmpich.rts -lmsglayer.rts -ldevices.rts -lrts.rts -ldevices.rts -lrts.rts\" "
	command=${command}" --prefix=${qmp_installdir}"
	echo "Configure command is ${command}"
	echo ${command} > ./configure_qmp.sh
	chmod u+x configure_qmp.sh
	./configure_qmp.sh
	gmake 
	gmake install
	popd
}


# Build BAGEL 
function build_bagel { 
	bageldir=$1
	bagel_install_dir=$2

	echo building BAGEL in ${bageldir}, installing in ${bagel_install_dir}
	rm -rf ./build_bagel
	rm -rf ${bagel_install_dir}
	mkdir -p build_bagel
	pushd ./build_bagel
	export PATH=/usr/bin:$PATH
	${bageldir}/configure --prefix=${bagel_install_dir}
	gmake 
	gmake install
	popd
}

function build_bagel_wilson_dslash {
	wilsondir=$1
	wilson_install_dir=$2
	bagel_install_dir=$3
	bagel_prec=$4
	bagel_comm=$5
	bagel_alloc=$6
	bagel_proc=$7
	host_sys=$8
	build_sys=$9

	echo Building BAGEL WilsonDslash in ${wilson_install_dir}

	if test "X${bagel_comm}X" == "XqmpX";
	then 
	    # qos=${10};
	    # qmp_cxxflags="-I${qos}/quser/include";
	    # qmp_ldflags="-L${qos}/quser/qcd_api";
	    # qmp_libs="-lqmp";
	    qmp_dir=${10}
	    qmp_cxxflags=`${qmp_dir}/bin/qmp-config --cflags`
	    qmp_ldflags=`${qmp_dir}/bin/qmp-config --ldflags`
	    qmp_libs=`${qmp_dir}/bin/qmp-config --libs`
        else
	    qmp_cxxflags="";
	    qmp_ldflags="";
	    qmp_libs="";
        fi
	rm -rf ${wilson_install_dir}
	rm -rf ./build_wilson_dslash_${bagel_prec}
	mkdir -p ./build_wilson_dslash_${bagel_prec}
	pushd ./build_wilson_dslash_${bagel_prec}
	export PATH=/bgl/BlueLight/DRV521_2004-050113/ppc/blrts-gnu/bin:$PATH
	export PATH=/home/y01/bjoo/bin:$PATH

	echo ${host_sys}
	echo ${build_sys}
	command="${wilsondir}/configure"
	command=${command}" --prefix=${wilson_install_dir}"
	command=${command}" --with-bagel=${bagel_install_dir}"
	command=${command}" --enable-precision=${bagel_prec}"
	command=${command}" --enable-comms=${bagel_comm}"
	command=${command}" --enable-allocator=${bagel_alloc}"
	command=${command}" --enable-target-cpu=${bagel_proc}"
	command=${command}" --host=${host_sys}"
	command=${command}" --build=${build_sys}"
	command=${command}" CXX=\"powerpc-bgl-blrts-gnu-g++\" "
	command=${command}" CXXFLAGS=\" ${qmp_cxxflags}\" "
	command=${command}" CFLAGS=\"${qmp_cxxflags}\" "
	command=${command}" ASFLAGS=\" \" "
	command=${command}" LDFLAGS=\" ${qmp_ldflags}\" "
	command=${command}" LIBS=\" ${qmp_libs}  -ldevices.rts -lrts.rts\" "
	command=${command}" RANLIB=\"powerpc-bgl-blrts-gnu-ranlib\" "
	echo ${command} > ./configure_wilson.sh
	chmod u+x ./configure_wilson.sh
	./configure_wilson.sh
	gmake
	gmake install
	popd
}

function build_qdp {
	qdpdir=$1
	qdp_install_dir=$2
	qmp_dir=$3
	libxml=$4
	precision=$5
	do_edram=$6
	do_blas=$7
	host_sys=$8
	build_sys=$9

        qdp_base=`basename ${qdp_install_dir}`
        qdp_builddir="./build_"${qdp_base}

	if test "X${do_edram}X" == "XyesX"; 
	then
	    qcdoc_edram="--enable-qcdoc-edram"; 
	else
	    qcdoc_edram="";
	fi

        qmp_cxxflags=`${qmp_dir}/bin/qmp-config --cflags`
        qmp_ldflags=`${qmp_dir}/bin/qmp-config --ldflags`
        qmp_libs=`${qmp_dir}/bin/qmp-config --libs`

        
	echo ${do_blas}
	if test "X${do_blas}X" == "XyesX";
	then 
	     qcdoc="--enable-qcdoc";
	     wilsondir=${10};
	     WILSON_CXXFLAGS=`${wilsondir}/bin/wfm-config --cxxflags`;
	     WILSON_LDFLAGS=`${wilsondir}/bin/wfm-config --ldflags`;
	     WILSON_LIBS=`${wilsondir}/bin/wfm-config --libs`;
        else
	     qcdoc="";
	     WILSON_CXXFLAGS="";
	     WILSON_LDFLAGS="";
	     WILSON_LIBS="";
        fi

	rm -rf ${qdp_install_dir}
	qdp_base=`basename ${qdp_install_dir}`
	qdp_builddir="./build_"${qdp_base}
	rm -rf ${qdp_builddir}
	mkdir -p ${qdp_builddir}
	pushd ${qdp_builddir}
	export PATH=/bgl/BlueLight/DRV521_2004-050113/ppc/blrts-gnu/bin:$PATH
	export PATH=/home/y01/bj/bin:$PATH

	command="${qdpdir}/configure --enable-parallel-arch=parscalar "
	command=${command}"  CXX=powerpc-bgl-blrts-gnu-g++ CC=powerpc-bgl-blrts-gnu-gcc RANLIB=powerpc-bgl-blrts-gnu-ranlib"
	command=${command}"  --with-qmp=${qmp_dir} "
	command=${command}"  CXXFLAGS=\"-O2 -finline-limit=50000 ${WILSON_CXXFLAGS}\" "
	command=${command}" CFLAGS=\"-O2\" "
	command=${command}" LDFLAGS=\"${WILSON_LDFLAGS}\" "
	command=${command}" LIBS=\"${WILSON_LIBS} -ldevices.rts -lrts.rts\" " 
	command=${command}" --enable-precision=${precision}"
	command=${command}" --disable-qmp-route"
#	command=${command}" --enable-slow-route"
	command=${command}" --with-libxml2=${libxml}"
	command=${command}" --host=powerpc-bgl-blrts-gnu"
	command=${command}" "${qcdoc_edram}" "${qcdoc}
	command=${command}" --build=none --prefix=${qdp_install_dir}"
	echo Configuring QDP++ with command:
	echo    ${command}
 	echo ${command} > configure_qdp.sh

	## Pick up right ar
	export PATH=/home/ed/bj/bin:$PATH
	chmod u+x configure_qdp.sh
	./configure_qdp.sh
	gmake 
	gmake install
	popd
}

function build_chroma { 
	chromadir=$1
	chroma_install_dir=$2
	qdp_dir=$3
	host_sys=$4
	build_sys=$5
	do_pab_dslash=$6

	echo Chromadir: ${chromadir}
	echo Chroma_Install_dir: ${chroma_install_dir}
	echo QDP_Install_dir: ${qdp_dir}
	echo HOST: ${host_sys}
	echo BUILD: ${build_sys}
	echo DO Dslash: ${do_pab_dslash}
	if test "X${do_pab_dslash}X" == "XyesX";
	then
	   pab_dslash="--enable-pab-wilson-dslash=noarch";
	else
	   pab_dslash="";
	fi

	rm -rf ${chroma_install_dir}
	install_base_name=`basename ${chroma_install_dir}`
	builddir=./build_${install_base_name}
	mm -rf ${builddir}

	mkdir -p ${chroma_install_dir}
	mkdir -p ${builddir}
	pushd ${builddir}

        export PATH=/bgl/BlueLight/DRV521_2004-050113/ppc/blrts-gnu/bin:$PATH
	export PATH=/home/y01/bj/bin:$PATH

	command="${chromadir}/configure CXXFLAGS=\"\" LDFLAGS=\"\" "
	command=${command}" CXX=powerpc-bgl-blrts-gnu-g++ "
	command=${command}" CXXFLAGS=\"-I/bgl/BlueLight/DRV521_2004-050113/ppc/bglsys/include\" "
	command=${command}" RANLIB=powerpc-bgl-blrts-gnu-ranlib"
	command=${command}" LDFLAGS=\"-L/bgl/BlueLight/DRV521_2004-050113/ppc/bglsys/lib\" "
	command=${command}" LIBS=\"-lmpich.rts -lmsglayer.rts -ldevices.rts -lrts.rts -ldevices.rts -lrts.rts\" "
	command=${command}" --with-qdp=${qdp_dir} "
	command=${command}" --prefix=${chroma_install_dir}"
#	command=${command}" --host=${host_sys}"
	command=${command}" --build=${build_sys}"
	command=${command}" ${pab_dslash} "	
	echo Configure command is:
	echo ${command}
	echo ${command} > ./configure_chroma.sh
	chmod u+x ./configure_chroma.sh
	source ./configure_chroma.sh	
	gmake
	gmake install
	popd
}

function build_libxml { 
    libxmldir=$1
    libxml_install_dir=$2
    host_sys=$3
    build_sys=$4
    
    echo "Libxml Source Dir: ${libxmldir}"
    echo "Libxml Install Dir: ${libxml_install_dir}"

    builddir="./build_libxml"
    rm -rf ${builddir}
    mkdir -p ${builddir}

    pushd ${builddir}
    
    command="${libxmldir}/configure --prefix=${libxml_install_dir}"
    command=${command}" CC=\"/bgl/BlueLight/DRV521_2004-050113/ppc/blrts-gnu/bin/powerpc-bgl-blrts-gnu-gcc\" "
    command=${command}" AR=\"/bgl/BlueLight/DRV521_2004-050113/ppc/blrts-gnu/bin/powerpc-bgl-blrts-gnu-ar\" "
    command=${command}" RANLIB=\"/bgl/BlueLight/DRV521_2004-050113/ppc/blrts-gnu/bin/powerpc-bgl-blrts-gnu-ranlib\" "
    command=${command}" CFLAGS=\"-O2\" "
#    command=${command}" --host=${host_sys}"
#    command=${command}" --build=${build_sys}"
    command=${command}" --disable-shared"
    command=${command}" --without-zlib"
    command=${command}" --without-python"
    command=${command}" --without-readline"
    command=${command}" --without-threads"
    command=${command}" --without-history"
    command=${command}" --with-output"
    command=${command}" --without-writer"
    command=${command}" --without-reader"
    command=${command}" --without-ftp"
    command=${command}" --without-http"
    command=${command}" --without-pattern"
    command=${command}" --without-catalog"
    command=${command}" --without-docbook"
    command=${command}" --without-iconv"
    command=${command}" --without-schemas"
    
    echo Configure command is:
    echo ${command}
    echo ${command} > ./configure_libxml.sh
    chmod u+x ./configure_libxml.sh
    source ./configure_libxml.sh	
    gmake
    gmake install
    popd
}
