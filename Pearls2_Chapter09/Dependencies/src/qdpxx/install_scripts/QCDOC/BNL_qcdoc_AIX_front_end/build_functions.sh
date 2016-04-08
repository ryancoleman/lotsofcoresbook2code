#!/usr/bin/bash
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

function build_bagel_qdp {
  bagelqdpdir=$1
  bagel_install_dir=$2
  bagelqdp_installdir=$3
  bagelqdp_prec=$4
  bagelqdp_proc=$5
  bagelqdp_hostsys=$6
  bagelqdp_buildsys=$7
    
  bagelqdp_builddir=build_bagel_qdp_${bagel_prec}

  echo Building BAGEL QDP for installing in ${bagelqdp_installdir}
  rm -rf ${bagelqdp_builddir}
  rm -rf ${bagelqdp_installdir}
  mkdir -p ${bagelqdp_builddir}
  pushd ${bagelqdp_builddir}
  command="${bagelqdpdir}/configure"
  command=${command}" --prefix=${bagelqdp_installdir}"
  command=${command}" --with-bagel=${bagel_install_dir}"
  command=${command}" --enable-precision=${bagel_prec}"
  command=${command}" --enable-target-cpu=${bagel_proc}"
  command=${command}" --host=${host_sys}"
  command=${command}" --build=${build_sys}"
  command=${command}" CXXFLAGS=\"-O2 \" "
  command=${command}" CFLAGS=\" -O2 \" "
  command=${command}" ASFLAGS=\" \" "
  command=${command}" LDFLAGS=\" \" "
  command=${command}" LIBS=\"  \" "
  echo ${command} > ./configure_bagel_qdp.sh
  chmod u+x ./configure_bagel_qdp.sh
  ./configure_bagel_qdp.sh
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
	    qos=${10};
	    qmp_cxxflags="-I${qos}/quser/include";
	    qmp_ldflags="-L${qos}/quser/qcd_api";
	    qmp_libs="-lqmp";
        else
	    qmp_cxxflags="";
	    qmp_ldflags="";
	    qmp_libs="";
        fi
	rm -rf ${wilson_install_dir}
	rm -rf ./build_wilson_dslash_${bagel_prec}
	mkdir -p ./build_wilson_dslash_${bagel_prec}
	pushd ./build_wilson_dslash_${bagel_prec}
	export PATH=/home/ed/bj/bin:$PATH

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
	command=${command}" CXXFLAGS=\"-O2 ${qmp_cxxflags}\" "
	command=${command}" CFLAGS=\" -O2 ${qmp_cxxflags}\" "
	command=${command}" ASFLAGS=\" \" "
	command=${command}" LDFLAGS=\" ${qmp_ldflags}\" "
	command=${command}" LIBS=\" ${qmp_libs} \" "
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
	qos=$3
	libxml=$4
	precision=$5
	do_edram=$6
	do_blas=$7
	host_sys=$8
	build_sys=$9
	bagel_qdp_dir=${10}

        qdp_base=`basename ${qdp_install_dir}`
        qdp_builddir="./build_"${qdp_base}

	if test "X${do_edram}X" == "XyesX"; 
	then
	    qcdoc_edram="--enable-qcdoc-edram"; 
	else
	    qcdoc_edram="";
	fi

	echo ${do_blas}
	if test "X${do_blas}X" == "XyesX";
	then 
	     qcdoc="--enable-qcdoc";
        else
	     qcdoc="";
        fi

	if test "X${bagel_qdp_dir}X" != "XX";
        then 
	     bagelqdp="--with-bagel-qdp=${bagel_qdp_dir}"
        else
	     bagelqdp=""
        fi

	rm -rf ${qdp_install_dir}
	qdp_base=`basename ${qdp_install_dir}`
	qdp_builddir="./build_"${qdp_base}
	rm -rf ${qdp_builddir}
	mkdir -p ${qdp_builddir}
	pushd ${qdp_builddir}

	command="${qdpdir}/configure --enable-parallel-arch=parscalar "
	command=${command}"  QMP_CFLAGS=\"-I$qos/quser/include\" "
	command=${command}"  QMP_LDFLAGS=\"-L$qos/quser/lib\" QMP_LIBS=\"-lqmp\" "
	command=${command}"  CXXFLAGS=\"-O2 -finline-limit=50000\" "
	command=${command}" CFLAGS=\"-O2 \" "
	command=${command}" --enable-precision=${precision}"
	command=${command}" --disable-qmp-route"
	command=${command}" --enable-slow-route"
	command=${command}" --with-libxml2=${libxml}"
	command=${command}" --host=powerpc-gnu-elf "
	command=${command}" "${qcdoc_edram}" "${qcdoc}" "${bagelqdp}
	command=${command}" --build=none --prefix=${qdp_install_dir}"
	echo Configuring QDP++ with command:
	echo    ${command}
 	echo ${command} > configure_qdp.sh

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
	pab_dslashdir=$7
	do_gmp=$8

	echo Chromadir: ${chromadir}
	echo Chroma_Install_dir: ${chroma_install_dir}
	echo QDP_Install_dir: ${qdp_dir}
	echo HOST: ${host_sys}
	echo BUILD: ${build_sys}
	echo DO Dslash: ${do_pab_dslash}
	echo DO GMP: ${do_gmp}
	if test "X${do_pab_dslash}X" == "XyesX";
	then
	   pab_dslash="--with-bagel-wilson-dslash=${pab_dslashdir}";
	else
	   pab_dslash="";
	fi

	if test "X${do_gmp}X" == "XyesX";
	then 
	   gmpdir=${9}
	   gmp="--with-gmp=${gmpdir}";
        else 
	   gmp=""
        fi

	rm -rf ${chroma_install_dir}
	install_base_name=`basename ${chroma_install_dir}`
	builddir=./build_${install_base_name}
	rm -rf ${builddir}

	mkdir -p ${chroma_install_dir}
	mkdir -p ${builddir}
	pushd ${builddir}
	export PATH=/home/ed/bj/bin:$PATH
	command="${chromadir}/configure CXXFLAGS=\"\" "
	command=${command}" --with-qdp=${qdp_dir} "
	command=${command}" --prefix=${chroma_install_dir}"
	command=${command}" --host=${host_sys}"
	command=${command}" --build=${build_sys}"
	command=${command}" ${gmp}"
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
    command=${command}" CFLAGS=\"\" "
    command=${command}" --host=${host_sys}"
    command=${command}" --build=${build_sys}"
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
