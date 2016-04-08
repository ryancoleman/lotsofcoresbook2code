#!/bin/bash

. /opt/intel/composerxe/bin/compilervars.sh intel64
export PATH=/opt/gcc-4.9.0/bin:/opt/binutils-2.24/bin:${PATH}
#export LD_LIBRARY_PATH=/opt/gcc-4.7.0/lib64:${LD_LIBRARY_PATH}
#export SINK_LD_LIBRARY_PATH=$MIC_LD_LIBRARY_PATH

#if [ "$1" == "knc" ]; then
#	export cmd="micnativeloadex ./gems"
#	export OUT_FILE=baseline_knc.csv
#        export plat=knc

#else
#	export cmd="./gems"
export OUT_FILE=result_xeon.csv
#fi

echo "Section,Compiler,Vectorized,Time" > ${OUT_FILE}
for section in c_simd_basic c_single_arch c_multi_arch c++_multi_arch fortran
do
	for comp in intel
	do
		for vect in yes no
		do
			pushd ${section} > /dev/null
			make clean build COMP=${comp} OMP=${vect} > /dev/null
			res=`./main`
                        popd > /dev/null
			echo $res
			time=`echo $res | cut -d ' ' -f 4`
			echo "$section,$comp,$vect,$time" >> $OUT_FILE
		done
	done
done

cat $OUT_FILE
