#!/bin/bash

. /opt/intel/composerxe/bin/compilervars.sh intel64
export SINK_LD_LIBRARY_PATH=$MIC_LD_LIBRARY_PATH

export OUT_FILE=result_knc.csv

echo "Section,Compiler,Vectorized,Time" > ${OUT_FILE}
for section in c_simd_basic fortran
do
	for vect in yes no
	do
		pushd ${section} > /dev/null
		make clean build OMP=${vect} PLAT=knc > /dev/null
		res=`micnativeloadex ./main`
                popd > /dev/null
		echo $res
		time=`echo $res | cut -d ' ' -f 4`
		echo "$section,intel,$vect,$time" >> $OUT_FILE
	done
done

cat $OUT_FILE
