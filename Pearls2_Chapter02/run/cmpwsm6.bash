#!/bin/bash

#set -x
case="real32level"
if (( $# == 1 )) ; then
  case=$1
fi

baselinedir="../../data/${case}"
dutdir="."
echo "comparing $dutdir with baseline in $baselinedir ..."
outfile="wsm6_output.dat"
echo "comparing $outfile ..."
cmp ${baselinedir}/${outfile} ${dutdir}/${outfile}
if (( $? == 0 )) ; then
  echo "OK"
else
  echo "FAIL"
fi
stdoutfile="stdout"
echo "diffing $stdoutfile ..."
xxdiff ${baselinedir}/${stdoutfile} ${dutdir}/${stdoutfile}

