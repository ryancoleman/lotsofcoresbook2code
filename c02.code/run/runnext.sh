#!/bin/sh
# 
# Usage:  runnext.sh previous_test_id

if [ "$#" -ne 1 ] ; then
  echo "Usage:  $0 previous_test_id"
  exit -1
fi

prevtest=$1

autodir=../autotests

if ! [ -d $autodir ] ; then
  echo "ERROR:  $autodir does not exist"
  exit -1
fi

mv -f stdout ${autodir}/stdout.${prevtest}
mv -f timing.summary ${autodir}/timing.summary.${prevtest}
mv -f timing.0 ${autodir}/timing.0.${prevtest}
rm -f wsm6_output.dat

./runscript

