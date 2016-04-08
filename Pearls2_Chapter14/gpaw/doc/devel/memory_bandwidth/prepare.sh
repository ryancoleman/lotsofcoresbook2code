#!/bin/sh

if test -z $NCORES;
then
   export NCORES=8
fi

if test -z $MACHINE;
then
   export MACHINE=TEST
fi

if [ ! -d "${MACHINE}" ]; then
    mkdir ${MACHINE}
    echo "${MACHINE} created"
    cd ${MACHINE}
fi
index=0
while [ "$index" -le "$NCORES" ];
do
  if [ "$index" -eq 0 ];
      then
      p=1
  else
      p=$index
  fi
  #
  if [ ! -d "${MACHINE}_py${p}_01" ]; then
      mkdir ${MACHINE}_py${p}_01
      echo "${MACHINE}_py${p}_01 created"
  fi
  index=`expr $index + 2`
done
cd ..
