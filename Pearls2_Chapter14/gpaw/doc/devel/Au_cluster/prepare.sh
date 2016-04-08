#!/bin/bash

FORMAT=1

setFORMAT () {
    # function setFORMAT takes integer as the argument $1
    # and returns integer in the format %05d or %d (printf' like format)
    # depending on the FORMAT variable (1 or 0)
    if [ ${FORMAT} -eq "1" ]; then
        integer_formatted=`echo $1 | awk '{if ($1<10) printf("0000%.0f", $1); else if ($1<100) printf("000%.0f", $1); else if ($1<1000) printf("00%.0f", $1); else if ($1<10000) printf("0%.0f", $1); else printf("%.0f", $1)}'`
    else
        integer_formatted=$1
    fi
    echo $integer_formatted
}

if test -z $PATTERN;
    then
    echo "Error: no directory pattern provided"
    exit
fi

for p in 256 512 1024 2048 4096
  do
  proc=`setFORMAT $p`
  dir="${PATTERN}_${proc}_"
  if [ ! -d "${dir}" ]; then
      mkdir ${dir}
      echo "${dir} created"
  fi
done
