#!/bin/bash

oldpwd=`pwd`

# clean GPAW
cd gpaw
./2-clean.sh
cd $oldpwd

# clean pyMIC
cd pyMIC
make realclean
cd $oldpwd
