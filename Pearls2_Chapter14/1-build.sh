#!/bin/bash

oldpwd=`pwd`

# build pyMIC module
cd pyMIC
make
cd $oldpwd

# build GPAW
cd gpaw
./1-build.sh
cd $oldpwd
