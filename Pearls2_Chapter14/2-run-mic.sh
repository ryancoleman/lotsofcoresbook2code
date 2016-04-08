#!/bin/bash

oldpwd=`pwd`
cd example
mpirun -np 12 -genv GPAW_OFFLOAD 1 -genv GPAW_PPN 12 ./wrapper.sh 12 gpaw-python C60.py
cd $oldpwd
