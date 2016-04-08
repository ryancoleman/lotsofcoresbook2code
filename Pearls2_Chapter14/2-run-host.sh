#!/bin/bash

oldpwd=`pwd`
cd example
mpirun -np 12 -genv GPAW_OFFLOAD 0 gpaw-python C60.py
cd $oldpwd
