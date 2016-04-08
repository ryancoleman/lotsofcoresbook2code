#!/bin/bash
#PBS -l nodes=1:ppn=4

. /home/camp/modulefiles.sh
module load openmpi/1.3.3-1.el5.fys.gfortran43.4.3.2
module load blacs-gfortran4364/1.1-24.el5.fys.gfortran43.4.3.2.openmpi.1.3.3
module load acml-gfortran4364/4.3.0-1.el5.fys
module load scalapack-gfortran4364/1.8.0-1.el5.fys.gfortran43.4.3.2.openmpi.1.3.3.acml.4.3.0.acml.4.3.0
mpirun -np 4 test.exe
