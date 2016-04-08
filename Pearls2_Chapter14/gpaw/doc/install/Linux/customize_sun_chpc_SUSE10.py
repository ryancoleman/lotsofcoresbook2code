compiler = 'gcc'
mpicompiler = '/opt/SUNWhpc/HPC8.2/gnu/bin/mpicc'
mpilinker = 'MPICH_CC=gcc mpicc -Xlinker --export-dynamic'

extra_compile_args = ['-O3', '-std=c99', '-fpic']

scalapack = False

mkl_dir = '/opt/gridware/intel/Compiler/11.1/056/mkl/lib/em64t/'

library_dirs += [mkl_dir]
libraries = ['mkl_intel_lp64' ,'mkl_sequential' ,'mkl_core',
             'mkl_lapack',
             ]

extra_link_args = ['-Wl,-rpath=' + mkl_dir]

