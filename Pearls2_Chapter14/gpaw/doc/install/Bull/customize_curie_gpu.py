import os

scalapack = True
hdf5 = False
# ld: /usr/local/phdf5-1.8.5/lib/libhdf5.a(H5.o): relocation R_X86_64_32 against `.rodata.str1.4' can not be used when making a shared object; recompile with -fPIC
compiler = 'icc'

mpi='/opt/mpi/bullxmpi/1.1.16.5'
mkl='/usr/local/Intel_compilers/c/composer_xe_2011_sp1.7.256/mkl/lib/intel64'
intel='/usr/local/Intel_compilers/c/composer_xe_2011_sp1.7.256/compiler/lib/intel64'
hdf='/usr/local/phdf5-1.8.5'
#
# cublasZdgmm does not exist in cuda 4.2
# /tmp/ipo_iccjq2M5h1.o: In function `cudgmm':
# ipo_out1.c:(.text.hot0001d+0x522b): undefined reference to `cublasZdgmm'
# strings /usr/local/cuda-4.2/lib64/libcublas.so | grep "cublasZdgmm"
cuda='/usr/local/cuda-4.2'  # comment out if no cuda

libraries =[
    'cublas', 'cufft', 'cuda',  # comment out if no cuda
    'cudart',  # comment out if no cuda
    #'mkl_def',
    'mkl_scalapack_lp64', 'mkl_intel_lp64', 'mkl_sequential',
    'mkl_core', 'mkl_blacs_openmpi_lp64',
    #'hdf5',
    'mpi',
    ]
library_dirs =[
    intel,
    os.path.join(mpi, 'lib'),
    mkl,
    os.path.join(cuda, 'lib64'),  # comment out if no cuda
    #os.path.join(hdf, 'lib'),
    ]
include_dirs +=[
    os.path.join(mpi, 'include'),
    os.path.join(cuda, 'include'),  # comment out if no cuda
    #os.path.join(hdf, 'include'),
    ]
extra_link_args =[
    '-Wl,-rpath=' + intel +
    ',-rpath=' + os.path.join(mpi, 'lib') +
    ',-rpath=' + os.path.join(cuda, 'lib64') +  # comment out if no cuda
    ',-rpath=' + mkl
    #',-rpath=' + os.path.join(hdf, 'lib')
    ]
extra_compile_args =['-xHOST', '-O3', '-ipo', '-std=c99', '-fPIC', '-Wall']
extra_objects += ['./c/cukernels.o']
define_macros += [('GPAW_NO_UNDERSCORE_CBLACS', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_CSCALAPACK', '1')]
mpicompiler = os.path.join(mpi, 'bin', 'mpicc')
mpilinker = mpicompiler

