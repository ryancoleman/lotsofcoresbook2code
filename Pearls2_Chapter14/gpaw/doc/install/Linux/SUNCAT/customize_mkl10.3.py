scalapack = False

compiler = 'icc'
libraries =['mkl_rt','pthread','m']

library_dirs = ['/nfs/slac/g/suncatfs/sw/external/intel11.1/openmpi/1.4.3/install/lib','/afs/slac/package/intel_tools/2011u8/mkl/lib/intel64/']

include_dirs += ['/nfs/slac/g/suncatfs/sw/external/numpy/1.4.1/install/lib64/python2.4/site-packages/numpy/core/include']

extra_link_args += ['-fPIC']

extra_compile_args = ['-I/afs/slac/package/intel_tools/2011u8/mkl/include','-xHOST','-O1','-ipo','-no-prec-div','-static','-std=c99','-fPIC']

mpicompiler = 'mpicc'
mpilinker = mpicompiler
