scalapack = False

compiler = 'icc'
libraries =['mkl_intel_lp64','mkl_sequential','mkl_cdft_core','mkl_core','pthread','m']

library_dirs = ['/nfs/slac/g/suncatfs/sw/external/intel11.1/openmpi/1.4.3/install/lib','/afs/slac/package/intel_tools/compiler11.1/mkl/lib/em64t/']

include_dirs += ['/nfs/slac/g/suncatfs/sw/external/numpy/1.4.1/install/lib64/python2.4/site-packages/numpy/core/include']

extra_link_args += ['-fPIC']

extra_compile_args = ['-I/afs/slac/package/intel_tools/compiler11.1/mkl/include','-xHOST','-O1','-ipo','-no-prec-div','-static','-std=c99','-fPIC']

define_macros =[('GPAW_NO_UNDERSCORE_CBLACS', '1'), ('GPAW_NO_UNDERSCORE_CSCALAPACK', '1')]

mpicompiler = 'mpicc'
mpilinker = mpicompiler
