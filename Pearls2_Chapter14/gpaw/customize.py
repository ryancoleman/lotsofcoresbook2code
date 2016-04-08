#User provided customizations for the gpaw setup
import os

#Here, one can override the default arguments, or append own
#arguments to default ones
#To override use the form
#     libraries = ['somelib','otherlib']
#To append use the form
#     libraries += ['somelib','otherlib']

compiler = 'mpiicc'
mklroot = os.environ['MKLROOT']

library_dirs += [mklroot + '/lib/intel64/']
#libraries = ['mkl_intel_lp64', 'mkl_sequential', 'mkl_core', 
#             'mkl_lapack95_lp64', 'mkl_scalapack_lp64', 'mkl_blacs_intelmpi_lp64', 
#			 'pthread'
#			]
libraries = ['mkl_intel_lp64', 'mkl_intel_thread', 'mkl_core', 
             'mkl_lapack95_lp64', 'mkl_scalapack_lp64', 'mkl_blacs_intelmpi_lp64', 
			 'pthread'
			]

# libxc
library_dirs += [os.environ['LIBXCDIR'] + '/lib']
include_dirs += [os.environ['LIBXCDIR'] + '/include']
libraries += ['xc']
			
# compiler settings for Intel Composer 2013 (ver 14.x)
extra_compile_args = ['-xHOST', '-openmp', '-g', '-O3', '-no-prec-div', '-static', '-std=c99', '-fPIC']
extra_compile_args += ['-offload-option,mic,compiler,"-openmp"']
extra_compile_args += ['-opt-report-phase=offload']
# extra_compile_args += ['-offload=none']

# for Intel Composer 2015 (ver 15.x)
# extra_compile_args += ['-qoffload-option,mic,compiler,"-qopenmp"']
# extra_compile_args += ['-qopt-report-phase=offload']
# extra_compile_args += ['-qoffload=none']

# linker settings for MKL on KNC
mic_mkl_lib = mklroot + '/lib/mic/'
extra_link_args += ['-openmp']
extra_link_args += ['-offload-option,mic,link,"-L' + mic_mkl_lib + ' -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread"']

# GPAW defines
define_macros += [('GPAW_NO_UNDERSCORE_CBLACS', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_CSCALAPACK', '1')]
define_macros += [("GPAW_ASYNC",1)]
define_macros += [("GPAW_MPI2",1)]

# use Intel MPI w/ Intel Composer
mpicompiler = 'mpiicc'
mpilinker = 'mpiicc'

platform_id = 'intel64'

#hdf5 = True

# Valid values for scalapack are False, or True:
# False (the default) - no ScaLapack compiled in
# True - ScaLapack compiled in
scalapack = True

