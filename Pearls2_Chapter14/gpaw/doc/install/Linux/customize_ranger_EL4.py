scalapack = True

extra_compile_args = ['-O3', '-std=c99', '-fpic']

compiler = 'gcc'
mpicompiler = '/opt/apps/gcc4_4/openmpi/1.3b/bin/mpicc'
mpilinker = mpicompiler

mkl_lib_path = '/opt/apps/intel/mkl/10.0.1.014/lib/em64t/'
ompi_lib_path = '/opt/apps/gcc4_4/openmpi/1.3b/lib'

libraries = []

# use static linking to avoid
# "cannot allocate memory for thread-local data: ABORT"
extra_link_args = [
mkl_lib_path+'libmkl_intel_lp64.a',
mkl_lib_path+'libmkl_sequential.a',
mkl_lib_path+'libmkl_core.a',
mkl_lib_path+'libmkl_blacs_openmpi_lp64.a',
mkl_lib_path+'libmkl_scalapack_lp64.a',
mkl_lib_path+'libmkl_blacs_openmpi_lp64.a',
mkl_lib_path+'libmkl_intel_lp64.a',
mkl_lib_path+'libmkl_sequential.a',
mkl_lib_path+'libmkl_core.a',
mkl_lib_path+'libmkl_intel_lp64.a',
mkl_lib_path+'libmkl_sequential.a',
mkl_lib_path+'libmkl_core.a',
]

extra_link_args += ['-Wl,-rpath='+ompi_lib_path]

define_macros += [('GPAW_NO_UNDERSCORE_CBLACS', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_CSCALAPACK', '1')]
