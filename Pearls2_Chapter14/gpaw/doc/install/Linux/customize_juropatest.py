#User provided customizations for the gpaw setup

#Here, one can override the default arguments, or append own
#arguments to default ones
#To override use the form
#     libraries = ['somelib','otherlib']
#To append use the form
#     libraries += ['somelib','otherlib']

compiler = 'mpicc'
libraries = ['mkl_intel_lp64' ,'mkl_sequential' ,'mkl_core',
             'mkl_lapack',
             'mkl_scalapack_lp64', 'mkl_blacs_intelmpi_lp64',
             'pthread', 
             'readline', 'termcap',
             'xc'
             ]
#libraries = []
#libraries += []

#library_dirs = []
#library_dirs += []

#include_dirs = []
#include_dirs += []

#extra_link_args = []
#extra_link_args += []

#extra_compile_args = []
#extra_compile_args += []

#runtime_library_dirs = []
#runtime_library_dirs += []

#extra_objects = []
#extra_objects += []

#define_macros = []
#define_macros += []

#mpicompiler = None
#mpilinker = None
#mpi_libraries = []
#mpi_libraries += []

#mpi_library_dirs = []
#mpi_library_dirs += []

#mpi_include_dirs = []
#mpi_include_dirs += []

#mpi_runtime_library_dirs = []
#mpi_runtime_library_dirs += []

#mpi_define_macros = []
#mpi_define_macros += []

#platform_id = ''

#hdf5 = True

# Valid values for scalapack are False, or True:
# False (the default) - no ScaLapack compiled in
# True - ScaLapack compiled in
# Warning! At least scalapack 2.0.1 is required!
# See https://trac.fysik.dtu.dk/projects/gpaw/ticket/230
scalapack = False

if scalapack:
    libraries += ['scalapack']
    library_dirs += []
    define_macros += [('GPAW_NO_UNDERSCORE_CBLACS', '1')]
    define_macros += [('GPAW_NO_UNDERSCORE_CSCALAPACK', '1')]

# In order to link libxc installed in a non-standard location
# (e.g.: configure --prefix=/home/user/libxc-2.0.1-1), use:
# - static linking:
#include_dirs += ['/home/user/libxc-2.0.1-1/include']
#extra_link_args += ['/home/user/libxc-2.0.1-1/lib/libxc.a']
#if 'xc' in libraries: libraries.remove('xc')
# - dynamic linking (requires also setting LD_LIBRARY_PATH at runtime):
#include_dirs += ['/home/user/libxc-2.0.1-1/include']
#library_dirs += ['/home/user/libxc-2.0.1-1/lib']
#if 'xc' not in libraries: libraries.append('xc')
