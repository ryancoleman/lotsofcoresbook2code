scalapack = True

compiler = 'cc'
mpicompiler = compiler
mpilinker = mpicompiler

extra_compile_args += [
    '-O3',
    '-funroll-all-loops',
    '-fPIC',
    ]

libraries= []

dir_base = '/autofs/na1_home/farberow/sw/xt5/'
acml_base = '/opt/acml/4.1.0/gfortran64/lib/'
numpy_base = dir_base+'numpy-1.2.1/build/lib.linux-x86_64-2.5/numpy/'


extra_link_args = [
    '-L/usr/lib64 -lreadline -lncurses',
    numpy_base+'core/multiarray.a',
    numpy_base+'core/_sort.a',
    numpy_base+'core/scalarmath.a',
    numpy_base+'core/umath.a',
    numpy_base+'lib/_compiled_base.a',
    numpy_base+'numarray/_capi.a',
    numpy_base+'fft/fftpack_lite.a',
    numpy_base+'linalg/lapack_lite.a',
    numpy_base+'random/mtrand.a',
    '-L'+dir_base+'zlib-1.2.3-1/lib -lz',
    '-L'+dir_base+'expat-2.0.1-1/lib -lexpat',
    '-L'+acml_base+' -lacml',
]

define_macros += [('GPAW_NO_UNDERSCORE_CBLACS', '1')]
define_macros += [('GPAW_NO_UNDERSCORE_CSCALAPACK', '1')]
define_macros += [('NO_SOCKET', '1')]
