scalapack = False

compiler = 'gcc'

extra_compile_args += [
    '-O3',
    '-funroll-all-loops',
    '-fPIC',
    ]

libraries = ['gfortran', 'util']

blas_lib_path = '/home/lv70174/gpaw/opt/acml-4.0.1/gfortran64/lib/'
lapack_lib_path = blas_lib_path

extra_link_args = [
    blas_lib_path+'libacml.a',
    lapack_lib_path+'libacml.a',
    ]
