/*  Copyright (C) 2003-2007  CAMP
 *  Copyright (C) 2007-2009  CAMd
 *  Copyright (C) 2007-2010  CSC - IT Center for Science Ltd.
 *  Please see the accompanying LICENSE file for further information. */

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL GPAW_ARRAY_API
#include <numpy/arrayobject.h>

#ifdef GPAW_WITH_HDF5
PyMODINIT_FUNC init_gpaw_hdf5(void);
#endif

#ifdef GPAW_HPM
PyObject* ibm_hpm_start(PyObject *self, PyObject *args);
PyObject* ibm_hpm_stop(PyObject *self, PyObject *args);
PyObject* ibm_mpi_start(PyObject *self);
PyObject* ibm_mpi_stop(PyObject *self);
#endif

#ifdef CRAYPAT
#include <pat_api.h>
PyObject* craypat_region_begin(PyObject *self, PyObject *args);
PyObject* craypat_region_end(PyObject *self, PyObject *args);
#endif


PyObject* symmetrize(PyObject *self, PyObject *args);
PyObject* symmetrize_ft(PyObject *self, PyObject *args);
PyObject* symmetrize_wavefunction(PyObject *self, PyObject *args);
PyObject* symmetrize_return_index(PyObject *self, PyObject *args);
PyObject* symmetrize_with_index(PyObject *self, PyObject *args);
PyObject* map_k_points(PyObject *self, PyObject *args);
PyObject* scal(PyObject *self, PyObject *args);
PyObject* mmm(PyObject *self, PyObject *args);
PyObject* gemm(PyObject *self, PyObject *args);
PyObject* gemv(PyObject *self, PyObject *args);
PyObject* axpy(PyObject *self, PyObject *args);
PyObject* czher(PyObject *self, PyObject *args);
PyObject* rk(PyObject *self, PyObject *args);
PyObject* r2k(PyObject *self, PyObject *args);
PyObject* dotc(PyObject *self, PyObject *args);
PyObject* dotu(PyObject *self, PyObject *args);
PyObject* multi_dotu(PyObject *self, PyObject *args);
PyObject* multi_axpy(PyObject *self, PyObject *args);
PyObject* diagonalize(PyObject *self, PyObject *args);
PyObject* diagonalize_mr3(PyObject *self, PyObject *args);
PyObject* general_diagonalize(PyObject *self, PyObject *args);
PyObject* inverse_cholesky(PyObject *self, PyObject *args);
PyObject* inverse_symmetric(PyObject *self, PyObject *args);
PyObject* inverse_general(PyObject *self, PyObject *args);
PyObject* linear_solve_band(PyObject *self, PyObject *args);
PyObject* linear_solve_tridiag(PyObject *self, PyObject *args);
PyObject* right_eigenvectors(PyObject *self, PyObject *args);
PyObject* NewLocalizedFunctionsObject(PyObject *self, PyObject *args);
PyObject* NewOperatorObject(PyObject *self, PyObject *args);
PyObject* NewSplineObject(PyObject *self, PyObject *args);
PyObject* NewTransformerObject(PyObject *self, PyObject *args);
PyObject* pc_potential(PyObject *self, PyObject *args);
PyObject* pc_potential_value(PyObject *self, PyObject *args);
PyObject* heap_mallinfo(PyObject *self);
PyObject* elementwise_multiply_add(PyObject *self, PyObject *args);
PyObject* utilities_gaussian_wave(PyObject *self, PyObject *args);
PyObject* utilities_vdot(PyObject *self, PyObject *args);
PyObject* utilities_vdot_self(PyObject *self, PyObject *args);
PyObject* errorfunction(PyObject *self, PyObject *args);
PyObject* cerf(PyObject *self, PyObject *args);
PyObject* pack(PyObject *self, PyObject *args);
PyObject* unpack(PyObject *self, PyObject *args);
PyObject* unpack_complex(PyObject *self, PyObject *args);
PyObject* hartree(PyObject *self, PyObject *args);
PyObject* localize(PyObject *self, PyObject *args);
PyObject* NewXCFunctionalObject(PyObject *self, PyObject *args);
PyObject* NewlxcXCFunctionalObject(PyObject *self, PyObject *args);
PyObject* lxcXCFuncNum(PyObject *self, PyObject *args);
PyObject* exterior_electron_density_region(PyObject *self, PyObject *args);
PyObject* plane_wave_grid(PyObject *self, PyObject *args);
PyObject* overlap(PyObject *self, PyObject *args);
PyObject* vdw(PyObject *self, PyObject *args);
PyObject* vdw2(PyObject *self, PyObject *args);
PyObject* spherical_harmonics(PyObject *self, PyObject *args);
PyObject* spline_to_grid(PyObject *self, PyObject *args);
PyObject* NewLFCObject(PyObject *self, PyObject *args);
#if defined(GPAW_WITH_SL) && defined(PARALLEL)
PyObject* new_blacs_context(PyObject *self, PyObject *args);
PyObject* get_blacs_gridinfo(PyObject* self, PyObject *args);
PyObject* get_blacs_local_shape(PyObject* self, PyObject *args);
PyObject* blacs_destroy(PyObject *self, PyObject *args);
PyObject* scalapack_set(PyObject *self, PyObject *args);
PyObject* scalapack_redist(PyObject *self, PyObject *args);
PyObject* scalapack_diagonalize_dc(PyObject *self, PyObject *args);
PyObject* scalapack_diagonalize_ex(PyObject *self, PyObject *args);
#ifdef GPAW_MR3
PyObject* scalapack_diagonalize_mr3(PyObject *self, PyObject *args);
#endif
PyObject* scalapack_general_diagonalize_dc(PyObject *self, PyObject *args);
PyObject* scalapack_general_diagonalize_ex(PyObject *self, PyObject *args);
#ifdef GPAW_MR3
PyObject* scalapack_general_diagonalize_mr3(PyObject *self, PyObject *args);
#endif
PyObject* scalapack_inverse_cholesky(PyObject *self, PyObject *args);
PyObject* scalapack_inverse(PyObject *self, PyObject *args);
PyObject* scalapack_solve(PyObject *self, PyObject *args);
PyObject* pblas_tran(PyObject *self, PyObject *args);
PyObject* pblas_gemm(PyObject *self, PyObject *args);
PyObject* pblas_hemm(PyObject *self, PyObject *args);
PyObject* pblas_gemv(PyObject *self, PyObject *args);
PyObject* pblas_r2k(PyObject *self, PyObject *args);
PyObject* pblas_rk(PyObject *self, PyObject *args);
#endif

PyObject* offload_report(PyObject *self, PyObject *args)
{
  int n;
  if (!PyArg_ParseTuple(args, "i", &n))
    return NULL;
  _Offload_report(n);
  Py_RETURN_NONE;
}

#ifdef GPAW_PAPI
PyObject* papi_mem_info(PyObject *self, PyObject *args);
#endif

// Moving least squares interpolation
PyObject* mlsqr(PyObject *self, PyObject *args);

static PyMethodDef functions[] = {
  {"symmetrize", symmetrize, METH_VARARGS, 0},
  {"symmetrize_ft", symmetrize_ft, METH_VARARGS, 0},
  {"symmetrize_wavefunction", symmetrize_wavefunction, METH_VARARGS, 0},
  {"symmetrize_return_index", symmetrize_return_index, METH_VARARGS, 0},
  {"symmetrize_with_index", symmetrize_with_index, METH_VARARGS, 0},
  {"map_k_points", map_k_points, METH_VARARGS, 0},
  {"scal", scal, METH_VARARGS, 0},
  {"mmm", mmm, METH_VARARGS, 0},
  {"gemm", gemm, METH_VARARGS, 0},
  {"gemv", gemv, METH_VARARGS, 0},
  {"axpy", axpy, METH_VARARGS, 0},
  {"czher", czher, METH_VARARGS, 0},
  {"rk",  rk,  METH_VARARGS, 0},
  {"r2k", r2k, METH_VARARGS, 0},
  {"dotc", dotc, METH_VARARGS, 0},
  {"dotu", dotu, METH_VARARGS, 0},
  {"multi_dotu", multi_dotu, METH_VARARGS, 0},
  {"multi_axpy", multi_axpy, METH_VARARGS, 0},
  {"diagonalize", diagonalize, METH_VARARGS, 0},
  {"diagonalize_mr3", diagonalize_mr3, METH_VARARGS, 0},
  {"general_diagonalize", general_diagonalize, METH_VARARGS, 0},
  {"inverse_cholesky", inverse_cholesky, METH_VARARGS, 0},
  {"inverse_symmetric", inverse_symmetric, METH_VARARGS, 0},
  {"inverse_general", inverse_general, METH_VARARGS, 0},
  {"linear_solve_band", linear_solve_band, METH_VARARGS, 0},
  {"linear_solve_tridiag", linear_solve_tridiag, METH_VARARGS, 0},
  {"right_eigenvectors", right_eigenvectors, METH_VARARGS, 0},
  {"LocalizedFunctions", NewLocalizedFunctionsObject, METH_VARARGS, 0},
  {"Operator", NewOperatorObject, METH_VARARGS, 0},
  {"Spline", NewSplineObject, METH_VARARGS, 0},
  {"Transformer", NewTransformerObject, METH_VARARGS, 0},
  {"heap_mallinfo", (PyCFunction) heap_mallinfo, METH_NOARGS, 0},
  {"elementwise_multiply_add", elementwise_multiply_add, METH_VARARGS, 0},
  {"utilities_gaussian_wave", utilities_gaussian_wave, METH_VARARGS, 0},
  {"utilities_vdot", utilities_vdot, METH_VARARGS, 0},
  {"utilities_vdot_self", utilities_vdot_self, METH_VARARGS, 0},
  {"eed_region", exterior_electron_density_region, METH_VARARGS, 0},
  {"plane_wave_grid", plane_wave_grid, METH_VARARGS, 0},
  {"erf",        errorfunction,        METH_VARARGS, 0},
  {"cerf",       cerf,        METH_VARARGS, 0},
  {"pack",       pack,           METH_VARARGS, 0},
  {"unpack",       unpack,           METH_VARARGS, 0},
  {"unpack_complex",       unpack_complex,           METH_VARARGS, 0},
  {"hartree",        hartree,        METH_VARARGS, 0},
  {"localize",       localize,        METH_VARARGS, 0},
  {"XCFunctional",    NewXCFunctionalObject,    METH_VARARGS, 0},
  /*  {"MGGAFunctional",    NewMGGAFunctionalObject,    METH_VARARGS, 0},*/
  {"lxcXCFunctional",    NewlxcXCFunctionalObject,    METH_VARARGS, 0},
  {"lxcXCFuncNum",    lxcXCFuncNum,    METH_VARARGS, 0},
  {"overlap",       overlap,        METH_VARARGS, 0},
  {"vdw", vdw, METH_VARARGS, 0},
  {"vdw2", vdw2, METH_VARARGS, 0},
  {"spherical_harmonics", spherical_harmonics, METH_VARARGS, 0},
  {"pc_potential", pc_potential, METH_VARARGS, 0},
  {"pc_potential_value", pc_potential_value, METH_VARARGS, 0},
  {"spline_to_grid", spline_to_grid, METH_VARARGS, 0},
  {"LFC", NewLFCObject, METH_VARARGS, 0},
  /*
  {"calculate_potential_matrix", calculate_potential_matrix, METH_VARARGS, 0},
  {"construct_density", construct_density, METH_VARARGS, 0},
  {"construct_density1", construct_density1, METH_VARARGS, 0},
  */
#if defined(GPAW_WITH_SL) && defined(PARALLEL)
  {"new_blacs_context", new_blacs_context, METH_VARARGS, NULL},
  {"get_blacs_gridinfo", get_blacs_gridinfo, METH_VARARGS, NULL},
  {"get_blacs_local_shape", get_blacs_local_shape, METH_VARARGS, NULL},
  {"blacs_destroy",     blacs_destroy,      METH_VARARGS, 0},
  {"scalapack_set", scalapack_set, METH_VARARGS, 0},
  {"scalapack_redist",      scalapack_redist,     METH_VARARGS, 0},
  {"scalapack_diagonalize_dc", scalapack_diagonalize_dc, METH_VARARGS, 0},
  {"scalapack_diagonalize_ex", scalapack_diagonalize_ex, METH_VARARGS, 0},
#ifdef GPAW_MR3
  {"scalapack_diagonalize_mr3", scalapack_diagonalize_mr3, METH_VARARGS, 0},
#endif // GPAW_MR3
  {"scalapack_general_diagonalize_dc",
   scalapack_general_diagonalize_dc, METH_VARARGS, 0},
  {"scalapack_general_diagonalize_ex",
   scalapack_general_diagonalize_ex, METH_VARARGS, 0},
#ifdef GPAW_MR3
  {"scalapack_general_diagonalize_mr3",
   scalapack_general_diagonalize_mr3, METH_VARARGS, 0},
#endif // GPAW_MR3
  {"scalapack_inverse_cholesky", scalapack_inverse_cholesky, METH_VARARGS, 0},
  {"scalapack_inverse", scalapack_inverse, METH_VARARGS, 0},
  {"scalapack_solve", scalapack_solve, METH_VARARGS, 0},
  {"pblas_tran", pblas_tran, METH_VARARGS, 0},
  {"pblas_gemm", pblas_gemm, METH_VARARGS, 0},
  {"pblas_hemm", pblas_hemm, METH_VARARGS, 0},
  {"pblas_gemv", pblas_gemv, METH_VARARGS, 0},
  {"pblas_r2k", pblas_r2k, METH_VARARGS, 0},
  {"pblas_rk", pblas_rk, METH_VARARGS, 0},
#endif // GPAW_WITH_SL && PARALLEL
#ifdef GPAW_HPM
  {"hpm_start", ibm_hpm_start, METH_VARARGS, 0},
  {"hpm_stop", ibm_hpm_stop, METH_VARARGS, 0},
  {"mpi_start", (PyCFunction) ibm_mpi_start, METH_NOARGS, 0},
  {"mpi_stop", (PyCFunction) ibm_mpi_stop, METH_NOARGS, 0},
#endif // GPAW_HPM
#ifdef CRAYPAT
  {"craypat_region_begin", craypat_region_begin, METH_VARARGS, 0},
  {"craypat_region_end", craypat_region_end, METH_VARARGS, 0},
#endif // CRAYPAT
#ifdef GPAW_PAPI
  {"papi_mem_info", papi_mem_info, METH_VARARGS, 0},
#endif // GPAW_PAPI
  {"mlsqr", mlsqr, METH_VARARGS, 0},
  {"offload_report", offload_report, METH_VARARGS, 0},
  {0, 0, 0, 0}
};

#ifdef PARALLEL
extern PyTypeObject MPIType;
extern PyTypeObject GPAW_MPI_Request_type;
#endif

extern PyTypeObject LFCType;
extern PyTypeObject LocalizedFunctionsType;
extern PyTypeObject OperatorType;
extern PyTypeObject SplineType;
extern PyTypeObject TransformerType;
extern PyTypeObject XCFunctionalType;
extern PyTypeObject lxcXCFunctionalType;

#ifndef GPAW_INTERPRETER
PyMODINIT_FUNC init_gpaw(void)
{
#ifdef PARALLEL
  if (PyType_Ready(&MPIType) < 0)
    return;
  if (PyType_Ready(&GPAW_MPI_Request_type) < 0)
    return;
#endif

  if (PyType_Ready(&LFCType) < 0)
    return;
  if (PyType_Ready(&LocalizedFunctionsType) < 0)
    return;
  if (PyType_Ready(&OperatorType) < 0)
    return;
  if (PyType_Ready(&SplineType) < 0)
    return;
  if (PyType_Ready(&TransformerType) < 0)
    return;
  if (PyType_Ready(&XCFunctionalType) < 0)
    return;
  if (PyType_Ready(&lxcXCFunctionalType) < 0)
    return;

  PyObject* m = Py_InitModule3("_gpaw", functions,
             "C-extension for GPAW\n\n...\n");
  if (m == NULL)
    return;

#ifdef PARALLEL
  Py_INCREF(&MPIType);
  Py_INCREF(&GPAW_MPI_Request_type);
  PyModule_AddObject(m, "Communicator", (PyObject *)&MPIType);
#endif

  Py_INCREF(&LFCType);
  Py_INCREF(&LocalizedFunctionsType);
  Py_INCREF(&OperatorType);
  Py_INCREF(&SplineType);
  Py_INCREF(&TransformerType);
  Py_INCREF(&XCFunctionalType);
  Py_INCREF(&lxcXCFunctionalType);

  import_array();
}
#endif


#ifdef GPAW_INTERPRETER
extern DL_EXPORT(int) Py_Main(int, char **);

// Performance measurement
int gpaw_perf_init();
void gpaw_perf_finalize();

#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

int gpaw_offload_enabled = 1;

__declspec(target(mic))
void init_openmp() {
#ifdef __MIC__
#pragma omp parallel
    {
        /* do nothing */
    }
#endif    
}

int
main(int argc, char **argv)
{
  int status;
  char* env = NULL;

  env = getenv("GPAW_OFFLOAD");
  if (env) {
      errno = 0;
      gpaw_offload_enabled = strtol(env, NULL, 10);
      if (errno) {
        fprintf(stderr, 
                "Wrong value for for GPAW_OFFLOAD.\nShould be either 0 or 1, but was %s\n",
                env);
      }
  }
  fprintf(stderr, "GPAW info: GPAW_OFFLOAD=%d\n", gpaw_offload_enabled);
  
#ifdef CRAYPAT
  PAT_region_begin(1, "C-Initializations");
#endif

#ifndef GPAW_OMP
  MPI_Init(&argc, &argv);
#else
  int granted;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &granted);
  if(granted != MPI_THREAD_MULTIPLE) exit(1);
#endif // GPAW_OMP

// Get initial timing
  double t0 = MPI_Wtime();

#ifdef GPAW_PERFORMANCE_REPORT
  gpaw_perf_init();
#endif

#ifdef GPAW_MPI_MAP
  int tag = 99;
  int myid, numprocs, i, procnamesize;
  char procname[MPI_MAX_PROCESSOR_NAME];
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs );
  MPI_Comm_rank(MPI_COMM_WORLD, &myid );
  MPI_Get_processor_name(procname, &procnamesize);
  if (myid > 0) {
      MPI_Send(&procnamesize, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
      MPI_Send(procname, procnamesize, MPI_CHAR, 0, tag, MPI_COMM_WORLD);
  }
  else {
      printf("MPI_COMM_SIZE is %d \n", numprocs);
      printf("%s \n", procname);
      
      for (i = 1; i < numprocs; ++i) {
          MPI_Recv(&procnamesize, 1, MPI_INT, i, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          MPI_Recv(procname, procnamesize, MPI_CHAR, i, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          printf("%s \n", procname);
      }
  }
#endif // GPAW_MPI_MAP

#ifdef GPAW_MPI_DEBUG
  // Default Errhandler is MPI_ERRORS_ARE_FATAL
  MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
#endif

  // Progname seems to be needed in some circumstances to resolve
  // correct default sys.path
  Py_SetProgramName(argv[0]);

  Py_Initialize();

#pragma offload target(mic) if(gpaw_offload_enabled)
    {
        init_openmp();
    }
  
  if (PyType_Ready(&MPIType) < 0)
    return -1;

  if (PyType_Ready(&LFCType) < 0)
    return -1;
  if (PyType_Ready(&LocalizedFunctionsType) < 0)
    return -1;
  if (PyType_Ready(&OperatorType) < 0)
    return -1;
  if (PyType_Ready(&SplineType) < 0)
    return -1;
  if (PyType_Ready(&TransformerType) < 0)
    return -1;
  if (PyType_Ready(&XCFunctionalType) < 0)
    return -1;
  if (PyType_Ready(&lxcXCFunctionalType) < 0)
    return -1;

  PyObject* m = Py_InitModule3("_gpaw", functions,
             "C-extension for GPAW\n\n...\n");
  if (m == NULL)
    return -1;

  Py_INCREF(&MPIType);
  PyModule_AddObject(m, "Communicator", (PyObject *)&MPIType);

  // Add initial time to _gpaw object
  PyModule_AddObject(m, "time0", PyFloat_FromDouble(t0));

  Py_INCREF(&LFCType);
  Py_INCREF(&LocalizedFunctionsType);
  Py_INCREF(&OperatorType);
  Py_INCREF(&SplineType);
  Py_INCREF(&TransformerType);
  Py_INCREF(&XCFunctionalType);
  Py_INCREF(&lxcXCFunctionalType);

#ifdef GPAW_WITH_HDF5
  init_gpaw_hdf5();
#endif
  import_array1(-1);
  MPI_Barrier(MPI_COMM_WORLD);
#ifdef CRAYPAT
  PAT_region_end(1);
  PAT_region_begin(2, "all other");
#endif
  status = Py_Main(argc, argv);
#ifdef CRAYPAT
  PAT_region_end(2);
#endif

#ifdef GPAW_PERFORMANCE_REPORT
  gpaw_perf_finalize();
#endif

  MPI_Finalize();
  return status;
}
#endif // GPAW_INTERPRETER
