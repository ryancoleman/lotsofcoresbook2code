/*
    Copyright 2005-2013 Intel Corporation.  All Rights Reserved.

    The source code contained or described herein and all documents related
    to the source code ("Material") are owned by Intel Corporation or its
    suppliers or licensors.  Title to the Material remains with Intel
    Corporation or its suppliers and licensors.  The Material is protected
    by worldwide copyright laws and treaty provisions.  No part of the
    Material may be used, copied, reproduced, modified, published, uploaded,
    posted, transmitted, distributed, or disclosed in any way without
    Intel's prior express written permission.

    No license under any patent, copyright, trade secret or other
    intellectual property right is granted to or conferred upon you by
    disclosure or delivery of the Materials, either expressly, by
    implication, inducement, estoppel or otherwise.  Any license under such
    intellectual property rights must be express and approved by Intel in
    writing.
*/
#include <mpi.h>
#include <iostream>
#include <sstream>
#include <immintrin.h>
#include <malloc.h>
#include <omp.h>
#include <cmath>
#include <vector>
#include "matrix.hpp"
using namespace std;

/** Test class to excecute dgemm in parallel using MPI and OMP threads
 *
 * Usage: mpirun -np MPI mkl_gemm_hybrid M N K
 *
 * OpenMP environments
 * export MKL_DYNAMIC=false
 * export MKL_NUM_THREADS=MKL
 * export OMP_NESTED=true
 * export OMP_NUM_THREADS=OMP,MKL
 * export OMP_PLACES=threads
 * export OMP_PROC_BIND=spread,close
 * export KMP_HOT_TEAMS_MODE=1
 * export KMP_HOT_TEAMS_MAX_LEVEL=2
 */
int main(int argc, char **argv)
{

  int mpi_mode=0;
  int rank;
  int np;
  MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&mpi_mode);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&np);

  int M=(argc>1)?atoi(argv[1]):256; 
  int N=(argc>2)?atoi(argv[2]):M; 
  int K=(argc>3)?atoi(argv[3]):M; 

  int ncores=omp_get_max_threads(),mkl_threads=1;
#pragma omp parallel
#pragma omp master
    mkl_threads=omp_get_max_threads();

  typedef Matrix<double,64> matrix64;

  std::vector<matrix64*> a;
  std::vector<matrix64*> b;
  std::vector<matrix64*> c;
  a.resize(ncores,nullptr);
  b.resize(ncores,nullptr);
  c.resize(ncores,nullptr);

//#pragma omp parallel for
//  for(int ip=0; ip<ncores; ++ip)
#pragma omp parallel 
  {
    int ip=omp_get_thread_num();
    a[ip]=new matrix64(M,K);
    b[ip]=new matrix64(K,N);
    c[ip]=new matrix64(M,N);
    a[ip]->init(3,-2,1);
    b[ip]->init(-2,1,3);
    c[ip]->init(0,0,0);
   
    multiply_mkl(*a[ip],*b[ip],*c[ip]);
  }

  double t0,t1;
  double tconst=2.0*M*N*K*ncores*np*1.0e-9;

  int niters= 32768/std::pow(M/32,3)+10;
  t0=dsecnd();
//#pragma omp parallel for
//  for(int ip=0; ip<ncores; ++ip)
#pragma omp parallel 
  {
    int ip=omp_get_thread_num();
    matrix64& aref(*a[ip]);
    matrix64& bref(*b[ip]);
    matrix64& cref(*c[ip]);

    for(int i=0; i<niters; i+=2)
    {
      multiply_mkl(aref,bref,cref);
      multiply_mkl(bref,cref,aref);
    }
  }
  t1=dsecnd()-t0;

  double tmax=0.0;
  MPI_Allreduce(&t1,&tmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  if(rank==0)
  {
    printf("DGEMM %4d %4d %4d %3d %3d %3d total(sec)= %.4f per-iter(msec) %.4f GFLOPS= %.4f\n",
        M,N,K,np,ncores,mkl_threads,tmax,tmax/niters*1000.,tconst*niters/tmax);
  }

  for(int ip=0; ip<ncores; ++ip) delete a[ip];
  for(int ip=0; ip<ncores; ++ip) delete b[ip];
  for(int ip=0; ip<ncores; ++ip) delete c[ip];

  MPI_Finalize();
  return 0;
}
