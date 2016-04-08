//////////////////////////////////////////////////////////////////
// (c) Copyright 2008-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file multidet.cpp
 * @brief Test codes for multidets
 */
#include <omp.h>
#include <mpi.h>
#include <complex>
#include <iostream>
#include <matrix.hpp>
#include <fft/fft.h>
#include <MKLRandom.h>
#include <getopt.h>
using namespace std;

template<typename T>
inline T abs_diff(std::complex<T> a, std::complex<T> b)
{
  return abs(a.real()-b.real())+abs(a.imag()-b.imag());
}

template<typename T>
inline T abs_diff(std::complex<T>& a, T b)
{
  return abs(a.real()-b);
}

template<typename T>
inline T verify(std::complex<T>* restrict in, std::complex<T>* restrict out,int n,T scale)
{
  T err=T();
  for(int i=0; i<n; ++i)
  {
    err += abs_diff(in[i],scale*out[i]);
  }
  return err;
}

using namespace fft;

int main(int argc, char** argv)
{

  int mpi_mode=0;
  int rank;
  int np;
  MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&mpi_mode);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&np);

  int howmany=1;
  int nsteps=100;
  int ngrid[]={64,64,64};
  int opt;
  bool debug=false;

  while((opt = getopt(argc, argv, "hdg:x:y:z:i:s:p:")) != -1)
  {
    switch(opt)
    {
    case 'h':
      printf("[-g grid| -x grid_x -y grid_y -z grid_z] -s states -i iterations\n");
      return 1;
    case 'd':
      debug=true;
      break;
    case 'g':
      ngrid[0]=ngrid[1]=ngrid[2]=atoi(optarg);
      break;
    case 'x':
      ngrid[0]=atoi(optarg);
      break;
    case 'y':
      ngrid[1]=atoi(optarg);
      break;
    case 'z':
      ngrid[2]=atoi(optarg);
      break;
    case 'i':
      nsteps=atoi(optarg);
      break;
    }
  }

  typedef double real_type;
  typedef complex<double> complex_type;
  const unsigned int fft_id=FFTMKL_ENG;
  double wtime=0.0;

  int fft_dim=ngrid[0]*ngrid[1]*ngrid[2];
  Matrix<complex_type,64> out(1,fft_dim);
  MKLRandom<double> myrand;
  myrand.generate_uniform(reinterpret_cast<double*>(out.data()),fft_dim*2);

  int ncores=omp_get_max_threads(),mkl_threads=1;
#pragma omp parallel
#pragma omp master
    mkl_threads=omp_get_max_threads();

  typedef Matrix<complex_type,64> matrix64;
  typedef fft_engine_base<double,3,FFTMKL_ENG> fft_eng_type;
  std::vector<fft_eng_type*> fft3d_engines(ncores);
  std::vector<matrix64*> inout(ncores);

  //initialization
#pragma omp parallel
  {
    int ip=omp_get_thread_num();
    fft_eng_type *fft3d=new fft_engine_base<double,3,FFTMKL_ENG>(1);
    matrix64* in=new matrix64(howmany,fft_dim);

    for(int i=0; i<howmany; ++i)
      std::copy(out.data(),out.data()+fft_dim,(*in)[i]);

    fft3d->create_plan(ngrid,howmany,in->data(),out.cols_max());
    fft3d->execute_fft(in->data());
    fft3d->execute_ifft(in->data());

    inout[ip]=in;
    fft3d_engines[ip]=fft3d;
  }

  double t0=dsecnd();
#pragma omp parallel 
  {
    int ip=omp_get_thread_num();

    matrix64& in=*inout[ip];
    fft_eng_type *fft3d=fft3d_engines[ip];

    for(int step=0; step<nsteps; ++step)
    {
      fft3d->execute_fft(in.data());
      fft3d->execute_ifft(in.data());
    }
  }
  double t1=dsecnd();
  double dt=t1-t0;

  double tsum=0.0;
  MPI_Allreduce(&dt,&wtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

  //print ngrd,howmany,omp,mkl,timing
  if(rank==0)
  {
    printf("FFT3D %5d %5d %4d %4d %4d total(sec)= %.3f per-call(sec) =%.3f Rate= %.3f\n",ngrid[0],howmany, np, ncores, mkl_threads,
        wtime,wtime/(nsteps*howmany),np*ncores*nsteps*howmany/wtime);
  }

  for(int i=0; i<ncores; ++i) delete inout[i];
  for(int i=0; i<ncores; ++i) delete fft3d_engines[i];
  return 0;
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1770 $   $Date: 2007-02-17 17:45:38 -0600 (Sat, 17 Feb 2007) $
 * $Id: OrbitalBase.h 1770 2007-02-17 23:45:38Z jnkim $
 ***************************************************************************/
