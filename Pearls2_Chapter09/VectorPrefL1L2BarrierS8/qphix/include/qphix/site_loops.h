#ifndef QPHIX_SITE_LOOPS_H
#define QPHIX_SITE_LOOPS_H

#include <qphix/geometry.h>
#include <qphix/comm.h>

#include <omp.h>
#include "qphix/print_utils.h"
#define MAXV 16

namespace QPhiX
{ 

#if defined (__GNUG__) && !defined (__INTEL_COMPILER)
  static double new_norm_array[240][MAXV] __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN)));
  static double new_iprod_array[240][2][MAXV]  __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
  __declspec(align(QPHIX_LLC_CACHE_ALIGN)) static double new_norm_array[240][MAXV];
  __declspec(align(QPHIX_LLC_CACHE_ALIGN)) static double new_iprod_array[240][2][MAXV];
#endif

  template<typename FT, int V, int S, bool compress,
	   typename SpinorFunctor >
  void siteLoopNoReduction( SpinorFunctor& theFunctor,
			    const Geometry<FT,V,S,compress>& geom, 
			    int n_blas_simt) 
  {
    
    
    const int n_simt = geom.getNSIMT();
    const int n_cores = geom.getNumCores();
    const int Nxh = geom.Nxh();
    const int Ny  = geom.Ny();
    const int Nz  = geom.Nz();
    const int Nt  = geom.Nt();
    const int Pxy = geom.getPxy();
    const int Pxyz = geom.getPxyz();
    const int Nvecs = geom.nVecs();
    
    // This is the total number of spinor vectors
    const int n_soavec = (Nxh*Ny*Nz*Nt)/S;
#pragma omp parallel
    {
      // Self ID
      int tid = omp_get_thread_num();
      int cid = tid / n_simt;
      int smtid = tid - n_simt*cid;
      int num_bthreads = n_cores*n_blas_simt;
      
      if( smtid < n_blas_simt ) { 
	int btid = smtid  + n_blas_simt*cid;
	int n_soavec_per_core = n_soavec / n_cores;
	if ( n_soavec % n_cores != 0) n_soavec_per_core++;
	
	int low = cid*n_soavec_per_core;
	int next =(cid +1)*n_soavec_per_core;
	int hi = n_soavec < next ? n_soavec : next;
	
	// Loop over spinors. Local 
	for(int spinor_i = low+smtid; spinor_i < hi;  spinor_i+= n_blas_simt ) {
	  
	  // Decompose spinor_i into vec, ybase, z and t coords
	  int t1 = spinor_i / Nvecs;
	  int vec = spinor_i - Nvecs*t1;
	  
	  int t2 = t1 / Ny;
	  int y = t1 - Ny*t2;
	  
	  int t = t2/Nz;
	  int z = t2 -Nz*t;
	  
	  // ztbase for padding
	  int ztbase = Pxy*z + Pxyz*t;
	  int block = vec + Nvecs*y + ztbase;
	  theFunctor.func(block);
	  
	}
      }
    }
  }

 template<typename FT, int V, int S, bool compress,
	   typename Reduce1Functor >
  void siteLoop1Reduction( Reduce1Functor& theFunctor,
			   double& reduction,
			   const Geometry<FT,V,S,compress>& geom, 
			   int n_blas_simt) 
  {
    
    reduction = (double)0;

    const int n_simt = geom.getNSIMT();
    const int n_cores = geom.getNumCores();
    const int Nxh = geom.Nxh();
    const int Ny  = geom.Ny();
    const int Nz  = geom.Nz();
    const int Nt  = geom.Nt();
    const int Pxy = geom.getPxy();
    const int Pxyz = geom.getPxyz();
    const int Nvecs = geom.nVecs();
    
    // This is the total number of spinor vectors
    const int n_soavec = (Nxh*Ny*Nz*Nt)/S;
#pragma omp parallel
    {
      // Self ID
      int tid = omp_get_thread_num();
      int cid = tid / n_simt;
      int smtid = tid - n_simt*cid;
      int num_bthreads = n_cores*n_blas_simt;
      
      if( smtid < n_blas_simt ) { 
	int btid = smtid  + n_blas_simt*cid;
	
	// Each thread zeroes
	for(int s=0; s < V; s++){ 
	  new_norm_array[btid][s] = (double)0;
	}

	int n_soavec_per_core = n_soavec / n_cores;
	if ( n_soavec % n_cores != 0) n_soavec_per_core++;
	
	int low = cid*n_soavec_per_core;
	int next =(cid +1)*n_soavec_per_core;
	int hi = n_soavec < next ? n_soavec : next;
	
	// Loop over spinors. Local 
	for(int spinor_i = low+smtid; spinor_i < hi;  spinor_i+= n_blas_simt ) {
	  
	  // Decompose spinor_i into vec, ybase, z and t coords
	  int t1 = spinor_i / Nvecs;
	  int vec = spinor_i - Nvecs*t1;
	  
	  int t2 = t1 / Ny;
	  int y = t1 - Ny*t2;
	  
	  int t = t2/Nz;
	  int z = t2 -Nz*t;
	  
	  // ztbase for padding
	  int ztbase = Pxy*z + Pxyz*t;
	  int block = vec + Nvecs*y + ztbase;
	  theFunctor.func(block, new_norm_array[btid]);
	  
	}
      }
    }

    // OK add up the norm array. Accumulate all the BLAS threads results onto 
    // BLAS Threads 0s.
    reduction=(double)0;
    for(int btid=0; btid < n_cores*n_blas_simt; btid++) { 
      for(int s=0; s<V; s++){ 
	reduction += new_norm_array[btid][s];
      }
    }
    CommsUtils::sumDouble(&reduction);
  }


 template<typename FT, int V, int S, bool compress,
	   typename Reduce2Functor >
  void siteLoop2Reductions( Reduce2Functor& theFunctor,
			    double reduction[2],
			    const Geometry<FT,V,S,compress>& geom, 
			    int n_blas_simt) 
	   {
    

    const int n_simt = geom.getNSIMT();
    const int n_cores = geom.getNumCores();
    const int Nxh = geom.Nxh();
    const int Ny  = geom.Ny();
    const int Nz  = geom.Nz();
    const int Nt  = geom.Nt();
    const int Pxy = geom.getPxy();
    const int Pxyz = geom.getPxyz();
    const int Nvecs = geom.nVecs();
    
    // This is the total number of spinor vectors
    const int n_soavec = (Nxh*Ny*Nz*Nt)/S;
#pragma omp parallel
    {
      // Self ID
      int tid = omp_get_thread_num();
      int cid = tid / n_simt;
      int smtid = tid - n_simt*cid;
      int num_bthreads = n_cores*n_blas_simt;
      
      if( smtid < n_blas_simt ) { 
	int btid = smtid  + n_blas_simt*cid;
	
	// Each thread zeroes
#pragma simd
	for(int s=0; s < V; s++){ 
	  new_iprod_array[btid][0][s] = (double)0;
	  new_iprod_array[btid][1][s]= (double)0;
	}

	int n_soavec_per_core = n_soavec / n_cores;
	if ( n_soavec % n_cores != 0) n_soavec_per_core++;
	
	int low = cid*n_soavec_per_core;
	int next =(cid +1)*n_soavec_per_core;
	int hi = n_soavec < next ? n_soavec : next;
	
	// Loop over spinors. Local 
	for(int spinor_i = low+smtid; spinor_i < hi;  spinor_i+= n_blas_simt ) {
	  
	  // Decompose spinor_i into vec, ybase, z and t coords
	  int t1 = spinor_i / Nvecs;
	  int vec = spinor_i - Nvecs*t1;
	  
	  int t2 = t1 / Ny;
	  int y = t1 - Ny*t2;
	  
	  int t = t2/Nz;
	  int z = t2 -Nz*t;
	  
	  // ztbase for padding
	  int ztbase = Pxy*z + Pxyz*t;
	  int block = vec + Nvecs*y + ztbase;
	  theFunctor.func(block, new_iprod_array[btid][0], new_iprod_array[btid][1]);
	  
	}
      }
    }

    // Horizontally sum the btid=0 result.
    reduction[0] = (double)0;
    reduction[1] = (double)0;

    for(int btid=0; btid < n_cores*n_blas_simt; btid++) { 
      for(int s=0; s < S; s++){ 
	reduction[0] += new_iprod_array[btid][0][s];
      }
      for(int s=0; s < S; s++){ 
	reduction[1] += new_iprod_array[btid][1][s];
      }

    }
        // DO A GLOBAL SUM HERE

    CommsUtils::sumDoubleArray(reduction,2);


  }



}; // Namespace

#endif
