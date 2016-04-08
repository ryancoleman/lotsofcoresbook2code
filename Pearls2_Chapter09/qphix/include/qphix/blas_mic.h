#ifndef QPHIX_BLAS_MIC_H
#define QPHIX_BLAS_MIC_H

#include <immintrin.h>

namespace QPhiX { 

template<>
inline
  void xmyNorm2Spinor<float,16>(float* restrict res, float* restrict x, float* restrict y, double& n2res, int n, int n_cores, int n_simt, int n_blas_simt) 
{
#if defined (__GNUG__) && !defined (__INTEL_COMPILER)
  double norm2res=0 __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
  __declspec(align(QPHIX_LLC_CACHE_ALIGN)) double norm2res=0;
#endif
  const int veclen = 16;
#pragma omp parallel shared(norm_array) 
 {
    // Self ID
    int tid = omp_get_thread_num();

    
    
    int cid = tid / n_simt;
    int smtid = tid - n_simt*cid;
    int bthreads = n_cores*n_blas_simt;
    int nvec = n/veclen;
    
    if ( smtid < n_blas_simt ) { 
      int btid = smtid  + n_blas_simt*cid;
      for(int i=0; i < 8;i++) norm_array[btid][i] = 0; // Every thread zeroes      
      // No of vectors per thread
      int n_per_thread= nvec / bthreads;
      if ( nvec % bthreads != 0 ) n_per_thread++;


      int low=btid*n_per_thread ;
      int next = (btid + 1)*n_per_thread;
      int hi= nvec  < next ? nvec  : next;
      
      low *=veclen;
      hi  *=veclen;

      // All sums in DP. Need to convert 16x floats into 2 8x doubles
      __m512d   local_norm = _mm512_setzero_pd();

      __m512d  lo_part = _mm512_setzero_pd();
      __m512d hi_part = _mm512_setzero_pd();

      __m512 xvec;
      __m512 yvec;
      __m512 resvec;
      __m512 permute;

      // L1dist is the distance in cache-lines to fetch ahead from L1
      //      const int L1dist = 1;
      const int L2dist = 3; 

      for(int i = low; i < hi; i+=16) {
	//	_mm_prefetch((const char *)( &x[i + L1dist*16]), _MM_HINT_T0);
	_mm_prefetch((const char *)( &x[i + L2dist*384]), _MM_HINT_T1);

       	//_mm_prefetch((const char *) (&y[i + L1dist*16]), _MM_HINT_T0);
	_mm_prefetch((const char *) (&y[i + L2dist*384]), _MM_HINT_T1);

	xvec = _mm512_load_ps(&x[i]);
	yvec = _mm512_load_ps(&y[i]);
	
	resvec = _mm512_sub_ps(xvec,yvec);

	permute = _mm512_permute4f128_ps(resvec, _MM_PERM_DCDC);
	_mm512_storenrngo_ps(&res[i], resvec);


	hi_part = _mm512_cvtpslo_pd(resvec);
	local_norm = _mm512_fmadd_pd( hi_part, hi_part, local_norm);

	lo_part = _mm512_cvtpslo_pd(permute);
	local_norm = _mm512_fmadd_pd( lo_part, lo_part, local_norm);
      }
      
      _mm512_store_pd( &norm_array[btid][0], local_norm );
    }

 }

 __m512d accum=_mm512_setzero_pd();
 for(int i=0; i < n_cores*n_blas_simt; i++) {
   __m512d tmp = _mm512_load_pd( &norm_array[i][0] );
   accum = _mm512_add_pd( accum, tmp);
 }
 norm2res = _mm512_reduce_add_pd(accum);

  CommsUtils::sumDouble(&norm2res);
  n2res = norm2res;
}

 template<>
inline
void rmammpNorm2rxpap<float,16>(float* restrict r, const float ar, float* restrict  mmp, double& cp, float* restrict  x, float* restrict p,int n, int n_cores, int n_simt, int n_blas_simt) 
{

#if defined (__GNUG__) && !defined (__INTEL_COMPILER)
  double norm2res=0 __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
  __declspec(align(QPHIX_LLC_CACHE_ALIGN)) double norm2res=0;
#endif

  const int veclen = 16;

#pragma omp parallel shared(norm_array) 
 {
    // Self ID
    int tid = omp_get_thread_num();

    
    
    int cid = tid / n_simt;
    int smtid = tid - n_simt*cid;
    int bthreads = n_cores*n_blas_simt;
    int nvec = n/veclen;
    
    if ( smtid < n_blas_simt ) { 
      int btid = smtid  + n_blas_simt*cid;
      for(int i=0; i < 8;i++) norm_array[btid][i] = 0; // Every thread zeroes      
      // No of vectors per thread
      int n_per_thread= nvec / bthreads;
      if ( nvec % bthreads != 0 ) n_per_thread++;


      int low=btid*n_per_thread ;
      int next = (btid + 1)*n_per_thread;
      int hi= nvec  < next ? nvec  : next;
      
      low *=veclen;
      hi  *=veclen;

      // All sums in DP. Need to convert 16x floats into 2 8x doubles
      __m512d   local_norm = _mm512_setzero_pd();

      __m512d  lo_part = _mm512_setzero_pd();
      __m512d hi_part = _mm512_setzero_pd();
      __m512 resvec;
      __m512 permute;

      __m512 arvec = _mm512_extload_ps(&ar, _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);
      __m512 rvec;
      __m512 mmpvec;
      __m512 xvec;
      __m512 pvec;


      // L1dist is the distance in cache-lines to fetch ahead from L1
      //      const int L1dist = 1;
      const int L2dist = 3; 

      for(int i = low; i < hi; i+=16) {
	_mm_prefetch((const char *)( &x[i + L2dist*384]), _MM_HINT_T1);
	_mm_prefetch((const char *) (&p[i + L2dist*384]), _MM_HINT_T1);


	// x[i] += ar * p[i]
	xvec = _mm512_load_ps(&x[i]);
	pvec = _mm512_load_ps(&p[i]);

	xvec = _mm512_fmadd_ps(arvec, pvec, xvec);
	_mm512_store_ps( &x[i], xvec);

	_mm_prefetch((const char *)( &r[i + L2dist*384]), _MM_HINT_T1);
	_mm_prefetch((const char *) (&mmp[i + L2dist*384]), _MM_HINT_T1);

	// r[i] -= ar* mmp[i]
	rvec = _mm512_load_ps(&r[i]);
	mmpvec = _mm512_load_ps(&mmp[i]);
	rvec = _mm512_fnmadd_ps( arvec, mmpvec, rvec);
	_mm512_store_ps( &r[i], rvec);


	permute = _mm512_permute4f128_ps(rvec, _MM_PERM_DCDC);
	lo_part = _mm512_cvtpslo_pd(rvec);
	local_norm = _mm512_fmadd_pd( lo_part, lo_part, local_norm);

	hi_part = _mm512_cvtpslo_pd(permute);
	local_norm = _mm512_fmadd_pd( hi_part, hi_part, local_norm);

      }
      
      _mm512_store_pd( &norm_array[btid][0], local_norm );
    }

 }

 __m512d accum=_mm512_setzero_pd();
 for(int i=0; i < n_cores*n_blas_simt; i++) {
   __m512d tmp = _mm512_load_pd( &norm_array[i][0] );
   accum = _mm512_add_pd( accum, tmp);
 }
 norm2res = _mm512_reduce_add_pd(accum);

 CommsUtils::sumDouble(&norm2res);
  cp = norm2res;
}


template<>
inline
double norm2Spinor<float,16>(float* restrict res, int n, int n_cores, int n_simt, int n_blas_simt)
{

#if defined (__GNUG__) && !defined (__INTEL_COMPILER)
  double norm2res=0 __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
  __declspec(align(QPHIX_LLC_CACHE_ALIGN)) double norm2res=0;
#endif

  const int veclen = 16;
#pragma omp parallel shared(norm_array)
  {
    // Self ID
    int tid = omp_get_thread_num();

    
    int cid = tid / n_simt;
    int smtid = tid - n_simt*cid;
    int bthreads = n_cores*n_blas_simt;
    int nvec = n/veclen;
    
    if ( smtid < n_blas_simt ) { 
      int btid = smtid  + n_blas_simt*cid;
      for(int i=0; i < 8; i++) norm_array[btid][i] = 0; // Every thread zeroes      
      // Divide N by the number cores
      int n_per_thread= nvec / bthreads;
      
      if ( nvec % bthreads != 0 ) n_per_thread++;
      
      int low=btid*n_per_thread ;
      int next = (btid + 1)*n_per_thread;
      int hi= nvec  < next ? nvec  : next;
      
      low *=veclen;
      hi  *=veclen;
      
      double lnorm=0;
      
      __m512 rvec;
      __m512 permute;
      __m512d lo_part;
      __m512d hi_part;
      __m512d local_norm;
      
      const int L1dist =4;
      const  int L2dist=3;
      
      for(int i = low; i < hi; i+=16) {
	// Prefetch
	//	_mm_prefetch((const char *)( &res[i + L1dist*16]), _MM_HINT_T0);
	_mm_prefetch((const char *)( &res[i + L2dist*384]), _MM_HINT_T1);

	// load
	rvec = _mm512_load_ps(&res[i]);
	
	// get high 8 into low 4 8 of permute
	permute = _mm512_permute4f128_ps(rvec, _MM_PERM_DCDC);
	
	lo_part = _mm512_cvtpslo_pd(rvec);
	local_norm = _mm512_fmadd_pd( lo_part, lo_part, local_norm);
	
	hi_part = _mm512_cvtpslo_pd(permute);
	local_norm = _mm512_fmadd_pd( hi_part, hi_part, local_norm);
	
      }
      _mm512_store_pd( &norm_array[btid][0], local_norm );
      
    }
    
  }
  
  __m512d accum=_mm512_setzero_pd();
  for(int i=0; i < n_cores*n_blas_simt; i++) {
    __m512d tmp = _mm512_load_pd( &norm_array[i][0] );
    accum = _mm512_add_pd( accum, tmp);
  }
  norm2res = _mm512_reduce_add_pd(accum);
  CommsUtils::sumDouble(&norm2res);
  return norm2res;
}



// Messed about with
template<>
inline
void aypx<float,16>(float a, float* restrict x, float* restrict y, int n, int n_cores, int threads_per_core, int blas_threads_per_core) 
{
  const int veclen = 16;
#pragma omp parallel
  {
    // Self ID
    int tid = omp_get_thread_num();
    int cid = tid / threads_per_core;
    int smtid = tid - threads_per_core*cid;
    int bthreads = n_cores*blas_threads_per_core;
    int nvec = n/veclen;
    
    if ( smtid < blas_threads_per_core ) { 
      __m512 avec = _mm512_extload_ps(&a, _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);

      int btid = smtid  + blas_threads_per_core*cid;
      
      // No of vectors per thread
      int n_per_thread= nvec / bthreads;
      
      // Add extra vector
      if ( nvec % bthreads != 0 ) n_per_thread++;
      
      // low vector
      int low=btid*n_per_thread ;
      
      // High vector 
      int next = (btid + 1)*n_per_thread;
      
      int hi= nvec  < next ? nvec  : next;
      
      low *=veclen;
      hi  *=veclen;
 
      __m512 yvec;
      __m512 xvec;

      const int L2dist = 3;
      for(int i = low; i < hi; i+=16) {
	_mm_prefetch((const char *)( &x[i + L2dist*384]), _MM_HINT_T1);
	_mm_prefetch((const char *)( &y[i + L2dist*384]), _MM_HINT_T1);
	yvec = _mm512_load_ps(&y[i]);
	xvec = _mm512_load_ps(&x[i]);
	yvec = _mm512_fmadd_ps(avec,yvec,xvec);
	_mm512_store_ps(&y[i], yvec);
      }
    }

  }
}

#if 1
// Messed about with
template<>
inline
void copySpinor<float,16>(float* restrict res, float* restrict src, int n, int n_cores, int threads_per_core, int blas_threads_per_core) 
{
  const int veclen = 16;
#pragma omp parallel
  {
    // Self ID
    int tid = omp_get_thread_num();
    int cid = tid / threads_per_core;
    int smtid = tid - threads_per_core*cid;
    int bthreads = n_cores*blas_threads_per_core;
    int nvec = n/veclen;
    
    if ( smtid < blas_threads_per_core ) { 

      int btid = smtid  + blas_threads_per_core*cid;
      
      // No of vectors per thread
      int n_per_thread= nvec / bthreads;
      
      // Add extra vector
      if ( nvec % bthreads != 0 ) n_per_thread++;
      
      // low vector
      int low=btid*n_per_thread ;
      
      // High vector 
      int next = (btid + 1)*n_per_thread;
      
      int hi= nvec  < next ? nvec  : next;
      
      low *=veclen;
      hi  *=veclen;
 
      __m512 vec;

      const int L2dist = 3;
      for(int i = low; i < hi; i+=16) {
	_mm_prefetch((const char *)( &src[i + L2dist*384]), _MM_HINT_T1);
	vec = _mm512_load_ps(&src[i]);
	_mm512_storenrngo_ps(&res[i],vec);
      }
    }

  }
}
#endif

} // Namespace

#endif

