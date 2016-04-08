#ifndef QPHIX_BLAS_C_H
#define QPHIX_BLAS_C_H

namespace QPhiX { 


template<typename FT, int veclen>
inline
void copySpinor(FT* restrict res, FT* restrict src, int n)
{

#pragma prefetch src
#pragma vector nontemporal (res)
#pragma omp parallel for
  for(int i=0; i < n; i++){ 
    res[i]=src[i];
  }
  
}

template<typename FT, int veclen>
inline
void xmyNorm2Spinor(FT* restrict res, FT* restrict x, FT* restrict y, double& n2res, int n) 
{
     
  double norm2res=0;

#if defined (__INTEL_COMPILER)  
#pragma prefetch x,y
#pragma unroll
#endif
#pragma omp parallel for reduction(+:norm2res)
  for(int i=0; i < n; i++){
      float tmpf = x[i]-y[i];
      res[i] = tmpf;
      double tmpd = (double)tmpf;
      norm2res += tmpd*tmpd;
  }
  
  CommsUtils::sumDouble(&norm2res);
  n2res = norm2res;
}

template<typename FT, int veclen>
inline
double norm2Spinor(FT* restrict res, int n)
{
  double norm2res=0;

#pragma prefetch res
#pragma omp parallel for reduction(+:norm2res)
  for(int i=0; i < n; i++){
    double tmpd=res[i];
    norm2res += tmpd*tmpd;
  }
  
  CommsUtils::sumDouble(&norm2res);
  return norm2res;
}

template<typename FT, int veclen>
inline
void rmammpNorm2rxpap(FT* restrict r, FT ar, FT* restrict  mmp, double& cp, FT* restrict  x, FT* restrict p,int n)
{
  double cp_internal = 0;


 #pragma prefetch x,p,r,mmp
#pragma omp parallel for reduction(+:cp_internal)
  for(int i=0; i < n; i++) { 
    x[i] += ar*p[i];
    float tmp = r[i] - ar*mmp[i];
    double tmpd=(double)tmp;
    cp_internal += tmpd * tmpd;
    r[i] = tmp;

  }
  
  
  CommsUtils::sumDouble(&cp_internal);
  cp= cp_internal;
}


// Original
template<typename FT, int veclen>
inline
void aypx(FT a, FT* restrict x, FT* restrict y, int n) 
{

#pragma prefetch x,y
#pragma omp parallel for
  for(int i=0; i < n; i++) {
    y[i] = a*y[i]+x[i];
  }
}





// Actually this is a MIC-ism (8 doubles=1 cacheline)
#if defined (__GNUG__) && !defined (__INTEL_COMPILER)
 static double norm_array[240][8] __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
 __declspec(align(QPHIX_LLC_CACHE_ALIGN)) static double norm_array[240][8];
#endif

template<typename FT, int veclen>
inline
  void xmyNorm2Spinor(FT* restrict res, FT* restrict x, FT* restrict y, double& n2res, int n, int n_cores, int n_simt, int n_blas_simt) 
{
 
#if defined (__GNUG__) && !defined (__INTEL_COMPILER)
  double norm2res __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN))) = 0;
#else
  __declspec(align(QPHIX_LLC_CACHE_ALIGN)) double norm2res = 0;
#endif

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

      double lnorm=0;
      

#pragma prefetch x,y,res
#pragma vector aligned(x,y,res)
      for(int i = low; i < hi; i++) {
	res[i]  = x[i] - y[i];
	double tmpd = (double)res[i];
	lnorm += (tmpd*tmpd);
      }
            
      norm_array[btid][0]=lnorm;

    }
 }
  
 for(int i=0; i < n_cores*n_blas_simt; i++) {
   norm2res += norm_array[i][0];
 }

  CommsUtils::sumDouble(&norm2res);
  n2res = norm2res;
}
 

template<typename FT, int veclen>
inline
double norm2Spinor(FT* restrict res, int n, int n_cores, int n_simt, int n_blas_simt)
{

#if defined (__GNUG__) && !defined (__INTEL_COMPILER)
  double norm2res __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN))) = 0;
#else
  __declspec(align(QPHIX_LLC_CACHE_ALIGN)) double norm2res = 0;
#endif

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

#pragma prefetch res
      for(int i = low; i < hi; i++) {
	double tmpd = (double) res[i]; 
	lnorm += (tmpd*tmpd);
      }
    
      norm_array[btid][0] = lnorm;
    }
   
  }
 
  for(int i=0; i < n_cores*n_blas_simt; i++) {
    norm2res += norm_array[i][0];
  }
  
  CommsUtils::sumDouble(&norm2res);
  return norm2res;
}

 
template<typename FT, int veclen>
inline
void rmammpNorm2rxpap(FT* restrict r, const FT ar, FT* restrict  mmp, double& cp, FT* restrict  x, FT* restrict p,int n, int n_cores, int n_simt, int n_blas_simt) 
{
  // 240 is max num threads...  -- fix this later 
 
#if defined (__GNUG__) && !defined (__INTEL_COMPILER)
  double norm2res __attribute__ ((aligned(QPHIX_LLC_CACHE_ALIGN))) = 0;
#else
  __declspec(align(QPHIX_LLC_CACHE_ALIGN)) double norm2res = 0;
#endif

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
      // Divide N by the number cores
      int n_per_thread= nvec / bthreads;
      
      if ( nvec % bthreads != 0 ) n_per_thread++;
      
      int low=btid*n_per_thread ;
      int next = (btid + 1)*n_per_thread;
      int hi= n  < next ? n  : next;
      
      low *=veclen;
      hi  *=veclen;

      double lnorm=0;

#pragma prefetch x,p,r,mmp
      for(int i = low; i < hi; i++) {
	x[i] += ar*p[i];
	float tmp = r[i] -ar*mmp[i];
	double tmpd=(double)tmp;
	r[i] = tmp;
	lnorm += tmpd * tmpd;
      }
      norm_array[btid][0]=lnorm;
    }  
  }

  for(int i=0; i < n_cores*n_blas_simt; i++) {
    norm2res += norm_array[i][0];
  }

  CommsUtils::sumDouble(&norm2res);
  cp = norm2res;
}


// Messed about with
template<typename FT, int veclen>
inline
  void copySpinor(FT* restrict res, FT* restrict src, int n, int n_cores, int threads_per_core, int blas_threads_per_core) 
{

#pragma omp parallel
  {
    int orig_nthreads=omp_get_num_threads();


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

#pragma prefetch src
#pragma vector aligned
#pragma vector nontemporal (res)
      for(int i = low; i < hi; i++) {
	res[i] = src[i];
      }
    }

  }
}

// Messed about with
template<typename FT, int veclen>
inline
  void aypx(FT a, FT* restrict x, FT* restrict y, int n, int n_cores, int threads_per_core, int blas_threads_per_core) 
{

#pragma omp parallel
  {
    int orig_nthreads=omp_get_num_threads();


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
      
#pragma simd
#pragma prefetch x,y
#pragma vector aligned
      for(int i = low; i < hi; i++) {
	y[i] = a * y[i] + x[i];
      }
    }

  }
}

};

#endif
