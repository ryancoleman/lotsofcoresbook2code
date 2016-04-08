// $Id: sse_blas_vaxpy4_double.cc,v 1.5 2009-07-14 20:08:42 bjoo Exp $

/*! @file
 *  @brief Generic Scalar VAXPY routine
 *
 */

#include "scalarsite_sse/sse_blas_vaxpy4_double.h"

namespace QDP {


#include <xmmintrin.h>
#include "scalarsite_sse/sse_prefetch.h"

#ifndef L2BY2
#define L2BY2 1365          /* L2 / 2 in SPINORS */
#endif


void vaxpy4(REAL64 *Out,REAL64 *scalep,REAL64 *InScale, int n_4spin)
{
  __m128d scalar;
  __m128d tmp1;
  __m128d tmp2;
  __m128d tmp3;
   __m128d in1;
  __m128d add1;
  __m128d in2;
  __m128d add2;
  __m128d in3;
  __m128d add3;
  __m128d in4;
  __m128d add4;
  __m128d out1;
  __m128d out2;
  __m128d out3;

  // Load the scalar into low bytes of scalar
  scalar = _mm_load_sd(scalep);
  
  // cross components into tmp 
  // Zero tmp
  tmp1 = _mm_setzero_pd();
  tmp1 = _mm_shuffle_pd(scalar, scalar, 0x1);
  scalar = _mm_add_pd(scalar, tmp1);

  double *in_p=InScale;
  double *out_p=Out;

  if( n_4spin < L2BY2 ) { 
    
    // Less than L2 size. 
    // Out p is sequential read/write: USE PREFETCHW
    // in_p is sequential read only: USE PREFETCH + HW PREFETCHER
    PREFETCH(((const char *)in_p)+16);

    for(int i=0; i < 3*n_4spin; i++) { 
      PREFETCHW(((const char *)out_p)+16);

      add1 = _mm_load_pd(out_p);
      in1  = _mm_load_pd(in_p);
      tmp1 = _mm_mul_pd(scalar, in1);
      out1 = _mm_add_pd(tmp1,add1);
      _mm_store_pd(out_p, out1);
      
      add2 = _mm_load_pd(out_p+2);
      in2  = _mm_load_pd(in_p+2);
      tmp2 = _mm_mul_pd(scalar, in2);
      out2 = _mm_add_pd(tmp2,add2);
      _mm_store_pd(out_p+2, out2);
      
      add3 = _mm_load_pd(out_p+4);
      in3 = _mm_load_pd(in_p+4);
      tmp3 = _mm_mul_pd(scalar, in3);
      out3 = _mm_add_pd(tmp3,add3);
      _mm_store_pd(out_p+4, out3);
      
      add1 = _mm_load_pd(out_p+6);
      in1  = _mm_load_pd(in_p+6);
      tmp1 = _mm_mul_pd(scalar, in1);
      out1 = _mm_add_pd(tmp1,add1);
      _mm_store_pd(out_p+6, out1);
      
      out_p+=8; in_p+=8;
      
    }
  }
  else {
    // > L2BY2 Size 
    // Out_p is sequential read/write: PREFETCHNTA
    // In_p is sequential read: PREFETCHNTA
    for(int i=0; i < 3*n_4spin; i++) { 

      PREFETCHNTA(((const char *)out_p)+56);
      PREFETCHNTA(((const char *)in_p)+56);
      add1 = _mm_load_pd(out_p);
      in1  = _mm_load_pd(in_p);
      tmp1 = _mm_mul_pd(scalar, in1);
      out1 = _mm_add_pd(tmp1,add1);
      _mm_store_pd(out_p, out1);
      
      add2 = _mm_load_pd(out_p+2);
      in2  = _mm_load_pd(in_p+2);
      tmp2 = _mm_mul_pd(scalar, in2);
      out2 = _mm_add_pd(tmp2,add2);
      _mm_store_pd(out_p+2, out2);
      
      add3 = _mm_load_pd(out_p+4);
      in3 = _mm_load_pd(in_p+4);
      tmp3 = _mm_mul_pd(scalar, in3);
      out3 = _mm_add_pd(tmp3,add3);
      _mm_store_pd(out_p+4, out3);
      
      add1 = _mm_load_pd(out_p+6);
      in1  = _mm_load_pd(in_p+6);
      tmp1 = _mm_mul_pd(scalar, in1);
      out1 = _mm_add_pd(tmp1,add1);
      _mm_store_pd(out_p+6, out1);
      
      out_p+=8; in_p+=8;
      
    }

  }
}



void vaxpyz4(REAL64 *Out,REAL64 *scalep,REAL64 *InScale, REAL64 *Add,int n_4vec)
{
 __m128d scalar;
  __m128d tmp1;
  __m128d tmp2;
  __m128d tmp3;
   __m128d in1;
  __m128d add1;
  __m128d in2;
  __m128d add2;
  __m128d in3;
  __m128d add3;
  __m128d in4;
  __m128d add4;
  __m128d out1;
  __m128d out2;
  __m128d out3; 

  // Load the scalar into low bytes of scalar
  scalar = _mm_load_sd(scalep);
  
  // cross components into tmp 
  // Zero tmp
  tmp1 = _mm_setzero_pd();

  tmp1 = _mm_shuffle_pd(scalar, scalar, 0x1);
  scalar = _mm_add_pd(scalar, tmp1);

  // Do n_3vec 3vectors.
  double *in_p=InScale;
  double *add_p=Add;
  double *out_p=Out;


  if( n_4vec < L2BY2 ) { 
    
    // Less than L2
    // in_p sequential read only   : 
    // add_p sequential read only  :
    // out_p sequential read/write :
    PREFETCH(((const char *)add_p)+16);
    PREFETCH(((const char *)in_p)+16);

    for(int i=0; i < 3*n_4vec; i++) { 
      PREFETCHW(((const char *)out_p)+16);
      add1 = _mm_load_pd(add_p);
      in1  = _mm_load_pd(in_p);
      tmp1 = _mm_mul_pd(scalar, in1);
      out1 = _mm_add_pd(tmp1,add1);
      _mm_store_pd(out_p, out1);
      
      add2 = _mm_load_pd(add_p+2);
      in2  = _mm_load_pd(in_p+2);
      tmp2 = _mm_mul_pd(scalar, in2);
      out2 = _mm_add_pd(tmp2,add2);
      _mm_store_pd(out_p+2, out2);
      
      add3 = _mm_load_pd(add_p+4);
      in3 = _mm_load_pd(in_p+4);
      tmp3 = _mm_mul_pd(scalar, in3);
      out3 = _mm_add_pd(tmp3,add3);
      _mm_store_pd(out_p+4, out3);
      
      add1 = _mm_load_pd(add_p+6);
      in1  = _mm_load_pd(in_p+6);
      tmp1 = _mm_mul_pd(scalar, in1);
      out1 = _mm_add_pd(tmp1,add1);
      _mm_store_pd(out_p+6, out1);
            
      out_p+=8; in_p+=8; add_p+=8;
      
    }
  }
  else {

    for(int i=0; i < 3*n_4vec; i++) { 

      PREFETCHNTA(((const char *)add_p)+56);
      PREFETCHNTA(((const char *)in_p)+56);

      add1 = _mm_load_pd(add_p);
      in1  = _mm_load_pd(in_p);
      tmp1 = _mm_mul_pd(scalar, in1);
      out1 = _mm_add_pd(tmp1,add1);
      _mm_stream_pd(out_p, out1);
      
      
      add2 = _mm_load_pd(add_p+2);
      in2  = _mm_load_pd(in_p+2);
      tmp2 = _mm_mul_pd(scalar, in2);
      out2 = _mm_add_pd(tmp2,add2);
      _mm_stream_pd(out_p+2, out2);
      
      add3 = _mm_load_pd(add_p+4);
      in3 = _mm_load_pd(in_p+4);
      tmp3 = _mm_mul_pd(scalar, in3);
      out3 = _mm_add_pd(tmp3,add3);
      _mm_stream_pd(out_p+4, out3);
      
      
      add1 = _mm_load_pd(add_p+6);
      in1  = _mm_load_pd(in_p+6);
      tmp1 = _mm_mul_pd(scalar, in1);
      out1 = _mm_add_pd(tmp1,add1);
      _mm_stream_pd(out_p+6, out1);
      
      out_p+=8; in_p+=8; add_p+=8;
      
    }
  }
}


} // namespace QDP;

