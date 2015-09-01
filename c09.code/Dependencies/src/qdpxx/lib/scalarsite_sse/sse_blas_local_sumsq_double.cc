// $Id: sse_blas_local_sumsq_double.cc,v 1.2 2008-06-19 13:50:20 bjoo Exp $

/*! @file
 *  @brief Generic Scalar VAXPY routine
 *
 */

#include "scalarsite_sse/sse_blas_local_sumsq_double.h"

namespace QDP {

#include <xmmintrin.h>


// (Vector) out = (Scalar) (*scalep) * (Vector) InScale + (Vector) Add
// #define DEBUG_VAXPY_DOUBLE
  void local_sumsq4(REAL64 *sum, REAL64 *vecptr, int n_4spin)
  {

    // Initialize the 4 sums to zero. Use _mm_setzero_pd() rather than explicit xor
    // Apparently we dont need volatile then.
    __m128d sum1 = _mm_setzero_pd();
    __m128d sum2 = _mm_setzero_pd();
    __m128d sum3 = _mm_setzero_pd();
    __m128d sum4 = _mm_setzero_pd();

  __m128d tmp1;
  __m128d tmp2;
  __m128d tmp3;
  __m128d tmp4;
  __m128d tmp5;
  __m128d tmp6;
  __m128d tmp7;
  __m128d tmp8;

  double *in=vecptr;

  for(int i=0; i < n_4spin; i++) { 

    tmp1=_mm_load_pd(in);
    tmp5=_mm_mul_pd(tmp1,tmp1);
    sum1=_mm_add_pd(sum1,tmp5);

    tmp2=_mm_load_pd(in+2);
    tmp6=_mm_mul_pd(tmp2,tmp2);
    sum2=_mm_add_pd(sum2,tmp6);

    tmp3=_mm_load_pd(in+4);
    tmp7=_mm_mul_pd(tmp3,tmp3);
    sum3=_mm_add_pd(sum3,tmp7);

    tmp4=_mm_load_pd(in+6);
    tmp8=_mm_mul_pd(tmp4,tmp4);
    sum4=_mm_add_pd(sum4,tmp8);

    tmp1=_mm_load_pd(in+8);
    tmp5=_mm_mul_pd(tmp1,tmp1);
    sum1=_mm_add_pd(sum1,tmp5);

    tmp2=_mm_load_pd(in+10);
    tmp6=_mm_mul_pd(tmp2,tmp2);
    sum2=_mm_add_pd(sum2,tmp6);

    tmp3=_mm_load_pd(in+12);
    tmp7=_mm_mul_pd(tmp3,tmp3);
    sum3=_mm_add_pd(sum3,tmp7);

    tmp4=_mm_load_pd(in+14);
    tmp8=_mm_mul_pd(tmp4,tmp4);
    sum4=_mm_add_pd(sum4,tmp8);

    tmp1=_mm_load_pd(in+16);
    tmp5=_mm_mul_pd(tmp1,tmp1);
    sum1=_mm_add_pd(sum1,tmp5);

    tmp2=_mm_load_pd(in+18);
    tmp6=_mm_mul_pd(tmp2,tmp2);
    sum2=_mm_add_pd(sum2,tmp6);

    tmp3=_mm_load_pd(in+20);
    tmp7=_mm_mul_pd(tmp3,tmp3);
    sum3=_mm_add_pd(sum3,tmp7);

    tmp4=_mm_load_pd(in+22);
    tmp8=_mm_mul_pd(tmp4,tmp4);
    sum4=_mm_add_pd(sum4,tmp8);

    in+=24;

  }
  
  sum1 = _mm_add_pd(sum1,sum2);
  sum3 = _mm_add_pd(sum3,sum4);
  sum1 = _mm_add_pd(sum1,sum3);

  tmp1 = _mm_shuffle_pd(sum1, sum1, 0x1);
  sum1 = _mm_add_pd(tmp1,sum1);
  _mm_storel_pd(sum,sum1);
  
}



} // namespace QDP;


