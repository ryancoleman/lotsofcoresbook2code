// $Id: sse_blas_vaypx4_double.cc,v 1.3 2009-07-14 20:08:42 bjoo Exp $

/*! @file
 *  @brief Generic Scalar VAXPY routine
 *
 */

#include "scalarsite_sse/sse_blas_vaypx4_double.h"

namespace QDP {

#include <xmmintrin.h>

void vaypx4(REAL64 *Out,REAL64 *scalep,REAL64 *InScale, int n_4spin)
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

 
  // We are dealing with n_4 spin 4spinor. 
  // 4 spinors 

  for(int i=0; i < 4*3*2*n_4spin; i+=2) { 
    add1 = _mm_load_pd(out_p+i);
    in1  = _mm_load_pd(in_p+i);
    tmp1 = _mm_mul_pd(scalar, add1);
    out1 = _mm_add_pd(tmp1,in1);
    _mm_store_pd(out_p+i, out1);

  }
}




} // namespace QDP;

