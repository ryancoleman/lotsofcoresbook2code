// $Id: sse_blas_local_vcdot_real_double.cc,v 1.2 2008-06-19 13:50:20 bjoo Exp $

/*! @file
 *  @brief Generic Scalar VAXPY routine
 *
 */

#include <qdp.h>
#include <xmmintrin.h>
#include "scalarsite_sse/sse_blas_local_vcdot_real_double.h"
#include <iostream>

namespace QDP {

  typedef union { 
    double c[2];
    __m128d vec;
  } VDU;

  // Re < y^\dag , x > 
  // =  sum ( y.re x.re + y.im x.im )
  // =  sum ( y.re x.re ) + sum( y.im x.im )
  //
  // Load in [ x.re | x.im ]
  //         [ y.re | y.im ]
  // Make    [ x.re y.re | x.im y.im ]
  // accumulate sum
  // 
  // At the end do a single crossing:
  //
  //       [ sum (x.re y.re) | sum(x.im y.im) ]
  //     + [ sum (x.im y.im) | sum(x.re y.re) ]
  //  =    [ innerProdReal   | innerProdReal  ]
  // 
  void local_vcdot_real4(REAL64 *sum, REAL64 *y, REAL64* x,int n_4spin)
{
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
  __m128d tmp9;
  __m128d tmp10;
  __m128d tmp11;
  __m128d tmp12;

  double *x_p=x;
  double *y_p=y;

  for(int i=0; i < n_4spin; i++) { 

    tmp1 = _mm_load_pd(x_p);
    tmp2 = _mm_load_pd(y_p);
    tmp3 = _mm_mul_pd(tmp1,tmp2);
    sum1 = _mm_add_pd(sum1,tmp3);


    tmp4 = _mm_load_pd(x_p+2);
    tmp5 = _mm_load_pd(y_p+2);
    tmp6 = _mm_mul_pd(tmp4,tmp5);
    sum2 = _mm_add_pd(sum2,tmp6);
    
    tmp7 = _mm_load_pd(x_p+4);
    tmp8 = _mm_load_pd(y_p+4);
    tmp9 = _mm_mul_pd(tmp7,tmp8);
    sum3 = _mm_add_pd(sum3,tmp9);

    tmp10 = _mm_load_pd(x_p+6);
    tmp11 = _mm_load_pd(y_p+6);
    tmp12 = _mm_mul_pd(tmp10,tmp11);
    sum4 = _mm_add_pd(sum4,tmp12);



    tmp1 = _mm_load_pd(x_p+8);
    tmp2 = _mm_load_pd(y_p+8);
    tmp3 = _mm_mul_pd(tmp1,tmp2);
    sum1 = _mm_add_pd(sum1,tmp3);

    tmp4 = _mm_load_pd(x_p+10);
    tmp5 = _mm_load_pd(y_p+10);
    tmp6 = _mm_mul_pd(tmp4,tmp5);
    sum2 = _mm_add_pd(sum2,tmp6);
    
    tmp7 = _mm_load_pd(x_p+12);
    tmp8 = _mm_load_pd(y_p+12);
    tmp9 = _mm_mul_pd(tmp7,tmp8);
    sum3 = _mm_add_pd(sum3,tmp9);

    tmp10 = _mm_load_pd(x_p+14);
    tmp11 = _mm_load_pd(y_p+14);
    tmp12 = _mm_mul_pd(tmp10,tmp11);
    sum4 = _mm_add_pd(sum4,tmp12);



    tmp1 = _mm_load_pd(x_p+16);
    tmp2 = _mm_load_pd(y_p+16);
    tmp3 = _mm_mul_pd(tmp1,tmp2);
    sum1 = _mm_add_pd(sum1,tmp3);

    tmp4 = _mm_load_pd(x_p+18);
    tmp5 = _mm_load_pd(y_p+18);
    tmp6 = _mm_mul_pd(tmp4,tmp5);
    sum2 = _mm_add_pd(sum2,tmp6);
    
    tmp7 = _mm_load_pd(x_p+20);
    tmp8 = _mm_load_pd(y_p+20);
    tmp9 = _mm_mul_pd(tmp7,tmp8);
    sum3 = _mm_add_pd(sum3,tmp9);

    tmp10 = _mm_load_pd(x_p+22);
    tmp11 = _mm_load_pd(y_p+22);
    tmp12 = _mm_mul_pd(tmp10,tmp11);
    sum4 = _mm_add_pd(sum4,tmp12);
    
    x_p += 24; y_p+=24;

  }

  // Accumulate into 2 vecs
  sum1 = _mm_add_pd(sum1,sum2);
  sum3 = _mm_add_pd(sum3,sum4);

  // Accumulate into 1 vec
  sum1 = _mm_add_pd(sum1,sum3);

  // Cross the vector and add
  tmp1 = _mm_shuffle_pd(sum1, sum1, 0x1);
  sum1 = _mm_add_pd(tmp1,sum1);

  // Store either half
  _mm_storeh_pd(sum,sum1);
  
}



} // namespace QDP;

