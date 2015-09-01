// $Id: sse_blas_local_vcdot_double.cc,v 1.2 2008-06-27 12:56:57 bjoo Exp $

/*! @file
 *  @brief Generic Scalar VAXPY routine
 *
 */

#include <xmmintrin.h>
#include "scalarsite_sse/sse_blas_local_vcdot_double.h"

namespace QDP {

#include "qdp_config.h"

#ifndef QDP_USE_SSE3

  /* SSE 2 */

#define CONJMUL(z,x,y)		\
  { \
    __m128d t1,t2,t3,t4; \
    t1 = _mm_mul_pd(x,y); \
    t2 = _mm_shuffle_pd(t1,t1,0x1); \
    t3 = _mm_shuffle_pd(y,y,0x1);\
    z = _mm_add_pd(t1,t2); \
    t2 = _mm_mul_pd(x,t3); \
    t3 = _mm_shuffle_pd(t2,t2,0x1); \
    t3 = _mm_sub_pd(t2,t3);	    \
    z= _mm_shuffle_pd(z,t3,0x2); \
  }

#define CONJMADD(z,x,y)				\
  { \
    __m128d t1,t2,t3,t4; \
    t1 = _mm_mul_pd(x,y); \
    t2 = _mm_shuffle_pd(t1,t1,0x1); \
    t3 = _mm_shuffle_pd(y,y,0x1);\
    t4 = _mm_add_pd(t1,t2); \
    t2 = _mm_mul_pd(x,t3); \
    t3 = _mm_shuffle_pd(t2,t2,0x1); \
    t3 = _mm_sub_pd(t2,t3); \
    t4= _mm_shuffle_pd(t4,t3,0x2); \
    z = _mm_add_pd(z,t4); \
  }

#else
#warning Using SSE3
  /* SSE 3 */
#include <pmmintrin.h>

#define CONJMUL(z,x,y)		\
  { \
    __m128d t1; \
    t1 = _mm_mul_pd((x),(y)); \
    (z) = _mm_hadd_pd(t1,t1);			\
    t1 = _mm_shuffle_pd((x),(x),0x1);\
    t1 = _mm_mul_pd((y),t1); \
    t1 = _mm_hsub_pd(t1,t1); \
    (z)= _mm_shuffle_pd((z),t1,0x2);		\
  }

#define CONJMADD(z,x,y)				\
  { \
    __m128d t1,t2; \
    t1 = _mm_mul_pd((x),(y)); \
    t1 = _mm_hadd_pd(t1,t1); \
    t2 = _mm_shuffle_pd((x),(x),0x1);\
    t2 = _mm_mul_pd((y),t2); \
    t2 = _mm_hsub_pd(t2,t2); \
    t1= _mm_shuffle_pd(t1,t2,0x2);		\
    (z) = _mm_add_pd((z),t1);			\
  }



#endif


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
  // then srore either half.
  void local_vcdot4(REAL64 *sum, REAL64 *y, REAL64* x,int n_4spin)
{
  // Use _mm_setzero_pd() to initialize the sums rather than xors
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

  
  double *x_p=x;
  double *y_p=y;


  for(int i=0; i < n_4spin; i++) { 
      
    tmp1 = _mm_load_pd(x_p);  // tmp1 = x
    tmp2 = _mm_load_pd(y_p);  // tmp2 = y

    CONJMADD(sum1,tmp1,tmp2);
    
    tmp3 = _mm_load_pd(x_p+2);  // tmp1 = x
    tmp4 = _mm_load_pd(y_p+2);  // tmp2 = y
    
    CONJMADD(sum2,tmp3,tmp4);

    tmp5 = _mm_load_pd(x_p+4);  // tmp1 = x
    tmp6 = _mm_load_pd(y_p+4);  // tmp2 = y
    
    CONJMADD(sum3,tmp5,tmp6);
    
    tmp7 = _mm_load_pd(x_p+6);  // tmp1 = x
    tmp8 = _mm_load_pd(y_p+6);  // tmp2 = y
    
    CONJMADD(sum4,tmp7,tmp8);
    
    tmp1 = _mm_load_pd(x_p+8);  // tmp1 = x
    tmp2 = _mm_load_pd(y_p+8);  // tmp2 = y
    
    CONJMADD(sum1,tmp1,tmp2);

    tmp3 = _mm_load_pd(x_p+10);  // tmp1 = x
    tmp4 = _mm_load_pd(y_p+10);  // tmp2 = y    
    
    CONJMADD(sum2,tmp3,tmp4);

    tmp5 = _mm_load_pd(x_p+12);  // tmp1 = x
    tmp6 = _mm_load_pd(y_p+12);  // tmp2 = y
    
    CONJMADD(sum3,tmp5,tmp6);

    tmp7 = _mm_load_pd(x_p+14);  // tmp1 = x
    tmp8 = _mm_load_pd(y_p+14);  // tmp2 = y
    
    CONJMADD(sum4,tmp7,tmp8);

    tmp1 = _mm_load_pd(x_p+16);  // tmp1 = x
    tmp2 = _mm_load_pd(y_p+16);  // tmp2 = y
    
    CONJMADD(sum1,tmp1,tmp2);

    tmp3 = _mm_load_pd(x_p+18);  // tmp1 = x
    tmp4 = _mm_load_pd(y_p+18);  // tmp2 = y
    
    CONJMADD(sum2,tmp3,tmp4);

    tmp5 = _mm_load_pd(x_p+20);  // tmp1 = x
    tmp6 = _mm_load_pd(y_p+20);  // tmp2 = y
    
    CONJMADD(sum3,tmp5,tmp6);

    tmp7 = _mm_load_pd(x_p+22);  // tmp1 = x
    tmp8 = _mm_load_pd(y_p+22);  // tmp2 = y
    
    CONJMADD(sum4,tmp7,tmp8);

    x_p+=24; y_p+=24;
  }


  // Collect the sums
  sum1 = _mm_add_pd(sum1,sum2);
  sum3 = _mm_add_pd(sum3,sum4);
  sum1 = _mm_add_pd(sum1,sum3);

  // Single store -- has to be unaligned in case
  // return value is not aligned. THe vectors should be aligned tho
  _mm_storeu_pd(sum,sum1);
  
}



} // namespace QDP;
