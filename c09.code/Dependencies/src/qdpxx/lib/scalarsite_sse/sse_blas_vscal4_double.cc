// $Id: sse_blas_vscal4_double.cc,v 1.1 2008-06-18 16:02:11 bjoo Exp $

/*! @file
 *  @brief Generic Scalar VAXPY routine
 *
 */

#include "scalarsite_sse/sse_blas_vscal4_double.h"

namespace QDP {

#include <xmmintrin.h>


void vscal4(REAL64 *z,REAL64 *a,REAL64 *x, int n_4spin)
{
  __m128d scalar;
  __m128d tmp1;

  __m128d z1;
  __m128d x1;

  __m128d z2;
  __m128d x2;

  __m128d z3;
  __m128d x3;

  __m128d z4;
  __m128d x4;

  __m128d z5;
  __m128d x5;

  __m128d z6;
  __m128d x6;



  // Load the scalar into low bytes of scalar
  scalar = _mm_load_sd(a);
  
  // cross components into tmp 
  // Zero tmp
  tmp1 = _mm_setzero_pd();
  tmp1 = _mm_shuffle_pd(scalar, scalar, 0x1);
  scalar = _mm_add_pd(scalar, tmp1);

  //QDPIO::cout << "In Balints Routine" << endl;

  double *x_p=x;
  double *z_p=z;


  for(int i=0; i < n_4spin; i++) { 
    x1 = _mm_load_pd(x_p);
    z1 = _mm_mul_pd(scalar,x1);
    
    x2 = _mm_load_pd(x_p+2);
    z2 = _mm_mul_pd(scalar,x2);
    
    x3 = _mm_load_pd(x_p+4);
    z3 = _mm_mul_pd(scalar,x3);
    
    x4 = _mm_load_pd(x_p+6);
    z4 = _mm_mul_pd(scalar,x4);
    
    x5 = _mm_load_pd(x_p+8);
    z5 = _mm_mul_pd(scalar,x5);
    
    x6 = _mm_load_pd(x_p+10);
    z6 = _mm_mul_pd(scalar,x6);
    
    _mm_store_pd(z_p, z1);
    _mm_store_pd(z_p+2, z2);
    _mm_store_pd(z_p+4, z3);
    _mm_store_pd(z_p+6, z4);
    _mm_store_pd(z_p+8, z5);
    _mm_store_pd(z_p+10, z6);
    
    x1 = _mm_load_pd(x_p+12);
    z1 = _mm_mul_pd(scalar,x1);
    
    x2 = _mm_load_pd(x_p+14);
    z2 = _mm_mul_pd(scalar,x2);
    
    x3 = _mm_load_pd(x_p+16);
    z3 = _mm_mul_pd(scalar,x3);
    
    x4 = _mm_load_pd(x_p+18);
    z4 = _mm_mul_pd(scalar,x4);
    
    x5 = _mm_load_pd(x_p+20);
    z5 = _mm_mul_pd(scalar,x5);
    
    x6 = _mm_load_pd(x_p+22);
    z6 = _mm_mul_pd(scalar,x6);
    
    _mm_store_pd(z_p+12, z1);
    _mm_store_pd(z_p+14, z2);
    _mm_store_pd(z_p+16, z3);
    _mm_store_pd(z_p+18, z4);
    _mm_store_pd(z_p+20, z5);
    _mm_store_pd(z_p+22, z6);

    x_p+=24; z_p+=24;
    
  }
  
}



} // namespace QDP;

