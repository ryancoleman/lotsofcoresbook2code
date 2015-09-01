// $Id: sse_linalg_m_eq_scal_m_double.cc,v 1.3 2008-09-23 15:23:46 bjoo Exp $

/*! @file
 *  @brief Generic Scalar VAXPY routine
 *
 */

#include "scalarsite_sse/sse_linalg_mm_su3_double.h"

namespace QDP {

#include <xmmintrin.h>

  /* M = a*M  a is scalar */
  void ssed_m_eq_scal_m(REAL64* m2, REAL64* a, REAL64 *m1, int n_mat)
  {
    __m128d scalar;
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
    __m128d tmp13;
    __m128d tmp14;
    __m128d tmp15;


    // Load the scalar into low bytes of scalar
    scalar = _mm_load_sd(a);
    
    // cross components into tmp 
    // Zero tmp
    tmp1 = _mm_setzero_pd();
    tmp1 = _mm_shuffle_pd(scalar, scalar, 0x1);
    scalar = _mm_add_pd(scalar, tmp1);

    REAL64* m1_p=m1;
    REAL64* m2_p=m2;

    for(int i=0; i < n_mat; i++) { 
      tmp1= _mm_loadu_pd(m1_p);
      tmp2= _mm_mul_pd(scalar, tmp1);
      _mm_storeu_pd(m2_p, tmp2);
      
      tmp3= _mm_loadu_pd(m1_p+2);
      tmp4= _mm_mul_pd(scalar, tmp3);
      _mm_storeu_pd(m2_p+2, tmp4);

      tmp5= _mm_loadu_pd(m1_p+4);
      tmp6= _mm_mul_pd(scalar, tmp5);
      _mm_storeu_pd(m2_p+4, tmp6);

      tmp7= _mm_loadu_pd(m1_p+6);
      tmp8= _mm_mul_pd(scalar, tmp7);
      _mm_storeu_pd(m2_p+6, tmp8);

      tmp9= _mm_loadu_pd(m1_p+8);
      tmp10= _mm_mul_pd(scalar, tmp9);
      _mm_storeu_pd(m2_p+8, tmp10);

      tmp11= _mm_loadu_pd(m1_p+10);
      tmp12= _mm_mul_pd(scalar, tmp11);
      _mm_storeu_pd(m2_p+10, tmp12);

      tmp13= _mm_loadu_pd(m1_p+12);
      tmp14= _mm_mul_pd(scalar, tmp13);
      _mm_storeu_pd(m2_p+12, tmp14);

      tmp15= _mm_loadu_pd(m1_p+14);
      tmp4= _mm_mul_pd(scalar, tmp15);
      _mm_storeu_pd(m2_p+14, tmp4);

      tmp5= _mm_loadu_pd(m1_p+16);
      tmp6= _mm_mul_pd(scalar, tmp5);
      _mm_storeu_pd(m2_p+16, tmp6);

      m1_p += 18; m2_p+=18;
    }

  }

  /* M *= a,  a is a scalar */
  void ssed_m_muleq_scal(REAL64* m, REAL64* a, int n_mat)
  {
    __m128d scalar;
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
    __m128d tmp13;
    __m128d tmp14;
    __m128d tmp15;


    // Load the scalar into low bytes of scalar
    scalar = _mm_load_sd(a);
    
    // cross components into tmp 
    // Zero tmp
    tmp1 = _mm_setzero_pd();
    tmp1 = _mm_shuffle_pd(scalar, scalar, 0x1);
    scalar = _mm_add_pd(scalar, tmp1);

    REAL64* m1_p=m;

    for(int i=0; i < n_mat; i++) { 
      tmp1= _mm_loadu_pd(m1_p);
      tmp2= _mm_mul_pd(scalar, tmp1);
      _mm_storeu_pd(m1_p, tmp2);
      
      tmp3= _mm_loadu_pd(m1_p+2);
      tmp4= _mm_mul_pd(scalar, tmp3);
      _mm_storeu_pd(m1_p+2, tmp4);

      tmp5= _mm_loadu_pd(m1_p+4);
      tmp6= _mm_mul_pd(scalar, tmp5);
      _mm_storeu_pd(m1_p+4, tmp6);

      tmp7= _mm_loadu_pd(m1_p+6);
      tmp8= _mm_mul_pd(scalar, tmp7);
      _mm_storeu_pd(m1_p+6, tmp8);

      tmp9= _mm_loadu_pd(m1_p+8);
      tmp10= _mm_mul_pd(scalar, tmp9);
      _mm_storeu_pd(m1_p+8, tmp10);

      tmp11= _mm_loadu_pd(m1_p+10);
      tmp12= _mm_mul_pd(scalar, tmp11);
      _mm_storeu_pd(m1_p+10, tmp12);

      tmp13= _mm_loadu_pd(m1_p+12);
      tmp14= _mm_mul_pd(scalar, tmp13);
      _mm_storeu_pd(m1_p+12, tmp14);

      tmp15= _mm_loadu_pd(m1_p+14);
      tmp4= _mm_mul_pd(scalar, tmp15);
      _mm_storeu_pd(m1_p+14, tmp4);

      tmp5= _mm_loadu_pd(m1_p+16);
      tmp6= _mm_mul_pd(scalar, tmp5);
      _mm_storeu_pd(m1_p+16, tmp6);

      m1_p += 18;
    }

  }


} // namespace QDP;

