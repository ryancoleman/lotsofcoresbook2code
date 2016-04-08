// $Id: sse_linalg_m_peq_m_double.cc,v 1.2 2008-09-23 15:23:46 bjoo Exp $

/*! @file
 *  @brief Generic Scalar VAXPY routine
 *
 */

#include "scalarsite_sse/sse_linalg_mm_su3_double.h"

namespace QDP {

#include <xmmintrin.h>

  typedef union { 
    double c[2];
    __m128d v;
  } VD;


  /* M2 += M1 */
  void ssed_m_peq_m(REAL64* m2, REAL64* m1, int n_mat)
  {

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


    REAL64* m1_p=m1;
    REAL64* m2_p=m2;

    for(int i=0; i < n_mat; i++) { 
      tmp1= _mm_loadu_pd(m1_p);
      tmp2= _mm_loadu_pd(m2_p);
      tmp3 = _mm_add_pd(tmp1,tmp2);
      _mm_storeu_pd(m2_p, tmp3);

      tmp4= _mm_loadu_pd(m1_p+2);
      tmp5= _mm_loadu_pd(m2_p+2);
      tmp6 = _mm_add_pd(tmp4,tmp5);
      _mm_storeu_pd(m2_p+2, tmp6);

      tmp7= _mm_loadu_pd(m1_p+4);
      tmp8= _mm_loadu_pd(m2_p+4);
      tmp9 = _mm_add_pd(tmp7,tmp8);
      _mm_storeu_pd(m2_p+4, tmp9);

      tmp10= _mm_loadu_pd(m1_p+6);
      tmp11= _mm_loadu_pd(m2_p+6);
      tmp12 = _mm_add_pd(tmp10,tmp11);
      _mm_storeu_pd(m2_p+6, tmp12);

      tmp13= _mm_loadu_pd(m1_p+8);
      tmp14= _mm_loadu_pd(m2_p+8);
      tmp15 = _mm_add_pd(tmp13,tmp14);
      _mm_storeu_pd(m2_p+8, tmp15);

      tmp1= _mm_loadu_pd(m1_p+10);
      tmp2= _mm_loadu_pd(m2_p+10);
      tmp3 = _mm_add_pd(tmp1,tmp2);
      _mm_storeu_pd(m2_p+10, tmp3);

      tmp4= _mm_loadu_pd(m1_p+12);
      tmp5= _mm_loadu_pd(m2_p+12);
      tmp6 = _mm_add_pd(tmp4,tmp5);
      _mm_storeu_pd(m2_p+12, tmp6);

      tmp7= _mm_loadu_pd(m1_p+14);
      tmp8= _mm_loadu_pd(m2_p+14);
      tmp9 = _mm_add_pd(tmp7,tmp8);
      _mm_storeu_pd(m2_p+14, tmp9);

      tmp10= _mm_loadu_pd(m1_p+16);
      tmp11= _mm_loadu_pd(m2_p+16);
      tmp12 = _mm_add_pd(tmp10,tmp11);
      _mm_storeu_pd(m2_p+16, tmp12);

      m1_p += 18; m2_p+=18;
    }

  }

  /* M2 -= M1 */
  void ssed_m_meq_m(REAL64* m2, REAL64* m1, int n_mat) 
  {
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


    REAL64* m1_p=m1;
    REAL64* m2_p=m2;

    for(int i=0; i < n_mat; i++) { 
      tmp1= _mm_loadu_pd(m1_p);
      tmp2= _mm_loadu_pd(m2_p);
      tmp3 = _mm_sub_pd(tmp2,tmp1);
      _mm_storeu_pd(m2_p, tmp3);

      tmp4= _mm_loadu_pd(m1_p+2);
      tmp5= _mm_loadu_pd(m2_p+2);
      tmp6 = _mm_sub_pd(tmp5,tmp4);
      _mm_storeu_pd(m2_p+2, tmp6);
      
      tmp7= _mm_loadu_pd(m1_p+4);
      tmp8= _mm_loadu_pd(m2_p+4);
      tmp9 = _mm_sub_pd(tmp8,tmp7);
      _mm_storeu_pd(m2_p+4, tmp9);
      
      tmp10= _mm_loadu_pd(m1_p+6);
      tmp11= _mm_loadu_pd(m2_p+6);
      tmp12 = _mm_sub_pd(tmp11,tmp10);
      _mm_storeu_pd(m2_p+6, tmp12);
      
      tmp13= _mm_loadu_pd(m1_p+8);
      tmp14= _mm_loadu_pd(m2_p+8);
      tmp15 = _mm_sub_pd(tmp14,tmp13);
      _mm_storeu_pd(m2_p+8, tmp15);
      
      tmp1= _mm_loadu_pd(m1_p+10);
      tmp2= _mm_loadu_pd(m2_p+10);
      tmp3 = _mm_sub_pd(tmp2,tmp1);
      _mm_storeu_pd(m2_p+10, tmp3);
      
      tmp4= _mm_loadu_pd(m1_p+12);
      tmp5= _mm_loadu_pd(m2_p+12);
      tmp6 = _mm_sub_pd(tmp5,tmp4);
      _mm_storeu_pd(m2_p+12, tmp6);
      
      tmp7= _mm_loadu_pd(m1_p+14);
      tmp8= _mm_loadu_pd(m2_p+14);
      tmp9 = _mm_sub_pd(tmp8,tmp7);
      _mm_storeu_pd(m2_p+14, tmp9);
      
      tmp10= _mm_loadu_pd(m1_p+16);
      tmp11= _mm_loadu_pd(m2_p+16);
      tmp12 = _mm_sub_pd(tmp11,tmp10);
      _mm_storeu_pd(m2_p+16, tmp12);

      m1_p += 18; m2_p+=18;
    }
  }


  /* M2 += adj(M1)*/
  void ssed_m_peq_h(REAL64* m2, REAL64* m1, int n_mat)
  {
    __m128d mfact = _mm_set_pd( (REAL64)(-1), (REAL64)(1) );

    __m128d m1_11;
    __m128d m1_12;
    __m128d m1_13;
    __m128d m1_21;
    __m128d m1_22;
    __m128d m1_23;
    __m128d m1_31;
    __m128d m1_32;
    __m128d m1_33;

    __m128d tmp1;
    __m128d tmp2;
    __m128d tmp3;
    __m128d tmp4;
    __m128d tmp5;
    __m128d tmp6;


    REAL64* m1_p=m1;
    REAL64* m2_p=m2;

    for(int i=0; i < n_mat; i++) { 
      // Stream in m1
      m1_11= _mm_loadu_pd(m1_p);
      m1_12= _mm_loadu_pd(m1_p+2);
      m1_13= _mm_loadu_pd(m1_p+4);
      m1_21= _mm_loadu_pd(m1_p+6);
      m1_22= _mm_loadu_pd(m1_p+8);
      m1_23= _mm_loadu_pd(m1_p+10);
      m1_31= _mm_loadu_pd(m1_p+12);
      m1_32= _mm_loadu_pd(m1_p+14);
      m1_33= _mm_loadu_pd(m1_p+16);
   
      tmp1 = _mm_loadu_pd(m2_p);
      tmp2 = _mm_mul_pd(mfact, m1_11);
      tmp3 = _mm_add_pd(tmp1, tmp2);
      _mm_storeu_pd(m2_p, tmp3);

      tmp4 = _mm_loadu_pd(m2_p+2);
      tmp5 = _mm_mul_pd(mfact, m1_21);
      tmp6 = _mm_add_pd(tmp4, tmp5);
      _mm_storeu_pd(m2_p+2,tmp6);

      tmp1 = _mm_loadu_pd(m2_p+4);
      tmp2 = _mm_mul_pd(mfact, m1_31);
      tmp3 = _mm_add_pd(tmp1, tmp2);
      _mm_storeu_pd(m2_p+4, tmp3);

      tmp4 = _mm_loadu_pd(m2_p+6);
      tmp5 = _mm_mul_pd(mfact, m1_12);
      tmp6 = _mm_add_pd(tmp4, tmp5);
      _mm_storeu_pd(m2_p+6, tmp6);


      tmp1 = _mm_loadu_pd(m2_p+8);
      tmp2 = _mm_mul_pd(mfact, m1_22);
      tmp3 = _mm_add_pd(tmp1, tmp2);
      _mm_storeu_pd(m2_p+8, tmp3);

      tmp4 = _mm_loadu_pd(m2_p+10);
      tmp5 = _mm_mul_pd(mfact, m1_32);
      tmp6 = _mm_add_pd(tmp4, tmp5);
      _mm_storeu_pd(m2_p+10,tmp6);

      tmp4 = _mm_loadu_pd(m2_p+12);
      tmp5 = _mm_mul_pd(mfact, m1_13);
      tmp6 = _mm_add_pd(tmp4, tmp5);
      _mm_storeu_pd(m2_p+12, tmp6);


      tmp1 = _mm_loadu_pd(m2_p+14);
      tmp2 = _mm_mul_pd(mfact, m1_23);
      tmp3 = _mm_add_pd(tmp1, tmp2);
      _mm_storeu_pd(m2_p+14, tmp3);

      tmp4 = _mm_loadu_pd(m2_p+16);
      tmp5 = _mm_mul_pd(mfact, m1_33);
      tmp6 = _mm_add_pd(tmp4, tmp5);
      _mm_storeu_pd(m2_p+16, tmp6);

      m1_p += 18; m2_p+=18;
    }

  }

  /* M2 -= adj(M1) */
  void ssed_m_meq_h(REAL64* m2, REAL64* m1, int n_mat) 
  {
    __m128d m1_11;
    __m128d m1_12;
    __m128d m1_13;
    __m128d m1_21;
    __m128d m1_22;
    __m128d m1_23;
    __m128d m1_31;
    __m128d m1_32;
    __m128d m1_33;

    __m128d tmp1;
    __m128d tmp2;
    __m128d tmp3;
    __m128d tmp4;
    __m128d tmp5;
    __m128d tmp6;


    __m128d  mfact = _mm_set_pd( (REAL64)(-1), (REAL64)(1)); 

    REAL64* m1_p=m1;
    REAL64* m2_p=m2;

    for(int i=0; i < n_mat; i++) {
      // Stream in m1
      m1_11= _mm_loadu_pd(m1_p);
      m1_12= _mm_loadu_pd(m1_p+2);
      m1_13= _mm_loadu_pd(m1_p+4);
      m1_21= _mm_loadu_pd(m1_p+6);
      m1_22= _mm_loadu_pd(m1_p+8);
      m1_23= _mm_loadu_pd(m1_p+10);
      m1_31= _mm_loadu_pd(m1_p+12);
      m1_32= _mm_loadu_pd(m1_p+14);
      m1_33= _mm_loadu_pd(m1_p+16);
   
      tmp1 = _mm_loadu_pd(m2_p);
      tmp2 = _mm_mul_pd(mfact, m1_11);
      tmp3 = _mm_sub_pd(tmp1, tmp2);
      _mm_storeu_pd(m2_p, tmp3);

      tmp4 = _mm_loadu_pd(m2_p+2);
      tmp5 = _mm_mul_pd(mfact, m1_21);
      tmp6 = _mm_sub_pd(tmp4, tmp5);
      _mm_storeu_pd(m2_p+2,tmp6);

      tmp1 = _mm_loadu_pd(m2_p+4);
      tmp2 = _mm_mul_pd(mfact, m1_31);
      tmp3 = _mm_sub_pd(tmp1, tmp2);
      _mm_storeu_pd(m2_p+4, tmp3);

      tmp4 = _mm_loadu_pd(m2_p+6);
      tmp5 = _mm_mul_pd(mfact, m1_12);
      tmp6 = _mm_sub_pd(tmp4, tmp5);
      _mm_storeu_pd(m2_p+6, tmp6);


      tmp1 = _mm_loadu_pd(m2_p+8);
      tmp2 = _mm_mul_pd(mfact, m1_22);
      tmp3 = _mm_sub_pd(tmp1, tmp2);
      _mm_storeu_pd(m2_p+8, tmp3);

      tmp4 = _mm_loadu_pd(m2_p+10);
      tmp5 = _mm_mul_pd(mfact, m1_32);
      tmp6 = _mm_sub_pd(tmp4, tmp5);
      _mm_storeu_pd(m2_p+10,tmp6);

      tmp4 = _mm_loadu_pd(m2_p+12);
      tmp5 = _mm_mul_pd(mfact, m1_13);
      tmp6 = _mm_sub_pd(tmp4, tmp5);
      _mm_storeu_pd(m2_p+12, tmp6);


      tmp1 = _mm_loadu_pd(m2_p+14);
      tmp2 = _mm_mul_pd(mfact, m1_23);
      tmp3 = _mm_sub_pd(tmp1, tmp2);
      _mm_storeu_pd(m2_p+14, tmp3);

      tmp4 = _mm_loadu_pd(m2_p+16);
      tmp5 = _mm_mul_pd(mfact, m1_33);
      tmp6 = _mm_sub_pd(tmp4, tmp5);
      _mm_storeu_pd(m2_p+16, tmp6);

      m1_p += 18; m2_p+=18;
    }
  }



} // namespace QDP;

