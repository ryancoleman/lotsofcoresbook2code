// $Id: sse_linalg_m_eq_mm_double.cc,v 1.6 2009-02-10 17:06:37 bjoo Exp $

/*! @file
 *  @brief Generic Scalar VAXPY routine
 *
 */

#include "scalarsite_sse/sse_linalg_mm_su3_double.h"

namespace QDP {

#include <xmmintrin.h>

#include "qdp_config.h"
#ifndef QDP_USE_SSE3

  // c = x*y;
#define CMUL(z,x,y)		\
  { \
    __m128d t1,t2,t3; \
    t1 = _mm_mul_pd((x),(y)); \
    t2 = _mm_shuffle_pd(t1,t1,0x1); \
    t3 = _mm_shuffle_pd((y),(y),0x1);\
    (z) = _mm_sub_pd(t1,t2); \
    t2 = _mm_mul_pd((x),t3); \
    t3 = _mm_shuffle_pd(t2,t2,0x1); \
    t3 = _mm_add_pd(t2,t3); \
    (z)= _mm_shuffle_pd((z),t3,0x2); \
  }

  // c+= x*y
#define CMADD(z,x,y)				\
  { \
    __m128d t1,t2,t3,t4; \
    t1 = _mm_mul_pd((x),(y)); \
    t2 = _mm_shuffle_pd(t1,t1,0x1); \
    t3 = _mm_shuffle_pd((y),(y),0x1);\
    t4 = _mm_sub_pd(t1,t2); \
    t2 = _mm_mul_pd((x),t3); \
    t3 = _mm_shuffle_pd(t2,t2,0x1); \
    t3 = _mm_add_pd(t2,t3); \
    t4= _mm_shuffle_pd(t4,t3,0x2); \
    (z) = _mm_add_pd((z),t4); \
  }



#else 
#warning Using SSE3
#include <pmmintrin.h>
  // Use SSE3

#define CMUL(z,x,y)		\
  { \
    __m128d t1; \
    t1 = _mm_mul_pd((x),(y)); \
    (z) = _mm_hsub_pd(t1,t1);			\
    t1 = _mm_shuffle_pd((y),(y),0x1);\
    t1 = _mm_mul_pd((x),t1); \
    t1 = _mm_hadd_pd(t1,t1); \
    (z)= _mm_shuffle_pd((z),t1,0x2);		\
  }

#define CMADD(z,x,y)				\
  { \
    __m128d t1,t2; \
    t1 = _mm_mul_pd((x),(y)); \
    t1 = _mm_hsub_pd(t1,t1); \
    t2 = _mm_shuffle_pd((y),(y),0x1);\
    t2 = _mm_mul_pd((x),t2); \
    t2 = _mm_hadd_pd(t2,t2); \
    t1= _mm_shuffle_pd(t1,t2,0x2);		\
    (z) = _mm_add_pd((z),t1);			\
  }
#endif

  /* M3 = M1*M2 */
  /* ALIGNED LOADS/STORES. Should be OK on evaluates over OLattices
     since these should always be aligned. */

  void ssed_m_eq_mm(REAL64* m3, REAL64* m2, REAL64* m1, int n_mat)
  {
    __m128d m1_1;
    __m128d m1_2;
    __m128d m1_3;

    __m128d m2_1;
    __m128d m2_2;
    
    __m128d m3_11;
    __m128d m3_12;
    __m128d m3_13;

    __m128d m3_21;
    __m128d m3_22;
    __m128d m3_23;

   

    REAL64* m1_p=m1;
    REAL64* m2_p=m2;
    REAL64* m3_p=m3;

    for(int i=0; i < n_mat; i++) { 
      m2_1 = _mm_load_pd(m2_p);
      m1_1 = _mm_load_pd(m1_p);      
      m1_2 = _mm_load_pd(m1_p+2);
      m1_3 = _mm_load_pd(m1_p+4);

      CMUL(m3_11, m2_1, m1_1);
      CMUL(m3_12, m2_1, m1_2);
      CMUL(m3_13, m2_1, m1_3);

      m2_1 = _mm_load_pd(m2_p+2);
      m1_1 = _mm_load_pd(m1_p+6);
      m1_2 = _mm_load_pd(m1_p+8);
      m1_3 = _mm_load_pd(m1_p+10);

      CMADD(m3_11, m2_1, m1_1);
      CMADD(m3_12, m2_1, m1_2);
      CMADD(m3_13, m2_1, m1_3);

      m2_1 = _mm_load_pd(m2_p+4);
      m1_1 = _mm_load_pd(m1_p+12);
      m1_2 = _mm_load_pd(m1_p+14);
      m1_3 = _mm_load_pd(m1_p+16);

      CMADD(m3_11, m2_1, m1_1);
      _mm_store_pd(m3_p, m3_11);

      CMADD(m3_12, m2_1, m1_2);
      _mm_store_pd(m3_p+2, m3_12);

      CMADD(m3_13, m2_1, m1_3);
      _mm_store_pd(m3_p+4, m3_13);


      m2_2 = _mm_load_pd(m2_p+6);
      m2_1 = _mm_load_pd(m2_p+12);

      m1_1 = _mm_load_pd(m1_p);      
      m1_2 = _mm_load_pd(m1_p+2);
      m1_3 = _mm_load_pd(m1_p+4);

      CMUL(m3_21, m2_2, m1_1);
      CMUL(m3_22, m2_2, m1_2);
      CMUL(m3_23, m2_2, m1_3);

      CMUL(m3_11, m2_1, m1_1);
      CMUL(m3_12, m2_1, m1_2);
      CMUL(m3_13, m2_1, m1_3);

      m2_2 = _mm_load_pd(m2_p+8);
      m2_1 = _mm_load_pd(m2_p+14);

      m1_1 = _mm_load_pd(m1_p+6);
      m1_2 = _mm_load_pd(m1_p+8);
      m1_3 = _mm_load_pd(m1_p+10);

      CMADD(m3_21, m2_2, m1_1);
      CMADD(m3_22, m2_2, m1_2);
      CMADD(m3_23, m2_2, m1_3);

      CMADD(m3_11, m2_1, m1_1);
      CMADD(m3_12, m2_1, m1_2);
      CMADD(m3_13, m2_1, m1_3);

      m2_2 = _mm_load_pd(m2_p+10);
      m2_1 = _mm_load_pd(m2_p+16);

      m1_1 = _mm_load_pd(m1_p+12);
      m1_2 = _mm_load_pd(m1_p+14);
      m1_3 = _mm_load_pd(m1_p+16);

      CMADD(m3_21, m2_2, m1_1);
      _mm_store_pd(m3_p+6, m3_21);

      CMADD(m3_22, m2_2, m1_2);
      _mm_store_pd(m3_p+8, m3_22);

      CMADD(m3_23, m2_2, m1_3);
      _mm_store_pd(m3_p+10, m3_23);

      CMADD(m3_11, m2_1, m1_1);
      _mm_store_pd(m3_p+12, m3_11);

      CMADD(m3_12, m2_1, m1_2);
      _mm_store_pd(m3_p+14, m3_12);

      CMADD(m3_13, m2_1, m1_3);
      _mm_store_pd(m3_p+16, m3_13);

      /* Next matrix */
      m1_p += 18; m2_p += 18; m3_p += 18;


    }

  }

  /* M3 += a M1*M2 */
  void ssed_m_peq_amm(REAL64* m3, REAL64* a, REAL64* m2, REAL64* m1, int n_mat)
  {
    __m128d m1_1;
    __m128d m1_2;
    __m128d m1_3;

    __m128d m2_1;
    __m128d m2_2;
    
    __m128d m3_11;
    __m128d m3_12;
    __m128d m3_13;

    __m128d m3_21;
    __m128d m3_22;
    __m128d m3_23;

    __m128d tmp1;
    __m128d scalar;
    
    scalar = _mm_load_sd(a);
  
    // cross components into tmp 
    // Zero tmp
    tmp1 = _mm_setzero_pd();
    tmp1 = _mm_shuffle_pd(scalar, scalar, 0x1);
    scalar = _mm_add_pd(scalar, tmp1);


    REAL64* m1_p=m1;
    REAL64* m2_p=m2;
    REAL64* m3_p=m3;


    // First row of M2 into all columns of M1
    for(int i=0; i < n_mat; i++) { 

      m3_11 = _mm_load_pd(m3_p);
      m3_12 = _mm_load_pd(m3_p+2);
      m3_13 = _mm_load_pd(m3_p+4);

      m2_1 = _mm_load_pd(m2_p);
      m2_1 = _mm_mul_pd(m2_1,scalar);

      m1_1 = _mm_load_pd(m1_p); 
      m1_2 = _mm_load_pd(m1_p+2);
      m1_3 = _mm_load_pd(m1_p+4);

      CMADD(m3_11, m2_1, m1_1);
      CMADD(m3_12, m2_1, m1_2);
      CMADD(m3_13, m2_1, m1_3);

      m2_1 = _mm_load_pd(m2_p+2);
      m2_1 = _mm_mul_pd(m2_1,scalar);

      m1_1 = _mm_load_pd(m1_p+6);
      m1_2 = _mm_load_pd(m1_p+8);
      m1_3 = _mm_load_pd(m1_p+10);

      CMADD(m3_11, m2_1, m1_1);
      CMADD(m3_12, m2_1, m1_2);
      CMADD(m3_13, m2_1, m1_3);

      m2_1 = _mm_load_pd(m2_p+4);
      m2_1 = _mm_mul_pd(m2_1,scalar);

      m1_1 = _mm_load_pd(m1_p+12);
      m1_2 = _mm_load_pd(m1_p+14);
      m1_3 = _mm_load_pd(m1_p+16);

      CMADD(m3_11, m2_1, m1_1);
      _mm_store_pd(m3_p, m3_11);

      CMADD(m3_12, m2_1, m1_2);
      _mm_store_pd(m3_p+2, m3_12);

      CMADD(m3_13, m2_1, m1_3);
      _mm_store_pd(m3_p+4, m3_13);


      m3_21 = _mm_load_pd(m3_p+6);
      m3_22 = _mm_load_pd(m3_p+8);
      m3_23 = _mm_load_pd(m3_p+10);
      m3_11 = _mm_load_pd(m3_p+12);
      m3_12 = _mm_load_pd(m3_p+14);
      m3_13 = _mm_load_pd(m3_p+16);

      m2_2 = _mm_load_pd(m2_p+6);
      m2_2 = _mm_mul_pd(m2_2,scalar);

      m2_1 = _mm_load_pd(m2_p+12);
      m2_1 = _mm_mul_pd(m2_1,scalar);

      m1_1 = _mm_load_pd(m1_p);      
      m1_2 = _mm_load_pd(m1_p+2);
      m1_3 = _mm_load_pd(m1_p+4);



      CMADD(m3_21, m2_2, m1_1)
      CMADD(m3_22, m2_2, m1_2);
      CMADD(m3_23, m2_2, m1_3);


      CMADD(m3_11, m2_1, m1_1);
      CMADD(m3_12, m2_1, m1_2);
      CMADD(m3_13, m2_1, m1_3);

      m2_2 = _mm_load_pd(m2_p+8);
      m2_2 = _mm_mul_pd(m2_2,scalar);

      m2_1 = _mm_load_pd(m2_p+14);
      m2_1 = _mm_mul_pd(m2_1,scalar);

      m1_1 = _mm_load_pd(m1_p+6);
      m1_2 = _mm_load_pd(m1_p+8);
      m1_3 = _mm_load_pd(m1_p+10);

      CMADD(m3_21, m2_2, m1_1);
      CMADD(m3_22, m2_2, m1_2);
      CMADD(m3_23, m2_2, m1_3);

      CMADD(m3_11, m2_1, m1_1);
      CMADD(m3_12, m2_1, m1_2);
      CMADD(m3_13, m2_1, m1_3);

      m2_2 = _mm_load_pd(m2_p+10);
      m2_2 = _mm_mul_pd(m2_2,scalar);

      m2_1 = _mm_load_pd(m2_p+16);
      m2_1 = _mm_mul_pd(m2_1,scalar);

      m1_1 = _mm_load_pd(m1_p+12);
      m1_2 = _mm_load_pd(m1_p+14);
      m1_3 = _mm_load_pd(m1_p+16);

      CMADD(m3_21, m2_2, m1_1);
      _mm_store_pd(m3_p+6, m3_21);

      CMADD(m3_22, m2_2, m1_2);
      _mm_store_pd(m3_p+8, m3_22);

      CMADD(m3_23, m2_2, m1_3);
      _mm_store_pd(m3_p+10, m3_23);

      CMADD(m3_11, m2_1, m1_1);
      _mm_store_pd(m3_p+12, m3_11);

      CMADD(m3_12, m2_1, m1_2);
      _mm_store_pd(m3_p+14, m3_12);

      CMADD(m3_13, m2_1, m1_3);
      _mm_store_pd(m3_p+16, m3_13);

      /* Next matrix */
      m1_p += 18; m2_p += 18; m3_p += 18;


    }


  }



  /* UNALIGNED LOADS/STORES. Should be OK when data is not in an olattice
     and the data may thus not be 16 byte aligned */

  /* M3 = M1*M2 */
  void ssed_m_eq_mm_u(REAL64* m3, REAL64* m2, REAL64* m1, int n_mat)
  {
    __m128d m1_1;
    __m128d m1_2;
    __m128d m1_3;

    __m128d m2_1;
    __m128d m2_2;
    
    __m128d m3_11;
    __m128d m3_12;
    __m128d m3_13;

    __m128d m3_21;
    __m128d m3_22;
    __m128d m3_23;

   

    REAL64* m1_p=m1;
    REAL64* m2_p=m2;
    REAL64* m3_p=m3;

    for(int i=0; i < n_mat; i++) { 
      m2_1 = _mm_loadu_pd(m2_p);
      m1_1 = _mm_loadu_pd(m1_p);      
      m1_2 = _mm_loadu_pd(m1_p+2);
      m1_3 = _mm_loadu_pd(m1_p+4);

      CMUL(m3_11, m2_1, m1_1);
      CMUL(m3_12, m2_1, m1_2);
      CMUL(m3_13, m2_1, m1_3);

      m2_1 = _mm_loadu_pd(m2_p+2);
      m1_1 = _mm_loadu_pd(m1_p+6);
      m1_2 = _mm_loadu_pd(m1_p+8);
      m1_3 = _mm_loadu_pd(m1_p+10);

      CMADD(m3_11, m2_1, m1_1);
      CMADD(m3_12, m2_1, m1_2);
      CMADD(m3_13, m2_1, m1_3);

      m2_1 = _mm_loadu_pd(m2_p+4);
      m1_1 = _mm_loadu_pd(m1_p+12);
      m1_2 = _mm_loadu_pd(m1_p+14);
      m1_3 = _mm_loadu_pd(m1_p+16);

      CMADD(m3_11, m2_1, m1_1);
      _mm_storeu_pd(m3_p, m3_11);

      CMADD(m3_12, m2_1, m1_2);
      _mm_storeu_pd(m3_p+2, m3_12);

      CMADD(m3_13, m2_1, m1_3);
      _mm_storeu_pd(m3_p+4, m3_13);


      m2_2 = _mm_loadu_pd(m2_p+6);
      m2_1 = _mm_loadu_pd(m2_p+12);

      m1_1 = _mm_loadu_pd(m1_p);      
      m1_2 = _mm_loadu_pd(m1_p+2);
      m1_3 = _mm_loadu_pd(m1_p+4);

      CMUL(m3_21, m2_2, m1_1);
      CMUL(m3_22, m2_2, m1_2);
      CMUL(m3_23, m2_2, m1_3);

      CMUL(m3_11, m2_1, m1_1);
      CMUL(m3_12, m2_1, m1_2);
      CMUL(m3_13, m2_1, m1_3);

      m2_2 = _mm_loadu_pd(m2_p+8);
      m2_1 = _mm_loadu_pd(m2_p+14);

      m1_1 = _mm_loadu_pd(m1_p+6);
      m1_2 = _mm_loadu_pd(m1_p+8);
      m1_3 = _mm_loadu_pd(m1_p+10);

      CMADD(m3_21, m2_2, m1_1);
      CMADD(m3_22, m2_2, m1_2);
      CMADD(m3_23, m2_2, m1_3);

      CMADD(m3_11, m2_1, m1_1);
      CMADD(m3_12, m2_1, m1_2);
      CMADD(m3_13, m2_1, m1_3);

      m2_2 = _mm_loadu_pd(m2_p+10);
      m2_1 = _mm_loadu_pd(m2_p+16);

      m1_1 = _mm_loadu_pd(m1_p+12);
      m1_2 = _mm_loadu_pd(m1_p+14);
      m1_3 = _mm_loadu_pd(m1_p+16);

      CMADD(m3_21, m2_2, m1_1);
      _mm_storeu_pd(m3_p+6, m3_21);

      CMADD(m3_22, m2_2, m1_2);
      _mm_storeu_pd(m3_p+8, m3_22);

      CMADD(m3_23, m2_2, m1_3);
      _mm_storeu_pd(m3_p+10, m3_23);

      CMADD(m3_11, m2_1, m1_1);
      _mm_storeu_pd(m3_p+12, m3_11);

      CMADD(m3_12, m2_1, m1_2);
      _mm_storeu_pd(m3_p+14, m3_12);

      CMADD(m3_13, m2_1, m1_3);
      _mm_storeu_pd(m3_p+16, m3_13);

      /* Next matrix */
      m1_p += 18; m2_p += 18; m3_p += 18;


    }

  }

  /* M3 += a M1*M2 */
  void ssed_m_peq_amm_u(REAL64* m3, REAL64* a, REAL64* m2, REAL64* m1, int n_mat)
  {
    __m128d m1_1;
    __m128d m1_2;
    __m128d m1_3;

    __m128d m2_1;
    __m128d m2_2;
    
    __m128d m3_11;
    __m128d m3_12;
    __m128d m3_13;

    __m128d m3_21;
    __m128d m3_22;
    __m128d m3_23;

    __m128d tmp1;
    __m128d scalar;
    
    scalar = _mm_load_sd(a);
  
    // cross components into tmp 
    // Zero tmp
    tmp1 = _mm_setzero_pd();
    tmp1 = _mm_shuffle_pd(scalar, scalar, 0x1);
    scalar = _mm_add_pd(scalar, tmp1);


    REAL64* m1_p=m1;
    REAL64* m2_p=m2;
    REAL64* m3_p=m3;


    // First row of M2 into all columns of M1
    for(int i=0; i < n_mat; i++) { 

      m3_11 = _mm_loadu_pd(m3_p);
      m3_12 = _mm_loadu_pd(m3_p+2);
      m3_13 = _mm_loadu_pd(m3_p+4);

      m2_1 = _mm_loadu_pd(m2_p);
      m2_1 = _mm_mul_pd(m2_1,scalar);

      m1_1 = _mm_loadu_pd(m1_p); 
      m1_2 = _mm_loadu_pd(m1_p+2);
      m1_3 = _mm_loadu_pd(m1_p+4);

      CMADD(m3_11, m2_1, m1_1);
      CMADD(m3_12, m2_1, m1_2);
      CMADD(m3_13, m2_1, m1_3);

      m2_1 = _mm_loadu_pd(m2_p+2);
      m2_1 = _mm_mul_pd(m2_1,scalar);

      m1_1 = _mm_loadu_pd(m1_p+6);
      m1_2 = _mm_loadu_pd(m1_p+8);
      m1_3 = _mm_loadu_pd(m1_p+10);

      CMADD(m3_11, m2_1, m1_1);
      CMADD(m3_12, m2_1, m1_2);
      CMADD(m3_13, m2_1, m1_3);

      m2_1 = _mm_loadu_pd(m2_p+4);
      m2_1 = _mm_mul_pd(m2_1,scalar);

      m1_1 = _mm_loadu_pd(m1_p+12);
      m1_2 = _mm_loadu_pd(m1_p+14);
      m1_3 = _mm_loadu_pd(m1_p+16);

      CMADD(m3_11, m2_1, m1_1);
      _mm_storeu_pd(m3_p, m3_11);

      CMADD(m3_12, m2_1, m1_2);
      _mm_storeu_pd(m3_p+2, m3_12);

      CMADD(m3_13, m2_1, m1_3);
      _mm_storeu_pd(m3_p+4, m3_13);


      m3_21 = _mm_loadu_pd(m3_p+6);
      m3_22 = _mm_loadu_pd(m3_p+8);
      m3_23 = _mm_loadu_pd(m3_p+10);
      m3_11 = _mm_loadu_pd(m3_p+12);
      m3_12 = _mm_loadu_pd(m3_p+14);
      m3_13 = _mm_loadu_pd(m3_p+16);

      m2_2 = _mm_loadu_pd(m2_p+6);
      m2_2 = _mm_mul_pd(m2_2,scalar);

      m2_1 = _mm_loadu_pd(m2_p+12);
      m2_1 = _mm_mul_pd(m2_1,scalar);

      m1_1 = _mm_loadu_pd(m1_p);      
      m1_2 = _mm_loadu_pd(m1_p+2);
      m1_3 = _mm_loadu_pd(m1_p+4);



      CMADD(m3_21, m2_2, m1_1)
      CMADD(m3_22, m2_2, m1_2);
      CMADD(m3_23, m2_2, m1_3);


      CMADD(m3_11, m2_1, m1_1);
      CMADD(m3_12, m2_1, m1_2);
      CMADD(m3_13, m2_1, m1_3);

      m2_2 = _mm_loadu_pd(m2_p+8);
      m2_2 = _mm_mul_pd(m2_2,scalar);

      m2_1 = _mm_loadu_pd(m2_p+14);
      m2_1 = _mm_mul_pd(m2_1,scalar);

      m1_1 = _mm_loadu_pd(m1_p+6);
      m1_2 = _mm_loadu_pd(m1_p+8);
      m1_3 = _mm_loadu_pd(m1_p+10);

      CMADD(m3_21, m2_2, m1_1);
      CMADD(m3_22, m2_2, m1_2);
      CMADD(m3_23, m2_2, m1_3);

      CMADD(m3_11, m2_1, m1_1);
      CMADD(m3_12, m2_1, m1_2);
      CMADD(m3_13, m2_1, m1_3);

      m2_2 = _mm_loadu_pd(m2_p+10);
      m2_2 = _mm_mul_pd(m2_2,scalar);

      m2_1 = _mm_loadu_pd(m2_p+16);
      m2_1 = _mm_mul_pd(m2_1,scalar);

      m1_1 = _mm_loadu_pd(m1_p+12);
      m1_2 = _mm_loadu_pd(m1_p+14);
      m1_3 = _mm_loadu_pd(m1_p+16);

      CMADD(m3_21, m2_2, m1_1);
      _mm_storeu_pd(m3_p+6, m3_21);

      CMADD(m3_22, m2_2, m1_2);
      _mm_storeu_pd(m3_p+8, m3_22);

      CMADD(m3_23, m2_2, m1_3);
      _mm_storeu_pd(m3_p+10, m3_23);

      CMADD(m3_11, m2_1, m1_1);
      _mm_storeu_pd(m3_p+12, m3_11);

      CMADD(m3_12, m2_1, m1_2);
      _mm_storeu_pd(m3_p+14, m3_12);

      CMADD(m3_13, m2_1, m1_3);
      _mm_storeu_pd(m3_p+16, m3_13);

      /* Next matrix */
      m1_p += 18; m2_p += 18; m3_p += 18;


    }


  }



} // namespace QDP;

