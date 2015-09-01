// $Id: sse_linalg_m_eq_hh_double.cc,v 1.5 2009-07-14 20:08:42 bjoo Exp $

/*! @file
 *  @brief Generic Scalar VAXPY routine
 *
 */

#include "scalarsite_sse/sse_linalg_mm_su3_double.h"

namespace QDP {

#include <xmmintrin.h>

typedef union {
  __m128d v;
  double  d[2];
} VD;

#include "qdp_config.h"

#ifndef QDP_USE_SSE3

  /* SSE 2 */
#define CCMUL(z,x,y)		\
  { \
    __m128d t1,t2,t3; \
    __m128d t4 = _mm_set_pd( (double)(-1),(double)1 );	\
    t1 = _mm_mul_pd(x,y); \
    t2 = _mm_shuffle_pd(t1,t1,0x1); \
    t3 = _mm_shuffle_pd(y,y,0x1);\
    z = _mm_sub_pd(t1,t2); \
    t2 = _mm_mul_pd(x,t3); \
    t3 = _mm_shuffle_pd(t2,t2,0x1); \
    t3 = _mm_add_pd(t2,t3); \
    z= _mm_shuffle_pd(z,t3,0x2); \
    z= _mm_mul_pd(z,t4); \
  }

#define CCMADD(z,x,y)				\
  { \
    __m128d t1,t2,t3,t4; \
    __m128d t5 = _mm_set_pd( (double)(-1),(double)1) ;	\
    t1 = _mm_mul_pd(x,y); \
    t2 = _mm_shuffle_pd(t1,t1,0x1); \
    t3 = _mm_shuffle_pd(y,y,0x1);\
    t4 = _mm_sub_pd(t1,t2); \
    t2 = _mm_mul_pd(x,t3); \
    t3 = _mm_shuffle_pd(t2,t2,0x1); \
    t3 = _mm_add_pd(t2,t3); \
    t4= _mm_shuffle_pd(t4,t3,0x2); \
    t4= _mm_mul_pd(t5, t4);	   \
    z = _mm_add_pd(z,t4); \
  }

#else
#warning Using SSE3
  /* SSE 3 */
#include <pmmintrin.h>
#define CCMUL(z,x,y)		\
  { \
    __m128d t1; \
    __m128d t2 = _mm_set_pd((double)(-1),(double)1);	\
    t1 = _mm_mul_pd((x),(y)); \
    (z) = _mm_hsub_pd(t1,t1);			\
    t1 = _mm_shuffle_pd((y),(y),0x1);\
    t1 = _mm_mul_pd((x),t1); \
    t1 = _mm_hadd_pd(t1,t1); \
    (z)= _mm_shuffle_pd((z),t1,0x2);		\
    (z)= _mm_mul_pd((z),t2); \
  }
#define CCMADD(z,x,y)				\
  { \
    __m128d t1,t2;	      \
    __m128d t3 = _mm_set_pd((double)(-1), (double)1);	\
    t1 = _mm_mul_pd((x),(y)); \
    t1 = _mm_hsub_pd(t1,t1); \
    t2 = _mm_shuffle_pd((y),(y),0x1);\
    t2 = _mm_mul_pd((x),t2); \
    t2 = _mm_hadd_pd(t2,t2); \
    t1= _mm_shuffle_pd(t1,t2,0x2);		\
    t1= _mm_mul_pd(t3,t1); \
    (z) = _mm_add_pd((z),t1);			\
  }

#endif

  // UNALIGNED

  /* M3 = M1*adj(M2) */
  void ssed_m_eq_hh_u(REAL64* m3, REAL64* m2, REAL64* m1, int n_mat)
  {
    __m128d m1_1;
    __m128d m1_2;
    __m128d m1_3;

    __m128d m2_1;
    __m128d m2_2;
    __m128d m2_3;

    __m128d m3_11;
    __m128d m3_12;
    __m128d m3_13;


    REAL64* m1_p=m1;
    REAL64* m2_p=m2;
    REAL64* m3_p=m3;


    for(int i=0; i < n_mat; i++) { 
      m1_1 = _mm_loadu_pd(m1_p);
      m1_2 = _mm_loadu_pd(m1_p+2);
      m1_3 = _mm_loadu_pd(m1_p+4);

      m2_1 = _mm_loadu_pd(m2_p);
      m2_2 = _mm_loadu_pd(m2_p+6);
      m2_3 = _mm_loadu_pd(m2_p+12);

    
      CCMUL(m3_11, m2_1, m1_1);
      m1_1 = _mm_loadu_pd(m1_p+6);
      CCMADD(m3_11, m2_2, m1_2);
      m1_2 = _mm_loadu_pd(m1_p+8);
      CCMADD(m3_11, m2_3, m1_3);
      m1_3 = _mm_loadu_pd(m1_p+10);


      CCMUL(m3_12, m2_1, m1_1);
      m1_1 = _mm_loadu_pd(m1_p+12);
      CCMADD(m3_12, m2_2, m1_2);
      m1_2 = _mm_loadu_pd(m1_p+14);
      CCMADD(m3_12, m2_3, m1_3);
      m1_3 = _mm_loadu_pd(m1_p+16);


      CCMUL(m3_13, m2_1, m1_1);
      m2_1 = _mm_loadu_pd(m2_p+2);
      CCMADD(m3_13, m2_2, m1_2)
      m2_2 = _mm_loadu_pd(m2_p+8);
      CCMADD(m3_13, m2_3, m1_3);
      m2_3 = _mm_loadu_pd(m2_p+14);

      _mm_storeu_pd(m3_p, m3_11);
      _mm_storeu_pd(m3_p+2, m3_12);
      _mm_storeu_pd(m3_p+4, m3_13);

      m1_1 = _mm_loadu_pd(m1_p);
      m1_2 = _mm_loadu_pd(m1_p+2);
      m1_3 = _mm_loadu_pd(m1_p+4);

      CCMUL(m3_11, m2_1, m1_1);
      m1_1 = _mm_loadu_pd(m1_p+6);
      CCMADD(m3_11, m2_2, m1_2);
      m1_2 = _mm_loadu_pd(m1_p+8);
      CCMADD(m3_11, m2_3, m1_3);
      m1_3 = _mm_loadu_pd(m1_p+10);


      CCMUL(m3_12, m2_1, m1_1);
      m1_1 = _mm_loadu_pd(m1_p+12);
      CCMADD(m3_12, m2_2, m1_2);
      m1_2 = _mm_loadu_pd(m1_p+14);
      CCMADD(m3_12, m2_3, m1_3);
      m1_3 = _mm_loadu_pd(m1_p+16);


      CCMUL(m3_13, m2_1, m1_1);
      m2_1 = _mm_loadu_pd(m2_p+4);
      CCMADD(m3_13, m2_2, m1_2)
      m2_2 = _mm_loadu_pd(m2_p+10);
      CCMADD(m3_13, m2_3, m1_3);
      m2_3 = _mm_loadu_pd(m2_p+16);

      _mm_storeu_pd(m3_p+6, m3_11);
      _mm_storeu_pd(m3_p+8, m3_12);
      _mm_storeu_pd(m3_p+10, m3_13);

      m1_1 = _mm_loadu_pd(m1_p);
      m1_2 = _mm_loadu_pd(m1_p+2);
      m1_3 = _mm_loadu_pd(m1_p+4);

      CCMUL(m3_11, m2_1, m1_1);
      m1_1 = _mm_loadu_pd(m1_p+6);
      CCMADD(m3_11, m2_2, m1_2);
      m1_2 = _mm_loadu_pd(m1_p+8);
      CCMADD(m3_11, m2_3, m1_3);
      m1_3 = _mm_loadu_pd(m1_p+10);

      CCMUL(m3_12, m2_1, m1_1);
      m1_1 = _mm_loadu_pd(m1_p+12);
      CCMADD(m3_12, m2_2, m1_2);
      m1_2 = _mm_loadu_pd(m1_p+14);
      CCMADD(m3_12, m2_3, m1_3);
      m1_3 = _mm_loadu_pd(m1_p+16);

      CCMUL(m3_13, m2_1, m1_1);
      CCMADD(m3_13, m2_2, m1_2);
      CCMADD(m3_13, m2_3, m1_3);
	
      _mm_storeu_pd(m3_p+12, m3_11);
      _mm_storeu_pd(m3_p+14, m3_12);
      _mm_storeu_pd(m3_p+16, m3_13);


      /* Next matrix */
      m1_p += 18; m2_p += 18; m3_p += 18;
    }

  }

  /* M3 += a M1*M2 */
  void ssed_m_peq_ahh_u(REAL64* m3, REAL64* a, REAL64* m2, REAL64* m1, int n_mat)
  {
    __m128d m1_1;
    __m128d m1_2;
    __m128d m1_3;

    __m128d m2_1;
    __m128d m2_2;
    __m128d m2_3;

    __m128d m3_11;
    __m128d m3_12;
    __m128d m3_13;

    __m128d res1,res2,res3;

    __m128d tmp1;
    __m128d scalar;

    // Load the scalar into low bytes of scalar
    scalar = _mm_load_sd(a);
  
    // cross components into tmp 
    // Zero tmp
    tmp1 = _mm_setzero_pd();
    tmp1 = _mm_shuffle_pd(scalar, scalar, 0x1);
    scalar = _mm_add_pd(scalar, tmp1);


    REAL64* m1_p=m1;
    REAL64* m2_p=m2;
    REAL64* m3_p=m3;

    
    for(int i =0; i < n_mat; i++) { 
      m1_1 = _mm_loadu_pd(m1_p);
      m1_2 = _mm_loadu_pd(m1_p+2);
      m1_3 = _mm_loadu_pd(m1_p+4);


      m2_1 = _mm_loadu_pd(m2_p);
      m2_2 = _mm_loadu_pd(m2_p+6);
      m2_3 = _mm_loadu_pd(m2_p+12);

      res1 = _mm_loadu_pd(m3_p);
      res2 = _mm_loadu_pd(m3_p+2);
      res3 = _mm_loadu_pd(m3_p+4);


      CCMUL(m3_11, m2_1, m1_1);
      m1_1 = _mm_loadu_pd(m1_p+6);
      CCMADD(m3_11, m2_2, m1_2);
      m1_2 = _mm_loadu_pd(m1_p+8);
      CCMADD(m3_11, m2_3, m1_3);
      m1_3 = _mm_loadu_pd(m1_p+10);

      m3_11 = _mm_mul_pd(scalar,m3_11);
      m3_11 = _mm_add_pd(res1, m3_11);
      

      CCMUL(m3_12, m2_1, m1_1);
      m1_1 = _mm_loadu_pd(m1_p+12);
      CCMADD(m3_12, m2_2, m1_2);
      m1_2 = _mm_loadu_pd(m1_p+14);
      CCMADD(m3_12, m2_3, m1_3);
      m1_3 = _mm_loadu_pd(m1_p+16);

      m3_12 = _mm_mul_pd(scalar,m3_12);
      m3_12 = _mm_add_pd(res2, m3_12);



      CCMUL(m3_13, m2_1, m1_1);
      m2_1 = _mm_loadu_pd(m2_p+2);
      CCMADD(m3_13, m2_2, m1_2)
      m2_2 = _mm_loadu_pd(m2_p+8);
      CCMADD(m3_13, m2_3, m1_3);
      m2_3 = _mm_loadu_pd(m2_p+14);

      m3_13 = _mm_mul_pd(scalar,m3_13);
      m3_13 = _mm_add_pd(res3, m3_13);

     


      _mm_storeu_pd(m3_p, m3_11);
      _mm_storeu_pd(m3_p+2, m3_12);
      _mm_storeu_pd(m3_p+4, m3_13);

      res1 = _mm_loadu_pd(m3_p+6);
      res2 = _mm_loadu_pd(m3_p+8);
      res3 = _mm_loadu_pd(m3_p+10);

      m1_1 = _mm_loadu_pd(m1_p);
      m1_2 = _mm_loadu_pd(m1_p+2);
      m1_3 = _mm_loadu_pd(m1_p+4);

      CCMUL(m3_11, m2_1, m1_1);
      m1_1 = _mm_loadu_pd(m1_p+6);
      CCMADD(m3_11, m2_2, m1_2);
      m1_2 = _mm_loadu_pd(m1_p+8);
      CCMADD(m3_11, m2_3, m1_3);
      m1_3 = _mm_loadu_pd(m1_p+10);      

      m3_11 = _mm_mul_pd(scalar,m3_11);
      m3_11 = _mm_add_pd(res1, m3_11);


      CCMUL(m3_12, m2_1, m1_1);
      m1_1 = _mm_loadu_pd(m1_p+12);
      CCMADD(m3_12, m2_2, m1_2);
      m1_2 = _mm_loadu_pd(m1_p+14);
      CCMADD(m3_12, m2_3, m1_3);
      m1_3 = _mm_loadu_pd(m1_p+16);

      m3_12 = _mm_mul_pd(scalar,m3_12);
      m3_12 = _mm_add_pd(res2, m3_12);



      CCMUL(m3_13, m2_1, m1_1);
      m2_1 = _mm_loadu_pd(m2_p+4);
      CCMADD(m3_13, m2_2, m1_2)
      m2_2 = _mm_loadu_pd(m2_p+10);
      CCMADD(m3_13, m2_3, m1_3);
      m2_3 = _mm_loadu_pd(m2_p+16);

      m3_13 = _mm_mul_pd(scalar,m3_13);
      m3_13 = _mm_add_pd(res3, m3_13);



      _mm_storeu_pd(m3_p+6, m3_11);
      _mm_storeu_pd(m3_p+8, m3_12);
      _mm_storeu_pd(m3_p+10, m3_13);
      res1 = _mm_loadu_pd(m3_p+12);
      res2 = _mm_loadu_pd(m3_p+14);
      res3 = _mm_loadu_pd(m3_p+16);


      m1_1 = _mm_loadu_pd(m1_p);
      m1_2 = _mm_loadu_pd(m1_p+2);
      m1_3 = _mm_loadu_pd(m1_p+4);

      CCMUL(m3_11, m2_1, m1_1);
      m1_1 = _mm_loadu_pd(m1_p+6);
      CCMADD(m3_11, m2_2, m1_2);
      m1_2 = _mm_loadu_pd(m1_p+8);
      CCMADD(m3_11, m2_3, m1_3);
      m1_3 = _mm_loadu_pd(m1_p+10);

      m3_11 = _mm_mul_pd(scalar,m3_11);
      m3_11 = _mm_add_pd(res1, m3_11);



      CCMUL(m3_12, m2_1, m1_1);
      m1_1 = _mm_loadu_pd(m1_p+12);
      CCMADD(m3_12, m2_2, m1_2);
      m1_2 = _mm_loadu_pd(m1_p+14);
      CCMADD(m3_12, m2_3, m1_3);
      m1_3 = _mm_loadu_pd(m1_p+16);

      m3_12 = _mm_mul_pd(scalar,m3_12);
      m3_12 = _mm_add_pd(res2, m3_12);

      CCMUL(m3_13, m2_1, m1_1);
      CCMADD(m3_13, m2_2, m1_2);
      CCMADD(m3_13, m2_3, m1_3);

      m3_13 = _mm_mul_pd(scalar,m3_13);
      m3_13 = _mm_add_pd(res3, m3_13);
	
      _mm_storeu_pd(m3_p+12, m3_11);
      _mm_storeu_pd(m3_p+14, m3_12);
      _mm_storeu_pd(m3_p+16, m3_13);

      m1_p += 18; m2_p += 18; m3_p += 18;
    }

  }


  /* ALIGNED */
  /* M3 = M1*adj(M2) */
  void ssed_m_eq_hh(REAL64* m3, REAL64* m2, REAL64* m1, int n_mat)
  {
    __m128d m1_1;
    __m128d m1_2;
    __m128d m1_3;

    __m128d m2_1;
    __m128d m2_2;
    __m128d m2_3;

    __m128d m3_11;
    __m128d m3_12;
    __m128d m3_13;


    REAL64* m1_p=m1;
    REAL64* m2_p=m2;
    REAL64* m3_p=m3;



    for(int i=0; i < n_mat; i++) { 
      m1_1 = _mm_load_pd(m1_p);
      m1_2 = _mm_load_pd(m1_p+2);
      m1_3 = _mm_load_pd(m1_p+4);

      m2_1 = _mm_load_pd(m2_p);
      m2_2 = _mm_load_pd(m2_p+6);
      m2_3 = _mm_load_pd(m2_p+12);

    
      CCMUL(m3_11, m2_1, m1_1);
      m1_1 = _mm_load_pd(m1_p+6);
      CCMADD(m3_11, m2_2, m1_2);
      m1_2 = _mm_load_pd(m1_p+8);
      CCMADD(m3_11, m2_3, m1_3);
      m1_3 = _mm_load_pd(m1_p+10);


      CCMUL(m3_12, m2_1, m1_1);
      m1_1 = _mm_load_pd(m1_p+12);
      CCMADD(m3_12, m2_2, m1_2);
      m1_2 = _mm_load_pd(m1_p+14);
      CCMADD(m3_12, m2_3, m1_3);
      m1_3 = _mm_load_pd(m1_p+16);


      CCMUL(m3_13, m2_1, m1_1);
      m2_1 = _mm_load_pd(m2_p+2);
      CCMADD(m3_13, m2_2, m1_2)
      m2_2 = _mm_load_pd(m2_p+8);
      CCMADD(m3_13, m2_3, m1_3);
      m2_3 = _mm_load_pd(m2_p+14);

      _mm_store_pd(m3_p, m3_11);
      _mm_store_pd(m3_p+2, m3_12);
      _mm_store_pd(m3_p+4, m3_13);

      m1_1 = _mm_load_pd(m1_p);
      m1_2 = _mm_load_pd(m1_p+2);
      m1_3 = _mm_load_pd(m1_p+4);

      CCMUL(m3_11, m2_1, m1_1);
      m1_1 = _mm_load_pd(m1_p+6);
      CCMADD(m3_11, m2_2, m1_2);
      m1_2 = _mm_load_pd(m1_p+8);
      CCMADD(m3_11, m2_3, m1_3);
      m1_3 = _mm_load_pd(m1_p+10);


      CCMUL(m3_12, m2_1, m1_1);
      m1_1 = _mm_load_pd(m1_p+12);
      CCMADD(m3_12, m2_2, m1_2);
      m1_2 = _mm_load_pd(m1_p+14);
      CCMADD(m3_12, m2_3, m1_3);
      m1_3 = _mm_load_pd(m1_p+16);


      CCMUL(m3_13, m2_1, m1_1);
      m2_1 = _mm_load_pd(m2_p+4);
      CCMADD(m3_13, m2_2, m1_2)
      m2_2 = _mm_load_pd(m2_p+10);
      CCMADD(m3_13, m2_3, m1_3);
      m2_3 = _mm_load_pd(m2_p+16);

      _mm_store_pd(m3_p+6, m3_11);
      _mm_store_pd(m3_p+8, m3_12);
      _mm_store_pd(m3_p+10, m3_13);

      m1_1 = _mm_load_pd(m1_p);
      m1_2 = _mm_load_pd(m1_p+2);
      m1_3 = _mm_load_pd(m1_p+4);

      CCMUL(m3_11, m2_1, m1_1);
      m1_1 = _mm_load_pd(m1_p+6);
      CCMADD(m3_11, m2_2, m1_2);
      m1_2 = _mm_load_pd(m1_p+8);
      CCMADD(m3_11, m2_3, m1_3);
      m1_3 = _mm_load_pd(m1_p+10);

      CCMUL(m3_12, m2_1, m1_1);
      m1_1 = _mm_load_pd(m1_p+12);
      CCMADD(m3_12, m2_2, m1_2);
      m1_2 = _mm_load_pd(m1_p+14);
      CCMADD(m3_12, m2_3, m1_3);
      m1_3 = _mm_load_pd(m1_p+16);

      CCMUL(m3_13, m2_1, m1_1);
      CCMADD(m3_13, m2_2, m1_2);
      CCMADD(m3_13, m2_3, m1_3);
	
      _mm_store_pd(m3_p+12, m3_11);
      _mm_store_pd(m3_p+14, m3_12);
      _mm_store_pd(m3_p+16, m3_13);


      /* Next matrix */
      m1_p += 18; m2_p += 18; m3_p += 18;
    }

  }

  /* M3 += a M1*M2 */
  void ssed_m_peq_ahh(REAL64* m3, REAL64* a, REAL64* m2, REAL64* m1, int n_mat)
  {
    __m128d m1_1;
    __m128d m1_2;
    __m128d m1_3;

    __m128d m2_1;
    __m128d m2_2;
    __m128d m2_3;

    __m128d m3_11;
    __m128d m3_12;
    __m128d m3_13;

    __m128d res1,res2,res3;

    __m128d tmp1;
    __m128d scalar;

    // Load the scalar into low bytes of scalar
    scalar = _mm_load_sd(a);
  
    // cross components into tmp 
    // Zero tmp
    tmp1 = _mm_setzero_pd();
    tmp1 = _mm_shuffle_pd(scalar, scalar, 0x1);
    scalar = _mm_add_pd(scalar, tmp1);


    REAL64* m1_p=m1;
    REAL64* m2_p=m2;
    REAL64* m3_p=m3;

    for(int i =0; i < n_mat; i++) { 
      m1_1 = _mm_load_pd(m1_p);
      m1_2 = _mm_load_pd(m1_p+2);
      m1_3 = _mm_load_pd(m1_p+4);


      m2_1 = _mm_load_pd(m2_p);
      m2_2 = _mm_load_pd(m2_p+6);
      m2_3 = _mm_load_pd(m2_p+12);

      res1 = _mm_load_pd(m3_p);
      res2 = _mm_load_pd(m3_p+2);
      res3 = _mm_load_pd(m3_p+4);


      CCMUL(m3_11, m2_1, m1_1);
      m1_1 = _mm_load_pd(m1_p+6);
      CCMADD(m3_11, m2_2, m1_2);
      m1_2 = _mm_load_pd(m1_p+8);
      CCMADD(m3_11, m2_3, m1_3);
      m1_3 = _mm_load_pd(m1_p+10);

      m3_11 = _mm_mul_pd(scalar,m3_11);
      m3_11 = _mm_add_pd(res1, m3_11);
      

      CCMUL(m3_12, m2_1, m1_1);
      m1_1 = _mm_load_pd(m1_p+12);
      CCMADD(m3_12, m2_2, m1_2);
      m1_2 = _mm_load_pd(m1_p+14);
      CCMADD(m3_12, m2_3, m1_3);
      m1_3 = _mm_load_pd(m1_p+16);

      m3_12 = _mm_mul_pd(scalar,m3_12);
      m3_12 = _mm_add_pd(res2, m3_12);



      CCMUL(m3_13, m2_1, m1_1);
      m2_1 = _mm_load_pd(m2_p+2);
      CCMADD(m3_13, m2_2, m1_2)
      m2_2 = _mm_load_pd(m2_p+8);
      CCMADD(m3_13, m2_3, m1_3);
      m2_3 = _mm_load_pd(m2_p+14);

      m3_13 = _mm_mul_pd(scalar,m3_13);
      m3_13 = _mm_add_pd(res3, m3_13);

     


      _mm_store_pd(m3_p, m3_11);
      _mm_store_pd(m3_p+2, m3_12);
      _mm_store_pd(m3_p+4, m3_13);

      res1 = _mm_load_pd(m3_p+6);
      res2 = _mm_load_pd(m3_p+8);
      res3 = _mm_load_pd(m3_p+10);

      m1_1 = _mm_load_pd(m1_p);
      m1_2 = _mm_load_pd(m1_p+2);
      m1_3 = _mm_load_pd(m1_p+4);

      CCMUL(m3_11, m2_1, m1_1);
      m1_1 = _mm_load_pd(m1_p+6);
      CCMADD(m3_11, m2_2, m1_2);
      m1_2 = _mm_load_pd(m1_p+8);
      CCMADD(m3_11, m2_3, m1_3);
      m1_3 = _mm_load_pd(m1_p+10);      

      m3_11 = _mm_mul_pd(scalar,m3_11);
      m3_11 = _mm_add_pd(res1, m3_11);


      CCMUL(m3_12, m2_1, m1_1);
      m1_1 = _mm_load_pd(m1_p+12);
      CCMADD(m3_12, m2_2, m1_2);
      m1_2 = _mm_load_pd(m1_p+14);
      CCMADD(m3_12, m2_3, m1_3);
      m1_3 = _mm_load_pd(m1_p+16);

      m3_12 = _mm_mul_pd(scalar,m3_12);
      m3_12 = _mm_add_pd(res2, m3_12);



      CCMUL(m3_13, m2_1, m1_1);
      m2_1 = _mm_load_pd(m2_p+4);
      CCMADD(m3_13, m2_2, m1_2)
      m2_2 = _mm_load_pd(m2_p+10);
      CCMADD(m3_13, m2_3, m1_3);
      m2_3 = _mm_load_pd(m2_p+16);

      m3_13 = _mm_mul_pd(scalar,m3_13);
      m3_13 = _mm_add_pd(res3, m3_13);



      _mm_store_pd(m3_p+6, m3_11);
      _mm_store_pd(m3_p+8, m3_12);
      _mm_store_pd(m3_p+10, m3_13);
      res1 = _mm_load_pd(m3_p+12);
      res2 = _mm_load_pd(m3_p+14);
      res3 = _mm_load_pd(m3_p+16);


      m1_1 = _mm_load_pd(m1_p);
      m1_2 = _mm_load_pd(m1_p+2);
      m1_3 = _mm_load_pd(m1_p+4);

      CCMUL(m3_11, m2_1, m1_1);
      m1_1 = _mm_load_pd(m1_p+6);
      CCMADD(m3_11, m2_2, m1_2);
      m1_2 = _mm_load_pd(m1_p+8);
      CCMADD(m3_11, m2_3, m1_3);
      m1_3 = _mm_load_pd(m1_p+10);

      m3_11 = _mm_mul_pd(scalar,m3_11);
      m3_11 = _mm_add_pd(res1, m3_11);



      CCMUL(m3_12, m2_1, m1_1);
      m1_1 = _mm_load_pd(m1_p+12);
      CCMADD(m3_12, m2_2, m1_2);
      m1_2 = _mm_load_pd(m1_p+14);
      CCMADD(m3_12, m2_3, m1_3);
      m1_3 = _mm_load_pd(m1_p+16);

      m3_12 = _mm_mul_pd(scalar,m3_12);
      m3_12 = _mm_add_pd(res2, m3_12);

      CCMUL(m3_13, m2_1, m1_1);
      CCMADD(m3_13, m2_2, m1_2);
      CCMADD(m3_13, m2_3, m1_3);

      m3_13 = _mm_mul_pd(scalar,m3_13);
      m3_13 = _mm_add_pd(res3, m3_13);
	
      _mm_store_pd(m3_p+12, m3_11);
      _mm_store_pd(m3_p+14, m3_12);
      _mm_store_pd(m3_p+16, m3_13);

      m1_p += 18; m2_p += 18; m3_p += 18;
    }

  }



} // namespace QDP;

