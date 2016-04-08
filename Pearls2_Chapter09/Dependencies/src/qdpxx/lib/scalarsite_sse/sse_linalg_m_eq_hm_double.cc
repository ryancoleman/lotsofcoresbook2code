// $Id: sse_linalg_m_eq_hm_double.cc,v 1.4 2009-02-10 17:06:37 bjoo Exp $

/*! @file
 *  @brief Generic Scalar VAXPY routine
 *
 */

#include "scalarsite_sse/sse_linalg_mm_su3_double.h"

namespace QDP {

#include <xmmintrin.h>
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




  /* M3 = adj(M2)*M1 */
  void ssed_m_eq_hm_u(REAL64* m3, REAL64* m2, REAL64* m1, int n_mat)
  {
    __m128d tmp1;
    __m128d tmp2;
    __m128d tmp3;

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
      m2_1 = _mm_loadu_pd(m2_p);    
      m2_2 = _mm_loadu_pd(m2_p+6);  
      m2_3 = _mm_loadu_pd(m2_p+12); 

      tmp1 = _mm_loadu_pd(m1_p);      
      tmp2 = _mm_loadu_pd(m1_p+2);
      tmp3 = _mm_loadu_pd(m1_p+4);

      CONJMUL(m3_11, tmp1, m2_1);   
      tmp1 = _mm_loadu_pd(m1_p+6);      
      CONJMUL(m3_12, tmp2, m2_1);  
      tmp2 = _mm_loadu_pd(m1_p+8);      
      CONJMUL(m3_13, tmp3, m2_1); 
      tmp3 = _mm_loadu_pd(m1_p+10);      

      CONJMADD(m3_11, tmp1, m2_2);   
      tmp1 = _mm_loadu_pd(m1_p+12);      
      CONJMADD(m3_12, tmp2, m2_2);  
      tmp2 = _mm_loadu_pd(m1_p+14);      
      CONJMADD(m3_13, tmp3, m2_2); 
      tmp3 = _mm_loadu_pd(m1_p+16);      

      CONJMADD(m3_11, tmp1, m2_3);  
      _mm_storeu_pd(m3_p, m3_11);
      CONJMADD(m3_12, tmp2, m2_3);  
      _mm_storeu_pd(m3_p+2, m3_12);
      CONJMADD(m3_13, tmp3, m2_3); 
      _mm_storeu_pd(m3_p+4, m3_13);

      m2_1 = _mm_loadu_pd(m2_p+2);    
      m2_2 = _mm_loadu_pd(m2_p+8);  
      m2_3 = _mm_loadu_pd(m2_p+14); 

      tmp1 = _mm_loadu_pd(m1_p);      
      tmp2 = _mm_loadu_pd(m1_p+2);
      tmp3 = _mm_loadu_pd(m1_p+4);

      CONJMUL(m3_11, tmp1, m2_1);   
      tmp1 = _mm_loadu_pd(m1_p+6);      
      CONJMUL(m3_12, tmp2, m2_1);  
      tmp2 = _mm_loadu_pd(m1_p+8);      
      CONJMUL(m3_13, tmp3, m2_1); 
      tmp3 = _mm_loadu_pd(m1_p+10);      

      CONJMADD(m3_11, tmp1, m2_2);   
      tmp1 = _mm_loadu_pd(m1_p+12);      
      CONJMADD(m3_12, tmp2, m2_2);  
      tmp2 = _mm_loadu_pd(m1_p+14);      
      CONJMADD(m3_13, tmp3, m2_2); 
      tmp3 = _mm_loadu_pd(m1_p+16);      

      CONJMADD(m3_11, tmp1, m2_3);  
      _mm_storeu_pd(m3_p+6, m3_11);
      CONJMADD(m3_12, tmp2, m2_3);  
      _mm_storeu_pd(m3_p+8, m3_12);
      CONJMADD(m3_13, tmp3, m2_3); 
      _mm_storeu_pd(m3_p+10, m3_13);

      m2_1 = _mm_loadu_pd(m2_p+4);    
      m2_2 = _mm_loadu_pd(m2_p+10);  
      m2_3 = _mm_loadu_pd(m2_p+16); 

      tmp1 = _mm_loadu_pd(m1_p);      
      tmp2 = _mm_loadu_pd(m1_p+2);
      tmp3 = _mm_loadu_pd(m1_p+4);

      CONJMUL(m3_11, tmp1, m2_1);   
      tmp1 = _mm_loadu_pd(m1_p+6);      
      CONJMUL(m3_12, tmp2, m2_1);  
      tmp2 = _mm_loadu_pd(m1_p+8);      
      CONJMUL(m3_13, tmp3, m2_1); 
      tmp3 = _mm_loadu_pd(m1_p+10);      

      CONJMADD(m3_11, tmp1, m2_2);   
      tmp1 = _mm_loadu_pd(m1_p+12);      
      CONJMADD(m3_12, tmp2, m2_2);  
      tmp2 = _mm_loadu_pd(m1_p+14);      
      CONJMADD(m3_13, tmp3, m2_2); 
      tmp3 = _mm_loadu_pd(m1_p+16);      

      CONJMADD(m3_11, tmp1, m2_3);  
      _mm_storeu_pd(m3_p+12, m3_11);
      CONJMADD(m3_12, tmp2, m2_3);  
      _mm_storeu_pd(m3_p+14, m3_12);
      CONJMADD(m3_13, tmp3, m2_3); 
      _mm_storeu_pd(m3_p+16, m3_13);

      /* Next matrix */
      m1_p += 18; m2_p += 18; m3_p += 18;
    }

  }

  /* M3 += a M1*M2 */
  void ssed_m_peq_ahm_u(REAL64* m3, REAL64* a, REAL64* m2, REAL64* m1, int n_mat)
  {
    __m128d m2_1;
    __m128d m2_2;
    __m128d m2_3;

    __m128d res1;
    __m128d res2;
    __m128d res3;

    __m128d m3_11;
    __m128d m3_12;
    __m128d m3_13;

    __m128d tmp1;
    __m128d tmp2;
    __m128d tmp3;

    __m128d scalar;
    
      // Zero tmp
    scalar = _mm_load_sd(a);
    tmp1 = _mm_setzero_pd();
    tmp1 = _mm_shuffle_pd(scalar, scalar, 0x1);
    scalar = _mm_add_pd(scalar, tmp1);

    REAL64* m1_p=m1;
    REAL64* m2_p=m2;
    REAL64* m3_p=m3;

    for(int i =0; i < n_mat; i++) { 
      res1 = _mm_loadu_pd(m3_p);
      res2 = _mm_loadu_pd(m3_p+2);
      res3 = _mm_loadu_pd(m3_p+4);

      m2_1 = _mm_loadu_pd(m2_p);    
      m2_2 = _mm_loadu_pd(m2_p+6);  
      m2_3 = _mm_loadu_pd(m2_p+12); 

      tmp1 = _mm_loadu_pd(m1_p);      
      tmp2 = _mm_loadu_pd(m1_p+2);
      tmp3 = _mm_loadu_pd(m1_p+4);

      CONJMUL(m3_11, tmp1, m2_1);   
      tmp1 = _mm_loadu_pd(m1_p+6);      
      CONJMUL(m3_12, tmp2, m2_1);  
      tmp2 = _mm_loadu_pd(m1_p+8);      
      CONJMUL(m3_13, tmp3, m2_1); 
      tmp3 = _mm_loadu_pd(m1_p+10);      

      CONJMADD(m3_11, tmp1, m2_2);   
      tmp1 = _mm_loadu_pd(m1_p+12);      
      CONJMADD(m3_12, tmp2, m2_2);  
      tmp2 = _mm_loadu_pd(m1_p+14);      
      CONJMADD(m3_13, tmp3, m2_2); 
      tmp3 = _mm_loadu_pd(m1_p+16);      

      CONJMADD(m3_11, tmp1, m2_3);  
      CONJMADD(m3_12, tmp2, m2_3);  
      CONJMADD(m3_13, tmp3, m2_3); 

      m3_11 = _mm_mul_pd(scalar, m3_11);
      m3_11 = _mm_add_pd(res1, m3_11);

      m3_12 = _mm_mul_pd(scalar, m3_12);
      m3_12 = _mm_add_pd(res2, m3_12);

      m3_13 = _mm_mul_pd(scalar, m3_13);
      m3_13 = _mm_add_pd(res3, m3_13);

      _mm_storeu_pd(m3_p, m3_11);
      _mm_storeu_pd(m3_p+2, m3_12);
      _mm_storeu_pd(m3_p+4, m3_13);

      res1=_mm_loadu_pd(m3_p+6);
      res2=_mm_loadu_pd(m3_p+8);
      res3=_mm_loadu_pd(m3_p+10);


      m2_1 = _mm_loadu_pd(m2_p+2);    
      m2_2 = _mm_loadu_pd(m2_p+8);  
      m2_3 = _mm_loadu_pd(m2_p+14); 

      tmp1 = _mm_loadu_pd(m1_p);      
      tmp2 = _mm_loadu_pd(m1_p+2);
      tmp3 = _mm_loadu_pd(m1_p+4);

      CONJMUL(m3_11, tmp1, m2_1);   
      tmp1 = _mm_loadu_pd(m1_p+6);      
      CONJMUL(m3_12, tmp2, m2_1);  
      tmp2 = _mm_loadu_pd(m1_p+8);      
      CONJMUL(m3_13, tmp3, m2_1); 
      tmp3 = _mm_loadu_pd(m1_p+10);      

      CONJMADD(m3_11, tmp1, m2_2);   
      tmp1 = _mm_loadu_pd(m1_p+12);      
      CONJMADD(m3_12, tmp2, m2_2);  
      tmp2 = _mm_loadu_pd(m1_p+14);      
      CONJMADD(m3_13, tmp3, m2_2); 
      tmp3 = _mm_loadu_pd(m1_p+16);      

      CONJMADD(m3_11, tmp1, m2_3);  
      CONJMADD(m3_12, tmp2, m2_3);  
      CONJMADD(m3_13, tmp3, m2_3);

      m3_11 = _mm_mul_pd(scalar, m3_11);
      m3_11 = _mm_add_pd(res1,m3_11);

      m3_12 = _mm_mul_pd(scalar, m3_12);
      m3_12 = _mm_add_pd(res2,m3_12);

      m3_13 = _mm_mul_pd(scalar, m3_13);
      m3_13 = _mm_add_pd(res3, m3_13);

      _mm_storeu_pd(m3_p+6, m3_11);
      _mm_storeu_pd(m3_p+8, m3_12);
      _mm_storeu_pd(m3_p+10, m3_13);

      res1=_mm_loadu_pd(m3_p+12);
      res2=_mm_loadu_pd(m3_p+14);
      res3=_mm_loadu_pd(m3_p+16);

      m2_1 = _mm_loadu_pd(m2_p+4);    
      m2_2 = _mm_loadu_pd(m2_p+10);  
      m2_3 = _mm_loadu_pd(m2_p+16); 

      tmp1 = _mm_loadu_pd(m1_p);      
      tmp2 = _mm_loadu_pd(m1_p+2);
      tmp3 = _mm_loadu_pd(m1_p+4);

      CONJMUL(m3_11, tmp1, m2_1);   
      tmp1 = _mm_loadu_pd(m1_p+6);      
      CONJMUL(m3_12, tmp2, m2_1);  
      tmp2 = _mm_loadu_pd(m1_p+8);      
      CONJMUL(m3_13, tmp3, m2_1); 
      tmp3 = _mm_loadu_pd(m1_p+10);      

      CONJMADD(m3_11, tmp1, m2_2);   
      tmp1 = _mm_loadu_pd(m1_p+12);      
      CONJMADD(m3_12, tmp2, m2_2);  
      tmp2 = _mm_loadu_pd(m1_p+14);      
      CONJMADD(m3_13, tmp3, m2_2); 
      tmp3 = _mm_loadu_pd(m1_p+16);      

      CONJMADD(m3_11, tmp1, m2_3);  
      CONJMADD(m3_12, tmp2, m2_3);  
      CONJMADD(m3_13, tmp3, m2_3); 

      m3_11 = _mm_mul_pd(scalar, m3_11);
      m3_11 = _mm_add_pd(res1,m3_11);

      m3_12 = _mm_mul_pd(scalar, m3_12);
      m3_12 = _mm_add_pd(res2,m3_12);

      m3_13 = _mm_mul_pd(scalar, m3_13);
      m3_13 = _mm_add_pd(res3, m3_13);

      _mm_storeu_pd(m3_p+12, m3_11);
      _mm_storeu_pd(m3_p+14, m3_12);
      _mm_storeu_pd(m3_p+16, m3_13);
      m1_p += 18; m2_p += 18; m3_p += 18;
    }

  }


  /* ALIGNED VERSIONS TO USE IN EVALUATES */

  /* M3 = adj(M2)*M1 */
  void ssed_m_eq_hm(REAL64* m3, REAL64* m2, REAL64* m1, int n_mat)
  {
    __m128d tmp1;
    __m128d tmp2;
    __m128d tmp3;

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
      m2_1 = _mm_load_pd(m2_p);    
      m2_2 = _mm_load_pd(m2_p+6);  
      m2_3 = _mm_load_pd(m2_p+12); 

      tmp1 = _mm_load_pd(m1_p);      
      tmp2 = _mm_load_pd(m1_p+2);
      tmp3 = _mm_load_pd(m1_p+4);

      CONJMUL(m3_11, tmp1, m2_1);   
      tmp1 = _mm_load_pd(m1_p+6);      
      CONJMUL(m3_12, tmp2, m2_1);  
      tmp2 = _mm_load_pd(m1_p+8);      
      CONJMUL(m3_13, tmp3, m2_1); 
      tmp3 = _mm_load_pd(m1_p+10);      

      CONJMADD(m3_11, tmp1, m2_2);   
      tmp1 = _mm_load_pd(m1_p+12);      
      CONJMADD(m3_12, tmp2, m2_2);  
      tmp2 = _mm_load_pd(m1_p+14);      
      CONJMADD(m3_13, tmp3, m2_2); 
      tmp3 = _mm_load_pd(m1_p+16);      

      CONJMADD(m3_11, tmp1, m2_3);  
      _mm_store_pd(m3_p, m3_11);
      CONJMADD(m3_12, tmp2, m2_3);  
      _mm_store_pd(m3_p+2, m3_12);
      CONJMADD(m3_13, tmp3, m2_3); 
      _mm_store_pd(m3_p+4, m3_13);

      m2_1 = _mm_load_pd(m2_p+2);    
      m2_2 = _mm_load_pd(m2_p+8);  
      m2_3 = _mm_load_pd(m2_p+14); 

      tmp1 = _mm_load_pd(m1_p);      
      tmp2 = _mm_load_pd(m1_p+2);
      tmp3 = _mm_load_pd(m1_p+4);

      CONJMUL(m3_11, tmp1, m2_1);   
      tmp1 = _mm_load_pd(m1_p+6);      
      CONJMUL(m3_12, tmp2, m2_1);  
      tmp2 = _mm_load_pd(m1_p+8);      
      CONJMUL(m3_13, tmp3, m2_1); 
      tmp3 = _mm_load_pd(m1_p+10);      

      CONJMADD(m3_11, tmp1, m2_2);   
      tmp1 = _mm_load_pd(m1_p+12);      
      CONJMADD(m3_12, tmp2, m2_2);  
      tmp2 = _mm_load_pd(m1_p+14);      
      CONJMADD(m3_13, tmp3, m2_2); 
      tmp3 = _mm_load_pd(m1_p+16);      

      CONJMADD(m3_11, tmp1, m2_3);  
      _mm_store_pd(m3_p+6, m3_11);
      CONJMADD(m3_12, tmp2, m2_3);  
      _mm_store_pd(m3_p+8, m3_12);
      CONJMADD(m3_13, tmp3, m2_3); 
      _mm_store_pd(m3_p+10, m3_13);

      m2_1 = _mm_load_pd(m2_p+4);    
      m2_2 = _mm_load_pd(m2_p+10);  
      m2_3 = _mm_load_pd(m2_p+16); 

      tmp1 = _mm_load_pd(m1_p);      
      tmp2 = _mm_load_pd(m1_p+2);
      tmp3 = _mm_load_pd(m1_p+4);

      CONJMUL(m3_11, tmp1, m2_1);   
      tmp1 = _mm_load_pd(m1_p+6);      
      CONJMUL(m3_12, tmp2, m2_1);  
      tmp2 = _mm_load_pd(m1_p+8);      
      CONJMUL(m3_13, tmp3, m2_1); 
      tmp3 = _mm_load_pd(m1_p+10);      

      CONJMADD(m3_11, tmp1, m2_2);   
      tmp1 = _mm_load_pd(m1_p+12);      
      CONJMADD(m3_12, tmp2, m2_2);  
      tmp2 = _mm_load_pd(m1_p+14);      
      CONJMADD(m3_13, tmp3, m2_2); 
      tmp3 = _mm_load_pd(m1_p+16);      

      CONJMADD(m3_11, tmp1, m2_3);  
      _mm_store_pd(m3_p+12, m3_11);
      CONJMADD(m3_12, tmp2, m2_3);  
      _mm_store_pd(m3_p+14, m3_12);
      CONJMADD(m3_13, tmp3, m2_3); 
      _mm_store_pd(m3_p+16, m3_13);

      /* Next matrix */
      m1_p += 18; m2_p += 18; m3_p += 18;
    }

  }

  /* M3 += a M1*M2 */
  void ssed_m_peq_ahm(REAL64* m3, REAL64* a, REAL64* m2, REAL64* m1, int n_mat)
  {
    __m128d m2_1;
    __m128d m2_2;
    __m128d m2_3;

    __m128d res1;
    __m128d res2;
    __m128d res3;

    __m128d m3_11;
    __m128d m3_12;
    __m128d m3_13;

    __m128d tmp1;
    __m128d tmp2;
    __m128d tmp3;

    __m128d scalar;
    
      // Zero tmp
    scalar = _mm_load_sd(a);
    tmp1 = _mm_setzero_pd();
    tmp1 = _mm_shuffle_pd(scalar, scalar, 0x1);
    scalar = _mm_add_pd(scalar, tmp1);

    REAL64* m1_p=m1;
    REAL64* m2_p=m2;
    REAL64* m3_p=m3;

    for(int i =0; i < n_mat; i++) { 
      res1 = _mm_load_pd(m3_p);
      res2 = _mm_load_pd(m3_p+2);
      res3 = _mm_load_pd(m3_p+4);

      m2_1 = _mm_load_pd(m2_p);    
      m2_2 = _mm_load_pd(m2_p+6);  
      m2_3 = _mm_load_pd(m2_p+12); 

      tmp1 = _mm_load_pd(m1_p);      
      tmp2 = _mm_load_pd(m1_p+2);
      tmp3 = _mm_load_pd(m1_p+4);

      CONJMUL(m3_11, tmp1, m2_1);   
      tmp1 = _mm_load_pd(m1_p+6);      
      CONJMUL(m3_12, tmp2, m2_1);  
      tmp2 = _mm_load_pd(m1_p+8);      
      CONJMUL(m3_13, tmp3, m2_1); 
      tmp3 = _mm_load_pd(m1_p+10);      

      CONJMADD(m3_11, tmp1, m2_2);   
      tmp1 = _mm_load_pd(m1_p+12);      
      CONJMADD(m3_12, tmp2, m2_2);  
      tmp2 = _mm_load_pd(m1_p+14);      
      CONJMADD(m3_13, tmp3, m2_2); 
      tmp3 = _mm_load_pd(m1_p+16);      

      CONJMADD(m3_11, tmp1, m2_3);  
      CONJMADD(m3_12, tmp2, m2_3);  
      CONJMADD(m3_13, tmp3, m2_3); 

      m3_11 = _mm_mul_pd(scalar, m3_11);
      m3_11 = _mm_add_pd(res1, m3_11);

      m3_12 = _mm_mul_pd(scalar, m3_12);
      m3_12 = _mm_add_pd(res2, m3_12);

      m3_13 = _mm_mul_pd(scalar, m3_13);
      m3_13 = _mm_add_pd(res3, m3_13);

      _mm_store_pd(m3_p, m3_11);
      _mm_store_pd(m3_p+2, m3_12);
      _mm_store_pd(m3_p+4, m3_13);

      res1=_mm_load_pd(m3_p+6);
      res2=_mm_load_pd(m3_p+8);
      res3=_mm_load_pd(m3_p+10);


      m2_1 = _mm_load_pd(m2_p+2);    
      m2_2 = _mm_load_pd(m2_p+8);  
      m2_3 = _mm_load_pd(m2_p+14); 

      tmp1 = _mm_load_pd(m1_p);      
      tmp2 = _mm_load_pd(m1_p+2);
      tmp3 = _mm_load_pd(m1_p+4);

      CONJMUL(m3_11, tmp1, m2_1);   
      tmp1 = _mm_load_pd(m1_p+6);      
      CONJMUL(m3_12, tmp2, m2_1);  
      tmp2 = _mm_load_pd(m1_p+8);      
      CONJMUL(m3_13, tmp3, m2_1); 
      tmp3 = _mm_load_pd(m1_p+10);      

      CONJMADD(m3_11, tmp1, m2_2);   
      tmp1 = _mm_load_pd(m1_p+12);      
      CONJMADD(m3_12, tmp2, m2_2);  
      tmp2 = _mm_load_pd(m1_p+14);      
      CONJMADD(m3_13, tmp3, m2_2); 
      tmp3 = _mm_load_pd(m1_p+16);      

      CONJMADD(m3_11, tmp1, m2_3);  
      CONJMADD(m3_12, tmp2, m2_3);  
      CONJMADD(m3_13, tmp3, m2_3);

      m3_11 = _mm_mul_pd(scalar, m3_11);
      m3_11 = _mm_add_pd(res1,m3_11);

      m3_12 = _mm_mul_pd(scalar, m3_12);
      m3_12 = _mm_add_pd(res2,m3_12);

      m3_13 = _mm_mul_pd(scalar, m3_13);
      m3_13 = _mm_add_pd(res3, m3_13);

      _mm_store_pd(m3_p+6, m3_11);
      _mm_store_pd(m3_p+8, m3_12);
      _mm_store_pd(m3_p+10, m3_13);

      res1=_mm_load_pd(m3_p+12);
      res2=_mm_load_pd(m3_p+14);
      res3=_mm_load_pd(m3_p+16);

      m2_1 = _mm_load_pd(m2_p+4);    
      m2_2 = _mm_load_pd(m2_p+10);  
      m2_3 = _mm_load_pd(m2_p+16); 

      tmp1 = _mm_load_pd(m1_p);      
      tmp2 = _mm_load_pd(m1_p+2);
      tmp3 = _mm_load_pd(m1_p+4);

      CONJMUL(m3_11, tmp1, m2_1);   
      tmp1 = _mm_load_pd(m1_p+6);      
      CONJMUL(m3_12, tmp2, m2_1);  
      tmp2 = _mm_load_pd(m1_p+8);      
      CONJMUL(m3_13, tmp3, m2_1); 
      tmp3 = _mm_load_pd(m1_p+10);      

      CONJMADD(m3_11, tmp1, m2_2);   
      tmp1 = _mm_load_pd(m1_p+12);      
      CONJMADD(m3_12, tmp2, m2_2);  
      tmp2 = _mm_load_pd(m1_p+14);      
      CONJMADD(m3_13, tmp3, m2_2); 
      tmp3 = _mm_load_pd(m1_p+16);      

      CONJMADD(m3_11, tmp1, m2_3);  
      CONJMADD(m3_12, tmp2, m2_3);  
      CONJMADD(m3_13, tmp3, m2_3); 

      m3_11 = _mm_mul_pd(scalar, m3_11);
      m3_11 = _mm_add_pd(res1,m3_11);

      m3_12 = _mm_mul_pd(scalar, m3_12);
      m3_12 = _mm_add_pd(res2,m3_12);

      m3_13 = _mm_mul_pd(scalar, m3_13);
      m3_13 = _mm_add_pd(res3, m3_13);

      _mm_store_pd(m3_p+12, m3_11);
      _mm_store_pd(m3_p+14, m3_12);
      _mm_store_pd(m3_p+16, m3_13);
      m1_p += 18; m2_p += 18; m3_p += 18;
    }

  }



} // namespace QDP;

