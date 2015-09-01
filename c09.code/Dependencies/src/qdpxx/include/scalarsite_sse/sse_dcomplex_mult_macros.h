#ifndef SSE_DCOMPLEX_MULT_MACROS
#define SSE_DCOMPLEX_MULT_MACROS

/* This selects the configuration */
#include "qdp_config.h"

/* SSE 2 Headers */
#include<xmmintrin.h>

/* A useful union type allows me to set values into the 
   vector  from code */

#ifndef QDP_USE_SSE3

/* SSE2 macros */
/* z = x*y    z, x, y are SSE registers containing complex numbers
              ordered with the real part in the low half, imag part 
	      in the upper half */
#define CMUL(z,x,y)		\
  { \
    __m128d t1,t2,t3,t4; \
    t1 = _mm_mul_pd(x,y); \
    t2 = _mm_shuffle_pd(t1,t1,0x1); \
    t3 = _mm_shuffle_pd(y,y,0x1);\
    z = _mm_sub_pd(t1,t2); \
    t2 = _mm_mul_pd(x,t3); \
    t3 = _mm_shuffle_pd(t2,t2,0x1); \
    t3 = _mm_add_pd(t2,t3); \
    z= _mm_shuffle_pd(z,t3,0x2); \
  }

/* z += x*y    z, x, y are SSE registers containing complex numbers
   ordered with the real part in the low half, imag part 
   in the upper half */
#define CMADD(z,x,y)				\
  { \
    __m128d t1,t2,t3,t4; \
    t1 = _mm_mul_pd(x,y); \
    t2 = _mm_shuffle_pd(t1,t1,0x1); \
    t3 = _mm_shuffle_pd(y,y,0x1);\
    t4 = _mm_sub_pd(t1,t2); \
    t2 = _mm_mul_pd(x,t3); \
    t3 = _mm_shuffle_pd(t2,t2,0x1); \
    t3 = _mm_add_pd(t2,t3); \
    t4= _mm_shuffle_pd(t4,t3,0x2); \
    z = _mm_add_pd(z,t4); \
  }

/* z = conj(x)*conj(y)    z, x, y are SSE registers containing complex numbers
   ordered with the real part in the low half, imag part 
   in the upper half */
#define CCMUL(z,x,y)		\
  { \
    __m128d t1,t2,t3; \
    __m128d t4 = _mm_set_pd( (double)(-1), (double)1 ); \
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

/* z += conj(x)*conj(y)    z, x, y are SSE registers containing complex numbers
   ordered with the real part in the low half, imag part 
   in the upper half */
#define CCMADD(z,x,y)				\
  { \
    __m128d t1,t2,t3,t4; \
    __m128d t5 = _mm_set_pd( (double)(-1), (double)1 ); \
    t1 = _mm_mul_pd(x,y); \
    t2 = _mm_shuffle_pd(t1,t1,0x1); \
    t3 = _mm_shuffle_pd(y,y,0x1);\
    t4 = _mm_sub_pd(t1,t2); \
    t2 = _mm_mul_pd(x,t3); \
    t3 = _mm_shuffle_pd(t2,t2,0x1); \
    t3 = _mm_add_pd(t2,t3); \
    t4= _mm_shuffle_pd(t4,t3,0x2); \
    t4= _mm_mul_pd(t5, t4); \
    z = _mm_add_pd(z,t4); \
  }


/* z = x*conj(y)    z, x, y are SSE registers containing complex numbers
   ordered with the real part in the low half, imag part 
   in the upper half */
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

/* z += x*conj(y)    z, x, y are SSE registers containing complex numbers
   ordered with the real part in the low half, imag part 
   in the upper half */
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
#warning "Using SSE3"
/* SSE 3 */
#include <pmmintrin.h>

/* z = x*y    z, x, y are SSE registers containing complex numbers
              ordered with the real part in the low half, imag part 
	      in the upper half */
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

/* z = conj(x)*conj(y)    z, x, y are SSE registers containing complex numbers
              ordered with the real part in the low half, imag part 
	      in the upper half */
#define CCMUL(z,x,y)		\
  { \
    __m128d t1; \
    __m128d t2 = _mm_set_pd((double)(-1), (double)1 ); \
    t1 = _mm_mul_pd((x),(y)); \
    (z) = _mm_hsub_pd(t1,t1);			\
    t1 = _mm_shuffle_pd((y),(y),0x1);\
    t1 = _mm_mul_pd((x),t1); \
    t1 = _mm_hadd_pd(t1,t1); \
    (z)= _mm_shuffle_pd((z),t1,0x2);		\
    (z)= _mm_mul_pd((z),t2); \
  }

/* z += x*y    z, x, y are SSE registers containing complex numbers
              ordered with the real part in the low half, imag part 
	      in the upper half */
#define CMADD(z,x,y)				\
  { \
    __m128d t1,t2;	      \
    t1 = _mm_mul_pd((x),(y)); \
    t1 = _mm_hsub_pd(t1,t1); \
    t2 = _mm_shuffle_pd((y),(y),0x1);\
    t2 = _mm_mul_pd((x),t2); \
    t2 = _mm_hadd_pd(t2,t2); \
    t1= _mm_shuffle_pd(t1,t2,0x2);		\
    (z) = _mm_add_pd((z),t1);			\
  }

/* z += conj(x)*conj(y)    z, x, y are SSE registers containing complex numbers
              ordered with the real part in the low half, imag part 
	      in the upper half */
#define CCMADD(z,x,y)				\
  { \
    __m128d t1,t2;	      \
    __m128d t3 = _mm_set_pd( (double)(-1), (double)1 ); \
    t1 = _mm_mul_pd((x),(y)); \
    t1 = _mm_hsub_pd(t1,t1); \
    t2 = _mm_shuffle_pd((y),(y),0x1);\
    t2 = _mm_mul_pd((x),t2); \
    t2 = _mm_hadd_pd(t2,t2); \
    t1= _mm_shuffle_pd(t1,t2,0x2);		\
    t1= _mm_mul_pd(t3,t1); \
    (z) = _mm_add_pd((z),t1);			\
  }

/* z = x*conj(y)    z, x, y are SSE registers containing complex numbers
   ordered with the real part in the low half, imag part 
   in the upper half */
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

/* z += x*conj(y)    z, x, y are SSE registers containing complex numbers
   ordered with the real part in the low half, imag part 
   in the upper half */
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
#endif
