#include "qdp.h"
#include "testCMulDouble.h"
#include "unittest.h"
#include <cmath>

using namespace QDP;
using namespace Assertions;




#include<xmmintrin.h>


typedef union {
  __m128d v;
  double  d[2];
} VD;

#include "qdp_config.h"
#ifndef QDP_USE_SSE3
// SSE2 

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

#define CCMUL(z,x,y)		\
  { \
    __m128d t1,t2,t3; \
    __m128d t4 = _mm_set_pd((double)(-1),(double)1);	\
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
    __m128d t5 = _mm_set_pd( (double)(-1),(double)1 );	\
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
#include <pmmintrin.h>

/* SSE 3? */
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

#define CCMUL(z,x,y)		\
  { \
    __m128d t1; \
    __m128d t2 = _mm_set_pd( (double)(-1),(double)1 ) ;	\
    t1 = _mm_mul_pd((x),(y)); \
    (z) = _mm_hsub_pd(t1,t1);			\
    t1 = _mm_shuffle_pd((y),(y),0x1);\
    t1 = _mm_mul_pd((x),t1); \
    t1 = _mm_hadd_pd(t1,t1); \
    (z)= _mm_shuffle_pd((z),t1,0x2);		\
    (z)= _mm_mul_pd((z),t2); \
  }

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

#define CCMADD(z,x,y)				\
  { \
    __m128d t1,t2;	      \
    __m128d t3= _mm_set_pd ( (double)(-1), (double)1 );	\
    t1 = _mm_mul_pd((x),(y)); \
    t1 = _mm_hsub_pd(t1,t1); \
    t2 = _mm_shuffle_pd((y),(y),0x1);\
    t2 = _mm_mul_pd((x),t2); \
    t2 = _mm_hadd_pd(t2,t2); \
    t1= _mm_shuffle_pd(t1,t2,0x2);		\
    t1= _mm_mul_pd(t3,t1); \
    (z) = _mm_add_pd((z),t1);			\
  }

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

void
testCMul::run()
{
  DComplex x=cmplx(Double(1.5),Double(4.0));
  DComplex y=cmplx(Double(2.5),Double(3.0));
  DComplex z1 = x*y;

  VD x_v;
  VD y_v;

  x_v.d[0] = toDouble(real(x));
  x_v.d[1] = toDouble(imag(x));

  y_v.d[0] = toDouble(real(y));
  y_v.d[1] = toDouble(imag(y));

  VD z_v;

  CMUL(z_v.v, x_v.v, y_v.v);

  DComplex z2=cmplx(Real(z_v.d[0]), Real(z_v.d[1]));

  double realdiff = fabs(toDouble( real(z2-z1) ));
  QDPIO::cout << std::endl << "Real diff = " << realdiff << std::endl;

  double imagdiff = fabs(toDouble( imag(z2-z1) ));
  QDPIO::cout << "Imag diff = " << imagdiff << std::endl;

  assertion( realdiff < 1.0e-14) ;
  assertion( imagdiff < 1.0e-14) ;
}

void
testCMadd::run()
{
  DComplex x=cmplx(Double(1.5),Double(4.0));
  DComplex y=cmplx(Double(2.5),Double(3.0));
  DComplex z1 = cmplx(Double(4), Double(-0.5));

 
  VD z_v;
  VD x_v;
  VD y_v;

  z_v.d[0] = toDouble(real(z1));
  z_v.d[1] = toDouble(imag(z1));

  z1 += x*y;

  x_v.d[0] = toDouble(real(x));
  x_v.d[1] = toDouble(imag(x));

  y_v.d[0] = toDouble(real(y));
  y_v.d[1] = toDouble(imag(y));


  CMADD(z_v.v, x_v.v, y_v.v);

  DComplex z2=cmplx(Real(z_v.d[0]), Real(z_v.d[1]));

  double realdiff = fabs(toDouble( real(z2-z1) ));
  QDPIO::cout << std::endl << "Real diff = " << realdiff << std::endl;

  double imagdiff = fabs(toDouble( imag(z2-z1) ));
  QDPIO::cout << "Imag diff = " << imagdiff << std::endl;

  assertion( realdiff < 1.0e-14) ;
  assertion( imagdiff < 1.0e-14) ;
}

void
testConjMul::run()
{
  DComplex x=cmplx(Double(1.5),Double(4.0));
  DComplex y=cmplx(Double(2.5),Double(3.0));
  DComplex z1 = x*conj(y);

  VD x_v;
  VD y_v;

  x_v.d[0] = toDouble(real(x));
  x_v.d[1] = toDouble(imag(x));

  y_v.d[0] = toDouble(real(y));
  y_v.d[1] = toDouble(imag(y));

  VD z_v;

  CONJMUL(z_v.v, x_v.v, y_v.v);

  DComplex z2=cmplx(Real(z_v.d[0]), Real(z_v.d[1]));

  double realdiff = fabs(toDouble( real(z2-z1) ));
  QDPIO::cout << std::endl << "Real diff = " << realdiff << std::endl;

  double imagdiff = fabs(toDouble( imag(z2-z1) ));
  QDPIO::cout << "Imag diff = " << imagdiff << std::endl;

  assertion( realdiff < 1.0e-14) ;
  assertion( imagdiff < 1.0e-14) ;
}

void
testConjMadd::run()
{
  DComplex x=cmplx(Double(1.5),Double(4.0));
  DComplex y=cmplx(Double(2.5),Double(3.0));
  DComplex z1 = cmplx(Double(4), Double(-0.5));

 
  VD z_v;
  VD x_v;
  VD y_v;

  z_v.d[0] = toDouble(real(z1));
  z_v.d[1] = toDouble(imag(z1));

  z1 += x*conj(y);

  x_v.d[0] = toDouble(real(x));
  x_v.d[1] = toDouble(imag(x));

  y_v.d[0] = toDouble(real(y));
  y_v.d[1] = toDouble(imag(y));


  CONJMADD(z_v.v, x_v.v, y_v.v);

  DComplex z2=cmplx(Real(z_v.d[0]), Real(z_v.d[1]));

  double realdiff = fabs(toDouble( real(z2-z1) ));
  QDPIO::cout << std::endl << "Real diff = " << realdiff << std::endl;

  double imagdiff = fabs(toDouble( imag(z2-z1) ));
  QDPIO::cout << "Imag diff = " << imagdiff << std::endl;

  assertion( realdiff < 1.0e-14) ;
  assertion( imagdiff < 1.0e-14) ;
}

void
testCCMul::run()
{
  DComplex x=cmplx(Double(1.5),Double(4.0));
  DComplex y=cmplx(Double(2.5),Double(3.0));
  DComplex z1 = conj(x)*conj(y);

  VD x_v;
  VD y_v;

  x_v.d[0] = toDouble(real(x));
  x_v.d[1] = toDouble(imag(x));

  y_v.d[0] = toDouble(real(y));
  y_v.d[1] = toDouble(imag(y));

  VD z_v;

  CCMUL(z_v.v, x_v.v, y_v.v);

  DComplex z2=cmplx(Real(z_v.d[0]), Real(z_v.d[1]));
  double realdiff = fabs(toDouble( real(z2-z1) ));
  QDPIO::cout << std::endl << "Real diff = " << realdiff << std::endl;

  double imagdiff = fabs(toDouble( imag(z2-z1) ));
  QDPIO::cout << "Imag diff = " << imagdiff << std::endl;

  assertion( realdiff < 1.0e-14) ;
  assertion( imagdiff < 1.0e-14) ;
}

void
testCCMadd::run()
{
  DComplex x=cmplx(Double(1.5),Double(4.0));
  DComplex y=cmplx(Double(2.5),Double(3.0));
  DComplex z1 = cmplx(Double(4), Double(-0.5));

 
  VD z_v;
  VD x_v;
  VD y_v;

  z_v.d[0] = toDouble(real(z1));
  z_v.d[1] = toDouble(imag(z1));

  z1 += conj(x)*conj(y);

  x_v.d[0] = toDouble(real(x));
  x_v.d[1] = toDouble(imag(x));

  y_v.d[0] = toDouble(real(y));
  y_v.d[1] = toDouble(imag(y));


  CCMADD(z_v.v, x_v.v, y_v.v);

  DComplex z2=cmplx(Real(z_v.d[0]), Real(z_v.d[1]));

  double realdiff = fabs(toDouble( real(z2-z1) ));
  QDPIO::cout << std::endl << "Real diff = " << realdiff << std::endl;

  double imagdiff = fabs(toDouble( imag(z2-z1) ));
  QDPIO::cout << "Imag diff = " << imagdiff << std::endl;

  assertion( realdiff < 1.0e-14) ;
  assertion( imagdiff < 1.0e-14) ;
}
