/*
    The MIT License (MIT)
    
    Copyright (c) 2015 OpenVec
    
    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
    
    
    Authors:
    Paulo Souza
    Leonardo Borges
    Cedric Andreolli
    Philippe Thierry

*/
/*
		DO NOT CALL DIRECTLY ANY FUNCTION WITH PREFIX ov_prv
		THESE FUNCTIONS ARE "PRIVATE" FOR INTERNAL USE ONLY.
*/
#ifndef _OPENVEC_
#define _OPENVEC_

/*******************************************************************************
                                C++ VS C 
*******************************************************************************/
#ifdef __cplusplus
#define OV_LANG "C++"
#define OV_INLINE static inline
#define OV_REF(a) &a
#define ov_restrict __restrict__

#else
#define OV_LANG "C"
#define OV_INLINE static __inline__
#define OV_REF(a) a

#ifdef __STDC_VERSION__
#if __STDC_VERSION__ >= 199900L
#define ov_restrict restrict
#endif
#endif

#endif


/*******************************************************************************
                          ARCHITECTURE AUTO-DETECTION
*******************************************************************************/
#ifndef _OV_NOAUTO
#if defined(__AVX512F__)
#define _OV_AVX512
#elif defined(__MIC__)
#define _OV_MIC
#elif defined(__AVX2__)
#define _OV_AVX2
#elif defined(__AVX__)
#define _OV_AVX
#elif defined(__SSE4_1__)
#define _OV_SSE4
#elif defined(__SSE2__)
#define _OV_SSE2
#elif defined(__SSE__)
#define _OV_SSE
#elif defined(__ARM_NEON__)
#define _OV_NEON
#elif defined(__aarch64__)
#define _OV_NEON
#elif defined(__VSX__)
#define _OV_ALTIVEC
#endif
#endif

#define OV_EPS 1e-23f


#define OV_SIZEOF_FLOAT 4

#define OV_ONES 0xffffffff

#define ov_sqrf(x) ov_mulf(x,x)
#define ov_prv_sqrf(x) ov_prv_mulf(x,x)



/*******************************************************************************
          Intrinsics for Intel(R)Many Integrated Core Architecture
*******************************************************************************/
#if defined _OV_MIC
#include"openvec.gccavx512.h"


/*******************************************************************************
          Intrinsics for Intel(R)Many Integrated Core Architecture
*******************************************************************************/
#elif defined _OV_AVX512
#include"openvec.gccavx512.h"


/*******************************************************************************
               Intrinsics for Intel(R) Advanced Vector Extensions
*******************************************************************************/
#elif defined _OV_AVX
#include"openvec.gccavx.h"

/*******************************************************************************
               Intrinsics for Intel(R) Advanced Vector Extensions 2
*******************************************************************************/
#elif defined _OV_AVX2
#include"openvec.gccavx.h"


/*******************************************************************************
         Intrinsics for Intel(R)Streaming SIMD Extensions (Intel(R)SSE)
*******************************************************************************/
#elif defined _OV_SSE
#include"openvec.gccsse.h"
#elif defined _OV_SSE2
#include"openvec.gccsse.h"
#elif defined _OV_SSE4
#include"openvec.gccsse.h"


/*******************************************************************************
                        Intrinsics for ARM NEON
*******************************************************************************/
#elif defined _OV_NEON
#define OV_CPP_OVERRIDE_TYPE
#include"openvec.gccneon.h"


/*******************************************************************************
                        Intrinsics for IBM Power 7/8
*******************************************************************************/
#elif defined _OV_ALTIVEC
#define OV_CPP_OVERRIDE_TYPE
#include"openvec.vsx.h"


/*******************************************************************************
                        Portable Scalar Intrinsics
*******************************************************************************/
#else
#include"openvec.scalar.h"
#endif



/*******************************************************************************
                            C++ Type Override 
*******************************************************************************/
#if defined(__cplusplus) && defined(OV_CPP_OVERRIDE_TYPE) && (OV_FLOAT_WIDTH>1)
struct __attribute__((aligned (OV_ALIGN))) ov_float 
{  
  ov_prv_float vreg;
};

#define OV_VREG(a) a.vreg
OV_INLINE ov_prv_float ov_prv_vreg(ov_float const OV_REF(a))
{
  return OV_VREG(a);
}

#else
#define ov_float ov_prv_float
#define OV_VREG(a) a
#define ov_prv_vreg(a) a
#endif



/*******************************************************************************
                             DOUBLE PRECISION
*******************************************************************************/
#ifdef ov_prv_double
#include"openvec.double.h"
#endif


/*******************************************************************************
                                   DEFAULTS
*******************************************************************************/








#ifndef ov_prefetch1
#define ov_prefetch1(a)
#endif

#ifndef ov_prefetch2
#define ov_prefetch2(a)
#endif

#ifndef ov_prefetch3
#define ov_prefetch3(a)
#endif

#ifndef ov_prefetch4
#define ov_prefetch4(a)
#endif

#ifndef ov_prefetch5
#define ov_prefetch5(a)
#endif



#if OV_FLOAT_WIDTH > 1 /* { */

#define ov_zerof ov_getzerof()

#ifdef __cplusplus /* { */
#ifdef OV_CPP_OVERRIDE_TYPE /* { */

static inline void ov_stf(float *addr, ov_float const OV_REF(a))
{
  ov_prv_stf(addr, OV_VREG(a));
}

static inline void ov_storef(float *addr, ov_float const OV_REF(a))
{
  ov_prv_storef(addr, OV_VREG(a));
}
static inline void ov_storeuf(float *addr, ov_float const OV_REF(a))
{
  ov_prv_storeuf(addr, OV_VREG(a));
}
static inline ov_float ov_addf(ov_float const OV_REF(a), ov_float const OV_REF(b))
{
  ov_float tmp;
  OV_VREG(tmp) = ov_prv_addf(OV_VREG(a), OV_VREG(b));
  return tmp;
}
static inline ov_float ov_subf(ov_float const OV_REF(a), ov_float const OV_REF(b))
{
  ov_float tmp;
  OV_VREG(tmp) = ov_prv_subf(OV_VREG(a), OV_VREG(b));
  return tmp;
}
static inline ov_float ov_mulf(ov_float const OV_REF(a), ov_float const OV_REF(b))
{
  ov_float tmp;
  OV_VREG(tmp) = ov_prv_mulf(OV_VREG(a), OV_VREG(b));
  return tmp;
}
static inline ov_float ov_divf(ov_float const OV_REF(a), ov_float const OV_REF(b))
{
  ov_float tmp;
  OV_VREG(tmp) = ov_prv_divf(OV_VREG(a), OV_VREG(b));
  return tmp;
}
static inline ov_float ov_getzerof()
{
  ov_float tmp;
  ov_prv_setzerof(OV_VREG(tmp));
  return tmp;
}
static inline void ov_setzerof(ov_float OV_REF(dest))
{
  ov_prv_setzerof(OV_VREG(dest));
}

static inline ov_float ov_uldf(float const *addr)
{
  ov_float tmp;
  ov_prv_uloadf(OV_VREG(tmp), addr);
  return tmp;
}
static inline ov_float ov_ldf(float const *b)
{
  ov_float tmp;
  OV_VREG(tmp) = ov_prv_ldf(b);
  return tmp;
}
static inline ov_float ov_loadf(float const *b)
{
  ov_float tmp;
  OV_VREG(tmp) = ov_prv_loadf(b);
  return tmp;
}
static inline void ov_uloadf(ov_float OV_REF(dest), float const *addr)
{
  ov_prv_uloadf(OV_VREG(dest), addr);
}

static inline void ov_stream_stf(float *addr, ov_float const OV_REF(a))
{
  ov_prv_stream_stf(addr, OV_VREG(a));
}

static inline void ov_nt_storef(float *addr, ov_float const OV_REF(a))
{
  ov_prv_nt_storef(addr, OV_VREG(a));
}

static inline ov_float ov_setf(float const OV_REF(a))
{
  ov_float tmp;
  OV_VREG(tmp) = ov_prv_setf(a);
  return tmp;
}


static inline ov_float ov_rsqrtf(ov_float const OV_REF(a))
{
  ov_float tmp;
  OV_VREG(tmp) = ov_prv_rsqrtf(OV_VREG(a));
  return tmp;
}

static inline ov_float ov_floorf(ov_float const OV_REF(a))
{
  ov_float tmp;
  OV_VREG(tmp) = ov_prv_floorf(OV_VREG(a));
  return tmp;
}

static inline ov_float ov_ceilf(ov_float const OV_REF(a))
{
  ov_float tmp;
  OV_VREG(tmp) = ov_prv_ceilf(OV_VREG(a));
  return tmp;
}

static inline ov_float ov_absf(ov_float const OV_REF(a))
{
  ov_float tmp;
  OV_VREG(tmp) = ov_prv_absf(OV_VREG(a));
  return tmp;
}

static inline ov_float ov_rcpf(ov_float const OV_REF(a))
{
  ov_float tmp;
  OV_VREG(tmp) = ov_prv_rcpf(OV_VREG(a));
  return tmp;
}

static inline ov_float ov_sqrtf(ov_float const OV_REF(a))
{
  ov_float tmp;
  OV_VREG(tmp) = ov_prv_sqrtf(OV_VREG(a));
  return tmp;
}


static inline ov_float ov_conditionalf(ov_maskf const OV_REF(mask), ov_float const OV_REF(a), ov_float const OV_REF(b))
{
  ov_float tmp;
  OV_VREG(tmp) = ov_prv_conditionalf(mask, OV_VREG(a), OV_VREG(b));
  return tmp;
}



static inline ov_float ov_maddf(ov_float const OV_REF(a), ov_float const OV_REF(b), ov_float const OV_REF(c))
{
  ov_float tmp;
  OV_VREG(tmp) = ov_prv_maddf(OV_VREG(a), OV_VREG(b), OV_VREG(c));
  return tmp;
}


static inline ov_float ov_msubf(ov_float const OV_REF(a), ov_float const OV_REF(b), ov_float const OV_REF(c))
{
  ov_float tmp;
  OV_VREG(tmp) = ov_prv_msubf(OV_VREG(a), OV_VREG(b), OV_VREG(c));
  return tmp;
}


static inline ov_maskf ov_eqf(ov_float const OV_REF(a), ov_float const OV_REF(b))
{
  return ov_prv_eqf(OV_VREG(a), OV_VREG(b));
}
static inline ov_maskf ov_gef(ov_float const OV_REF(a), ov_float const OV_REF(b))
{
  return ov_prv_gef(OV_VREG(a), OV_VREG(b));
}
static inline ov_maskf ov_gtf(ov_float const OV_REF(a), ov_float const OV_REF(b))
{
  return ov_prv_gtf(OV_VREG(a), OV_VREG(b));
}
static inline ov_maskf ov_lef(ov_float const OV_REF(a), ov_float const OV_REF(b))
{
  return ov_prv_lef(OV_VREG(a), OV_VREG(b));
}
static inline ov_maskf ov_ltf(ov_float const OV_REF(a), ov_float const OV_REF(b))
{
  return ov_prv_ltf(OV_VREG(a), OV_VREG(b));
}
/*
#define ov_gef ov_prv_gef
#define ov_lef ov_prv_lef
#define ov_gtf ov_prv_gtf
#define ov_ltf ov_prv_ltf
*/

#else /* }{ */

#define ov_addf ov_prv_addf
#define ov_subf ov_prv_subf
#define ov_mulf ov_prv_mulf
#define ov_divf ov_prv_divf

#define ov_setf ov_prv_setf
#define ov_sqrtf ov_prv_sqrtf
#define ov_rsqrtf ov_prv_rsqrtf
#define ov_rcpf ov_prv_rcpf
#define ov_storef ov_prv_storef
#define ov_storeuf ov_prv_storeuf
#define ov_stream_stf ov_prv_stream_stf
#define ov_nt_storef ov_prv_nt_storef


#define ov_msubf ov_prv_msubf
#define ov_maddf ov_prv_maddf
#define ov_conditionalf ov_prv_conditionalf

#define ov_floorf ov_prv_floorf
#define ov_ceilf ov_prv_ceilf
#define ov_absf ov_prv_absf
#define ov_getzerof ov_prv_getzerof
#define ov_setzerof ov_prv_setzerof

#define ov_loadf ov_prv_loadf
#define ov_ldf ov_prv_ldf
#define ov_uldf ov_prv_uldf
#define ov_stf ov_prv_stf

/* Masked comparisom */
#define ov_gef ov_prv_gef
#define ov_lef ov_prv_lef
#define ov_gtf ov_prv_gtf
#define ov_ltf ov_prv_ltf
#define ov_eqf ov_prv_eqf

#endif /*}*/

/* C++ common code */

static inline ov_float sqrtf(ov_float const OV_REF(a))
{
  ov_float tmp;
  OV_VREG(tmp) = ov_prv_sqrtf(OV_VREG(a));
  return tmp;
}

static inline ov_float fabsf(ov_float const OV_REF(a))
{
  ov_float tmp;
  OV_VREG(tmp) = ov_prv_absf(OV_VREG(a));
  return tmp;
}
static inline ov_float floorf(ov_float const OV_REF(a))
{
  ov_float tmp;
  OV_VREG(tmp) = ov_prv_floorf(OV_VREG(a));
  return tmp;
}
static inline ov_float ceilf(ov_float const OV_REF(a))
{
  ov_float tmp;
  OV_VREG(tmp) = ov_prv_ceilf(OV_VREG(a));
  return tmp;
}

/* MAXIMUM */
static inline ov_float ov_maxf(ov_float const OV_REF(a), ov_float const OV_REF(b))
{
  ov_float tmp;
  OV_VREG(tmp) = ov_prv_maxf(OV_VREG(a), OV_VREG(b));
  return tmp;
}
static inline ov_float ov_maxf(ov_float const OV_REF(a), float const OV_REF(b))
{
  ov_float tmp;
  OV_VREG(tmp) = ov_prv_maxf(OV_VREG(a), ov_prv_setf(b));
  return tmp;
}
static inline ov_float ov_maxf(float const OV_REF(a), ov_float const OV_REF(b))
{
  ov_float tmp;
  OV_VREG(tmp) = ov_prv_maxf(ov_prv_setf(a), OV_VREG(b));
  return tmp;
}

/* MINIMUM */
static inline ov_float ov_minf(ov_float const OV_REF(a), ov_float const OV_REF(b))
{
  ov_float tmp;
  OV_VREG(tmp) = ov_prv_minf(OV_VREG(a), OV_VREG(b));
  return tmp;
}
static inline ov_float ov_minf(ov_float const OV_REF(a), float const OV_REF(b))
{
  ov_float tmp;
  OV_VREG(tmp) = ov_prv_minf(OV_VREG(a), ov_prv_setf(b));
  return tmp;
}
static inline ov_float ov_minf(float const OV_REF(a), ov_float const OV_REF(b))
{
  ov_float tmp;
  OV_VREG(tmp) = ov_prv_minf(ov_prv_setf(a), OV_VREG(b));
  return tmp;
}
/* END MINIMUM */


/*
   C++ operators
*/
/* Conditional vector vector */
static inline ov_maskf operator==(ov_float const &a, ov_float const &b)
{
   return ov_prv_eqf(OV_VREG(a),OV_VREG(b));
}
static inline ov_maskf operator<(ov_float const &a, ov_float const &b)
{
   return ov_prv_ltf(OV_VREG(a),OV_VREG(b));
}
static inline ov_maskf operator>(ov_float const &a, ov_float const &b)
{
   return ov_prv_gtf(OV_VREG(a),OV_VREG(b));
}
static inline ov_maskf operator<=(ov_float const &a, ov_float const &b)
{
   return ov_prv_lef(OV_VREG(a),OV_VREG(b));
}
static inline ov_maskf operator>=(ov_float const &a, ov_float const &b)
{
   return ov_prv_gef(OV_VREG(a),OV_VREG(b));
}



/* Conditional float vector */
static inline ov_maskf operator==(float const &a, ov_float const &b)
{
   return ov_prv_eqf(ov_prv_setf(a),OV_VREG(b));
}
static inline ov_maskf operator<(float const &a, ov_float const &b)
{
   return ov_prv_ltf(ov_prv_setf(a),OV_VREG(b));
}
static inline ov_maskf operator>(float const &a, ov_float const &b)
{
   return ov_prv_gtf(ov_prv_setf(a),OV_VREG(b));
}
static inline ov_maskf operator<=(float const &a, ov_float const &b)
{
   return ov_prv_lef(ov_prv_setf(a),OV_VREG(b));
}
static inline ov_maskf operator>=(float const &a, ov_float const &b)
{
   return ov_prv_gef(ov_prv_setf(a),OV_VREG(b));
}



/* Conditional vector float */
static inline ov_maskf operator==(ov_float const &a, float const &b)
{
   return ov_prv_eqf(OV_VREG(a),ov_prv_setf(b));
}
static inline ov_maskf operator<(ov_float const &a, float const &b)
{
   return ov_prv_ltf(OV_VREG(a),ov_prv_setf(b));
}
static inline ov_maskf operator>(ov_float const &a, float const &b)
{
   return ov_prv_gtf(OV_VREG(a),ov_prv_setf(b));
}
static inline ov_maskf operator<=(ov_float const &a, float const &b)
{
   return ov_prv_lef(OV_VREG(a),ov_prv_setf(b));
}
static inline ov_maskf operator>=(ov_float const &a, float const &b)
{
   return ov_prv_gef(OV_VREG(a),ov_prv_setf(b));
}


/* Add */
static inline ov_float operator+(ov_float const &a, ov_float const &b)
{
   ov_float tmp;
   OV_VREG(tmp) = ov_prv_addf(OV_VREG(a),OV_VREG(b));
   return tmp;
}
static inline ov_float operator+(ov_float const &a, float const &b)
{
   ov_float tmp;
   OV_VREG(tmp) = ov_prv_addf(ov_prv_setf(b),OV_VREG(a));
   return tmp;
}
static inline ov_float operator+(float const &a, ov_float const &b)
{
   ov_float tmp;
   OV_VREG(tmp) = ov_prv_addf(ov_prv_setf(a),OV_VREG(b));
   return tmp;
}


/* add += */
void operator+=(ov_float &a, ov_float const b)
{
   OV_VREG(a) = ov_prv_addf(OV_VREG(a), OV_VREG(b));
}
void operator+=(ov_float &a, float const b)
{
   OV_VREG(a) = ov_prv_addf(OV_VREG(a), ov_prv_setf(b));
}


/* subtraction */
static inline ov_float operator-(ov_float const &a, ov_float const &b)
{
   ov_float tmp;
   OV_VREG(tmp) = ov_prv_subf(OV_VREG(a),OV_VREG(b));
   return tmp;
}
static inline ov_float operator-(ov_float const &a, float const &b)
{
   ov_float tmp;
   OV_VREG(tmp) = ov_prv_subf(OV_VREG(a),ov_prv_setf(b));
   return tmp;
}
static inline ov_float operator-(float const &a, ov_float const &b)
{
   ov_float tmp;
   OV_VREG(tmp) = ov_prv_subf(ov_prv_setf(a), OV_VREG(b));
   return tmp;
}


/* sub -= */
void operator-=(ov_float &a, ov_float const b)
{
   OV_VREG(a) = ov_prv_subf(OV_VREG(a), OV_VREG(b));
}
void operator-=(ov_float &a, float const b)
{
   OV_VREG(a) = ov_prv_subf(OV_VREG(a), ov_prv_setf(b));
}


/* multiplication */
static inline ov_float operator*(ov_float const &a, ov_float const &b)
{
   ov_float tmp;
   OV_VREG(tmp) = ov_prv_mulf(OV_VREG(a),OV_VREG(b));
   return tmp;
}
static inline ov_float operator*(ov_float const &a, float const &b)
{
   ov_float tmp;
   OV_VREG(tmp) = ov_prv_mulf(OV_VREG(a),ov_prv_setf(b));
   return tmp;
}
static inline ov_float operator*(float const &a, ov_float const &b)
{
   ov_float tmp;
   OV_VREG(tmp) = ov_prv_mulf(ov_prv_setf(a),OV_VREG(b));
   return tmp;
}


/* Division */
static inline ov_float operator/(ov_float const &a, ov_float const &b)
{
   ov_float tmp;
   OV_VREG(tmp) = ov_prv_divf(OV_VREG(a),OV_VREG(b));
   return tmp;
}
static inline ov_float operator/(ov_float const &a, float const &b)
{
   ov_float tmp;
   OV_VREG(tmp) = ov_prv_divf(OV_VREG(a),ov_prv_setf(b));
   return tmp;
}
static inline ov_float operator/(float const &a, ov_float const &b)
{
   ov_float tmp;
   OV_VREG(tmp) = ov_prv_divf(ov_prv_setf(a), OV_VREG(b));
   return tmp;
}


#else /* }{ NOTcplusplus */

/* C zero overhead one to one map */

#define ov_addf ov_prv_addf
#define ov_subf ov_prv_subf
#define ov_mulf ov_prv_mulf
#define ov_divf ov_prv_divf
#define ov_maddf ov_prv_maddf
#define ov_conditionalf ov_prv_conditionalf
#define ov_getzerof ov_prv_getzerof
#define ov_uloadf ov_prv_uloadf
#define ov_setf ov_prv_setf
#define ov_loadf ov_prv_loadf
#define ov_storef ov_prv_storef
#define ov_storeuf ov_prv_storeuf
#define ov_nt_storef ov_prv_nt_storef
#define ov_msubf ov_prv_msubf
#define ov_stream_stf ov_prv_stream_stf
#define ov_setzerof ov_prv_setzerof

#define ov_rsqrtf ov_prv_rsqrtf
#define ov_rcpf ov_prv_rcpf

#define ov_sqrtf ov_prv_sqrtf

#define ov_ldf ov_prv_ldf
#define ov_uldf ov_prv_uldf
#define ov_stf ov_prv_stf
#define ov_maxf ov_prv_maxf
#define ov_minf ov_prv_minf
#define ov_floorf ov_prv_floorf
#define ov_ceilf ov_prv_ceilf
#define ov_absf ov_prv_absf

/* Masked comparisom */
#define ov_gef ov_prv_gef
#define ov_lef ov_prv_lef
#define ov_gtf ov_prv_gtf
#define ov_ltf ov_prv_ltf
#define ov_eqf ov_prv_eqf

#endif /* }} */


int posix_memalign(void **memptr, size_t alignment, size_t size);


static void *ov_malloc(size_t size)
{
  void *memptr;

  /*printf("Allocating %ld MB\n", (long)size/(1024*1024));*/

  int ierr = posix_memalign(&memptr, (size_t)OV_ALIGN, size+OV_ALIGN);

  if (ierr) return NULL;
  else      return memptr;

/*
  char *pt, *pt_ini;
  pt_ini=(char*)malloc(size+2*OV_ALIGN+sizeof(char*));
  pt=pt_ini + sizeof(char*);
  size_t addr=(size_t)pt;
  while(addr%OV_ALIGN)
  {
    pt++;
    addr=(size_t)pt;
  }

  void **pini=(void**)pt;
  pini--;
  *pini=pt_ini;
  return (void*)pt;
*/
}

static void *ov_calloc(size_t nmemb, size_t size)
{
  size_t const sizebytes = nmemb*size;
  size_t i;

  ov_float *array;
  array = (ov_float*)ov_malloc(sizebytes);

  if (array != NULL)
  {
    /*printf("Making zeroooooo %ld MB\n", (long)sizebytes/(1024*1024));*/
    size_t const n = sizebytes/OV_ALIGN + 1;
    for(i=0; i<n; i++) array[i]=ov_zerof;
  }

  return (void*)array;
}

  

/*
static void ov_free(void *ptr)
{
  free(ptr);
*/
/*
  void **pini=(void**)ptr;
  pini--;
  void *pt=*pini;
  free(pt);
*/
/*
}
*/
#define malloc ov_malloc
#define calloc ov_calloc
#define ov_free free
/* #define free ov_free */

OV_INLINE float ov_all_sumf(ov_float const OV_REF(a))
{
  float const *fp = (float*)&a;
  int i;

  float result=fp[0];

  for (i=1; i<OV_FLOAT_WIDTH; i++) result += fp[i];
  return result;
}
OV_INLINE float ov_all_prodf(ov_float const OV_REF(a))
{
  float const *fp = (float*)&a;
  int i;

  float result=fp[0];

  for (i=1; i<OV_FLOAT_WIDTH; i++) result *= fp[i];
  return result;
}
OV_INLINE float ov_all_maxf(ov_float const OV_REF(a))
{
  float const *fp = (float*)&a;
  int i;

  float result=fp[0];

  for (i=1; i<OV_FLOAT_WIDTH; i++) if (fp[i]>result) result=fp[i];
  return result;
}
OV_INLINE float ov_all_minf(ov_float const OV_REF(a))
{
  float const *fp = (float*)&a;
  int i;

  float result=fp[0];

  for (i=1; i<OV_FLOAT_WIDTH; i++) if (fp[i]<result) result=fp[i];
  return result;
}

#else /* OV_FLOAT_WIDTH > 1 */

#define ov_all_sumf(x) (x)
#define ov_all_prodf(x) (x)
#define ov_all_minf(x) (x)
#define ov_all_maxf(x) (x)

#endif /* } OV_FLOAT_WIDTH > 1 */

#ifndef ov_fastsqrtf
#define ov_fastsqrtf(x) ov_mulf(x, ov_rsqrtf(ov_maxf(x, ov_setf(OV_EPS))))
#endif

#define ov_ustf ov_storeuf

#define OV_FLOAT_TAIL ((OV_FLOAT_WIDTH)-1)





#ifndef ov_restrict
#define ov_restrict
#endif

#endif /* _OPENVEC_ */
