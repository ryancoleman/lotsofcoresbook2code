#ifndef _INLINE_SSE_H_
#define _INLINE_SSE_H_


#include <xmmintrin.h>

#ifdef __cplusplus
extern "C" { 
#endif

typedef struct {
   float real;		   
   float imag;
} complexf;

typedef struct { complexf e[3][3]; } su3_matrixf ;
typedef struct { complexf c[3]; } su3_vectorf;
typedef struct { su3_vectorf d[4]; } wilson_vectorf;
typedef struct { su3_vectorf h[2]; } half_wilson_vectorf;

typedef union
{
  unsigned int c1[4];
  __m128 xmm;
} sse_mask  __attribute__ ((aligned (16)));

static sse_mask _sse_sgn13  __attribute__ ((unused)) ={{0x80000000, 0x00000000, 0x80000000, 0x00000000}};
static sse_mask _sse_sgn24  __attribute__ ((unused)) ={{0x00000000, 0x80000000, 0x00000000, 0x80000000}};
static sse_mask _sse_sgn3   __attribute__ ((unused)) ={{0x00000000, 0x00000000, 0x80000000, 0x00000000}};
static sse_mask _sse_sgn4   __attribute__ ((unused)) ={{0x00000000, 0x00000000, 0x00000000, 0x80000000}};


#ifdef __cplusplus
};
#endif

#endif
