// $Id: jlab_sse.h,v 1.3 2003-08-05 21:15:02 edwards Exp $

#ifndef JLAB_SSE_H
#define JLAB_SSE_H

#if 0
typedef struct
{
   unsigned int c1,c2,c3,c4;
} sse_mask __attribute__ ((aligned (16)));

static sse_mask _sse_sgn13 __attribute__ ((unused)) ={0x80000000, 0x00000000, 0x80000000, 0x00000000};
static sse_mask _sse_sgn24 __attribute__ ((unused)) ={0x00000000, 0x80000000, 0x00000000, 0x80000000};
static sse_mask _sse_sgn3 __attribute__  ((unused)) ={0x00000000, 0x00000000, 0x80000000, 0x00000000};
static sse_mask _sse_sgn4 __attribute__  ((unused)) ={0x00000000, 0x00000000, 0x00000000, 0x80000000};
#endif


typedef struct
{
   float re,im;
} complex;

typedef struct 
{
   complex c11,c12,c13,c21,c22,c23,c31,c32,c33;
} su3;

#include "sse32.h"
#include "sse_align.h"

#define _c11 .c11
#define _c12 .c12
#define _c13 .c13
#define _c21 .c21
#define _c22 .c22
#define _c23 .c23
#define _c31 .c31
#define _c32 .c32
#define _c33 .c33

#define _c11re .c11.re
#define _c12re .c12.re
#define _c13re .c13.re
#define _c21re .c21.re
#define _c22re .c22.re
#define _c23re .c23.re
#define _c31re .c31.re
#define _c32re .c32.re
#define _c33re .c33.re

#define _c11im .c11.im
#define _c12im .c12.im
#define _c13im .c13.im
#define _c21im .c21.im
#define _c22im .c22.im
#define _c23im .c23.im
#define _c31im .c31.im
#define _c32im .c32.im
#define _c33im .c33.im



#define _sse_pair_store_up_c1_c2(r) \
__asm__ __volatile__ ("movlps %%xmm3, %0 \n\t" \
                      "movlps %%xmm4, %1 \n\t" \
                      "movlps %%xmm5, %2 \n\t" \
                      "movhps %%xmm3, %3 \n\t" \
                      "movhps %%xmm4, %4 \n\t" \
                      "movhps %%xmm5, %5" \
                      : \
                      "=m" ((r)_c11), \
                      "=m" ((r)_c21), \
                      "=m" ((r)_c31), \
                      "=m" ((r)_c12), \
                      "=m" ((r)_c22), \
                      "=m" ((r)_c32))

#define _sse_pair_store_up_c2_c3(r) \
__asm__ __volatile__ ("movlps %%xmm3, %0 \n\t" \
                      "movlps %%xmm4, %1 \n\t" \
                      "movlps %%xmm5, %2 \n\t" \
                      "movhps %%xmm3, %3 \n\t" \
                      "movhps %%xmm4, %4 \n\t" \
                      "movhps %%xmm5, %5" \
                      : \
                      "=m" ((r)_c12), \
                      "=m" ((r)_c22), \
                      "=m" ((r)_c32), \
                      "=m" ((r)_c13), \
                      "=m" ((r)_c23), \
                      "=m" ((r)_c33))

#define _sse_pair_store_up_c3_c1(r1,r2) \
__asm__ __volatile__ ("movlps %%xmm3, %0 \n\t" \
                      "movlps %%xmm4, %1 \n\t" \
                      "movlps %%xmm5, %2 \n\t" \
                      "movhps %%xmm3, %3 \n\t" \
                      "movhps %%xmm4, %4 \n\t" \
                      "movhps %%xmm5, %5" \
                      : \
                      "=m" ((r1)_c13), \
                      "=m" ((r1)_c23), \
                      "=m" ((r1)_c33), \
                      "=m" ((r2)_c11), \
                      "=m" ((r2)_c21), \
                      "=m" ((r2)_c31))


/* SSE or SSE2 specific routines  */
#if defined(SSE2)



#define _sse_pair_load_c1_c2(s) \
__asm__ __volatile__ ("movsd %0, %%xmm0 \n\t" \
                      "movsd %1, %%xmm1 \n\t" \
                      "movsd %2, %%xmm2 \n\t" \
                      "movhps %3, %%xmm0 \n\t" \
                      "movhps %4, %%xmm1 \n\t" \
                      "movhps %5, %%xmm2" \
                      : \
                      : \
                      "m" ((s)_c11), \
                      "m" ((s)_c21), \
                      "m" ((s)_c31), \
                      "m" ((s)_c12), \
                      "m" ((s)_c22), \
                      "m" ((s)_c32))


#define _sse_pair_load_c3_c4(s) \
__asm__ __volatile__ ("movsd %0, %%xmm0 \n\t" \
                      "movsd %1, %%xmm1 \n\t" \
                      "movsd %2, %%xmm2 \n\t" \
                      "movhps %3, %%xmm0 \n\t" \
                      "movhps %4, %%xmm1 \n\t" \
                      "movhps %5, %%xmm2" \
                      : \
                      : \
                      "m" ((s)_c13), \
                      "m" ((s)_c23), \
                      "m" ((s)_c33), \
                      "m" ((s)_c14), \
                      "m" ((s)_c24), \
                      "m" ((s)_c34))

#define _sse_pair_load_c1_c2_2sites(s1,s2) \
__asm__ __volatile__ ("movsd %0, %%xmm0 \n\t" \
                      "movsd %1, %%xmm1 \n\t" \
                      "movsd %2, %%xmm2 \n\t" \
                      "movhps %3, %%xmm0 \n\t" \
                      "movhps %4, %%xmm1 \n\t" \
                      "movhps %5, %%xmm2" \
                      : \
                      : \
                      "m" ((s1)_c11), \
                      "m" ((s1)_c21), \
                      "m" ((s1)_c31), \
                      "m" ((s2)_c11), \
                      "m" ((s2)_c21), \
                      "m" ((s2)_c31))


#define _sse_pair_load_c2_c3(s) \
__asm__ __volatile__ ("movsd %0, %%xmm0 \n\t" \
                      "movsd %1, %%xmm1 \n\t" \
                      "movsd %2, %%xmm2 \n\t" \
                      "movhps %3, %%xmm0 \n\t" \
                      "movhps %4, %%xmm1 \n\t" \
                      "movhps %5, %%xmm2" \
                      : \
                      : \
                      "m" ((s)_c12), \
                      "m" ((s)_c22), \
                      "m" ((s)_c32), \
                      "m" ((s)_c13), \
                      "m" ((s)_c23), \
                      "m" ((s)_c33))

#define _sse_pair_load_c3_c1(s1,s2) \
__asm__ __volatile__ ("movsd %0, %%xmm0 \n\t" \
                      "movsd %1, %%xmm1 \n\t" \
                      "movsd %2, %%xmm2 \n\t" \
                      "movhps %3, %%xmm0 \n\t" \
                      "movhps %4, %%xmm1 \n\t" \
                      "movhps %5, %%xmm2" \
                      : \
                      : \
                      "m" ((s1)_c13), \
                      "m" ((s1)_c23), \
                      "m" ((s1)_c33), \
                      "m" ((s2)_c11), \
                      "m" ((s2)_c21), \
                      "m" ((s2)_c31))



#define _sse_su3_multiply_3x1_2sites(u,u1) \
__asm__ __volatile__ ("movss %0, %%xmm3 \n\t" \
                      "movss %1, %%xmm6 \n\t"   /*move initial scalars in */ \
                      "movss %2, %%xmm4 \n\t" \
                      "movss %3, %%xmm7 \n\t" \
                      "movss %4, %%xmm5 \n\t" \
  		      "movhpd  %5, %%xmm3 \n\t"     /*now the low quadword has the current site, high quadword next site */ \
		      "movhpd %6, %%xmm6 \n\t"    /*added */  \
		      "movhpd %7, %%xmm4 \n\t"    /* is the cast to a double required? */ \
                      "shufps $0xa0, %%xmm3, %%xmm3 \n\t"    /* for now assume most significant bit is in the second to highest register */ \
                      "shufps $0xa0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0xa0, %%xmm4, %%xmm4 \n\t"     /*now the top quadword has the current site, low quadword next site */ \
		      "movhpd %8, %%xmm7 \n\t"       /*moves next sites' u to low quadword of xmm7 */ \
		      "movhpd %9, %%xmm5 \n\t"       /* how to overlap these loads and ops ? */ \
		      "mulps %%xmm0, %%xmm3 \n\t" \
		      "shufps $0xa0, %%xmm7, %%xmm7 \n\t" \
		      "mulps %%xmm1, %%xmm6 \n\t" \
		      "shufps $0xa0, %%xmm5, %%xmm5 \n\t" \
		      /* break here this asm call */		\
                      : \
                      : \
                      "m" ((u)_c11re), \
                      "m" ((u)_c12re), \
                      "m" ((u)_c21re), \
                      "m" ((u)_c23re), \
                      "m" ((u)_c31re), \
		      "m" ((u1)_c11re),  /* typecast to a double not needed and would be incorrect */ \
                      "m" (((u1)_c12re)), \
                      "m" (((u1)_c21re)), \
                      "m" (((u1)_c23re)), \
                      "m" (((u1)_c31re))); \
__asm__ __volatile__ ("mulps %%xmm0, %%xmm4 \n\t" \
                      "addps %%xmm6, %%xmm3 \n\t" \
                      "mulps %%xmm2, %%xmm7 \n\t" \
                      "mulps %%xmm0, %%xmm5 \n\t" \
                      "addps %%xmm7, %%xmm4 \n\t" \
		      "movss %0, %%xmm6 \n\t" \
                      "movss %1, %%xmm7 \n\t" \
		      "movhpd %4, %%xmm6 \n\t"                  /*move next sites' u in to low quad of xmm6 */ \
		      "movhpd %5, %%xmm7 \n\t" \
                      "shufps $0xa0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0xa0, %%xmm7, %%xmm7 \n\t" \
		      "mulps %%xmm1, %%xmm6 \n\t" \
                      "mulps %%xmm2, %%xmm7 \n\t" \
                      "addps %%xmm6, %%xmm5 \n\t" \
                      "addps %%xmm7, %%xmm3 \n\t" \
		      "movss %2, %%xmm6 \n\t" \
                      "movss %3, %%xmm7 \n\t" \
		      "movhpd %6, %%xmm6 \n\t"           /*next sites' worth */ \
		      "movhpd %7, %%xmm7 \n\t" \
                      "shufps $0xa0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0xa0, %%xmm7, %%xmm7 \n\t" \
		      "mulps %%xmm1, %%xmm6 \n\t" \
                      "mulps %%xmm2, %%xmm7 \n\t" \
                      "addps %%xmm6, %%xmm4 \n\t" \
                      "addps %%xmm7, %%xmm5" \
		      : \
		      : \
  		      "m" ((u)_c32re), \
                      "m" ((u)_c13re), \
                      "m" ((u)_c22re), \
                      "m" ((u)_c33re), \
                      "m" ((u1)_c32re), \
                      "m" ((u1)_c13re), \
                      "m" ((u1)_c22re), \
                      "m" ((u1)_c33re)); \
__asm__ __volatile__ ("movss %0, %%xmm6 \n\t" \
                      "movss %1, %%xmm7 \n\t" \
                      "movhpd %4, %%xmm6 \n\t" \
                      "movhpd %5, %%xmm7 \n\t" \
                      "shufps $0xb1, %%xmm0, %%xmm0 \n\t" /* optimize later */ \
                      "shufps $0xb1, %%xmm1, %%xmm1 \n\t" \
                      "shufps $0xb1, %%xmm2, %%xmm2 \n\t" \
                      "shufps $0xa0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0xa0, %%xmm7, %%xmm7 \n\t" \
                      "mulps %8, %%xmm0 \n\t" \
                      "mulps %8, %%xmm1 \n\t" \
                      "mulps %8, %%xmm2 \n\t" \
                      "mulps %%xmm0, %%xmm6 \n\t" \
                      "mulps %%xmm1, %%xmm7 \n\t" \
                      "addps %%xmm6, %%xmm3 \n\t" \
                      "addps %%xmm7, %%xmm4 \n\t" \
                      "movss %2, %%xmm6 \n\t" \
                      "movss %3, %%xmm7 \n\t" \
                      "movhpd %6, %%xmm6 \n\t" \
                      "movhpd %7, %%xmm7 \n\t" \
                      "shufps $0xa0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0xa0, %%xmm7, %%xmm7 \n\t" \
                      : \
                      : \
                      "m" ((u)_c11im), \
                      "m" ((u)_c22im), \
                      "m" ((u)_c33im), \
                      "m" ((u)_c21im), \
                      "m" ((u1)_c11im), \
                      "m" ((u1)_c22im), \
                      "m" ((u1)_c33im), \
                      "m" ((u1)_c21im), \
                      "m" (_sse_sgn13)); \
__asm__ __volatile__ ("mulps %%xmm2, %%xmm6 \n\t" \
                      "mulps %%xmm0, %%xmm7 \n\t" \
                      "addps %%xmm6, %%xmm5 \n\t" \
                      "addps %%xmm7, %%xmm4 \n\t" \
                      "movss %0, %%xmm6 \n\t" \
                      "movss %1, %%xmm7 \n\t" \
                      "movhpd %5, %%xmm6 \n\t" \
                      "movhpd %6, %%xmm7 \n\t" \
                      "shufps $0xa0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0xa0, %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm1, %%xmm6 \n\t" \
                      "mulps %%xmm0, %%xmm7 \n\t" \
                      "addps %%xmm6, %%xmm3 \n\t" \
                      "addps %%xmm7, %%xmm5 \n\t" \
                      "movss %2, %%xmm0 \n\t" \
                      "movss %3, %%xmm6 \n\t" \
                      "movss %4, %%xmm7 \n\t" \
                      "movhpd %7, %%xmm0 \n\t" \
                      "movhpd %8, %%xmm6 \n\t" \
                      "movhpd %9, %%xmm7 \n\t" \
                      "shufps $0xa0, %%xmm0, %%xmm0 \n\t" \
                      "shufps $0xa0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0xa0, %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm2, %%xmm0 \n\t" \
                      "mulps %%xmm1, %%xmm6 \n\t" \
                      "mulps %%xmm2, %%xmm7 \n\t" \
                      "addps %%xmm0, %%xmm3 \n\t" \
                      "addps %%xmm6, %%xmm5 \n\t" \
                      "addps %%xmm7, %%xmm4" \
                      : \
                      : \
                      "m" ((u)_c12im), \
                      "m" ((u)_c31im), \
                      "m" ((u)_c13im), \
                      "m" ((u)_c32im), \
                      "m" ((u)_c23im), \
                      "m" ((u1)_c12im), \
                      "m" ((u1)_c31im), \
                      "m" ((u1)_c13im), \
                      "m" ((u1)_c32im), \
                      "m" ((u1)_c23im)); 
                     

#else  /* ! SSE2 */

#define _sse_pair_load_c1_c2(s) \
__asm__ __volatile__ ("movlps %0, %%xmm0 \n\t" \
                      "movlps %1, %%xmm1 \n\t" \
                      "movlps %2, %%xmm2 \n\t" \
                      "movhps %3, %%xmm0 \n\t" \
                      "movhps %4, %%xmm1 \n\t" \
                      "movhps %5, %%xmm2" \
                      : \
                      : \
                      "m" ((s)_c11), \
                      "m" ((s)_c21), \
                      "m" ((s)_c31), \
                      "m" ((s)_c12), \
                      "m" ((s)_c22), \
                      "m" ((s)_c32))


#define _sse_pair_load_c3_c4(s) \
__asm__ __volatile__ ("movlps %0, %%xmm0 \n\t" \
                      "movlps %1, %%xmm1 \n\t" \
                      "movlps %2, %%xmm2 \n\t" \
                      "movhps %3, %%xmm0 \n\t" \
                      "movhps %4, %%xmm1 \n\t" \
                      "movhps %5, %%xmm2" \
                      : \
                      : \
                      "m" ((s)_c13), \
                      "m" ((s)_c23), \
                      "m" ((s)_c33), \
                      "m" ((s)_c14), \
                      "m" ((s)_c24), \
                      "m" ((s)_c34))

#define _sse_pair_load_c1_c2_2sites(s1,s2) \
__asm__ __volatile__ ("movlps %0, %%xmm0 \n\t" \
                      "movlps %1, %%xmm1 \n\t" \
                      "movlps %2, %%xmm2 \n\t" \
                      "movhps %3, %%xmm0 \n\t" \
                      "movhps %4, %%xmm1 \n\t" \
                      "movhps %5, %%xmm2" \
                      : \
                      : \
                      "m" ((s1)_c11), \
                      "m" ((s1)_c21), \
                      "m" ((s1)_c31), \
                      "m" ((s2)_c11), \
                      "m" ((s2)_c21), \
                      "m" ((s2)_c31))


#define _sse_pair_load_c2_c3(s) \
__asm__ __volatile__ ("movlps %0, %%xmm0 \n\t" \
                      "movlps %1, %%xmm1 \n\t" \
                      "movlps %2, %%xmm2 \n\t" \
                      "movhps %3, %%xmm0 \n\t" \
                      "movhps %4, %%xmm1 \n\t" \
                      "movhps %5, %%xmm2" \
                      : \
                      : \
                      "m" ((s)_c12), \
                      "m" ((s)_c22), \
                      "m" ((s)_c32), \
                      "m" ((s)_c13), \
                      "m" ((s)_c23), \
                      "m" ((s)_c33))

#define _sse_pair_load_c3_c1(s1,s2) \
__asm__ __volatile__ ("movlps %0, %%xmm0 \n\t" \
                      "movlps %1, %%xmm1 \n\t" \
                      "movlps %2, %%xmm2 \n\t" \
                      "movhps %3, %%xmm0 \n\t" \
                      "movhps %4, %%xmm1 \n\t" \
                      "movhps %5, %%xmm2" \
                      : \
                      : \
                      "m" ((s1)_c13), \
                      "m" ((s1)_c23), \
                      "m" ((s1)_c33), \
                      "m" ((s2)_c11), \
                      "m" ((s2)_c21), \
                      "m" ((s2)_c31))

#define _sse_su3_multiply_3x1_2sites(u,u1) \
__asm__ __volatile__ ("movss %0, %%xmm3 \n\t" \
                      "movss %1, %%xmm6 \n\t"   /*move initial scalars in */ \
                      "movss %2, %%xmm4 \n\t" \
                      "movss %3, %%xmm7 \n\t" \
                      "movss %4, %%xmm5 \n\t" \
  		      "movhps  %5, %%xmm3 \n\t"     /*now the low quadword has the current site, high quadword next site */ \
		      "movhps %6, %%xmm6 \n\t"    /*added */  \
		      "movhps %7, %%xmm4 \n\t"    /* is the cast to a double required? */ \
                      "shufps $0xa0, %%xmm3, %%xmm3 \n\t"    /* for now assume most significant bit is in the second to highest register */ \
                      "shufps $0xa0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0xa0, %%xmm4, %%xmm4 \n\t"     /*now the top quadword has the current site, low quadword next site */ \
		      "movhps %8, %%xmm7 \n\t"       /*moves next sites' u to low quadword of xmm7 */ \
		      "movhps %9, %%xmm5 \n\t"       /* how to overlap these loads and ops ? */ \
		      "mulps %%xmm0, %%xmm3 \n\t" \
		      "shufps $0xa0, %%xmm7, %%xmm7 \n\t" \
		      "mulps %%xmm1, %%xmm6 \n\t" \
		      "shufps $0xa0, %%xmm5, %%xmm5 \n\t" \
		      /* break here this asm call */		\
                      : \
                      : \
                      "m" ((u)_c11re), \
                      "m" ((u)_c12re), \
                      "m" ((u)_c21re), \
                      "m" ((u)_c23re), \
                      "m" ((u)_c31re), \
		      "m" ((u1)_c11re),  /* typecast to a double not needed and would be incorrect */ \
                      "m" (((u1)_c12re)), \
                      "m" (((u1)_c21re)), \
                      "m" (((u1)_c23re)), \
                      "m" (((u1)_c31re))); \
__asm__ __volatile__ ("mulps %%xmm0, %%xmm4 \n\t" \
                      "addps %%xmm6, %%xmm3 \n\t" \
                      "mulps %%xmm2, %%xmm7 \n\t" \
                      "mulps %%xmm0, %%xmm5 \n\t" \
                      "addps %%xmm7, %%xmm4 \n\t" \
		      "movss %0, %%xmm6 \n\t" \
                      "movss %1, %%xmm7 \n\t" \
		      "movhps %4, %%xmm6 \n\t"                  /*move next sites' u in to low quad of xmm6 */ \
		      "movhps %5, %%xmm7 \n\t" \
                      "shufps $0xa0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0xa0, %%xmm7, %%xmm7 \n\t" \
		      "mulps %%xmm1, %%xmm6 \n\t" \
                      "mulps %%xmm2, %%xmm7 \n\t" \
                      "addps %%xmm6, %%xmm5 \n\t" \
                      "addps %%xmm7, %%xmm3 \n\t" \
		      "movss %2, %%xmm6 \n\t" \
                      "movss %3, %%xmm7 \n\t" \
		      "movhps %6, %%xmm6 \n\t"           /*next sites' worth */ \
		      "movhps %7, %%xmm7 \n\t" \
                      "shufps $0xa0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0xa0, %%xmm7, %%xmm7 \n\t" \
		      "mulps %%xmm1, %%xmm6 \n\t" \
                      "mulps %%xmm2, %%xmm7 \n\t" \
                      "addps %%xmm6, %%xmm4 \n\t" \
                      "addps %%xmm7, %%xmm5" \
		      : \
		      : \
  		      "m" ((u)_c32re), \
                      "m" ((u)_c13re), \
                      "m" ((u)_c22re), \
                      "m" ((u)_c33re), \
                      "m" ((u1)_c32re), \
                      "m" ((u1)_c13re), \
                      "m" ((u1)_c22re), \
                      "m" ((u1)_c33re)); \
__asm__ __volatile__ ("movss %0, %%xmm6 \n\t" \
                      "movss %1, %%xmm7 \n\t" \
                      "movhps %4, %%xmm6 \n\t" \
                      "movhps %5, %%xmm7 \n\t" \
                      "shufps $0xb1, %%xmm0, %%xmm0 \n\t" /* optimize later */ \
                      "shufps $0xb1, %%xmm1, %%xmm1 \n\t" \
                      "shufps $0xb1, %%xmm2, %%xmm2 \n\t" \
                      "shufps $0xa0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0xa0, %%xmm7, %%xmm7 \n\t" \
                      "mulps %8, %%xmm0 \n\t" \
                      "mulps %8, %%xmm1 \n\t" \
                      "mulps %8, %%xmm2 \n\t" \
                      "mulps %%xmm0, %%xmm6 \n\t" \
                      "mulps %%xmm1, %%xmm7 \n\t" \
                      "addps %%xmm6, %%xmm3 \n\t" \
                      "addps %%xmm7, %%xmm4 \n\t" \
                      "movss %2, %%xmm6 \n\t" \
                      "movss %3, %%xmm7 \n\t" \
                      "movhps %6, %%xmm6 \n\t" \
                      "movhps %7, %%xmm7 \n\t" \
                      "shufps $0xa0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0xa0, %%xmm7, %%xmm7 \n\t" \
                      : \
                      : \
                      "m" ((u)_c11im), \
                      "m" ((u)_c22im), \
                      "m" ((u)_c33im), \
                      "m" ((u)_c21im), \
                      "m" ((u1)_c11im), \
                      "m" ((u1)_c22im), \
                      "m" ((u1)_c33im), \
                      "m" ((u1)_c21im), \
                      "m" (_sse_sgn13)); \
__asm__ __volatile__ ("mulps %%xmm2, %%xmm6 \n\t" \
                      "mulps %%xmm0, %%xmm7 \n\t" \
                      "addps %%xmm6, %%xmm5 \n\t" \
                      "addps %%xmm7, %%xmm4 \n\t" \
                      "movss %0, %%xmm6 \n\t" \
                      "movss %1, %%xmm7 \n\t" \
                      "movhps %5, %%xmm6 \n\t" \
                      "movhps %6, %%xmm7 \n\t" \
                      "shufps $0xa0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0xa0, %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm1, %%xmm6 \n\t" \
                      "mulps %%xmm0, %%xmm7 \n\t" \
                      "addps %%xmm6, %%xmm3 \n\t" \
                      "addps %%xmm7, %%xmm5 \n\t" \
                      "movss %2, %%xmm0 \n\t" \
                      "movss %3, %%xmm6 \n\t" \
                      "movss %4, %%xmm7 \n\t" \
                      "movhps %7, %%xmm0 \n\t" \
                      "movhps %8, %%xmm6 \n\t" \
                      "movhps %9, %%xmm7 \n\t" \
                      "shufps $0xa0, %%xmm0, %%xmm0 \n\t" \
                      "shufps $0xa0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0xa0, %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm2, %%xmm0 \n\t" \
                      "mulps %%xmm1, %%xmm6 \n\t" \
                      "mulps %%xmm2, %%xmm7 \n\t" \
                      "addps %%xmm0, %%xmm3 \n\t" \
                      "addps %%xmm6, %%xmm5 \n\t" \
                      "addps %%xmm7, %%xmm4" \
                      : \
                      : \
                      "m" ((u)_c12im), \
                      "m" ((u)_c31im), \
                      "m" ((u)_c13im), \
                      "m" ((u)_c32im), \
                      "m" ((u)_c23im), \
                      "m" ((u1)_c12im), \
                      "m" ((u1)_c31im), \
                      "m" ((u1)_c13im), \
                      "m" ((u1)_c32im), \
                      "m" ((u1)_c23im)); 

#endif

#endif
