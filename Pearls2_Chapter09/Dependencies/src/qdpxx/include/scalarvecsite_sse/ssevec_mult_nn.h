#if 1
#define PREFETCH(addr)  __asm__ __volatile__("prefetcht0 %0"::"m"(*(addr)))
#else
#define PREFETCH(addr)
#endif

#define _inline_ssevec_mult_su3_nn(cc,aa,bb,j) \
{ \
__asm__ __volatile__ (                    \
              "movaps %0,%%xmm0\n\t"      \
              "movaps %%xmm0,%%xmm1\n\t"  \
              "mulps  %2,%%xmm1\n\t"      \
              "movaps %%xmm0,%%xmm2\n\t"  \
              "mulps  %3,%%xmm2\n\t"      \
              "movaps %%xmm0,%%xmm3\n\t"  \
              "mulps  %4,%%xmm3\n\t"      \
              "movaps %%xmm0,%%xmm4\n\t"  \
              "mulps  %5,%%xmm4\n\t"      \
              "movaps %%xmm0,%%xmm5\n\t"  \
              "mulps  %6,%%xmm5\n\t"      \
              "mulps  %1,%%xmm0\n\t"      \
	      :                           \
	      : "m" (*(bb+0+8*j)),     \
		"m" (*(aa+0)),     \
		"m" (*(aa+4)),     \
		"m" (*(aa+24)),     \
		"m" (*(aa+28)),     \
		"m" (*(aa+48)),     \
		"m" (*(aa+52)));    \
__asm__ __volatile__ (                    \
              "movaps %0,%%xmm6\n\t"      \
              "movaps %%xmm6,%%xmm7\n\t"  \
              "mulps  %2,%%xmm7\n\t"      \
              "subps  %%xmm7,%%xmm0\n\t"  \
              "movaps %%xmm6,%%xmm7\n\t"  \
              "mulps  %1,%%xmm7\n\t"      \
              "addps  %%xmm7,%%xmm1\n\t"  \
              "movaps %%xmm6,%%xmm7\n\t"  \
              "mulps  %4,%%xmm7\n\t"      \
              "subps  %%xmm7,%%xmm2\n\t"  \
              "movaps %%xmm6,%%xmm7\n\t"  \
              "mulps  %3,%%xmm7\n\t"      \
              "addps  %%xmm7,%%xmm3\n\t"  \
              "movaps %%xmm6,%%xmm7\n\t"  \
              "mulps  %6,%%xmm7\n\t"      \
              "subps  %%xmm7,%%xmm4\n\t"  \
              "mulps  %5,%%xmm6\n\t"      \
              "addps  %%xmm6,%%xmm5\n\t"  \
	      :                           \
	      : "m" (*(bb+4+8*j)),     \
		"m" (*(aa+0)),     \
		"m" (*(aa+4)),     \
		"m" (*(aa+24)),     \
		"m" (*(aa+28)),     \
		"m" (*(aa+48)),     \
		"m" (*(aa+52)));    \
__asm__ __volatile__ (                    \
              "movaps %0,%%xmm6\n\t"      \
              "movaps %%xmm6,%%xmm7\n\t"  \
              "mulps  %1,%%xmm7\n\t"      \
              "addps  %%xmm7,%%xmm0\n\t"  \
              "movaps %%xmm6,%%xmm7\n\t"  \
              "mulps  %2,%%xmm7\n\t"      \
              "addps  %%xmm7,%%xmm1\n\t"  \
              "movaps %%xmm6,%%xmm7\n\t"  \
              "mulps  %3,%%xmm7\n\t"      \
              "addps  %%xmm7,%%xmm2\n\t"  \
              "movaps %%xmm6,%%xmm7\n\t"  \
              "mulps  %4,%%xmm7\n\t"      \
              "addps  %%xmm7,%%xmm3\n\t"  \
              "movaps %%xmm6,%%xmm7\n\t"  \
              "mulps  %5,%%xmm7\n\t"      \
              "addps  %%xmm7,%%xmm4\n\t"  \
              "mulps  %6,%%xmm6\n\t"      \
              "addps  %%xmm6,%%xmm5\n\t"  \
	      :                           \
	      : "m" (*(bb+24+8*j)),     \
		"m" (*(aa+8)),     \
		"m" (*(aa+12)),     \
		"m" (*(aa+32)),     \
		"m" (*(aa+36)),     \
		"m" (*(aa+56)),     \
		"m" (*(aa+60)));    \
__asm__ __volatile__ (                    \
              "movaps %0,%%xmm6\n\t"      \
              "movaps %%xmm6,%%xmm7\n\t"  \
              "mulps  %2,%%xmm7\n\t"      \
              "subps  %%xmm7,%%xmm0\n\t"  \
              "movaps %%xmm6,%%xmm7\n\t"  \
              "mulps  %1,%%xmm7\n\t"      \
              "addps  %%xmm7,%%xmm1\n\t"  \
              "movaps %%xmm6,%%xmm7\n\t"  \
              "mulps  %4,%%xmm7\n\t"      \
              "subps  %%xmm7,%%xmm2\n\t"  \
              "movaps %%xmm6,%%xmm7\n\t"  \
              "mulps  %3,%%xmm7\n\t"      \
              "addps  %%xmm7,%%xmm3\n\t"  \
              "movaps %%xmm6,%%xmm7\n\t"  \
              "mulps  %6,%%xmm7\n\t"      \
              "subps  %%xmm7,%%xmm4\n\t"  \
              "mulps  %5,%%xmm6\n\t"      \
              "addps  %%xmm6,%%xmm5\n\t"  \
	      :                           \
	      : "m" (*(bb+28+8*j)),     \
		"m" (*(aa+8)),     \
		"m" (*(aa+12)),     \
		"m" (*(aa+32)),     \
		"m" (*(aa+36)),     \
		"m" (*(aa+56)),     \
		"m" (*(aa+60)));    \
__asm__ __volatile__ (                    \
              "movaps %0,%%xmm6\n\t"      \
              "movaps %%xmm6,%%xmm7\n\t"  \
              "mulps  %1,%%xmm7\n\t"      \
              "addps  %%xmm7,%%xmm0\n\t"  \
              "movaps %%xmm6,%%xmm7\n\t"  \
              "mulps  %2,%%xmm7\n\t"      \
              "addps  %%xmm7,%%xmm1\n\t"  \
              "movaps %%xmm6,%%xmm7\n\t"  \
              "mulps  %3,%%xmm7\n\t"      \
              "addps  %%xmm7,%%xmm2\n\t"  \
              "movaps %%xmm6,%%xmm7\n\t"  \
              "mulps  %4,%%xmm7\n\t"      \
              "addps  %%xmm7,%%xmm3\n\t"  \
              "movaps %%xmm6,%%xmm7\n\t"  \
              "mulps  %5,%%xmm7\n\t"      \
              "addps  %%xmm7,%%xmm4\n\t"  \
              "mulps  %6,%%xmm6\n\t"      \
              "addps  %%xmm6,%%xmm5\n\t"  \
	      :                           \
	      : "m" (*(bb+48+8*j)),     \
		"m" (*(aa+16)),     \
		"m" (*(aa+20)),     \
		"m" (*(aa+40)),     \
		"m" (*(aa+44)),     \
		"m" (*(aa+64)),     \
		"m" (*(aa+68)));    \
__asm__ __volatile__ (                    \
              "movaps %0,%%xmm6\n\t"      \
              "movaps %%xmm6,%%xmm7\n\t"  \
              "mulps  %2,%%xmm7\n\t"      \
              "subps  %%xmm7,%%xmm0\n\t"  \
              "movaps %%xmm6,%%xmm7\n\t"  \
              "mulps  %1,%%xmm7\n\t"      \
              "addps  %%xmm7,%%xmm1\n\t"  \
              "movaps %%xmm6,%%xmm7\n\t"  \
              "mulps  %4,%%xmm7\n\t"      \
              "subps  %%xmm7,%%xmm2\n\t"  \
              "movaps %%xmm6,%%xmm7\n\t"  \
              "mulps  %3,%%xmm7\n\t"      \
              "addps  %%xmm7,%%xmm3\n\t"  \
              "movaps %%xmm6,%%xmm7\n\t"  \
              "mulps  %6,%%xmm7\n\t"      \
              "subps  %%xmm7,%%xmm4\n\t"  \
              "mulps  %5,%%xmm6\n\t"      \
              "addps  %%xmm6,%%xmm5\n\t"  \
	      :                           \
	      : "m" (*(bb+52+8*j)),     \
		"m" (*(aa+16)),     \
		"m" (*(aa+20)),     \
		"m" (*(aa+40)),     \
		"m" (*(aa+44)),     \
		"m" (*(aa+64)),     \
		"m" (*(aa+68)));    \
__asm__ __volatile__ (                    \
              "movaps %%xmm0,%0\n\t"      \
              "movaps %%xmm1,%1\n\t"      \
              "movaps %%xmm2,%2\n\t"      \
              "movaps %%xmm3,%3\n\t"      \
              "movaps %%xmm4,%4\n\t"      \
              "movaps %%xmm5,%5\n\t"      \
	      : "=m" (*(cc+0+8*j)),    \
		"=m" (*(cc+4+8*j)),    \
		"=m" (*(cc+24+8*j)),    \
		"=m" (*(cc+28+8*j)),    \
		"=m" (*(cc+48+8*j)),    \
		"=m" (*(cc+52+8*j)));   \
}

