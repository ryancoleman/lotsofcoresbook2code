#include "sse_sub_four_su3_vecs.h"
#include <xmmintrin.h>

#ifdef __cplusplus
extern "C" { 
#endif

void
intrin_sse_sub_four_su3_vecs(su3_vectorf* aa, su3_vectorf* bb0, su3_vectorf* bb1, su3_vectorf* bb2, su3_vectorf* bb3)
{

	 /* XMM Variables */
	 __m128 xmm2, xmm3, xmm0, xmm1;

	xmm0 = _mm_loadu_ps((float *)&((aa)->c[0]) );
	xmm1 = _mm_loadl_pi(xmm1, (__m64 *)&((aa)->c[2]) );
	xmm1 = _mm_shuffle_ps( xmm1, xmm1, 0x44 );
	xmm2 = _mm_loadu_ps((float *)&((bb0)->c[0]) );
	xmm3 = _mm_loadl_pi(xmm3, (__m64 *)&((bb0)->c[2]) );
	xmm3 = _mm_shuffle_ps( xmm3, xmm3, 0x44 );
	xmm0 = _mm_sub_ps( xmm0, xmm2 );
	xmm1 = _mm_sub_ps( xmm1, xmm3 );
	xmm2 = _mm_loadu_ps((float *)&((bb1)->c[0]) );
	xmm3 = _mm_loadl_pi(xmm3, (__m64 *)&((bb1)->c[2]) );
	xmm3 = _mm_shuffle_ps( xmm3, xmm3, 0x44 );
	xmm0 = _mm_sub_ps( xmm0, xmm2 );
	xmm1 = _mm_sub_ps( xmm1, xmm3 );
	xmm2 = _mm_loadu_ps((float *)&((bb2)->c[0]) );
	xmm3 = _mm_loadl_pi(xmm3, (__m64 *)&((bb2)->c[2]) );
	xmm3 = _mm_shuffle_ps( xmm3, xmm3, 0x44 );
	xmm0 = _mm_sub_ps( xmm0, xmm2 );
	xmm1 = _mm_sub_ps( xmm1, xmm3 );
	xmm2 = _mm_loadu_ps((float *)&((bb3)->c[0]) );
	xmm3 = _mm_loadl_pi(xmm3, (__m64 *)&((bb3)->c[2]) );
	xmm3 = _mm_shuffle_ps( xmm3, xmm3, 0x44 );
	xmm0 = _mm_sub_ps( xmm0, xmm2 );
	xmm1 = _mm_sub_ps( xmm1, xmm3 );
	_mm_storeu_ps((float *)&((aa)->c[0]), xmm0 );
	_mm_storel_pi((__m64 *)&((aa)->c[2]), xmm1 );
}

#ifdef __cplusplus
}; 
#endif

