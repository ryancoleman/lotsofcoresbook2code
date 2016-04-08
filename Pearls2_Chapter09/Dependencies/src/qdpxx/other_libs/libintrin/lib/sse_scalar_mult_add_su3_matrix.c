#include "sse_scalar_mult_add_su3_matrix.h"
#include <xmmintrin.h>

#ifdef __cplusplus
extern "C" { 
#endif

void
intrin_sse_scalar_mult_add_su3_matrix(su3_matrixf* aa, su3_matrixf* bb, float cc, su3_matrixf* dd)
{

	 /* XMM Variables */
	 __m128 xmm0, xmm1, xmm4;

	xmm4 = _mm_load_ss((float *)&((cc)) );
	xmm4 = _mm_shuffle_ps( xmm4, xmm4, 0x00 );
	xmm0 = _mm_loadu_ps((float *)&((aa)->e[0][0]) );
	xmm1 = _mm_loadu_ps((float *)&((bb)->e[0][0]) );
	xmm1 = _mm_mul_ps( xmm1, xmm4 );
	xmm0 = _mm_add_ps( xmm0, xmm1 );
	_mm_storeu_ps((float *)&((dd)->e[0][0]), xmm0 );
	xmm0 = _mm_loadu_ps((float *)&((aa)->e[0][2]) );
	xmm1 = _mm_loadu_ps((float *)&((bb)->e[0][2]) );
	xmm1 = _mm_mul_ps( xmm1, xmm4 );
	xmm0 = _mm_add_ps( xmm0, xmm1 );
	_mm_storeu_ps((float *)&((dd)->e[0][2]), xmm0 );
	xmm0 = _mm_loadu_ps((float *)&((aa)->e[1][1]) );
	xmm1 = _mm_loadu_ps((float *)&((bb)->e[1][1]) );
	xmm1 = _mm_mul_ps( xmm1, xmm4 );
	xmm0 = _mm_add_ps( xmm0, xmm1 );
	_mm_storeu_ps((float *)&((dd)->e[1][1]), xmm0 );
	xmm0 = _mm_loadu_ps((float *)&((aa)->e[2][0]) );
	xmm1 = _mm_loadu_ps((float *)&((bb)->e[2][0]) );
	xmm1 = _mm_mul_ps( xmm1, xmm4 );
	xmm0 = _mm_add_ps( xmm0, xmm1 );
	_mm_storeu_ps((float *)&((dd)->e[2][0]), xmm0 );
	xmm0 = _mm_loadl_pi(xmm0, (__m64 *)&((aa)->e[2][2]) );
	xmm1 = _mm_loadl_pi(xmm1, (__m64 *)&((bb)->e[2][2]) );
	xmm1 = _mm_mul_ps( xmm1, xmm4 );
	xmm0 = _mm_add_ps( xmm0, xmm1 );
	_mm_storel_pi((__m64 *)&((dd)->e[2][2]), xmm0 );
}

#ifdef __cplusplus
}; 
#endif

