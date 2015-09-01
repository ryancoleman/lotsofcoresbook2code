#include "sse_scalar_mult_add_su3_vector.h"
#include <xmmintrin.h>

#ifdef __cplusplus
extern "C" { 
#endif

void
intrin_sse_scalar_mult_add_su3_vector(su3_vectorf* aa, su3_vectorf* bb, float cc, su3_vectorf* dd)
{

	 /* XMM Variables */
	 __m128 xmm2, xmm3, xmm0, xmm1, xmm4;

	xmm4 = _mm_load_ss((float *)&((cc)) );
	xmm4 = _mm_shuffle_ps( xmm4, xmm4, 0x00 );
	xmm0 = _mm_loadu_ps((float *)&((aa)->c[0]) );
	xmm1 = _mm_loadl_pi(xmm1, (__m64 *)&((aa)->c[2]) );
	xmm1 = _mm_shuffle_ps( xmm1, xmm1, 0x44 );
	xmm2 = _mm_loadu_ps((float *)&((bb)->c[0]) );
	xmm3 = _mm_loadl_pi(xmm3, (__m64 *)&((bb)->c[2]) );
	xmm3 = _mm_shuffle_ps( xmm3, xmm3, 0x44 );
	xmm2 = _mm_mul_ps( xmm2, xmm4 );
	xmm3 = _mm_mul_ps( xmm3, xmm4 );
	xmm0 = _mm_add_ps( xmm0, xmm2 );
	xmm1 = _mm_add_ps( xmm1, xmm3 );
	_mm_storeu_ps((float *)&((dd)->c[0]), xmm0 );
	_mm_storel_pi((__m64 *)&((dd)->c[2]), xmm1 );
}

#ifdef __cplusplus
}; 
#endif

