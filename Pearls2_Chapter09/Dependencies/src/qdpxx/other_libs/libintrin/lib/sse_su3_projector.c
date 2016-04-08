#include "sse_su3_projector.h"
#include <xmmintrin.h>

#ifdef __cplusplus
extern "C" { 
#endif

void
intrin_sse_su3_projector(su3_vectorf* aa, su3_vectorf* bb, su3_matrixf* cc)
{

	 /* XMM Variables */
	 __m128 xmm2, xmm3, xmm0, xmm1;

	xmm0 = _mm_loadl_pi(xmm0, (__m64 *)&((bb)->c[0]) );
	xmm0 = _mm_loadh_pi(xmm0, (__m64 *)&((bb)->c[1]) );
	xmm1 = xmm0 ; 
	xmm1 = _mm_shuffle_ps( xmm1, xmm1, 0xb1 );
	 	 xmm1 = _mm_xor_ps( xmm1, _sse_sgn24.xmm );
	xmm2 = _mm_load_ss((float *)&((aa)->c[0].real) );
	xmm2 = _mm_shuffle_ps( xmm2, xmm2, 0x00 );
	xmm3 = _mm_load_ss((float *)&((aa)->c[0].imag) );
	xmm3 = _mm_shuffle_ps( xmm3, xmm3, 0x00 );
	xmm2 = _mm_mul_ps( xmm2, xmm0 );
	xmm3 = _mm_mul_ps( xmm3, xmm1 );
	xmm2 = _mm_add_ps( xmm2, xmm3 );
	 	 xmm2 = _mm_xor_ps( xmm2, _sse_sgn24.xmm );
	_mm_storeu_ps((float *)&((cc)->e[0][0]), xmm2 );
	xmm2 = _mm_load_ss((float *)&((aa)->c[1].real) );
	xmm2 = _mm_shuffle_ps( xmm2, xmm2, 0x00 );
	xmm3 = _mm_load_ss((float *)&((aa)->c[1].imag) );
	xmm3 = _mm_shuffle_ps( xmm3, xmm3, 0x00 );
	xmm2 = _mm_mul_ps( xmm2, xmm0 );
	xmm3 = _mm_mul_ps( xmm3, xmm1 );
	xmm2 = _mm_add_ps( xmm2, xmm3 );
	 	 xmm2 = _mm_xor_ps( xmm2, _sse_sgn24.xmm );
	_mm_storeu_ps((float *)&((cc)->e[1][0]), xmm2 );
	xmm2 = _mm_load_ss((float *)&((aa)->c[2].real) );
	xmm2 = _mm_shuffle_ps( xmm2, xmm2, 0x00 );
	xmm3 = _mm_load_ss((float *)&((aa)->c[2].imag) );
	xmm3 = _mm_shuffle_ps( xmm3, xmm3, 0x00 );
	xmm2 = _mm_mul_ps( xmm2, xmm0 );
	xmm3 = _mm_mul_ps( xmm3, xmm1 );
	xmm2 = _mm_add_ps( xmm2, xmm3 );
	 	 xmm2 = _mm_xor_ps( xmm2, _sse_sgn24.xmm );
	_mm_storeu_ps((float *)&((cc)->e[2][0]), xmm2 );
	xmm0 = _mm_loadl_pi(xmm0, (__m64 *)&((aa)->c[0]) );
	xmm0 = _mm_loadh_pi(xmm0, (__m64 *)&((aa)->c[1]) );
	xmm1 = xmm0 ; 
	xmm1 = _mm_shuffle_ps( xmm1, xmm1, 0xb1 );
	 	 xmm1 = _mm_xor_ps( xmm1, _sse_sgn24.xmm );
	xmm2 = _mm_load_ss((float *)&((bb)->c[2].real) );
	xmm2 = _mm_shuffle_ps( xmm2, xmm2, 0x00 );
	xmm3 = _mm_load_ss((float *)&((bb)->c[2].imag) );
	xmm3 = _mm_shuffle_ps( xmm3, xmm3, 0x00 );
	xmm2 = _mm_mul_ps( xmm2, xmm0 );
	xmm3 = _mm_mul_ps( xmm3, xmm1 );
	xmm2 = _mm_add_ps( xmm2, xmm3 );
	_mm_storel_pi((__m64 *)&((cc)->e[0][2]), xmm2 );
	_mm_storeh_pi((__m64 *)&((cc)->e[1][2]), xmm2 );
	xmm0 = _mm_loadl_pi(xmm0, (__m64 *)&((aa)->c[2]) );
	xmm0 = _mm_shuffle_ps( xmm0, xmm0, 0x14 );
	xmm2 = _mm_loadl_pi(xmm2, (__m64 *)&((bb)->c[2]) );
	xmm2 = _mm_shuffle_ps( xmm2, xmm2, 0x44 );
	 	 xmm2 = _mm_xor_ps( xmm2, _sse_sgn4.xmm );
	xmm2 = _mm_mul_ps( xmm2, xmm0 );
	xmm1 = xmm2 ; 
	xmm1 = _mm_shuffle_ps( xmm1, xmm1, 0xd4 );
	xmm2 = _mm_shuffle_ps( xmm2, xmm2, 0x8c );
	xmm2 = _mm_add_ps( xmm2, xmm1 );
	_mm_storeh_pi((__m64 *)&((cc)->e[2][2]), xmm2 );
}

#ifdef __cplusplus
}; 
#endif

