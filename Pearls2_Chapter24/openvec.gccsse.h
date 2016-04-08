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
/*******************************************************************************
               Intrinsics for Streaming SIMD Extensions (SSE)
*******************************************************************************/
#include"immintrin.h"
#define OV_FLOAT_WIDTH 4
#define OV_ALIGN (OV_FLOAT_WIDTH*OV_SIZEOF_FLOAT)


#ifdef __INTEL_COMPILER
#define OV_COMPILER "Intel"
#define OV_CPP_OVERRIDE_TYPE
#else
#define OV_COMPILER "GCC"
#define OV_CPP_OVERRIDE_TYPE
#endif


#if defined(_OV_SSE4)
#define OV_PLATFORM OV_COMPILER " COMPILER SSE 4"
#elif defined(_OV_SSE2)
#define OV_PLATFORM OV_COMPILER " COMPILER SSE 2"
#else
#define OV_PLATFORM OV_COMPILER " COMPILER SSE"
#endif


#define ov_prv_float __m128
#define ov_maskf ov_prv_float
#define OV_ALL_TRUE 15


#define OV_MASK_ABSF 0x7FFFFFFF

#define ov_prv_addf _mm_add_ps
#define ov_prv_subf _mm_sub_ps
#define ov_prv_mulf _mm_mul_ps
#define ov_prv_divf _mm_div_ps
#define ov_prv_loadf _mm_load_ps
#define ov_prv_ldf _mm_load_ps
#define ov_prv_uloadf(dest,addr) dest=_mm_loadu_ps(addr)
#define ov_prv_uldf _mm_loadu_ps

#define ov_prv_storef _mm_store_ps
#define ov_prv_stf _mm_store_ps
#define ov_prv_storeuf _mm_storeu_ps
#define ov_prv_stream_stf _mm_stream_ps
#define ov_prv_nt_storef _mm_stream_ps
#define ov_prv_setzerof(a) a=_mm_setzero_ps()
#define ov_prv_getzerof _mm_setzero_ps

#define ov_prv_setf _mm_set1_ps
#define ov_prv_seti _mm_set1_epi32

#define ov_casti_f _mm_castsi128_ps
#define ov_prv_casti_f _mm_castsi128_ps
#define ov_float_to_int _mm_cvtps_epi32
#define ov_prv_maxf _mm_max_ps
#define ov_prv_minf _mm_min_ps

#define ov_prv_rsqrtf _mm_rsqrt_ps

#define ov_prv_sqrtf _mm_sqrt_ps
#define ov_prv_rcpf _mm_rcp_ps



#define ov_prv_absf(x) _mm_and_ps(x, _mm_castsi128_ps(ov_prv_seti(OV_MASK_ABSF)))

#define ov_prv_maddf(a,b,c) (ov_prv_addf(c, ov_prv_mulf(a,b)))
#define ov_prv_msubf(a,b,c) (ov_prv_subf(ov_prv_mulf(a,b),c))

#define ov_prefetch1(a) _mm_prefetch((char*)a, _MM_HINT_T0)
#define ov_prefetch2(a) _mm_prefetch((char*)a, _MM_HINT_T1)
#define ov_prefetch3(a) _mm_prefetch((char*)a, _MM_HINT_T2)

#define ov_prv_eqf _mm_cmpeq_ps
#define ov_prv_gtf _mm_cmpgt_ps
#define ov_prv_gef _mm_cmpge_ps
#define ov_prv_ltf _mm_cmplt_ps
#define ov_prv_lef _mm_cmple_ps

#define ov_prv_andf _mm_and_ps
#define ov_andmaskf _mm_and_ps
#define ov_ormaskf _mm_or_ps

#define ov_prv_conditionalf(mask, val1, val2)  _mm_or_ps(ov_prv_andf(mask, val1), \
                                                           _mm_andnot_ps(mask, val2))


#define ov_allf(a)  (_mm_movemask_ps(a)==OV_ALL_TRUE)
#define ov_anyf _mm_movemask_ps




/* FLOOR and CEIL SSE 4 */
#if defined(_OV_SSE4)

#define ov_prv_floorf(a) _mm_floor_ps(a)
#define ov_prv_ceilf(a)  _mm_ceil_ps(a)

/* FLOOR and CEIL SSE */
#else
#define OV_ONEF ov_prv_setf(1.0f)
#define ov_prv_round(a) _mm_cvtepi32_ps(_mm_cvtps_epi32(a))

#define ov_prv_fastceilf(a) ov_prv_conditionalf(ov_prv_ltf(ov_prv_round(a), a), ov_prv_addf(ov_prv_round(a), OV_ONEF), ov_prv_round(a))
#define ov_prv_fastfloorf(a) ov_prv_conditionalf(ov_prv_gtf(ov_prv_round(a), a), ov_prv_subf(ov_prv_round(a), OV_ONEF), ov_prv_round(a))

#define OV_PRV_BIGINT ov_prv_setf(2145336164.0f)

#define ov_prv_ceilf(a) ov_prv_conditionalf(ov_prv_ltf(ov_prv_absf(a), OV_PRV_BIGINT), ov_prv_fastceilf(a), a);
#define ov_prv_floorf(a) ov_prv_conditionalf(ov_prv_ltf(ov_prv_absf(a), OV_PRV_BIGINT), ov_prv_fastfloorf(a), a);

#endif



#define ov_any_lt_0f(a) _mm_movemask_ps(ov_prv_vreg(a))
#define ov_all_ge_0f(a) (!ov_any_lt_0f(a))
#define ov_any_ge_0f(a) (_mm_movemask_ps(ov_prv_vreg(a))!=OV_ALL_TRUE)
#define ov_all_lt_0f(a) (_mm_movemask_ps(ov_prv_vreg(a))==OV_ALL_TRUE)
                                                                          

/*
    DOUBLE PRECISION
*/
#if defined(_OV_SSE4) || defined(_OV_SSE2)
#define OV_DOUBLE_WIDTH 2
#define ov_prv_double __m128d
#define ov_maskd ov_prv_double
#define OV_ALL_TRUED 3



#define ov_prv_addd _mm_add_pd
#define ov_prv_subd _mm_sub_pd
#define ov_prv_muld _mm_mul_pd
#define ov_prv_divd _mm_div_pd
#define ov_prv_loadd _mm_load_pd
#define ov_prv_ldd _mm_load_pd
#define ov_prv_uloadd(dest,addr) dest=_mm_loadu_pd(addr)
#define ov_prv_uldd _mm_loadu_pd

#define ov_prv_stored _mm_store_pd
#define ov_prv_std _mm_store_pd
#define ov_prv_storeud _mm_storeu_pd
#define ov_prv_stream_std _mm_stream_pd
#define ov_prv_nt_stored _mm_stream_pd
#define ov_prv_setzerod(a) a=_mm_setzero_pd()
#define ov_prv_getzerod _mm_setzero_pd

#define ov_prv_setd _mm_set1_pd

#define ov_casti_d _mm_castsi128_pd
#define ov_prv_casti_d _mm_castsi128_pd
#define ov_prv_maxd _mm_max_pd
#define ov_prv_mind _mm_min_pd

#define ov_prv_sqrtd _mm_sqrt_pd

#define ov_prv_rcpd(x) ov_prv_divd(ov_prv_setd(1.0),x)
#define ov_prv_rsqrtd(x) ov_prv_rcpd(ov_prv_sqrtd(x))

#ifdef __LP64__
static long const OV_MASK_ABS64=0x7FFFFFFFFFFFFFFFl;
#else
static long long const OV_MASK_ABS64=0x7FFFFFFFFFFFFFFFll;
#endif
#define ov_prv_absd(x) _mm_and_pd(x, ov_prv_setd(*(double*)&OV_MASK_ABS64))

#define ov_prv_maddd(a,b,c) (ov_prv_addd(c, ov_prv_muld(a,b)))
#define ov_prv_msubd(a,b,c) (ov_prv_subd(ov_prv_muld(a,b),c))


#define ov_prv_eqd _mm_cmpeq_pd
#define ov_prv_gtd _mm_cmpgt_pd
#define ov_prv_ged _mm_cmpge_pd
#define ov_prv_ltd _mm_cmplt_pd
#define ov_prv_led _mm_cmple_pd

#define ov_prv_andd _mm_and_pd

#define ov_andmaskd _mm_and_pd
#define ov_ormaskd _mm_or_pd

#define ov_prv_conditionald(mask, val1, val2)  _mm_or_pd(ov_prv_andd(mask, val1), \
                                                           _mm_andnot_pd(mask, val2))


#define ov_alld(a)  (_mm_movemask_pd(a)==OV_ALL_TRUED)
#define ov_anyd _mm_movemask_pd




/* FLOOR and CEIL SSE 4 */
#if defined(_OV_SSE4)

#define ov_prv_floord(a) _mm_floor_pd(a)
#define ov_prv_ceild(a)  _mm_ceil_pd(a)

/* FLOOR and CEIL SSE */
#else

#define OV_ONEd ov_prv_setd(1.0)

#define ov_prv_roundlow(a) _mm_cvtsi64_sd(a, (_mm_cvtsd_si64(a)))
#define ov_prv_roundd(a) _mm_unpacklo_pd(ov_prv_roundlow(a), ov_prv_roundlow(_mm_unpackhi_pd(a,a)))

#define ov_prv_fastceild(a) ov_prv_conditionald(ov_prv_ltd(ov_prv_roundd(a), a), ov_prv_addd(ov_prv_roundd(a), OV_ONEd), ov_prv_roundd(a))
#define ov_prv_fastfloord(a) ov_prv_conditionald(ov_prv_gtd(ov_prv_roundd(a), a), ov_prv_subd(ov_prv_roundd(a), OV_ONEd), ov_prv_roundd(a))

#define OV_PRV_BIGINT64 ov_prv_setd(9.e18)

#define ov_prv_ceild(a) ov_prv_conditionald(ov_prv_ltd(ov_prv_absd(a), OV_PRV_BIGINT64), ov_prv_fastceild(a), a);
#define ov_prv_floord(a) ov_prv_conditionald(ov_prv_ltd(ov_prv_absd(a), OV_PRV_BIGINT64), ov_prv_fastfloord(a), a);

#endif



#define ov_any_lt_0d(a) _mm_movemask_pd(ov_prv_vreg(a))
#define ov_all_ge_0d(a) (!ov_any_lt_0d(a))
#define ov_any_ge_0d(a) (_mm_movemask_pd(ov_prv_vreg(a))!=OV_ALL_TRUED)
#define ov_all_lt_0d(a) (_mm_movemask_pd(ov_prv_vreg(a))==OV_ALL_TRUED)

#else
#include"openvec.double_scalar.h"                                                                           
#endif /* defined(_OV_SSE4) || defined(_OV_SSE2) */
