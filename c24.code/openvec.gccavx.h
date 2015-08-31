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
               Intrinsics for Advanced Vector Extensions
*******************************************************************************/
#include"immintrin.h"
#define OV_FLOAT_WIDTH 8
#define OV_ALIGN (OV_FLOAT_WIDTH*OV_SIZEOF_FLOAT)

#ifdef __INTEL_COMPILER
#define OV_COMPILER "Intel"
#define OV_CPP_OVERRIDE_TYPE
#else
#define OV_COMPILER "GCC"
#define OV_CPP_OVERRIDE_TYPE
#endif

#if defined(_OV_AVX2)
#define OV_PLATFORM OV_COMPILER " COMPILER AVX 2"
#else
#define OV_PLATFORM OV_COMPILER " COMPILER AVX"
#endif


#define ov_prv_float __m256
#define ov_maskf ov_prv_float
#define OV_ALL_TRUE 255


#define ov_prv_addf _mm256_add_ps
#define ov_prv_subf _mm256_sub_ps
#define ov_prv_mulf _mm256_mul_ps
#define ov_prv_divf _mm256_div_ps
#define ov_prv_loadf _mm256_load_ps
#define ov_prv_uloadf(dest,addr) dest=_mm256_loadu_ps(addr)

#define ov_prv_ldf _mm256_load_ps
#define ov_prv_uldf _mm256_loadu_ps

#define ov_prv_storef _mm256_store_ps
#define ov_prv_stf _mm256_store_ps

#define ov_prv_storeuf _mm256_storeu_ps
#define ov_prv_stream_stf _mm256_stream_ps
#define ov_prv_nt_storef _mm256_stream_ps
#define ov_prv_setzerof(a) a=_mm256_setzero_ps()
#define ov_prv_getzerof _mm256_setzero_ps

#define ov_prv_setf _mm256_set1_ps
#define ov_prv_seti _mm256_set1_epi32

#define ov_casti_f _mm256_castsi256_ps
#define ov_float_to_int _mm256_cvtps_epi32
#define ov_prv_maxf _mm256_max_ps
#define ov_prv_minf _mm256_min_ps

#define ov_prv_rsqrtf _mm256_rsqrt_ps

#define ov_prv_sqrtf _mm256_sqrt_ps
#define ov_prv_rcpf _mm256_rcp_ps

#define ov_prv_floorf(a) _mm256_round_ps((a), 0x09)
#define ov_prv_ceilf(a) _mm256_round_ps((a), 0x0A)

#define ov_andmaskf _mm256_and_ps
#define ov_ormaskf _mm256_or_ps

#define OV_MASK_ABSF 0x7FFFFFFF
#define ov_prv_absf(x) _mm256_and_ps(x, _mm256_castsi256_ps(ov_prv_seti(OV_MASK_ABSF)))

#if defined(_OV_AVX2)
#define ov_prv_maddf _mm256_fmadd_ps
#define ov_prv_msubf _mm256_fmsub_ps
#else
#define ov_prv_maddf(a,b,c) ov_prv_addf(c, ov_prv_mulf(a,b))
#define ov_prv_msubf(a,b,c) ov_prv_subf(ov_prv_mulf(a,b),c)

#endif


#define ov_prefetch1(a) _mm_prefetch((char*)a, _MM_HINT_T0)
#define ov_prefetch2(a) _mm_prefetch((char*)a, _MM_HINT_T1)
#define ov_prefetch3(a) _mm_prefetch((char*)a, _MM_HINT_T2)


#define ov_prv_eqf(a,b) _mm256_cmp_ps(a, b, _CMP_EQ_OS)
#define ov_prv_gtf(a,b) _mm256_cmp_ps(a, b, _CMP_GT_OS)
#define ov_prv_gef(a,b) _mm256_cmp_ps(a, b, _CMP_GE_OS)
#define ov_prv_ltf(a,b) _mm256_cmp_ps(a, b, _CMP_LT_OS)
#define ov_prv_lef(a,b) _mm256_cmp_ps(a, b, _CMP_LE_OS)


#define ov_allf(mask) (_mm256_movemask_ps(mask) == OV_ALL_TRUE)
#define ov_anyf _mm256_movemask_ps


#define ov_prv_conditionalf(mask, val1, val2) _mm256_blendv_ps(val2, val1, mask)
  

#define ov_any_lt_0f(a) _mm256_movemask_ps(ov_prv_vreg(a))
#define ov_all_ge_0f(a) (!ov_any_lt_0f(a))
#define ov_any_ge_0f(a) (_mm256_movemask_ps(ov_prv_vreg(a))!=OV_ALL_TRUE)
#define ov_all_lt_0f(a) (_mm256_movemask_ps(ov_prv_vreg(a))==OV_ALL_TRUE)

/* 
    DOUBLE PRECISION
*/
#define OV_DOUBLE_WIDTH 4

#define ov_prv_double __m256d
#define ov_maskd ov_prv_double
#define OV_ALL_TRUED 15


#define ov_prv_addd _mm256_add_pd
#define ov_prv_subd _mm256_sub_pd
#define ov_prv_muld _mm256_mul_pd
#define ov_prv_divd _mm256_div_pd
#define ov_prv_loadd _mm256_load_pd
#define ov_prv_uloadd(dest,addr) dest=_mm256_loadu_pd(addr)

#define ov_prv_ldd _mm256_load_pd
#define ov_prv_uldd _mm256_loadu_pd

#define ov_prv_stored _mm256_store_pd
#define ov_prv_std _mm256_store_pd

#define ov_prv_storeud _mm256_storeu_pd
#define ov_prv_stream_std _mm256_stream_pd
#define ov_prv_nt_stored _mm256_stream_pd
#define ov_prv_setzerod(a) a=_mm256_setzero_pd()
#define ov_prv_getzerod _mm256_setzero_pd

#define ov_prv_setd _mm256_set1_pd
#define ov_prv_seti _mm256_set1_epi32

#define ov_casti_d _mm256_castsi256_pd
#define ov_float_to_int _mm256_cvtps_epi32
#define ov_prv_maxd _mm256_max_pd
#define ov_prv_mind _mm256_min_pd


#define ov_prv_sqrtd _mm256_sqrt_pd
#define ov_prv_rcpd(x) ov_prv_divd(ov_prv_setd(1.0),x)
#define ov_prv_rsqrtd(x) ov_prv_rcpd(ov_prv_sqrtd(x))

#define ov_prv_floord(a) _mm256_round_pd((a), 0x09)
#define ov_prv_ceild(a) _mm256_round_pd((a), 0x0A)

#define ov_andmaskd _mm256_and_pd
#define ov_ormaskd _mm256_or_pd

#ifdef __LP64__
static long const OV_MASK_ABS64=0x7FFFFFFFFFFFFFFFl;
#else
static long long const OV_MASK_ABS64=0x7FFFFFFFFFFFFFFFll;
#endif
#define ov_prv_absd(x) _mm256_and_pd(x, ov_prv_setd(*(double*)&OV_MASK_ABS64))

#if defined(_OV_AVX2)
#define ov_prv_maddd _mm256_fmadd_pd
#define ov_prv_msubd _mm256_fmsub_pd
#else
#define ov_prv_maddd(a,b,c) ov_prv_addd(c, ov_prv_muld(a,b))
#define ov_prv_msubd(a,b,c) ov_prv_subd(ov_prv_muld(a,b),c)

#endif


#define ov_prefetch1(a) _mm_prefetch((char*)a, _MM_HINT_T0)
#define ov_prefetch2(a) _mm_prefetch((char*)a, _MM_HINT_T1)
#define ov_prefetch3(a) _mm_prefetch((char*)a, _MM_HINT_T2)


#define ov_prv_eqd(a,b) _mm256_cmp_pd(a, b, _CMP_EQ_OS)
#define ov_prv_gtd(a,b) _mm256_cmp_pd(a, b, _CMP_GT_OS)
#define ov_prv_ged(a,b) _mm256_cmp_pd(a, b, _CMP_GE_OS)
#define ov_prv_ltd(a,b) _mm256_cmp_pd(a, b, _CMP_LT_OS)
#define ov_prv_led(a,b) _mm256_cmp_pd(a, b, _CMP_LE_OS)


#define ov_alld(mask) (_mm256_movemask_pd(mask) == OV_ALL_TRUED)
#define ov_anyd _mm256_movemask_pd


#define ov_prv_conditionald(mask, val1, val2) _mm256_blendv_pd(val2, val1, mask)
  

#define ov_any_lt_0d(a) _mm256_movemask_pd(ov_prv_vreg(a))
#define ov_all_ge_0d(a) (!ov_any_lt_0d(a))
#define ov_any_ge_0d(a) (_mm256_movemask_pd(ov_prv_vreg(a))!=OV_ALL_TRUED)
#define ov_all_lt_0d(a) (_mm256_movemask_pd(ov_prv_vreg(a))==OV_ALL_TRUED)
