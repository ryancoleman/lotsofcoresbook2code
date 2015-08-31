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
          Intrinsics for GCC Advanced Vector Extensions 512 and MIC
*******************************************************************************/
#include"immintrin.h"
#define OV_FLOAT_WIDTH 16
#define OV_ALIGN (OV_FLOAT_WIDTH*OV_SIZEOF_FLOAT)

#ifdef __INTEL_COMPILER
#define OV_COMPILER "Intel"
#define OV_CPP_OVERRIDE_TYPE
#else
#define OV_COMPILER "GCC"
#define OV_CPP_OVERRIDE_TYPE
#endif

#ifdef _OV_MIC
#define OV_PLATFORM OV_COMPILER " COMPILER MIC"
#else
#define OV_PLATFORM OV_COMPILER " COMPILER AVX 512"
#endif


#define ov_prv_float __m512
#define ov_maskf __mmask16
#define OV_ALL_TRUE 0xFFFF

#define ov_andmaskf(a,b) ((a)&(b))
#define ov_ormaskf(a,b) ((a)|(b))



#define ov_prv_addf _mm512_add_ps
#define ov_prv_subf _mm512_sub_ps
#define ov_prv_mulf _mm512_mul_ps
#define ov_prv_divf _mm512_div_ps
#define ov_prv_loadf _mm512_load_ps

#define ov_prv_ldf _mm512_load_ps
#define ov_prv_storef _mm512_store_ps
#define ov_prv_stf _mm512_store_ps


#ifdef _OV_MIC
/* MIC KNC INTRINSICS */

#define ov_prv_uloadf(dest, addr) {\
dest = _mm512_loadunpacklo_ps(dest, addr);\
dest = _mm512_loadunpackhi_ps(dest, (addr)+OV_FLOAT_WIDTH);}

OV_INLINE ov_prv_float ov_prv_uldf(float const *addr)
{
  ov_prv_float tmp;
  ov_prv_uloadf(tmp, addr);
  return tmp;
}

#define ov_prv_storeuf(addr, a) {\
  _mm512_packstorelo_ps(addr, a);\
  _mm512_packstorehi_ps((addr)+OV_FLOAT_WIDTH, a);}
/* Need to improve nt stores for MIC */
#define ov_prv_stream_stf _mm512_storenr_ps
#define ov_prv_nt_storef ov_prv_stf
#else
/* AVX-512 INTRINSICS */
#define ov_prv_uldf _mm512_loadu_ps
#define ov_prv_uloadf(dest,addr) dest=_mm512_loadu_ps(addr)
#define ov_prv_storeuf _mm512_storeu_ps
#define ov_prv_stream_stf _mm512_stream_ps
#define ov_prv_nt_storef _mm512_stream_ps
#endif



#define ov_prv_setzerof(a) a=_mm512_setzero_ps()
#define ov_prv_getzerof _mm512_setzero_ps

#define ov_prv_setf _mm512_set1_ps
#define ov_prv_seti _mm512_set1_epi32

#define ov_casti_f _mm512_castsi512_ps
#define ov_castf_i _mm512_castps_si512
#define ov_float_to_int _mm512_cvtps_epi32
#define ov_prv_maxf _mm512_max_ps
#define ov_prv_minf _mm512_min_ps

#define ov_prv_rsqrtf _mm512_invsqrt_ps

#define ov_prv_sqrtf _mm512_sqrt_ps
#define ov_prv_rcpf(x) ov_prv_divf(ov_prv_setf(1.0f), x)

#define ov_prv_floorf _mm512_floor_ps
#define ov_prv_ceilf _mm512_ceil_ps

#define ov_prv_absf _mm512_abs_ps

#define ov_prv_maddf _mm512_fmadd_ps
#define ov_prv_msubf _mm512_fmsub_ps


#define ov_prefetch1(a) _mm_prefetch((char*)a, _MM_HINT_T0)
#define ov_prefetch2(a) _mm_prefetch((char*)a, _MM_HINT_T1)
#define ov_prefetch3(a) _mm_prefetch((char*)a, _MM_HINT_T2)


#define ov_prv_eqf(a,b) _mm512_cmp_ps_mask(a, b, _MM_CMPINT_EQ)
#define ov_prv_gtf(a,b) _mm512_cmp_ps_mask(a, b, _MM_CMPINT_GT)
#define ov_prv_gef(a,b) _mm512_cmp_ps_mask(a, b, _MM_CMPINT_GE)
#define ov_prv_ltf(a,b) _mm512_cmp_ps_mask(a, b, _MM_CMPINT_LT)
#define ov_prv_lef(a,b) _mm512_cmp_ps_mask(a, b, _MM_CMPINT_LE)


#define ov_allf(mask) ((mask)==OV_ALL_TRUE)
#define ov_anyf(mask) (mask)


#define ov_prv_conditionalf(mask,a,b) _mm512_mask_blend_ps(mask, b, a)


#define ov_any_lt_0f(a) ov_anyf(ov_prv_ltf(ov_prv_vreg(a), ov_prv_getzerof()))
#define ov_all_ge_0f(a) ov_allf(ov_prv_gef(ov_prv_vreg(a), ov_prv_getzerof()))
#define ov_any_ge_0f(a) ov_anyf(ov_prv_gef(ov_prv_vreg(a), ov_prv_getzerof()))
#define ov_all_lt_0f(a) ov_allf(ov_prv_ltf(ov_prv_vreg(a), ov_prv_getzerof()))



/*
  DOUBLE PRECISION
*/
#define OV_DOUBLE_WIDTH 8
#define ov_prv_double __m512d
#define ov_maskd __mmask8
#define OV_ALL_TRUED 255


#define ov_andmaskd(a,b) ((a)&(b))
#define ov_ormaskd(a,b) ((a)|(b))


#define ov_prv_addd _mm512_add_pd
#define ov_prv_subd _mm512_sub_pd
#define ov_prv_muld _mm512_mul_pd
#define ov_prv_divd _mm512_div_pd
#define ov_prv_loadd _mm512_load_pd

#define ov_prv_ldd _mm512_load_pd
#define ov_prv_stored _mm512_store_pd
#define ov_prv_std _mm512_store_pd


#ifdef _OV_MIC
/* MIC KNC INTRINSICS */

#define ov_prv_uloadd(dest, addr) {\
dest = _mm512_loadunpacklo_pd(dest, addr);\
dest = _mm512_loadunpackhi_pd(dest, (addr)+OV_DOUBLE_WIDTH);}

OV_INLINE ov_prv_double ov_prv_uldd(double const *addr)
{
  ov_prv_double tmp;
  ov_prv_uloadd(tmp, addr);
  return tmp;
}

#define ov_prv_storeud(addr, a) {\
  _mm512_packstorelo_pd(addr, a);\
  _mm512_packstorehi_pd((addr)+OV_DOUBLE_WIDTH, a);}
/* Need to improve nt stores for MIC */
#define ov_prv_stream_std _mm512_storenr_pd
#define ov_prv_nt_stored ov_prv_std
#else
/* AVX-512 INTRINSICS */
#define ov_prv_uldd _mm512_loadu_pd
#define ov_prv_uloadd(dest,addr) dest=_mm512_loadu_pd(addr)
#define ov_prv_storeud _mm512_storeu_pd
#define ov_prv_stream_std _mm512_stream_pd
#define ov_prv_nt_stored _mm512_stream_pd
#endif



#define ov_prv_setzerod(a) a=_mm512_setzero_pd()
#define ov_prv_getzerod _mm512_setzero_pd

#define ov_prv_setd _mm512_set1_pd
#define ov_prv_seti _mm512_set1_epi32

#define ov_casti_d _mm512_castsi512_pd
#define ov_castf_i _mm512_castps_si512
#define ov_double_to_int _mm512_cvtps_epi32
#define ov_prv_maxd _mm512_max_pd
#define ov_prv_mind _mm512_min_pd


#define ov_prv_sqrtd _mm512_sqrt_pd
#define ov_prv_rcpd(x) ov_prv_divd(ov_prv_setd(1.0),x)
#define ov_prv_rsqrtd(x) ov_prv_rcpd(ov_prv_sqrtd(x))


#define ov_prv_floord _mm512_floor_pd
#define ov_prv_ceild _mm512_ceil_pd

#define ov_prv_absd _mm512_abs_pd

#define ov_prv_maddd _mm512_fmadd_pd
#define ov_prv_msubd _mm512_fmsub_pd


#define ov_prefetch1(a) _mm_prefetch((char*)a, _MM_HINT_T0)
#define ov_prefetch2(a) _mm_prefetch((char*)a, _MM_HINT_T1)
#define ov_prefetch3(a) _mm_prefetch((char*)a, _MM_HINT_T2)


#define ov_prv_eqd(a,b) _mm512_cmp_pd_mask(a, b, _MM_CMPINT_EQ)
#define ov_prv_gtd(a,b) _mm512_cmp_pd_mask(a, b, _MM_CMPINT_GT)
#define ov_prv_ged(a,b) _mm512_cmp_pd_mask(a, b, _MM_CMPINT_GE)
#define ov_prv_ltd(a,b) _mm512_cmp_pd_mask(a, b, _MM_CMPINT_LT)
#define ov_prv_led(a,b) _mm512_cmp_pd_mask(a, b, _MM_CMPINT_LE)


#define ov_alld(mask) ((mask)==OV_ALL_TRUED)
#define ov_anyd(mask) (mask)


#define ov_prv_conditionald(mask,a,b) _mm512_mask_blend_pd(mask, b, a)


#define ov_any_lt_0d(a) ov_anyd(ov_prv_ltd(ov_prv_vreg(a), ov_prv_getzerod()))
#define ov_all_ge_0d(a) ov_alld(ov_prv_ged(ov_prv_vreg(a), ov_prv_getzerod()))
#define ov_any_ge_0d(a) ov_anyd(ov_prv_ged(ov_prv_vreg(a), ov_prv_getzerod()))
#define ov_all_lt_0d(a) ov_alld(ov_prv_ltd(ov_prv_vreg(a), ov_prv_getzerod()))


