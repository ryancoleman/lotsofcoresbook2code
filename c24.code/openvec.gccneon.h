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
               Intrinsics for ARM NEON

Compile flags -O3 -mfpu=neon

*******************************************************************************/

#define __STDC_CONSTANT_MACROS
#include"arm_neon.h"
#define OV_FLOAT_WIDTH 4
#define OV_ALIGN (OV_FLOAT_WIDTH*OV_SIZEOF_FLOAT)

#define OV_PLATFORM "GCC COMPILER ARM/NEON"

#define OV_CPP_OVERRIDE_TYPE

#define ov_prv_float float32x4_t
#define ov_maskf ov_prv_float

#define ov_prv_int int32x4_t



#define ov_prv_addf vaddq_f32
#define ov_prv_subf vsubq_f32
#define ov_prv_mulf vmulq_f32

/* DIVISION */
#define ov_prv_rcpf__ vrecpeq_f32
#define ov_prv_rcpf_(x) ov_prv_mulf(vrecpsq_f32(x, ov_prv_rcpf__(x)), ov_prv_rcpf__(x))
#define ov_prv_rcpf(x) ov_prv_mulf(vrecpsq_f32(x, ov_prv_rcpf_(x)), ov_prv_rcpf_(x))

#define ov_prv_divf__(a,b) ov_prv_mulf(a, ov_prv_rcpf__(b))
#define ov_prv_divf_(a,b)  ov_prv_mulf(a, ov_prv_rcpf_(b))
#define ov_prv_divf(a,b)  ov_prv_mulf(a, ov_prv_rcpf(b))


#define ov_prv_loadf vld1q_f32
#define ov_prv_ldf vld1q_f32
#define ov_prv_uloadf(dest,addr) dest=vld1q_f32(addr)
#define ov_prv_uldf vld1q_f32

#define ov_prv_storef vst1q_f32
#define ov_prv_stf ov_prv_storef
#define ov_prv_storeuf ov_prv_storef
#define ov_prv_stream_stf ov_prv_storef
#define ov_prv_nt_storef ov_prv_storef

#define ov_prv_setf vdupq_n_f32
#define ov_prv_setzerof(a) a=ov_prv_setf(0.0f)
#define ov_prv_getzerof() ov_prv_setf(0.0f)


#define ov_prv_seti vdupq_n_s32

#define ov_prv_casti_f vreinterpretq_s32_f32 
#define ov_prv_castf_i vreinterpretq_f32_s32
#define ov_prv_castf_ui vreinterpretq_f32_u32

#define ov_prv_float_to_int vcvtq_s32_f32
#define ov_prv_int_to_float vcvtq_f32_s32

#define ov_prv_maxf vmaxq_f32
#define ov_prv_minf vminq_f32

#define ov_prv_rsqrtf vrsqrteq_f32
#define ov_prv_xhalf(x) ov_prv_mulf(ov_prv_setf(0.5f), x)
#define ov_prv_prec_rsqrtf(x) ov_prv_mulf(ov_prv_rsqrtf(x), ov_prv_subf(ov_prv_setf(1.5f), ov_prv_mulf(ov_prv_xhalf(x), ov_prv_sqrf(ov_prv_rsqrtf(x)))))
#define ov_prv_sqrtf(y) ov_prv_mulf(y, ov_prv_prec_rsqrtf(ov_prv_maxf(y, ov_prv_setf(OV_EPS))))


#define ov_prv_maddf(a,b,c) vmlaq_f32(c,b,a)

//#define ov_prv_msubf(a,b,c) ov_prv_mulf(ov_setf(-1.0f), vmlsq_f32(c,b,a))
#define ov_prv_msubf(a,b,c) ov_prv_subf(ov_prv_mulf(a,b),c)



#define ov_prv_absf vabsq_f32

 


#define ov_prefetch1(a) __builtin_prefetch(a,0,3)
#define ov_prefetch2(a) __builtin_prefetch(a,0,2)
#define ov_prefetch3(a) __builtin_prefetch(a,0,1)

#define ov_bool int64x1_t


#define ov_prv_eqf(a,b) ov_prv_castf_ui(vceqq_f32(a,b))
#define ov_prv_gtf(a,b) ov_prv_castf_ui(vcgtq_f32(a,b))
#define ov_prv_gef(a,b) ov_prv_castf_ui(vcgeq_f32(a,b))
#define ov_prv_ltf(a,b) ov_prv_castf_ui(vcltq_f32(a,b))
#define ov_prv_lef(a,b) ov_prv_castf_ui(vcleq_f32(a,b))


#define ov_prv_xorf ov_prv_castf_i(veorq_s32(ov_prv_casti_f(a), ov_prv_casti_f(b)))
#define ov_prv_andf(a,b) ov_prv_castf_i(vandq_s32(ov_prv_casti_f(a), ov_prv_casti_f(b)))
#define ov_prv_andnf(b,a) ov_prv_castf_i(vbicq_s32(ov_prv_casti_f(a), ov_prv_casti_f(b)))
#define ov_prv_invmaskf(x) ov_prv_castf_i(vmvnq_s32(ov_prv_casti_f(x)))

#define ov_andmaskf(a,b) ov_prv_castf_i(vandq_s32(ov_prv_casti_f(a), ov_prv_casti_f(b)))
#define ov_ormaskf(a,b) ov_prv_castf_i(vorrq_s32(ov_prv_casti_f(a), ov_prv_casti_f(b)))


#define ov_prv_conditionalf(mask, val1, val2) vbslq_f32(vreinterpretq_u32_f32(mask), val1, val2)
/*
#define ov_prv_conditionalf(mask, val1, val2) ov_prv_addf(ov_prv_andf(                mask , val1),\
                                                   ov_prv_andnf(mask, val2)) */



#define OV_ONEF ov_prv_setf(1.0f)
#define ov_prv_round(a) ov_prv_int_to_float(ov_prv_float_to_int(a))


#define ov_prv_fastceilf(a) ov_prv_conditionalf(ov_prv_ltf(ov_prv_round(a), a), ov_prv_addf(ov_prv_round(a), OV_ONEF), ov_prv_round(a))
#define ov_prv_fastfloorf(a) ov_prv_conditionalf(ov_prv_gtf(ov_prv_round(a), a), ov_prv_subf(ov_prv_round(a), OV_ONEF), ov_prv_round(a))

#define OV_PRV_BIGINT ov_prv_setf(2145336164.0f)

#define ov_prv_ceilf(a) ov_prv_conditionalf(ov_prv_ltf(ov_prv_absf(a), OV_PRV_BIGINT), ov_prv_fastceilf(a), a);
#define ov_prv_floorf(a) ov_prv_conditionalf(ov_prv_ltf(ov_prv_absf(a), OV_PRV_BIGINT), ov_prv_fastfloorf(a), a);


#define ov_anyf(mask) vreinterpret_s64_s32(vand_s32(vget_low_s32(ov_prv_casti_f(mask)), vget_high_s32(ov_prv_casti_f(mask))))

#ifdef __LP64__
#define ov_allf(mask) (ov_anyf(mask) == (ov_bool)0xFFFFFFFFFFFFFFFFL)
#else
#define ov_allf(mask) (ov_anyf(mask) == (ov_bool)0xFFFFFFFFFFFFFFFFLL)
#endif


#define ov_any_lt_0f(a) ov_anyf(ov_prv_ltf(ov_prv_vreg(a), ov_prv_getzerof()))
#define ov_all_ge_0f(a) ov_allf(ov_prv_gef(ov_prv_vreg(a), ov_prv_getzerof()))
#define ov_any_ge_0f(a) ov_anyf(ov_prv_gef(ov_prv_vreg(a), ov_prv_getzerof()))
#define ov_all_lt_0f(a) ov_allf(ov_prv_ltf(ov_prv_vreg(a), ov_prv_getzerof()))

/*
   DOUBLE PRECISION
*/
#include"openvec.double_scalar.h"
