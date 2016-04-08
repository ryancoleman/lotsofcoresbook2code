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
/*
*** NEED TO FIX: unaligned loads and stores *****

UNALIGNED LOAD: ov_prv_uldf with vec_xld2 intrinsic 

UNALIGNED STORE: ov_prv_storeuf with vec_xstd2 intrinsic 

*** END NEED TO FIX *****

Compiler Options
–qarch=pwr7 –qtune=pwr7 –O3 –qhot –qsimd

gcc -std=c99 -Wall -O3 -mvsx -maltivec -mtune=power7 -mcpu=power7

https://gcc.gnu.org/onlinedocs/gcc/PowerPC-AltiVec_002fVSX-Built-in-Functions.html]

*/
/*#include<builtins.h>*/
#include<altivec.h>
#define OV_FLOAT_WIDTH 4
#define OV_ALIGN (OV_FLOAT_WIDTH*OV_SIZEOF_FLOAT)

#define OV_PLATFORM "XLC COMPILER VSX"
#define OV_CPP_OVERRIDE_TYPE

#define ov_prv_float vector float
#define ov_maskf ov_prv_float




#define ov_prv_maddf vec_madd
#define ov_prv_msubf vec_msub
#define ov_prv_addf vec_add
#define ov_prv_subf vec_sub
#define ov_prv_mulf vec_mul
#define ov_prv_divf vec_div
#define ov_prv_loadf(addr) vec_ld(0, addr)
//#define ov_prv_loadf(addr) vec_xld2(0, addr)
/*#define ov_prv_loadf(addr) vec_vsx_ld(0, addr) */

#define ov_prv_ldf ov_prv_loadf
/*
#define ov_prv_uldf ov_prv_ldf
#define ov_prv_uloadf(dest,addr) dest=ov_prv_ldf(addr)
*/
#define ov_prv_uldf(a) (*((ov_prv_float *)(a)))
#define ov_prv_uloadf(dest,addr) dest=ov_prv_uldf(addr)

#define ov_prv_storef(addr,val) vec_st(val, 0, addr)
//#define ov_prv_storef(addr,val) vec_xstd2(val, 0, addr)
//#define ov_prv_storeuf ov_prv_storef
#define ov_prv_storeuf(addr,val) *(ov_prv_float*)addr=val
#define ov_prv_stf ov_prv_storef

#define ov_prv_stream_stf(addr,val) ov_prv_storef(addr,val)
#define ov_prv_nt_storef(addr,val) ov_prv_storef(addr,val)

#define ov_prv_getzerof() ((ov_prv_float)vec_splat_u32(0))
#define ov_prv_setzerof(a) a=ov_prv_getzerof()

#define ov_prv_setf vec_splats

#define ov_prv_maxf vec_max
#define ov_prv_minf vec_min

#define ov_prv_rsqrtf vec_rsqrte
#define ov_prv_rcpf vec_re

#define ov_prv_floorf vec_floor
#define ov_prv_ceilf vec_ceil


#define ov_andmaskf vec_and
#define ov_ormaskf vec_or


#define ov_prv_absf vec_abs
#define ov_prv_andf vec_and
#define ov_prv_xorf vec_xor

#define ov_prv_sqrtf(x) ov_prv_mulf(x, ov_prv_rsqrtf(ov_prv_maxf(x, ov_prv_setf(OV_EPS))))

#define ov_prv_gtf (ov_prv_float)vec_cmpgt
#define ov_prv_gef (ov_prv_float)vec_cmpge
#define ov_prv_ltf (ov_prv_float)vec_cmplt
#define ov_prv_lef (ov_prv_float)vec_cmple
#define ov_prv_eqf (ov_prv_float)vec_cmpeq


#define ov_allf(mask) vec_all_ne(mask, ov_prv_getzerof())
#define ov_anyf(mask) vec_any_ne(mask, ov_prv_getzerof())


#define ov_prv_conditionalf(mask, val1, val2) vec_sel(val2, val1, mask)


#ifdef __XLC__
#define ov_prefetch1 __prefetch_by_load
#define ov_prefetch2 __prefetch_by_load
#define ov_prefetch3 __prefetch_by_load
#else
#define ov_prefetch1 __builtin_prefetch
#define ov_prefetch2 __builtin_prefetch
#define ov_prefetch3 __builtin_prefetch
#endif


#define ov_any_lt_0f(a) ov_anyf(ov_prv_ltf(ov_prv_vreg(a), ov_prv_getzerof()))
#define ov_all_ge_0f(a) ov_allf(ov_prv_gef(ov_prv_vreg(a), ov_prv_getzerof()))
#define ov_any_ge_0f(a) ov_anyf(ov_prv_gef(ov_prv_vreg(a), ov_prv_getzerof()))
#define ov_all_lt_0f(a) ov_allf(ov_prv_ltf(ov_prv_vreg(a), ov_prv_getzerof()))

/*
    DOUBLE PRECISION
*/
#include"openvec.double_scalar.h"
