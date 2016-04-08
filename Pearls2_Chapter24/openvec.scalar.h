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
#include<math.h>

#define OV_FLOAT_WIDTH 1
#define OV_PLATFORM "SCALAR"
#define OV_ALIGN (OV_FLOAT_WIDTH*OV_SIZEOF_FLOAT)
#define ov_malloc malloc
#define ov_calloc calloc
#define ov_free free

#define ov_prv_float float
#define ov_addf(a,b) ((a)+(b))
#define ov_subf(a,b) ((a)-(b))
#define ov_mulf(a,b) ((a)*(b))
#define ov_divf(a,b) ((a)/(b))
#define ov_loadf(a) (*(a))
#define ov_ldf(a) (*(a))
#define ov_uloadf(dest,addr) dest=ov_ldf(addr)
#define ov_uldf ov_ldf
#define ov_storef(a,b) *(a)=b
#define ov_storeuf ov_storef
#define ov_stf ov_storeuf
#define ov_setzerof(a) a=0.0f
#define ov_setf(a) (a)


#define ov_andmaskf(a,b) ((a)&(b))
#define ov_ormaskf(a,b) ((a)|(b))



#define ov_minf(X,Y) ((X) < (Y) ? (X) : (Y))
#define ov_maxf(X,Y) ((X) > (Y) ? (X) : (Y))

#define ov_rsqrtf(x) (1.0f/ov_sqrtf(x))
#define ov_rcpf(x) (1.0f/(x))
#define ov_sqrtf sqrtf
#define ov_fastsqrtf sqrtf
#define ov_absf fabsf

/*float ov_casti_f(int i){return *(float*)&i;}*/
#define ov_casti_f (float &)

#define ov_float_to_int(f) ((int)(f))

#define ov_any_gt(a,b) ((a)>(b))
#define ov_all_gt(a,b) ((a)>(b))
#define ov_any_gt_zero(a) ((a)>0.0f)
#define ov_any_ge_zero(a) ((a)>=0.0f)

#define ov_any_lt(a,b) ((a)<(b))
#define ov_all_lt(a,b) ((a)<(b))
#define ov_any_lt_zero(a) ((a)<0.0f)

#define ov_maddf(a,b,c) ((c)+((a)*(b)))
#define ov_msubf(a,b,c) (ov_subf(ov_mulf(a,b),c))

#define ov_ceilf ceilf
#define ov_floorf floorf

#define ov_zerof (0.0f)
#define ov_getzerof() (0.0f)

#define ov_gtf(X,Y) ((X) >  (Y))
#define ov_gef(X,Y) ((X) >= (Y))

#define ov_ltf(X,Y) ((X) <  (Y))
#define ov_lef(X,Y) ((X) <= (Y))

#define ov_eqf(X,Y) ((X) == (Y))

#define ov_conditionalf(mask, val1, val2)  ((mask) ? (val1) : (val2))

#define ov_allf(expr) (expr)
#define ov_anyf(expr) (expr)

#define ov_any_lt_0f(a) ((a)<0.0f)
#define ov_all_ge_0f(a) ((a)>=0.0f)
#define ov_any_ge_0f ov_all_ge_0f
#define ov_all_lt_0f ov_any_lt_0f

/*
    DOUBLE PRECISION
*/
#include"openvec.double_scalar.h"
