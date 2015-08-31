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
    DOUBLE PRECISION
*/
#define OV_DOUBLE_WIDTH 1

#define ov_prv_double double
#define ov_addd(a,b) ((a)+(b))
#define ov_subd(a,b) ((a)-(b))
#define ov_muld(a,b) ((a)*(b))
#define ov_divd(a,b) ((a)/(b))
#define ov_loadd(a) (*(a))
#define ov_ldd(a) (*(double*)(a))
#define ov_uloadd(dest,addr) dest=ov_ldd(addr)
#define ov_uldd ov_ldd
#define ov_stored(a,b) *(a)=b
#define ov_storeud ov_stored
#define ov_std ov_storeud
#define ov_setzerod(a) a=0.0
#define ov_setd(a) (a)


#define ov_andmaskd(a,b) ((a)&(b))
#define ov_ormaskd(a,b) ((a)|(b))



#define ov_mind(X,Y) ((X) < (Y) ? (X) : (Y))
#define ov_maxd(X,Y) ((X) > (Y) ? (X) : (Y))

#define ov_sqrtd sqrt
#define ov_rcpd(x) (1.0/(x))
#define ov_rsqrtd(x) (1.0/ov_sqrtd(x))
#define ov_fastsqrtd sqrt
#define ov_absd abs

/*double ov_casti_d(int i){return *(double*)&i;}*/
#define ov_casti_d (double &)

#define ov_double_to_int(d) ((int)(d))


#define ov_maddd(a,b,c) ((c)+((a)*(b)))
#define ov_msubd(a,b,c) (ov_subd(ov_muld(a,b),c))

#define ov_ceild ceil
#define ov_floord floor

#define ov_zerod (0.0)
#define ov_getzerod() (0.0)

#define ov_gtd(X,Y) ((X) >  (Y))
#define ov_ged(X,Y) ((X) >= (Y))

#define ov_ltd(X,Y) ((X) <  (Y))
#define ov_led(X,Y) ((X) <= (Y))

#define ov_eqd(X,Y) ((X) == (Y))

#define ov_conditionald(mask, val1, val2)  ((mask) ? (val1) : (val2))

#define ov_alld(expr) (expr)
#define ov_anyd(expr) (expr)

#define ov_any_lt_0d(a) ((a)<0.0)
#define ov_all_ge_0d(a) ((a)>=0.0)
#define ov_any_ge_0d ov_all_ge_0d
#define ov_all_lt_0d ov_any_lt_0d
