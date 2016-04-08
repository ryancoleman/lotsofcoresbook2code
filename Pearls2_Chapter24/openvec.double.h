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
		DO NOT CALL DIRECTLY ANY FUNCTION WITH PREFIX ov_prv
		THESE FUNCTIONS ARE "PRIVATE" FOR INTERNAL USE ONLY.
*/
#ifndef _OPENVEC_DOUBLE_
#define _OPENVEC_DOUBLE_


#define OV_EPSD 1.0e-40


#define OV_SIZEOF_DOUBLE 8


#define ov_sqrd(x) ov_muld(x,x)
#define ov_prv_sqrd(x) ov_prv_muld(x,x)






/*******************************************************************************
                            C++ Type Override 
*******************************************************************************/
#if defined(__cplusplus) && defined(OV_CPP_OVERRIDE_TYPE) && (OV_DOUBLE_WIDTH>1)
struct __attribute__((aligned (OV_ALIGN))) ov_double 
{  
  ov_prv_double vreg;
};
OV_INLINE ov_prv_double ov_prv_vreg(ov_double const OV_REF(a))
{
  return OV_VREG(a);
}

#else
#define ov_double ov_prv_double
#endif



/*******************************************************************************
                                   DEFAULTS
*******************************************************************************/











#if OV_DOUBLE_WIDTH > 1 /* { */

#define ov_zerod ov_getzerod()

#ifdef __cplusplus /* { */
#ifdef OV_CPP_OVERRIDE_TYPE /* { */

static inline void ov_std(double *addr, ov_double const OV_REF(a))
{
  ov_prv_std(addr, OV_VREG(a));
}

static inline void ov_stored(double *addr, ov_double const OV_REF(a))
{
  ov_prv_stored(addr, OV_VREG(a));
}
static inline void ov_storeud(double *addr, ov_double const OV_REF(a))
{
  ov_prv_storeud(addr, OV_VREG(a));
}
static inline ov_double ov_addd(ov_double const OV_REF(a), ov_double const OV_REF(b))
{
  ov_double tmp;
  OV_VREG(tmp) = ov_prv_addd(OV_VREG(a), OV_VREG(b));
  return tmp;
}
static inline ov_double ov_subd(ov_double const OV_REF(a), ov_double const OV_REF(b))
{
  ov_double tmp;
  OV_VREG(tmp) = ov_prv_subd(OV_VREG(a), OV_VREG(b));
  return tmp;
}
static inline ov_double ov_muld(ov_double const OV_REF(a), ov_double const OV_REF(b))
{
  ov_double tmp;
  OV_VREG(tmp) = ov_prv_muld(OV_VREG(a), OV_VREG(b));
  return tmp;
}
static inline ov_double ov_divd(ov_double const OV_REF(a), ov_double const OV_REF(b))
{
  ov_double tmp;
  OV_VREG(tmp) = ov_prv_divd(OV_VREG(a), OV_VREG(b));
  return tmp;
}
static inline ov_double ov_getzerod()
{
  ov_double tmp;
  ov_prv_setzerod(OV_VREG(tmp));
  return tmp;
}
static inline void ov_setzerod(ov_double OV_REF(dest))
{
  ov_prv_setzerod(OV_VREG(dest));
}

static inline ov_double ov_uldd(double const *addr)
{
  ov_double tmp;
  ov_prv_uloadd(OV_VREG(tmp), addr);
  return tmp;
}
static inline ov_double ov_ldd(double const *b)
{
  ov_double tmp;
  OV_VREG(tmp) = ov_prv_ldd(b);
  return tmp;
}
static inline ov_double ov_loadd(double const *b)
{
  ov_double tmp;
  OV_VREG(tmp) = ov_prv_loadd(b);
  return tmp;
}
static inline void ov_uloadd(ov_double OV_REF(dest), double const *addr)
{
  ov_prv_uloadd(OV_VREG(dest), addr);
}

static inline void ov_stream_std(double *addr, ov_double const OV_REF(a))
{
  ov_prv_stream_std(addr, OV_VREG(a));
}

static inline void ov_nt_stored(double *addr, ov_double const OV_REF(a))
{
  ov_prv_nt_stored(addr, OV_VREG(a));
}

static inline ov_double ov_setd(double const OV_REF(a))
{
  ov_double tmp;
  OV_VREG(tmp) = ov_prv_setd(a);
  return tmp;
}


static inline ov_double ov_rsqrtd(ov_double const OV_REF(a))
{
  ov_double tmp;
  OV_VREG(tmp) = ov_prv_rsqrtd(OV_VREG(a));
  return tmp;
}

static inline ov_double ov_floord(ov_double const OV_REF(a))
{
  ov_double tmp;
  OV_VREG(tmp) = ov_prv_floord(OV_VREG(a));
  return tmp;
}

static inline ov_double ov_ceild(ov_double const OV_REF(a))
{
  ov_double tmp;
  OV_VREG(tmp) = ov_prv_ceild(OV_VREG(a));
  return tmp;
}

static inline ov_double ov_absd(ov_double const OV_REF(a))
{
  ov_double tmp;
  OV_VREG(tmp) = ov_prv_absd(OV_VREG(a));
  return tmp;
}

static inline ov_double ov_rcpd(ov_double const OV_REF(a))
{
  ov_double tmp;
  OV_VREG(tmp) = ov_prv_rcpd(OV_VREG(a));
  return tmp;
}

static inline ov_double ov_sqrtd(ov_double const OV_REF(a))
{
  ov_double tmp;
  OV_VREG(tmp) = ov_prv_sqrtd(OV_VREG(a));
  return tmp;
}


static inline ov_double ov_conditionald(ov_maskd const OV_REF(mask), ov_double const OV_REF(a), ov_double const OV_REF(b))
{
  ov_double tmp;
  OV_VREG(tmp) = ov_prv_conditionald(mask, OV_VREG(a), OV_VREG(b));
  return tmp;
}



static inline ov_double ov_maddd(ov_double const OV_REF(a), ov_double const OV_REF(b), ov_double const OV_REF(c))
{
  ov_double tmp;
  OV_VREG(tmp) = ov_prv_maddd(OV_VREG(a), OV_VREG(b), OV_VREG(c));
  return tmp;
}


static inline ov_double ov_msubd(ov_double const OV_REF(a), ov_double const OV_REF(b), ov_double const OV_REF(c))
{
  ov_double tmp;
  OV_VREG(tmp) = ov_prv_msubd(OV_VREG(a), OV_VREG(b), OV_VREG(c));
  return tmp;
}


static inline ov_maskd ov_eqd(ov_double const OV_REF(a), ov_double const OV_REF(b))
{
  return ov_prv_eqd(OV_VREG(a), OV_VREG(b));
}
static inline ov_maskd ov_ged(ov_double const OV_REF(a), ov_double const OV_REF(b))
{
  return ov_prv_ged(OV_VREG(a), OV_VREG(b));
}
static inline ov_maskd ov_gtd(ov_double const OV_REF(a), ov_double const OV_REF(b))
{
  return ov_prv_gtd(OV_VREG(a), OV_VREG(b));
}
static inline ov_maskd ov_led(ov_double const OV_REF(a), ov_double const OV_REF(b))
{
  return ov_prv_led(OV_VREG(a), OV_VREG(b));
}
static inline ov_maskd ov_ltd(ov_double const OV_REF(a), ov_double const OV_REF(b))
{
  return ov_prv_ltd(OV_VREG(a), OV_VREG(b));
}
/*
#define ov_ged ov_prv_ged
#define ov_led ov_prv_led
#define ov_gtd ov_prv_gtd
#define ov_ltd ov_prv_ltd
*/

#else /* }{ */

#define ov_addd ov_prv_addd
#define ov_subd ov_prv_subd
#define ov_muld ov_prv_muld
#define ov_divd ov_prv_divd

#define ov_setd ov_prv_setd
#define ov_sqrtd ov_prv_sqrtd
#define ov_rsqrtd ov_prv_rsqrtd
#define ov_rcpd ov_prv_rcpd
#define ov_stored ov_prv_stored
#define ov_storeud ov_prv_storeud
#define ov_stream_std ov_prv_stream_std
#define ov_nt_stored ov_prv_nt_stored


#define ov_msubd ov_prv_msubd
#define ov_maddd ov_prv_maddd
#define ov_conditionald ov_prv_conditionald

#define ov_floord ov_prv_floord
#define ov_ceild ov_prv_ceild
#define ov_absd ov_prv_absd
#define ov_getzerod ov_prv_getzerod
#define ov_setzerod ov_prv_setzerod

#define ov_loadd ov_prv_loadd
#define ov_ldd ov_prv_ldd
#define ov_uldd ov_prv_uldd
#define ov_std ov_prv_std

/* Masked comparisom */
#define ov_ged ov_prv_ged
#define ov_led ov_prv_led
#define ov_gtd ov_prv_gtd
#define ov_ltd ov_prv_ltd
#define ov_eqd ov_prv_eqd

#endif /*}*/

/* C++ common code */

static inline ov_double sqrtd(ov_double const OV_REF(a))
{
  ov_double tmp;
  OV_VREG(tmp) = ov_prv_sqrtd(OV_VREG(a));
  return tmp;
}

static inline ov_double fabs(ov_double const OV_REF(a))
{
  ov_double tmp;
  OV_VREG(tmp) = ov_prv_absd(OV_VREG(a));
  return tmp;
}
static inline ov_double floor(ov_double const OV_REF(a))
{
  ov_double tmp;
  OV_VREG(tmp) = ov_prv_floord(OV_VREG(a));
  return tmp;
}
static inline ov_double ceil(ov_double const OV_REF(a))
{
  ov_double tmp;
  OV_VREG(tmp) = ov_prv_ceild(OV_VREG(a));
  return tmp;
}

/* MAXIMUM */
static inline ov_double ov_maxd(ov_double const OV_REF(a), ov_double const OV_REF(b))
{
  ov_double tmp;
  OV_VREG(tmp) = ov_prv_maxd(OV_VREG(a), OV_VREG(b));
  return tmp;
}
static inline ov_double ov_maxd(ov_double const OV_REF(a), double const OV_REF(b))
{
  ov_double tmp;
  OV_VREG(tmp) = ov_prv_maxd(OV_VREG(a), ov_prv_setd(b));
  return tmp;
}
static inline ov_double ov_maxd(double const OV_REF(a), ov_double const OV_REF(b))
{
  ov_double tmp;
  OV_VREG(tmp) = ov_prv_maxd(ov_prv_setd(a), OV_VREG(b));
  return tmp;
}

/* MINIMUM */
static inline ov_double ov_mind(ov_double const OV_REF(a), ov_double const OV_REF(b))
{
  ov_double tmp;
  OV_VREG(tmp) = ov_prv_mind(OV_VREG(a), OV_VREG(b));
  return tmp;
}
static inline ov_double ov_mind(ov_double const OV_REF(a), double const OV_REF(b))
{
  ov_double tmp;
  OV_VREG(tmp) = ov_prv_mind(OV_VREG(a), ov_prv_setd(b));
  return tmp;
}
static inline ov_double ov_mind(double const OV_REF(a), ov_double const OV_REF(b))
{
  ov_double tmp;
  OV_VREG(tmp) = ov_prv_mind(ov_prv_setd(a), OV_VREG(b));
  return tmp;
}
/* END MINIMUM */


/*
   C++ operators
*/
/* Conditional vector vector */
static inline ov_maskd operator==(ov_double const &a, ov_double const &b)
{
   return ov_prv_eqd(OV_VREG(a),OV_VREG(b));
}
static inline ov_maskd operator<(ov_double const &a, ov_double const &b)
{
   return ov_prv_ltd(OV_VREG(a),OV_VREG(b));
}
static inline ov_maskd operator>(ov_double const &a, ov_double const &b)
{
   return ov_prv_gtd(OV_VREG(a),OV_VREG(b));
}
static inline ov_maskd operator<=(ov_double const &a, ov_double const &b)
{
   return ov_prv_led(OV_VREG(a),OV_VREG(b));
}
static inline ov_maskd operator>=(ov_double const &a, ov_double const &b)
{
   return ov_prv_ged(OV_VREG(a),OV_VREG(b));
}



/* Conditional double vector */
static inline ov_maskd operator==(double const &a, ov_double const &b)
{
   return ov_prv_eqd(ov_prv_setd(a),OV_VREG(b));
}
static inline ov_maskd operator<(double const &a, ov_double const &b)
{
   return ov_prv_ltd(ov_prv_setd(a),OV_VREG(b));
}
static inline ov_maskd operator>(double const &a, ov_double const &b)
{
   return ov_prv_gtd(ov_prv_setd(a),OV_VREG(b));
}
static inline ov_maskd operator<=(double const &a, ov_double const &b)
{
   return ov_prv_led(ov_prv_setd(a),OV_VREG(b));
}
static inline ov_maskd operator>=(double const &a, ov_double const &b)
{
   return ov_prv_ged(ov_prv_setd(a),OV_VREG(b));
}



/* Conditional vector double */
static inline ov_maskd operator==(ov_double const &a, double const &b)
{
   return ov_prv_eqd(OV_VREG(a),ov_prv_setd(b));
}
static inline ov_maskd operator<(ov_double const &a, double const &b)
{
   return ov_prv_ltd(OV_VREG(a),ov_prv_setd(b));
}
static inline ov_maskd operator>(ov_double const &a, double const &b)
{
   return ov_prv_gtd(OV_VREG(a),ov_prv_setd(b));
}
static inline ov_maskd operator<=(ov_double const &a, double const &b)
{
   return ov_prv_led(OV_VREG(a),ov_prv_setd(b));
}
static inline ov_maskd operator>=(ov_double const &a, double const &b)
{
   return ov_prv_ged(OV_VREG(a),ov_prv_setd(b));
}


/* Add */
static inline ov_double operator+(ov_double const &a, ov_double const &b)
{
   ov_double tmp;
   OV_VREG(tmp) = ov_prv_addd(OV_VREG(a),OV_VREG(b));
   return tmp;
}
static inline ov_double operator+(ov_double const &a, double const &b)
{
   ov_double tmp;
   OV_VREG(tmp) = ov_prv_addd(ov_prv_setd(b),OV_VREG(a));
   return tmp;
}
static inline ov_double operator+(double const &a, ov_double const &b)
{
   ov_double tmp;
   OV_VREG(tmp) = ov_prv_addd(ov_prv_setd(a),OV_VREG(b));
   return tmp;
}


/* add += */
void operator+=(ov_double &a, ov_double const b)
{
   OV_VREG(a) = ov_prv_addd(OV_VREG(a), OV_VREG(b));
}
void operator+=(ov_double &a, double const b)
{
   OV_VREG(a) = ov_prv_addd(OV_VREG(a), ov_prv_setd(b));
}


/* subtraction */
static inline ov_double operator-(ov_double const &a, ov_double const &b)
{
   ov_double tmp;
   OV_VREG(tmp) = ov_prv_subd(OV_VREG(a),OV_VREG(b));
   return tmp;
}
static inline ov_double operator-(ov_double const &a, double const &b)
{
   ov_double tmp;
   OV_VREG(tmp) = ov_prv_subd(OV_VREG(a),ov_prv_setd(b));
   return tmp;
}
static inline ov_double operator-(double const &a, ov_double const &b)
{
   ov_double tmp;
   OV_VREG(tmp) = ov_prv_subd(ov_prv_setd(a), OV_VREG(b));
   return tmp;
}


/* sub -= */
void operator-=(ov_double &a, ov_double const b)
{
   OV_VREG(a) = ov_prv_subd(OV_VREG(a), OV_VREG(b));
}
void operator-=(ov_double &a, double const b)
{
   OV_VREG(a) = ov_prv_subd(OV_VREG(a), ov_prv_setd(b));
}


/* multiplication */
static inline ov_double operator*(ov_double const &a, ov_double const &b)
{
   ov_double tmp;
   OV_VREG(tmp) = ov_prv_muld(OV_VREG(a),OV_VREG(b));
   return tmp;
}
static inline ov_double operator*(ov_double const &a, double const &b)
{
   ov_double tmp;
   OV_VREG(tmp) = ov_prv_muld(OV_VREG(a),ov_prv_setd(b));
   return tmp;
}
static inline ov_double operator*(double const &a, ov_double const &b)
{
   ov_double tmp;
   OV_VREG(tmp) = ov_prv_muld(ov_prv_setd(a),OV_VREG(b));
   return tmp;
}


/* Division */
static inline ov_double operator/(ov_double const &a, ov_double const &b)
{
   ov_double tmp;
   OV_VREG(tmp) = ov_prv_divd(OV_VREG(a),OV_VREG(b));
   return tmp;
}
static inline ov_double operator/(ov_double const &a, double const &b)
{
   ov_double tmp;
   OV_VREG(tmp) = ov_prv_divd(OV_VREG(a),ov_prv_setd(b));
   return tmp;
}
static inline ov_double operator/(double const &a, ov_double const &b)
{
   ov_double tmp;
   OV_VREG(tmp) = ov_prv_divd(ov_prv_setd(a), OV_VREG(b));
   return tmp;
}


#else /* }{ NOTcplusplus */

/* C zero overhead one to one map */

#define ov_addd ov_prv_addd
#define ov_subd ov_prv_subd
#define ov_muld ov_prv_muld
#define ov_divd ov_prv_divd
#define ov_maddd ov_prv_maddd
#define ov_conditionald ov_prv_conditionald
#define ov_getzerod ov_prv_getzerod
#define ov_uloadd ov_prv_uloadd
#define ov_setd ov_prv_setd
#define ov_loadd ov_prv_loadd
#define ov_stored ov_prv_stored
#define ov_storeud ov_prv_storeud
#define ov_nt_stored ov_prv_nt_stored
#define ov_msubd ov_prv_msubd
#define ov_stream_std ov_prv_stream_std
#define ov_setzerod ov_prv_setzerod

#define ov_rsqrtd ov_prv_rsqrtd
#define ov_rcpd ov_prv_rcpd

#define ov_sqrtd ov_prv_sqrtd

#define ov_ldd ov_prv_ldd
#define ov_uldd ov_prv_uldd
#define ov_std ov_prv_std
#define ov_maxd ov_prv_maxd
#define ov_mind ov_prv_mind
#define ov_floord ov_prv_floord
#define ov_ceild ov_prv_ceild
#define ov_absd ov_prv_absd

/* Masked comparisom */
#define ov_ged ov_prv_ged
#define ov_led ov_prv_led
#define ov_gtd ov_prv_gtd
#define ov_ltd ov_prv_ltd
#define ov_eqd ov_prv_eqd

#endif /* }} */

OV_INLINE double ov_all_sumd(ov_double const OV_REF(a))
{
  double const *fp = (double*)&a;
  int i;

  double result=fp[0];

  for (i=1; i<OV_DOUBLE_WIDTH; i++) result += fp[i];
  return result;
}
OV_INLINE double ov_all_prodd(ov_double const OV_REF(a))
{
  double const *fp = (double*)&a;
  int i;

  double result=fp[0];

  for (i=1; i<OV_DOUBLE_WIDTH; i++) result *= fp[i];
  return result;
}
OV_INLINE double ov_all_maxd(ov_double const OV_REF(a))
{
  double const *fp = (double*)&a;
  int i;

  double result=fp[0];

  for (i=1; i<OV_DOUBLE_WIDTH; i++) if (fp[i]>result) result=fp[i];
  return result;
}
OV_INLINE double ov_all_mind(ov_double const OV_REF(a))
{
  double const *fp = (double*)&a;
  int i;

  double result=fp[0];

  for (i=1; i<OV_DOUBLE_WIDTH; i++) if (fp[i]<result) result=fp[i];
  return result;
}

#else

#define ov_all_sumd(x) (x)
#define ov_all_prodd(x) (x)
#define ov_all_maxd(x) (x)
#define ov_all_mind(x) (x)

#endif /* } OV_DOUBLE_WIDTH > 1 */


#ifndef ov_fastsqrtd
#define ov_fastsqrtd(x) ov_muld(x, ov_rsqrtd(ov_maxd(x, ov_setd(OV_EPS))))
#endif

#define ov_ustd ov_storeud

#define OV_DOUBLE_TAIL ((OV_DOUBLE_WIDTH)-1)





#endif /* _OPENVEC_DOUBLE_ */
