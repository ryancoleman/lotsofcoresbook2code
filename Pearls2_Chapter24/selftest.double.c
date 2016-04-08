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
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include"openvec.h"
#include"test_utils.h"


#define FMT_IN "%4.1lf"
#define FMT_OUT "%7.3lf"


#define DO_TEST(init, vectorfunc, scalarfunc, funcname) do{\
  for (i=0; i<n; i++) x[i]=(double)init; \
  for (i=0; i<n; i+=OV_DOUBLE_WIDTH) \
  { \
    ov_double vx = ov_ldd(&x[i]); \
    ov_double vy = vectorfunc(vx); \
    ov_stored(&y[i], vy); \
  } \
  for (i=0; i<n; i++) r[i]=scalarfunc; \
\
  for (i=0; i<n; i++) printf(funcname "(" FMT_IN ") = " FMT_OUT " " FMT_OUT"\n", \
                                 x[i], y[i], r[i]); \
  for (i=1; i<n; i++) TEST(r[i], y[i], funcname); \
 \
 }while(0)



#define DO_COND_TEST(op, op_symb,  funcname) do{\
  for (i=0; i<n; i++) x[i]=i-n/2;\
  for (i=0; i<n; i++) y[i]=1+i/10;\
  for (i=0; i<n; i+=OV_DOUBLE_WIDTH) \
  { \
    ov_double vx=ov_ldd(&x[i]);\
    ov_double vy=ov_ldd(&y[i]);\
    ov_double vw = ov_conditionald(op(vx, vy), vx, ov_addd(vy,vx));\
    ov_std(&w[i], vw);\
  } \
\
  for (i=0; i<n; i++)\
  {\
    printf("x=%lf y=%lf\n", x[i],y[i]);\
    if (x[i] op_symb y[i]) r[i]=x[i];\
    else           r[i]=y[i]+x[i];\
  }\
\
  for (i=0; i<n; i++) printf(funcname " MASKED VEC " FMT_OUT\
                                 " SCALAR " FMT_OUT "\n",\
                                 w[i], r[i]); \
\
  for (i=1; i<n; i++) TEST(w[i], r[i], funcname); \
}while(0)


#define DO_COND_TESTCPP(op_symb,  funcname) do{\
  for (i=0; i<n; i++) x[i]=i-n/2;\
  for (i=0; i<n; i++) y[i]=1+i/10;\
  for (i=0; i<n; i+=OV_DOUBLE_WIDTH) \
  { \
    ov_double vx=ov_ldd(&x[i]);\
    ov_double vy=ov_ldd(&y[i]);\
    ov_double vw = ov_conditionald((vx) op_symb (vy), vx, ov_addd(vy,vx));\
    ov_std(&w[i], vw);\
  } \
\
  for (i=0; i<n; i++)\
  {\
    printf("x=%lf y=%lf\n", x[i],y[i]);\
    if (x[i] op_symb y[i]) r[i]=x[i];\
    else           r[i]=y[i]+x[i];\
  }\
\
  for (i=0; i<n; i++) printf(funcname " MASKED VEC " FMT_OUT\
                                 " SCALAR " FMT_OUT "\n",\
                                 w[i], r[i]); \
\
  for (i=1; i<n; i++) TEST(w[i], r[i], funcname); \
}while(0)


#define DO_COND_EXPR(ov_expr, sca_expr, funcname) do{\
  for (i=0; i<n; i++) x[i]=i-n/2;\
  for (i=0; i<n; i++) y[i]=1+i/10;\
  for (i=0; i<n; i+=OV_DOUBLE_WIDTH) \
  { \
    ov_double vx=ov_ldd(&x[i]);\
    ov_double vy=ov_ldd(&y[i]);\
    ov_double vw = ov_conditionald(ov_expr, vx, ov_addd(vy,vx));\
    ov_std(&w[i], vw);\
  } \
\
  for (i=0; i<n; i++)\
  {\
    if (sca_expr) r[i]=x[i];\
    else           r[i]=y[i]+x[i];\
    printf(funcname ": x=%f y=%f -> r=%f\n", x[i],y[i],r[i]);\
  }\
\
  for (i=0; i<n; i++) printf(funcname " MASKED VEC " FMT_OUT\
                                 " SCALAR " FMT_OUT "\n",\
                                 w[i], r[i]); \
\
  for (i=1; i<n; i++) TEST(w[i], r[i], funcname); \
}while(0)



  /*
     Coditional Reduction Functions
  */
#define DO_COND_REDUCTION_TEST(init, red_func, op, op_symb, funcname) do{\
  for (i=0; i<n+OV_DOUBLE_TAIL; i++) x[i]=init;\
  for (i=0; i<n+OV_DOUBLE_TAIL; i++) y[i]=i;\
  for (i=0; i<n+OV_DOUBLE_TAIL; i++) r[i]=-1.0;\
  for (i=0; i<n+OV_DOUBLE_TAIL; i++) w[i]=-1.0;\
\
  for (i=0; i<n; i++) \
  {\
    double vx=x[i];\
    double vy=y[i];\
    double vw=-2.0;\
    if (vx op_symb vy) vw=vx;\
    r[i]=vw;\
  }\
\
  for (i=0; i<n; i+=OV_DOUBLE_WIDTH) \
  {\
    ov_double vx=ov_ldd(&x[i]);\
    ov_double vy=ov_ldd(&y[i]);\
    ov_double vw=ov_setd(-2.0);\
    if (red_func(op(vx, vy))) vw=vx;\
    ov_std(&w[i], vw); \
  }\
\
  for (i=0; i<n; i++) printf("REDUCTION " funcname " VEC " FMT_OUT\
                                 " SCALAR " FMT_OUT "\n",\
                                 w[i], r[i]); \
\
  for (i=1; i<n; i++) TEST(w[i], r[i], funcname " REDUCTION");}while(0)

int main(int argc, char *argv[])
{
  int const n=40;
  ov_double mysizeof;
  double *x, *y, *z, *w, *r;
	int i;

  printf("OV_PLATFORM: %s\n", OV_PLATFORM);
  printf("OV_DOUBLE_WIDTH %d\n", OV_DOUBLE_WIDTH);

  if(sizeof(mysizeof) != (OV_DOUBLE_WIDTH*sizeof(double)))
  {
    printf("sizeof error.\n");
    exit(1);
  }

  x=(double*)malloc(sizeof(double)*(n));
  y=(double*)malloc(sizeof(double)*(n));
  z=(double*)malloc(sizeof(double)*(n));
  w=(double*)malloc(sizeof(double)*(n));
  r=(double*)malloc(sizeof(double)*(n));


  /*
     Square Root
  */
  DO_TEST(i, ov_sqrtd, sqrt(x[i]), "SQRT");

  /*
     Fast square Root
  */
  DO_TEST(i, ov_fastsqrtd, sqrt(x[i]), "FAST SQRT");

  /*
     Reciprocal Square Root
  */
  DO_TEST(i, ov_rsqrtd, 1.0/sqrt(x[i]), "RECIPROCAL SQRT");

  /*
     Reciprocal
  */
  DO_TEST(i, ov_rcpd, 1.0/x[i], "RECIPROCAL");

  /*
     Floor
  */
  DO_TEST((2*i -n)/8.0, ov_floord, floor(x[i]), "FLOOR");

  /*
     Ceil
  */
  DO_TEST(i/10.0 -n/20.0, ov_ceild, ceil(x[i]), "CEIL");

  /*
     Abs
  */
  DO_TEST(i-n/2, ov_absd, abs(x[i]), "ABS");


  /*
     GET ZERO
  */
  for (i=0; i<n; i++) y[i]=1.0;
  for (i=0; i<n; i+=OV_DOUBLE_WIDTH) 
  { 
    ov_double vy = ov_getzerod();
    ov_std(&y[i], vy); 
  } 
  for (i=1; i<n; i++) TEST(y[i], 0.0, "OV_GETZERO"); 


  /*
     Masked Ops
  */
  DO_COND_TEST(ov_ltd, <,  "LT");
  DO_COND_TEST(ov_led, <=, "LE");
  DO_COND_TEST(ov_gtd, >,  "GT");
  DO_COND_TEST(ov_ged, >=, "GE");
  DO_COND_TEST(ov_eqd, ==, "GE");

  /*
     Masked Ops with C++ operators
  */
#ifdef __cplusplus
  DO_COND_TESTCPP(<,  "LT");
  DO_COND_TESTCPP(<=, "LE");
  DO_COND_TESTCPP(>,  "GT");
  DO_COND_TESTCPP(>=, "GE");
  DO_COND_TESTCPP(==, "GE");

 /*
     Logical operators on MASKS
  */
  DO_COND_EXPR(ov_andmaskd((vx>vy), (vx> 10.0)), (x[i]>y[i]) && (x[i]> 10.0), "OV_AND");
  DO_COND_EXPR( ov_ormaskd((vx>vy), (vx>-10.0)), (x[i]>y[i]) || (x[i]>-10.0), "OV_OR");

#endif

  /*
     Logical operators on MASKS
  */
  DO_COND_EXPR(ov_andmaskd(ov_gtd(vx,vy), ov_gtd(vx,ov_setd(10.0))), (x[i]>y[i]) && (x[i]> 10.0), "OV_AND");
  DO_COND_EXPR( ov_ormaskd(ov_gtd(vx,vy), ov_gtd(vx,ov_setd(-10.0))), (x[i]>y[i]) || (x[i]>-10.0), "OV_OR");

  /*
     Coditional Reduction Functions
  */
  DO_COND_REDUCTION_TEST(i+1,   ov_alld, ov_gtd, >,  "ALL GT");
  DO_COND_REDUCTION_TEST(i+i%2, ov_alld, ov_ged, >=, "ALL GE");
  DO_COND_REDUCTION_TEST(i  ,   ov_alld, ov_eqd, ==, "ALL EQ");
  DO_COND_REDUCTION_TEST(i-1,   ov_alld, ov_ltd, <,  "ALL LT");
  DO_COND_REDUCTION_TEST(i-i%2, ov_alld, ov_led, <=, "ALL LE");

  DO_COND_REDUCTION_TEST(i+1,   ov_anyd, ov_gtd, >,  "ANY GT");
  DO_COND_REDUCTION_TEST(i+i%2, ov_anyd, ov_ged, >=, "ANY GE");
  DO_COND_REDUCTION_TEST(i  ,   ov_anyd, ov_eqd, ==, "ANY EQ");
  DO_COND_REDUCTION_TEST(i-1,   ov_anyd, ov_ltd, <,  "ANY LT");
  DO_COND_REDUCTION_TEST(i-i%2, ov_anyd, ov_led, <=, "ANY LE");


  /*
     MATH
  */
  int const ihig=n-2;
  for (i=0; i<n; i++) x[i]=i-n/2;
  for (i=0; i<n; i++) y[i]=10.0*(i+1);
  for (i=0; i<n; i++) z[i]=i*0.5;
  for (i=0; i<ihig; i++) 
  {
    double vx=x[i];
    double vy=y[i+1];
    double vz=z[i+2];
    vz=(vx*vy)+vz;
    double vw=(vz+vx)/vy;
    vw += 100.0;
    vw = vw - vx;
    vw = vw * vx;
    vw = vx * vy -vw;
    vw = MAX(vw, -250.0);
    vw = MIN(vw, 5000.0);
    r[i]=vw;
  }

  for (i=0; i<n; i++) x[i]=i-n/2;
  for (i=0; i<n; i++) y[i]=10.0*(i+1);
  for (i=0; i<n; i++) z[i]=i*0.5;

#ifdef __cplusplus
  for (i=0; i<n; i++) w[i]=9999.9;
  for (i=0; i<ihig; i+=OV_DOUBLE_WIDTH) 
  {
    ov_double vx=ov_ldd(&x[i]);
    ov_double vy=ov_uldd(&y[i+1]);
    ov_double vz=ov_uldd(&z[i+2]);
    vz=(vx*vy)+vz;
    ov_double vw=(vz+vx)/vy;
    vw += 100.0;
    vw = 1.0*vw - vx*1.0;
    vw = vw * (vx/1.0);
    vw = vx * vy -vw;
    vw = ov_maxd(vw, -250.0);
    vw = ov_mind(vw, 5000.0);
   
    ov_std(&w[i], vw); 
  }

  for (i=0; i<ihig; i++) printf("MATH C++ VEC " FMT_OUT
                                 " SCALAR " FMT_OUT "\n",
                                 w[i], r[i]); \

  for (i=1; i<ihig; i++) TEST(w[i], r[i], "MATH C++");

#endif

  for (i=0; i<n; i++) w[i]=9999.9;
  for (i=0; i<ihig; i+=OV_DOUBLE_WIDTH) 
  {
    ov_double vx=ov_ldd(&x[i]);
    ov_double vy=ov_uldd(&y[i+1]);
    ov_double vz=ov_uldd(&z[i+2]);
    vz=ov_maddd(vx,vy,vz);
    ov_double vw=ov_divd(ov_addd(vz, vx), vy);
    vw = ov_addd(vw, ov_setd(100.0));
    vw = ov_subd(vw, vx);
    vw = ov_muld(vw, vx);
    vw = ov_msubd(vx, vy, vw);
    vw = ov_maxd(vw, ov_setd(-250.0));
    vw = ov_mind(vw, ov_setd(5000.0));
    ov_std(&w[i], vw); 
  }


  for (i=0; i<ihig; i++) printf("MATH VEC " FMT_OUT
                                 " SCALAR " FMT_OUT "\n",
                                 w[i], r[i]); \

  for (i=1; i<ihig; i++) TEST(w[i], r[i], "MATH");


  free(x);
  free(y);
  free(z);
  free(w);
  free(r);

  printf("Passed all tests.\n");

  return 0;
}
