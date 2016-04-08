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


#define FMT_IN "%4.1f"
#define FMT_OUT "%7.3f"


#define DO_TEST(init, vectorfunc, scalarfunc, funcname) do{\
  for (i=0; i<n; i++) x[i]=init; \
  for (i=0; i<n; i+=OV_FLOAT_WIDTH) \
  { \
    ov_float vx = ov_ldf(&x[i]); \
    ov_float vy = vectorfunc(vx); \
    ov_storef(&y[i], vy); \
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
  for (i=0; i<n; i+=OV_FLOAT_WIDTH) \
  { \
    ov_float vx=ov_ldf(&x[i]);\
    ov_float vy=ov_ldf(&y[i]);\
    ov_float vw = ov_conditionalf(op(vx, vy), vx, ov_addf(vy,vx));\
    ov_stf(&w[i], vw);\
  } \
\
  for (i=0; i<n; i++)\
  {\
    printf("x=%f y=%f\n", x[i],y[i]);\
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
  for (i=0; i<n; i+=OV_FLOAT_WIDTH) \
  { \
    ov_float vx=ov_ldf(&x[i]);\
    ov_float vy=ov_ldf(&y[i]);\
    ov_float vw = ov_conditionalf((vx) op_symb (vy), vx, ov_addf(vy,vx));\
    ov_stf(&w[i], vw);\
  } \
\
  for (i=0; i<n; i++)\
  {\
    printf("x=%f y=%f\n", x[i],y[i]);\
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
  for (i=0; i<n; i+=OV_FLOAT_WIDTH) \
  { \
    ov_float vx=ov_ldf(&x[i]);\
    ov_float vy=ov_ldf(&y[i]);\
    ov_float vw = ov_conditionalf(ov_expr, vx, ov_addf(vy,vx));\
    ov_stf(&w[i], vw);\
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
  for (i=0; i<n+OV_FLOAT_TAIL; i++) x[i]=init;\
  for (i=0; i<n+OV_FLOAT_TAIL; i++) y[i]=i;\
  for (i=0; i<n+OV_FLOAT_TAIL; i++) r[i]=-1.0f;\
  for (i=0; i<n+OV_FLOAT_TAIL; i++) w[i]=-1.0f;\
\
  for (i=0; i<n; i++) \
  {\
    float vx=x[i];\
    float vy=y[i];\
    float vw=-2.0f;\
    if (vx op_symb vy) vw=vx;\
    r[i]=vw;\
  }\
\
  for (i=0; i<n; i+=OV_FLOAT_WIDTH) \
  {\
    ov_float vx=ov_ldf(&x[i]);\
    ov_float vy=ov_ldf(&y[i]);\
    ov_float vw=ov_setf(-2.0f);\
    if (red_func(op(vx, vy))) vw=vx;\
    ov_stf(&w[i], vw); \
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
  ov_float mysizeof;
  float *x, *y, *z, *w, *r;
	int i;

  printf("OV_PLATFORM: %s\n", OV_PLATFORM);
  printf("OV_FLOAT_WIDTH %d\n", OV_FLOAT_WIDTH);

  if(sizeof(mysizeof) != (OV_FLOAT_WIDTH*sizeof(float)))
  {
    printf("Sizeof error.\n");
    exit(1);
  }

  x=(float*)malloc(sizeof(float)*(n));
  y=(float*)malloc(sizeof(float)*(n));
  z=(float*)malloc(sizeof(float)*(n));
  w=(float*)malloc(sizeof(float)*(n));
  r=(float*)malloc(sizeof(float)*(n));


  /*
     Square Root
  */
  DO_TEST(i, ov_sqrtf, sqrtf(x[i]), "SQRT");

  /*
     Fast square Root
  */
  DO_TEST(i, ov_fastsqrtf, sqrtf(x[i]), "FAST SQRT");

  /*
     Reciprocal Square Root
  */
  DO_TEST(i, ov_rsqrtf, 1.0f/sqrtf(x[i]), "RECIPROCAL SQRT");

  /*
     Reciprocal
  */
  DO_TEST(i, ov_rcpf, 1.0f/x[i], "RECIPROCAL");

  /*
     Floor
  */
  DO_TEST((2*i -n)/8.0f, ov_floorf, floorf(x[i]), "FLOOR");

  /*
     Ceil
  */
  DO_TEST(i/10.0f -n/20.0f, ov_ceilf, ceilf(x[i]), "CEIL");

  /*
     Abs
  */
  DO_TEST(i-n/2, ov_absf, fabs(x[i]), "ABS");


  /*
     GET ZERO
  */
  for (i=0; i<n; i++) y[i]=1.0f;
  for (i=0; i<n; i+=OV_FLOAT_WIDTH) 
  { 
    ov_float vy = ov_getzerof();
    ov_stf(&y[i], vy); 
  } 
  for (i=1; i<n; i++) TEST(y[i], 0.0f, "OV_GETZERO"); 


  /*
     Masked Ops
  */
  DO_COND_TEST(ov_ltf, <,  "LT");
  DO_COND_TEST(ov_lef, <=, "LE");
  DO_COND_TEST(ov_gtf, >,  "GT");
  DO_COND_TEST(ov_gef, >=, "GE");
  DO_COND_TEST(ov_eqf, ==, "GE");

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
  DO_COND_EXPR(ov_andmaskf((vx>vy), (vx> 10.0f)), (x[i]>y[i]) && (x[i]> 10.0f), "OV_AND");
  DO_COND_EXPR( ov_ormaskf((vx>vy), (vx>-10.0f)), (x[i]>y[i]) || (x[i]>-10.0f), "OV_OR");

#endif

  /*
     Logical operators on MASKS
  */
  DO_COND_EXPR(ov_andmaskf(ov_gtf(vx,vy), ov_gtf(vx,ov_setf(10.0f))), (x[i]>y[i]) && (x[i]> 10.0f), "OV_AND");
  DO_COND_EXPR( ov_ormaskf(ov_gtf(vx,vy), ov_gtf(vx,ov_setf(-10.0f))), (x[i]>y[i]) || (x[i]>-10.0f), "OV_OR");

  /*
     Coditional Reduction Functions
  */
  DO_COND_REDUCTION_TEST(i+1,   ov_allf, ov_gtf, >,  "ALL GT");
  DO_COND_REDUCTION_TEST(i+i%2, ov_allf, ov_gef, >=, "ALL GE");
  DO_COND_REDUCTION_TEST(i  ,   ov_allf, ov_eqf, ==, "ALL EQ");
  DO_COND_REDUCTION_TEST(i-1,   ov_allf, ov_ltf, <,  "ALL LT");
  DO_COND_REDUCTION_TEST(i-i%2, ov_allf, ov_lef, <=, "ALL LE");

  DO_COND_REDUCTION_TEST(i+1,   ov_anyf, ov_gtf, >,  "ANY GT");
  DO_COND_REDUCTION_TEST(i+i%2, ov_anyf, ov_gef, >=, "ANY GE");
  DO_COND_REDUCTION_TEST(i  ,   ov_anyf, ov_eqf, ==, "ANY EQ");
  DO_COND_REDUCTION_TEST(i-1,   ov_anyf, ov_ltf, <,  "ANY LT");
  DO_COND_REDUCTION_TEST(i-i%2, ov_anyf, ov_lef, <=, "ANY LE");


  /*
     MATH
  */
  int const ihig=n-2;
  for (i=0; i<n; i++) x[i]=i-n/2;
  for (i=0; i<n; i++) y[i]=10.0f*(i+1);
  for (i=0; i<n; i++) z[i]=i*0.5f;
  for (i=0; i<ihig; i++) 
  {
    float vx=x[i];
    float vy=y[i+1];
    float vz=z[i+2];
    vz=(vx*vy)+vz;
    float vw=(vz+vx)/vy;
    vw += 100.0f;
    vw = vw - vx;
    vw = vw * vx;
    vw = vx * vy -vw;
    vw = MAX(vw, -250.0f);
    vw = MIN(vw, 5000.0f);
    r[i]=vw;
  }

  for (i=0; i<n; i++) x[i]=i-n/2;
  for (i=0; i<n; i++) y[i]=10.0f*(i+1);
  for (i=0; i<n; i++) z[i]=i*0.5f;

#ifdef __cplusplus
  for (i=0; i<n; i++) w[i]=9999.9f;
  for (i=0; i<ihig; i+=OV_FLOAT_WIDTH) 
  {
    ov_float vx=ov_ldf(&x[i]);
    ov_float vy=ov_uldf(&y[i+1]);
    ov_float vz=ov_uldf(&z[i+2]);
    vz=(vx*vy)+vz;
    ov_float vw=(vz+vx)/vy;
    vw += 100.0f;
    vw = 1.0f*vw - vx*1.0f;
    vw = vw * (vx/1.0f);
    vw = vx * vy -vw;
    vw = ov_maxf(vw, -250.0f);
    vw = ov_minf(vw, 5000.0f);
   
    ov_stf(&w[i], vw); 
  }

  for (i=0; i<ihig; i++) printf("MATH C++ VEC " FMT_OUT
                                 " SCALAR " FMT_OUT "\n",
                                 w[i], r[i]); \

  for (i=1; i<ihig; i++) TEST(w[i], r[i], "MATH C++");

#endif

  for (i=0; i<n; i++) w[i]=9999.9f;
  for (i=0; i<ihig; i+=OV_FLOAT_WIDTH) 
  {
    ov_float vx=ov_ldf(&x[i]);
    ov_float vy=ov_uldf(&y[i+1]);
    ov_float vz=ov_uldf(&z[i+2]);
    vz=ov_maddf(vx,vy,vz);
    ov_float vw=ov_divf(ov_addf(vz, vx), vy);
    vw = ov_addf(vw, ov_setf(100.0f));
    vw = ov_subf(vw, vx);
    vw = ov_mulf(vw, vx);
    vw = ov_msubf(vx, vy, vw);
    vw = ov_maxf(vw, ov_setf(-250.0f));
    vw = ov_minf(vw, ov_setf(5000.0f));
    ov_stf(&w[i], vw); 
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
