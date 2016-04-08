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
#include<string.h>

#include"math.h"
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

void memcopy_float(float *dst, float const *src, int const upperLimit)
{
  /* upperLimit is not always multiple of OV_FLOAT_WIDTH
     and we can't fill more than upperLimit elements
     The first loop will do multiple of OV_FLOAT_WIDTH elements
     The second loop will do the rest up to upperLimit (tail)
  */
  

  int i=0; /* first element */

#if OV_FLOAT_WIDTH > 1 /* Loop below won't be compiled in the scalar backend */
  for (; i<upperLimit-OV_FLOAT_TAIL; i+=OV_FLOAT_WIDTH) 
  { 
    ov_float vx = ov_ldf(&src[i]);
    ov_stf(&dst[i], vx); 
  } 
#endif 

  /* Second loop goes up to upper limit 
     and is simple plain scalar code */
  for (; i<upperLimit; i++) dst[i]=src[i];
  /* End vector version */
}

int main(int argc, char *argv[])
{
  int const n=40;
  float *x, *y, *z, *w, *r;
	int i;

  printf("OV_PLATFORM: %s\n", OV_PLATFORM);
  printf("OV_FLOAT_WIDTH %d\n", OV_FLOAT_WIDTH);

  x=(float*)malloc(sizeof(float)*(n));
  y=(float*)malloc(sizeof(float)*(n));
  z=(float*)malloc(sizeof(float)*(n));
  w=(float*)malloc(sizeof(float)*(n));
  r=(float*)malloc(sizeof(float)*(n));


  /*
     GET ZERO
  */
  for (i=0; i<n; i++) x[i]=i;
  for (i=0; i<n; i++) y[i]=0;
  for (i=0; i<n; i++) r[i]=0;

  int const upperLimit = n-7;

  /* Scalar version */
  memcpy(r, x, upperLimit*sizeof(float));

  /* Vector version */
  memcopy_float(y, x, upperLimit);


  /* Test for correctness */
  for (i=1; i<n; i++) TEST(r[i], y[i], "OV_GETZERO"); 

  free(x);
  free(y);
  free(z);
  free(w);
  free(r);

  printf("Passed all tests.\n");

  return 0;
}
