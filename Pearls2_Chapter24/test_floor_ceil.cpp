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

#include"math.h"
#include"openvec.h"


#include"test_utils.h"

#define DO_TEST(ov_round, round) do{\
  for (i=0; i<n; i+=OV_FLOAT_WIDTH) \
  {\
    ov_float vx = ov_ldf(&x[i]);  \
    \
    ov_float vy=ov_round(vx);\
    ov_stf(&y[i], vy);              \
  }\
  for (i=1; i<n-1; i++) TEST(round(x[i]), y[i], "ROUND");\
}while(0)

int main(int argc, char *argv[])
{
  int const n=1024;
  float *x, *y, *z;
  int i;

  printf("OV_PLATFORM: %s\n", OV_PLATFORM);
  printf("OV_FLOAT_WIDTH %d\n", OV_FLOAT_WIDTH);

  x=(float*)malloc(sizeof(float)*(n));
  y=(float*)calloc(n, sizeof(float));
  z=(float*)calloc(n, sizeof(float));


  /*
     Init
  */
//  float f=8388608.0f;
//  float f=-2147483647.0f;
  float f=-1.0e20;
  int k=*(int*)&f;
  printf("k=%d, f=%f\n", k, f);
  while(fabs(f)>0.4f)
  {
  //printf("k=%d, f=%f\n", k, f);
  for (i=0; i<n; i++)
  {
  //  printf("k=%d, f=%f\n", k, f);
    k--;
    f=*(float*)&k;
    x[i]=f;
  }
  
  //for (i=0; i<n; i++) printf("x[%d]=%f\n",i,x[i]);

  DO_TEST(ov_floorf, floorf);
  DO_TEST(ov_ceilf, ceilf);

  if ((f<0.0f) && (f>-.5))
  {
    f=2147483647.0f;
    k=*(int*)&f;
  }
  }
  free(x);
  free(y);
  free(z);

  printf("Passed all tests.\n");

  return 0;
}
