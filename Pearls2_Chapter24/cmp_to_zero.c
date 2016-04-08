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

int main(int argc, char *argv[])
{
  int const n=77;
  float *x, *y;
  int i;

  printf("OV_PLATFORM: %s\n", OV_PLATFORM);
  printf("OV_FLOAT_WIDTH %d\n", OV_FLOAT_WIDTH);

  x=(float*)malloc(sizeof(float)*(n));
  y=(float*)malloc(sizeof(float)*(n));


  /*
     Init
  */
  for (i=0; i<n; i++) x[i]=i/OV_FLOAT_WIDTH -n/(2*OV_FLOAT_WIDTH);
  for (i=0; i<n; i++)
  {
    if (i%2==0) y[i]=-(i+1);
    else y[i]=i+1;
  /*  printf("x[%d]=%f y[%d]=%f\n",i, x[i], i, y[i]); */
  }


  for (i=0; i<n-OV_FLOAT_TAIL; i+=OV_FLOAT_WIDTH) 
  { 
    ov_float vx = ov_ldf(&x[i]);
    ov_float vy = ov_ldf(&y[i]);
    if (ov_all_lt_0f(vx) && (x[i] >= 0)) FAIL(x[i], x[i], "all_lt_0f");
    if (ov_all_ge_0f(vx) && (x[i] <  0)) FAIL(x[i], x[i], "all_ge_0f");
    if (ov_any_lt_0f(vy) && (y[i] >= 0)) FAIL(y[i], y[i], "any_lt_0f");
#if OV_FLOAT_WIDTH >1
    if (ov_any_ge_0f(vy) && (y[i+1] <  0)) FAIL(y[i], y[i+1], "any_ge_0f");
    if (!ov_any_ge_0f(vy)) FAIL(y[i], y[i+1], "any_ge_0f");
    if (ov_all_ge_0f(vy)) FAIL(y[i], y[i+1], "any_ge_0f");
    if (ov_all_lt_0f(vy)) FAIL(y[i], y[i+1], "any_ge_0f");
#endif
    if (ov_any_ge_0f(vy) && ((y[i+1] <  0) && (y[i]<0))) FAIL(y[i], y[i+1], "any_ge_0f");
    if (ov_all_lt_0f(vx) && ov_any_ge_0f(vx)) FAIL(x[i], x[i], "test 5");
    if (ov_all_lt_0f(vx) && !ov_any_lt_0f(vx)) FAIL(x[i], x[i], "test 6");
    if (ov_all_ge_0f(vx) && ov_any_lt_0f(vx)) FAIL(x[i], x[i], "test 7");
    if (ov_all_ge_0f(vx) && !ov_any_ge_0f(vx)) FAIL(x[i], x[i], "test 8");

    if (ov_any_lt_0f(vy) && ov_all_ge_0f(vy)) FAIL(x[i], x[i], "test 9");
    if (ov_any_ge_0f(vy) && ov_all_lt_0f(vy)) FAIL(x[i], x[i], "test 11");
  }

  free(x);
  free(y);

  printf("Passed all tests.\n");

  return 0;
}
