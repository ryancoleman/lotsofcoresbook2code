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


#include <sys/time.h>

#include"openvec.h"

#define UNROLL 4

double wallt()
{
  struct timeval tv;
  /*
  struct timezone tz;
  gettimeofday(&tv, &tz);
  */
  gettimeofday(&tv, (struct timezone*)0);
  return tv.tv_sec + tv.tv_usec*1e-6;
}



void saxpy(int n, float a, float *x, float *y)
{
  for (int i=0; i<n; i+=OV_FLOAT_WIDTH)
  {
    ov_float const va = ov_setf(a);
    ov_float vy = ov_loadf(&y[i]);
    ov_float vx = ov_loadf(&x[i]);
    vy = ov_maddf(va,vx,vy);
    ov_storef(&y[i], vy);
  }
}


void vsaxpy(int nv, float a, const ov_float* __restrict__ x, ov_float* __restrict__ y)
{
  ov_float const va = ov_setf(a);
//  for (int i=0; i<nv; i++)
  for (int i=0; i<nv; i+=UNROLL)
  {
#ifdef __cplusplus
    y[i] = va*x[i] + y[i];
    y[i+1] = va*x[i+1] + y[i+1];
    y[i+2] = va*x[i+2] + y[i+2];
    y[i+3] = va*x[i+3] + y[i+3];
#else
    y[i] = ov_maddf(va, x[i], y[i]);
    y[i+1] = ov_maddf(va, x[i+1], y[i+1]);
    y[i+2] = ov_maddf(va, x[i+2], y[i+2]);
    y[i+3] = ov_maddf(va, x[i+3], y[i+3]);
#endif
  }
}

void do_test(int n, const float* x, float* y)
{
  int const nv = ((n-1)/OV_FLOAT_WIDTH) + 1;
 // saxpy(n, 1.0f, x, y); 
  int const niter=512*(100000000/n);
//  int const niter=10000;
  double t0=wallt();
  for (int i=0; i<niter; i++) 
  {
    vsaxpy(nv, 1.0f, (ov_float*)x, (ov_float*)y); 
  }
  double t1=wallt();
  double samples=niter*(n*2.0)/(t1-t0);

  printf("CSV, " OV_LANG " " OV_PLATFORM ", %d, %.2lf\n", n, samples/1000000000.0);
}
int main(int argc, char *argv[])
{
  int const n=1<<19;
//  int const n=1<<18;
  float *x, *y;

  printf("OV_PLATFORM: " OV_LANG " " OV_PLATFORM "\n");
  printf("OV_FLOAT_WIDTH %d\n", OV_FLOAT_WIDTH);

  x=(float*)malloc(sizeof(float)*n);
  y=(float*)malloc(sizeof(float)*n);

  for (int i=0; i<n; i++)
  {
    y[i]=i;
    x[i]=i*100.0f;
  }

  printf("CSV, plat, size, GFlops/s\n");
  do_test(1<<9, x, y);
  do_test(1<<14, x, y);
  do_test(n, x, y);

  free(x);
  free(y);

  return 0;
}
