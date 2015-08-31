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

#include"openvec.h"


int main(int argc, char *argv[])
{
  int const n=40;
  float *x, *y, *z;

  printf("OV_PLATFORM: %s\n", OV_PLATFORM);
  printf("OV_FLOAT_WIDTH %d\n", OV_FLOAT_WIDTH);

  x=(float*)malloc(sizeof(float)*(n));
  y=(float*)malloc(sizeof(float)*(n));
  z=(float*)malloc(sizeof(float)*(n));

  for (int i=0; i<n; i++)
  {
    y[i]=0.0f;
    if (i<n/2) x[i]=-i-i/10.0f;
    else       x[i]= i+i/10.0f;
  }

  for (int i=0; i<n; i+=OV_FLOAT_WIDTH)
  {
    ov_float vx = ov_loadf(&x[i]);
    ov_float vy = floorf(vx);
    ov_storef(&y[i], vy);
    ov_float vz = ceilf(vx);
    ov_storef(&z[i], vz);
  }
  for (int i=0; i<n; i++) printf("FLOOR(%6.1f) = %6.1f "
                                 "CEIL(%6.1f) = %6.1f\n",
                                 x[i], y[i], x[i], z[i]);

  free(x);
  free(y);
  free(z);

  return 0;
}
