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
  float *x;
  int i;

  printf("OV_PLATFORM: %s\n", OV_PLATFORM);
  printf("OV_FLOAT_WIDTH %d\n", OV_FLOAT_WIDTH);

  x=(float*)malloc(sizeof(float)*(n));


  /*
     Init
  */
  for (i=0; i<n; i++) x[i]=i;


  /* Scalar version just one line */
  float sum=0.0;
  for (i=0; i<n; i++) sum += x[i];
  printf("sum = %.1f\n", sum);


  /* n is not always multiple of OV_FLOAT_WIDTH
     and we can't sum more than n elements
     The first loop will do multiple of OV_FLOAT_WIDTH elements
     The second loop will do the rest up to upperLimit
  */
  

  /* Vector version */
  float vsum=0.0f;
  i=0; /* first element */
#if OV_FLOAT_WIDTH > 1 /* In the scalar mode this LOOP won`t be compiled */
  ov_float tmp=ov_zerof;
  for (; i<n-OV_FLOAT_TAIL; i+=OV_FLOAT_WIDTH) 
  { 
    ov_float vx = ov_ldf(&x[i]);
    tmp = ov_addf(tmp, vx);
  }
  vsum = ov_all_sumf(tmp);
  printf("next i=%d\n", i);
#endif 
  
  /* Second loop goes up to upper limit 
     and is simple plain scalar code */
  for (; i<n; i++) vsum += x[i];
  printf("vsum = %.1f\n", vsum);
  /* End vector version */

  TEST(sum, vsum, "REDUCTION");

  free(x);

  printf("Passed all tests.\n");

  return 0;
}
