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
  for (i=0; i<n; i++) x[i]=i;
  

  /* Moving average of 3 points */
  for (i=1; i<n-1; i+=OV_FLOAT_WIDTH) 
  {
	  /* Since loop starts in i=1 we can't guarantee aligment on indexes i and i+1 */
	  ov_float left  = ov_ldf(&x[i-1]);  // ALIGNED
	  ov_float curr  = ov_uldf(&x[i]);   // UNALIGNED
	  ov_float right = ov_uldf(&x[i+1]); // UNALIGNED
    
    curr = (left + curr + right)/3.0f;

    ov_ustf(&y[i], curr);              // UNALIGNED
		ov_stream_stf(&z[i-1], curr);      // ALIGNED

  }
  
  for (i=1; i<n-1; i++) printf("AVG[%d]=%f\n", i, y[i]);

  for (i=1; i<n-1; i++) TEST(x[i], y[i], "UNALIGNED LD/ST");
  for (i=1; i<n-1; i++) TEST(x[i], z[i-1], "STREAM ST");

  free(x);
  free(y);
  free(z);

  printf("Passed all tests.\n");

  return 0;
}
