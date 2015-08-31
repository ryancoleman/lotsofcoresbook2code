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


#define FMT_IN "%4.1f"
#define FMT_OUT "%7.3f"


int main(int argc, char *argv[])
{
  int const n=177;
  float *x, *y, *simd_result1, *simd_result2, *scalar_result;
	int i;

  printf("OV_PLATFORM: %s\n", OV_PLATFORM);
  printf("OV_FLOAT_WIDTH %d\n", OV_FLOAT_WIDTH);

  /*
     Data allocation
  */
  x=(float*)malloc(sizeof(float)*(n));
  y=(float*)malloc(sizeof(float)*(n));
  simd_result1=(float*)malloc(sizeof(float)*(n));
  simd_result2=(float*)malloc(sizeof(float)*(n));
  scalar_result=(float*)malloc(sizeof(float)*(n));

  /*
     Data initialization
  */
  for (i=0; i<n; i++) x[i]=2*i;
  for (i=0; i<n; i++) y[i]=i+1;


  /*
    Simple if else scalar construct
  */
  for (i=0; i<n; i++)
  {
    if (x[i] > y[i]) scalar_result[i] = x[i]-y[i];
    else             scalar_result[i] = y[i]+x[i];
  }



  /*
    Masked equivalent if else SIMD construct
  */
  for (i=0; i<n; i+=OV_FLOAT_WIDTH) /* for every SIMD element */
  { 
    /* Load x and y vectors */
    ov_float const vx=ov_ldf(&x[i]);
    ov_float const vy=ov_ldf(&y[i]);
     

    /* Merge two vector depending on a comparisson expression 
		 *ov_conditionalf(vx > vy, ov_subf(vx,vy), ov_addf(vy,vx));
		 */
    ov_float resul = ov_conditionalf(ov_gtf(vx,vy), ov_subf(vx,vy), ov_addf(vy,vx));

    /* Store the result */
    ov_stf(&simd_result1[i], resul);
  } 




  /*
    Optimized masked equivalent if else SIMD construct
    Apply the mask only on the corner cases
    In this code sample, only in the first two floats x > y is not true
    (you must know the corner cases of your code)
    (you can also split the loop trip count to avoid masked instructions)
  */
  for (i=0; i<n; i+=OV_FLOAT_WIDTH) /* for every SIMD element */
  { 
    /* Load x and y vectors */
    ov_float const vx=ov_ldf(&x[i]);
    ov_float const vy=ov_ldf(&y[i]);
 
    ov_float result;

    /*
       Optimization: if on all SIMD elements x[] > y[]
       Then there's no divergent paths and no need for MASKED operations
    */
    if (ov_allf(ov_gtf(vx,vy))) result = ov_subf(vx,vy); /* No MASK involved */
    else                result = ov_conditionalf(ov_gtf(vx,vy), ov_subf(vx,vy), ov_addf(vy,vx));

    ov_stf(&simd_result2[i], result);
  } 


  for (i=0; i<n; i++) printf(" MASKED VEC " FMT_OUT
                                 " SCALAR " FMT_OUT "\n",
                                 simd_result2[i], scalar_result[i]); 

  /* Check for correctness */
  for (i=1; i<n; i++) TEST(simd_result1[i], scalar_result[i], "Greater Than"); 
  for (i=1; i<n; i++) TEST(simd_result2[i], scalar_result[i], "Optimized greater Than"); 

  free(x);
  free(y);
  free(simd_result1);
  free(simd_result2);
  free(scalar_result);

  printf("Passed all tests.\n");

  return 0;
}
