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
#include <stdio.h>
#include <stdlib.h>

#include <sys/time.h>

#include "openvec.h"
#include "stencil_common.h"


/*
 *
 VECTOR STENCIL
 *
 */
#define DER_BLOCK(e) \
  u0 = ov_maddf(ov_setf(coef[e]), ov_addf(ov_addf(ov_addf(ov_uldf(&U0[i+e      ]), ov_uldf(&U0[i-e      ])), \
                                  ov_addf( ov_ldf(&U0[i+e*nx   ]),  ov_ldf(&U0[i-e*nx   ]))),\
                                  ov_addf( ov_ldf(&U0[i+e*nx_ny]),  ov_ldf(&U0[i-e*nx_ny]))), u0)

void Kernel_vector(size_t const nx,
                                                        size_t const ny,
                                                        size_t const nz,
                                                        const float* P,
                                                        const float* U0,
                                                              float* U1)
{
  float const coef [] = { W0, W1, W2, W3, W4, W5, W6 , W7, W8};

  size_t const nx_ny = nx*ny;

#pragma omp parallel for collapse(2)
  for(size_t ibz=NPOP; ibz<(nz-NPOP); ibz+=BSIZEZ)
  for(size_t iby=0;    iby<ny;        iby+=BSIZEY)
  {
  
  for(size_t iz=ibz; iz<MIN(ibz+BSIZEZ, nz-NPOP); iz++)
  {
    size_t ind = iz*nx_ny;
    for(size_t iy=iby; iy<MIN(iby+BSIZEY, ny); iy++)
    {
      size_t i=ind+iy*nx;
      for(size_t ix=0; ix<nx; ix+=OV_FLOAT_WIDTH)
      {

        ov_float curr = ov_ldf(&U0[i]);

        ov_float u0 = ov_zerof;
   
        DER_BLOCK(8);
        DER_BLOCK(7);
        DER_BLOCK(6);
        DER_BLOCK(5);
        DER_BLOCK(4);
        DER_BLOCK(3);
        DER_BLOCK(3);
        DER_BLOCK(2);
        DER_BLOCK(1);

        u0 = ov_maddf(ov_setf(coef[0]), curr, u0);

        ov_float u1 = ov_ldf(&U1[i]);
        ov_float p  = ov_ldf(&P[i]);

        u1 = ov_addf(ov_mulf(ov_mulf(p,p),u0), ov_msubf(ov_setf(2.0f), curr, u1));

        ov_stf(&U1[i], u1);

        i+=OV_FLOAT_WIDTH;
      }
    }
  }
  }

}



int main(int argc, char *argv[])
{
  printf("OV_PLATFORM: %s\n", OV_PLATFORM);
  printf("OV_FLOAT_WIDTH %d\n", OV_FLOAT_WIDTH);
  printf("OV_ALIGN %d\n", OV_ALIGN);

  ARRAY_SIZES;

  size_t const size = n1*n2*n3;
  size_t const sizebytes = size*sizeof(float);

  /*
     malloc will use ov_malloc to allocate aligned memory on OV_ALIGN boundary
  */
  float *P  = (float*) malloc(sizebytes);
  float *U0 = (float*) malloc(sizebytes);
  float *U1 = (float*) malloc(sizebytes);

  for (size_t i=0; i<size; i++)
  {
    P[i]  = 1500.0f;
    U0[i] = 0.0f;
    U1[i] = 0.0f;
  }

  /* WARMUP */
  Kernel_vector(n1, n2, n3, P, U0, U1);
  Kernel_vector(n1, n2, n3, P, U1, U0);

  double t0=wallt();
  int const niter=1;
  for(int iter=0; iter<niter; iter++)
  {
    Kernel_vector(n1, n2, n3, P, U0, U1);
    Kernel_vector(n1, n2, n3, P, U1, U0);
  }
  double t1=wallt();
  double msamples = niter*2l*(size/1000000.0);

  printf("%.1lf Msamples/s\n", msamples/(t1-t0));

  free(U1);
  free(U0);
  free(P);

  printf("Success.\n");
  return 0;
}
