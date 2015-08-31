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
#include <stdlib.h>
//#include <stdexcept>
#include <stdio.h>
#include <sys/time.h>

#include "openvec.h"

#include "stencil_common.h"


#define LOOP_BODY(OFFSET) \
        curr = ov_ldf(&U0[i+OFFSET]);\
\
        u0 = W8 * (ov_uldf(&U0[i+OFFSET+8      ]) + ov_uldf(&U0[i+OFFSET-8      ])+ \
                             ov_ldf(&U0[i+OFFSET+8*nx   ]) +  ov_ldf(&U0[i+OFFSET-8*nx   ])+ \
                             ov_ldf(&U0[i+OFFSET+8*nx_ny]) +  ov_ldf(&U0[i+OFFSET-8*nx_ny]))+\
\
                      W7 * (ov_uldf(&U0[i+OFFSET+7      ]) + ov_uldf(&U0[i+OFFSET-7      ])+ \
                             ov_ldf(&U0[i+OFFSET+7*nx   ]) +  ov_ldf(&U0[i+OFFSET-7*nx   ])+ \
                             ov_ldf(&U0[i+OFFSET+7*nx_ny]) +  ov_ldf(&U0[i+OFFSET-7*nx_ny]))+\
\
                      W6 * (ov_uldf(&U0[i+OFFSET+6      ]) + ov_uldf(&U0[i+OFFSET-6      ])+ \
                             ov_ldf(&U0[i+OFFSET+6*nx   ]) +  ov_ldf(&U0[i+OFFSET-6*nx   ])+ \
                             ov_ldf(&U0[i+OFFSET+6*nx_ny]) +  ov_ldf(&U0[i+OFFSET-6*nx_ny]))+\
\
                      W5 * (ov_uldf(&U0[i+OFFSET+5      ]) + ov_uldf(&U0[i+OFFSET-5      ])+ \
                             ov_ldf(&U0[i+OFFSET+5*nx   ]) +  ov_ldf(&U0[i+OFFSET-5*nx   ])+ \
                             ov_ldf(&U0[i+OFFSET+5*nx_ny]) +  ov_ldf(&U0[i+OFFSET-5*nx_ny]))+\
\
                      W4 * (ov_uldf(&U0[i+OFFSET+4      ]) + ov_uldf(&U0[i+OFFSET-4      ])+ \
                             ov_ldf(&U0[i+OFFSET+4*nx   ]) +  ov_ldf(&U0[i+OFFSET-4*nx   ])+ \
                             ov_ldf(&U0[i+OFFSET+4*nx_ny]) +  ov_ldf(&U0[i+OFFSET-4*nx_ny]))+\
\
                      W3 * (ov_uldf(&U0[i+OFFSET+3      ]) + ov_uldf(&U0[i+OFFSET-3      ])+ \
                             ov_ldf(&U0[i+OFFSET+3*nx   ]) +  ov_ldf(&U0[i+OFFSET-3*nx   ])+ \
                             ov_ldf(&U0[i+OFFSET+3*nx_ny]) +  ov_ldf(&U0[i+OFFSET-3*nx_ny]))+\
\
                      W2 * (ov_uldf(&U0[i+OFFSET+2      ]) + ov_uldf(&U0[i+OFFSET-2      ])+ \
                             ov_ldf(&U0[i+OFFSET+2*nx   ]) +  ov_ldf(&U0[i+OFFSET-2*nx   ])+ \
                             ov_ldf(&U0[i+OFFSET+2*nx_ny]) +  ov_ldf(&U0[i+OFFSET-2*nx_ny]))+\
\
                      W1 * (ov_uldf(&U0[i+OFFSET+1      ]) + ov_uldf(&U0[i+OFFSET-1      ])+ \
                             ov_ldf(&U0[i+OFFSET+1*nx   ]) +  ov_ldf(&U0[i+OFFSET-1*nx   ])+ \
                             ov_ldf(&U0[i+OFFSET+1*nx_ny]) +  ov_ldf(&U0[i+OFFSET-1*nx_ny]))+\
\
                      W0 * curr;\
\
        u1 = ov_ldf(&U1[i+OFFSET]);\
        p  = ov_ldf(&P[i+OFFSET]);\
\
        u1 = p*p*u0 + 2.0f*curr- u1;\
\
        ov_stf(&U1[i+OFFSET], u1)

//
//
// vector stencil
//
//
void Kernel_vector(size_t const nx,
                                                        size_t const ny,
                                                        size_t const nz,
                                                        const float* P,
                                                        const float* U0,
                                                              float* U1)
{
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
      for(size_t ix=0; ix<nx; ix+=OV_FLOAT_WIDTH*UNROLL)
      {

        ov_float curr, u0, u1, p;

        LOOP_BODY(0);
        LOOP_BODY(OV_FLOAT_WIDTH);
        LOOP_BODY(1*OV_FLOAT_WIDTH);
        LOOP_BODY(2*OV_FLOAT_WIDTH);

        i+=OV_FLOAT_WIDTH*UNROLL;
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

  //
  // malloc will use ov_malloc to allocate aligned memory on OV_ALIGN boundary
  //
  float *P  = (float*) malloc(sizebytes);
  float *U0 = (float*) malloc(sizebytes);
  float *U1 = (float*) malloc(sizebytes);

  for (size_t i=0; i<size; i++)
  {
    P[i]  = 1500.0f;
    U0[i] = 0.0f;
    U1[i] = 0.0f;
  }

  //WARMUP
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
