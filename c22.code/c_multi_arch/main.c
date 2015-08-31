/*
  This file is provided under a BSD license.

  Copyright (c) 2015 Intel Corporation. All rights reserved.
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the
      distribution.
    * Neither the name of Intel Corporation nor the names of its
      contributors may be used to endorse or promote products derived
      from this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define FUNC(x) exp(x)

#include <get_time.h>

#include "interpolate.h"

const int steps = 512;
const int ARRAY_SIZE = 2048;

int main(int argc, char* argv[])
{
    point* vals = malloc(sizeof(point)*(steps+1));
    double src[ARRAY_SIZE] __attribute__((aligned(128)));
    double dst[ARRAY_SIZE] __attribute__((aligned(128)));
	
	double prev_val = 0.;
	int ind = 0;
	int i = 0;
	
    /* Set interpolation array */
    double delta = 1. / steps;
    double delta_inv = 1. / delta;

    vals[0].c0 = 0.;
    vals[0].c1 = FUNC(0.);

    prev_val = vals[0].c1;
    /* Fill interpolation array */
    for(ind=1; ind <= steps; ++ind) {
        double x = ind * delta;
        double val = FUNC(x);
        double c0 = (val - prev_val)*delta_inv;
        double c1 = val - c0*x;
        vals[ind].c0 = c0;
        vals[ind].c1 = c1;
        prev_val = val;
    }

    /* Initialize input array */
    for(i=0;i<ARRAY_SIZE;++i) {
        src[i]=((double)rand())/RAND_MAX;
    }

    START_LOOP_TIME()
    /*Critical loop requires vectorization */   
#ifdef _OPENMP
    #pragma omp simd
#endif
    for(i=0; i<ARRAY_SIZE;++i) {
        dst[i] = Interpolate(src[i],vals);
    }
    END_LOOP_TIME()

    /* Test results */   
    for(i=0;i<ARRAY_SIZE;++i) {
        double ref_val = FUNC(src[i]);

        if ( (ref_val-dst[i])/ref_val > 0.01 ) {
            printf("Error at %i, val %lf, got %lf, expected %lf\n", i, src[i], dst[i], ref_val);
	    return -1;
        }
    }

    free(vals);
    return 0;
}

