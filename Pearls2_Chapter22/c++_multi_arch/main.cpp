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

#include<math.h>
#include<stdlib.h>
#include<vector>
#include<iostream>

#ifdef __INTEL_COMPILER
#include<tbb/cache_aligned_allocator.h>
typedef std::vector<double, tbb::cache_aligned_allocator<double> > data_array_t;
#else
typedef std::vector<double> data_array_t;
#endif

#include <get_time.h>

#define FUNC(x) exp(x)

static int ARRAY_SIZE = 2048;

class Interpolator
{
public:
    Interpolator(int steps) {
        // Set interpolation array
        vals.resize(steps+1);

        delta = 1.0 / steps;
        delta_inv = 1./delta;

        vals[0].c0 = 0.;
        vals[0].c1 = FUNC(0.);

        double prev_val = vals[0].c1;
        // Fill interpolation array
        for(int ind=1; ind <= steps; ++ind) {
            double x = ind * delta;
            double val = FUNC(x);
            double c0 = (val - prev_val)*delta_inv;
            double c1 = val - c0*x;
            vals[ind].c0 = c0;
            vals[ind].c1 = c1;
            prev_val = val;
        }
        val_ptr = vals.data();
    }

#ifdef _OPENMP
#ifndef __MIC__
    #pragma omp declare simd simdlen(4) processor(core_4th_gen_avx)
    #pragma omp declare simd simdlen(2) processor(core_i7_sse4_2)
#ifdef __INTEL_COMPILER	
    #pragma omp declare simd simdlen(4) uniform(this) processor(core_4th_gen_avx)
    #pragma omp declare simd simdlen(2) uniform(this) processor(core_i7_sse4_2)
#endif
#else
    #pragma omp declare simd simdlen(8)
#endif
#endif
    int FindPosition(double x) const {
        return (int)(log(exp(x*delta_inv)));
    }

#ifdef _OPENMP
#ifndef __MIC__
    #pragma omp declare simd processor(core_4th_gen_avx)
    #pragma omp declare simd processor(core_i7_sse4_2)
#ifdef __INTEL_COMPILER	
    #pragma omp declare simd uniform(this) processor(core_4th_gen_avx)
    #pragma omp declare simd uniform(this) processor(core_i7_sse4_2)
#endif
#else
    #pragma omp declare simd uniform(this)
#endif
#endif
    double Interpolate(double x) const {
        int ind = FindPosition(x);
        const point& pnt = val_ptr[ind];
        double res = log(exp(pnt.c0*x+pnt.c1));
        return res;
    }

private:

    struct point {
        double c0;
        double c1;
    };
    
    double delta;
    double delta_inv;
    std::vector<point> vals;
    const point* val_ptr;
};

int main(int argc, char* aargv[])
{
    Interpolator interpolator(512);

    data_array_t src(ARRAY_SIZE);
    data_array_t dst(ARRAY_SIZE);

    // Create initial values
    for(int i=0;i<ARRAY_SIZE;++i) {
        src[i]=((double)rand())/RAND_MAX;
    }

    const double* src_ptr = src.data();
    double* dst_ptr = dst.data();

    START_LOOP_TIME()
#ifdef _OPENMP
#ifdef __INTEL_COMPILER
    #pragma omp simd aligned(src_ptr:64, dst_ptr:64)
#else
    #pragma omp simd
#endif
#endif
    for(int i=0; i<ARRAY_SIZE; ++i) {
        dst_ptr[i] = interpolator.Interpolate(src_ptr[i]);
    }
    END_LOOP_TIME()


    /* Test results */
    for(int i=0;i<ARRAY_SIZE;++i) {
        double ref_val = FUNC(src[i]);

        if ( (ref_val-dst[i])/ref_val > 0.01 ) {
            printf("Error at %i, val %lf, got %lf, expected %lf\n", i, src[i], dst[i], ref_val);
	    return -1;
        }
    }

    return 0;
}
