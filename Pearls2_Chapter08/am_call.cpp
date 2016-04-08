//Copyright (c) 2015 Intel Corporation
//All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright
//   notice, this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
// * Neither the name of Intel Corporation
//   nor the names of its contributors may be used to endorse or promote
//   products derived from this software without specific prior written
//   permission.
//
//   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//   ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//   WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//   DISCLAIMED. IN NO EVENT SHALL INTEL CORPORATION BE LIABLE FOR ANY
//   DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//   (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
//   ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <sys/time.h>
#include <time.h>
#include <tbb/scalable_allocator.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
#ifndef PI 
#define PI 3.141592653589793238462643f
#endif

const float ACCURACY = 1.0e-6;
const int  ALIGNMENT = 1024;
static int     OPT_N = 61*3*3*7*1024*32;
const float RISKFREE = 0.02f;
double sTime, eTime;

double second()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec + (double)tv.tv_usec / 1000000.0;
}


inline float RandFloat(float low, float high, unsigned int *seed){
    float t = (float)rand_r(seed) / (float)RAND_MAX;
    return (1.0f - t) * low + t * high;
}
__forceinline
float cnd_opt(float d){
 return 0.5f +0.5f*erff(float(M_SQRT1_2)*d);
}

__forceinline
float n_opt(const float z) {
    return (1.0f/sqrtf(2.0f*PI))*expf(-0.5f*z*z);
};

__forceinline
float european_call_opt( const float S, 
		const float X,
		const float r,
		const float q,
		const float sigma,
		const float time)
{
    float sigma_sqr = sigma*sigma;
    float time_sqrt = sqrtf(time);
    float d1 = (logf(S/X) + (r-q + 0.5f*sigma_sqr)*time)/(sigma*time_sqrt);
    float d2 = d1-(sigma*time_sqrt);
    float call_price = S * expf(-q*time)* cnd_opt(d1) - X * expf(-r*time) * cnd_opt(d2);
    return call_price;
};
 __forceinline
float baw_scalaropt( const float S,
		 const float X, 
		 const float r,
		 const float b,
		 const float sigma,
		 const float time)
{
    float sigma_sqr = sigma*sigma;
    float time_sqrt = sqrtf(time);
    float nn_1 = 2.0f*b/sigma_sqr-1; 
    float m = 2.0f*r/sigma_sqr;  
    float K = 1.0f-expf(-r*time); 
    float rq2 = 1/((-(nn_1)+sqrtf((nn_1)*(nn_1) +(4.f*m/K)))*0.5f);

    float rq2_inf = 1/(0.5f * ( -(nn_1) + sqrtf(nn_1*nn_1+4.0f*m)));
    float S_star_inf = X / (1.0f - rq2_inf);
    float h2 = -(b*time+2.0f*sigma*time_sqrt)*(X/(S_star_inf-X));
    float S_seed = X + (S_star_inf-X)*(1.0f-expf(h2));
    float cndd1 = 0; 
    float Si=S_seed;         
    float g=1.f;
    float gprime=1.0f;
    float expbr=expf((b-r)*time);
    for (  int no_iterations =0; no_iterations<100; no_iterations++) {
	float c  = european_call_opt(Si,X,r,b,sigma,time);
	float d1 = (logf(Si/X)+(b+0.5f*sigma_sqr)*time)/(sigma*time_sqrt);
    	float cndd1=cnd_opt(d1);
	g=(1.0f-rq2)*Si-X-c+rq2*Si*expbr*cndd1;
	gprime=( 1.0f-rq2)*(1.0f-expbr*cndd1)+rq2*expbr*n_opt(d1)*(1.0f/(sigma*time_sqrt));
	Si=Si-(g/gprime); 
    };
    float S_star = 0;
    if (fabs(g)>ACCURACY) { S_star = S_seed; }
    else { S_star = Si; };
    float C=0;
    float c  = european_call_opt(S,X,r,b,sigma,time);
    if (S>=S_star) {
	C=S-X;
    } 
    else {
	float d1 = (logf(S_star/X)+(b+0.5f*sigma_sqr)*time)/(sigma*time_sqrt);
	float A2 =  (1.0f-expbr*cnd_opt(d1))* (S_star*rq2); 
	C=c+A2*powf((S/S_star),1/rq2);
    };
    return (C>c)?C:c;
};

void main()
{
    unsigned long long start_cyc;
    unsigned long long end_cyc;

    float S = 100;   float X = 100;     float sigma = 0.20;
    float r = 0.08;  float b = -0.04;   float time = 0.25;
    start_cyc = _rdtsc();
    float res = baw_scalaropt(S,X,r,b,sigma,time); 
    end_cyc = _rdtsc();
    cout << " Call price using Barone-Adesi Whaley approximation Optimized = " 
	 << std::setprecision(7) << res << endl << "  cycles consumed is " << end_cyc - start_cyc  << endl;
#ifdef _OPENMP
    	kmp_set_defaults("KMP_AFFINITY=scatter,granularity=fine");
	int ThreadNum = omp_get_max_threads();
	omp_set_num_threads(ThreadNum); 
#else
	int ThreadNum = 1; 
#endif
	int OptPerThread = OPT_N / ThreadNum;
	int mem_size = sizeof(float) * OptPerThread; 
	setlocale(LC_ALL,"");
	printf("Pricing American Options using BAW Approximation in %d threads, Number of options = %d.\n", ThreadNum, OPT_N);
	int threadID = 0; 
#pragma omp parallel 
{
#ifdef _OPENMP
	threadID = omp_get_thread_num();
#else
	threadID = 0; 
#endif

	float *CallResult = (float *)scalable_aligned_malloc(mem_size, ALIGNMENT);
	float *CallResult2 = (float *)scalable_aligned_malloc(mem_size, ALIGNMENT);
	float *StockPrice = (float *)scalable_aligned_malloc(mem_size, ALIGNMENT);
	float *OptionStrike = (float *)scalable_aligned_malloc(mem_size, ALIGNMENT);
	float *OptionYears = (float *)scalable_aligned_malloc(mem_size, ALIGNMENT);
	float *CostofCarry = (float *)scalable_aligned_malloc(mem_size, ALIGNMENT);
	float *Volatility = (float *)scalable_aligned_malloc(mem_size, ALIGNMENT);

        unsigned int seed = 123;
	unsigned int thread_seed = seed + threadID;
	for(int i = OptPerThread-1; i > -1; i--)
	{
		CallResult[i] = 0.0f;
		StockPrice[i]    = RandFloat(5.0f, 30.0f, &thread_seed);
		OptionStrike[i]  = RandFloat(1.0f, 100.0f, &thread_seed);
		OptionYears[i]   = RandFloat(0.25f, 2.0f, &thread_seed);
		CostofCarry[i]   = RandFloat(0.02f, 0.05f, &thread_seed);
		Volatility[i]   = RandFloat(0.10, 0.60, &thread_seed);
	}
#pragma omp barrier
#pragma omp master
	sTime = second();
#ifdef SCALAR
#else
#pragma vector nontemporal (CallResult)
#pragma simd
#pragma vector aligned
#endif
	for (int opt = 0; opt < OptPerThread; opt++)
	{

		float T = OptionYears[opt];
		float S = StockPrice[opt];
		float X = OptionStrike[opt];
		float b = CostofCarry[opt];
		float v = Volatility[opt];

		CallResult[opt] = baw_scalaropt(S,X,RISKFREE,b, v, T);
	}
#pragma omp barrier
#pragma omp master
	{
		eTime = second();
        	printf("Completed pricing %7.5f million options in %7.5f seconds:\n", OPT_N/1e6, eTime-sTime);
        	printf("Scalar %s version runs at %7.5f Thousand options per second.\n",(ThreadNum > 1)?"parallel":"serial", OPT_N/(1e3*(eTime-sTime)));
	}
	scalable_aligned_free(CallResult);
	scalable_aligned_free(CallResult2);
	scalable_aligned_free(StockPrice);
	scalable_aligned_free(OptionStrike);
	scalable_aligned_free(OptionYears);
	scalable_aligned_free(CostofCarry);
	scalable_aligned_free(Volatility);
}
};
