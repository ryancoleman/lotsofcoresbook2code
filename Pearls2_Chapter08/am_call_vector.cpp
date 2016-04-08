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
#include <dvec.h>
#include <sys/time.h>
#include <time.h>
#include <tbb/scalable_allocator.h>

using namespace std;
#ifndef PI 
#define PI 3.141592653589793238462643
#endif

const float ACCURACY=1.0e-6;

const int ALIGNMENT = 64;

static int      OPT_N = 61*3*3*7*1024*32;
const float  RISKFREE = 0.02f;
const	F32vec8 R = F32vec8(RISKFREE);
double sTime, eTime;

double second()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec + (double)tv.tv_usec / 1000000.0;
}
inline float RandFloat(float low, float high){
    float t = (float)rand() / (float)RAND_MAX;
    return (1.0f - t) * low + t * high;
}
const F32vec8 twof = F32vec8(2.0f);
const F32vec8 onef = F32vec8(1.0f);
const F32vec8 halff = F32vec8(0.5f);
const F32vec8 fourf = F32vec8(4.0f);
const F32vec8 rsqrtf2(float(M_SQRT1_2));


__forceinline
F32vec8 vn(F32vec8 z) {  // normal distribution function    
   F32vec8 rsqrt2pi(1.0/sqrt(2.0*PI));
    return rsqrt2pi*exp(-halff*z*z);
};
__forceinline
F32vec8 N(F32vec8 x)
{
   return halff + halff*erf(rsqrtf2*x); 
}
__forceinline
F32vec8 bs_call( F32vec8 S, // spot price
	 F32vec8 X,         // Strike (exercise) price,
	 F32vec8 r,         // interest rate
	 F32vec8 q,         // yield on underlying
	 F32vec8 sigma,     // volatility
	 F32vec8 time)      // time to maturity
{
    F32vec8 sigma_sqr = sigma*sigma;
    F32vec8 time_sqrt = sqrt(time);
    F32vec8 d1 = (log(S*rcp_nr(X)) + (r-q + halff*sigma_sqr)*time)*rcp_nr(sigma*time_sqrt);
    F32vec8 d2 = d1-(sigma*time_sqrt);
    F32vec8 call_price = S * exp(-q*time)* N(d1) - X * exp(-r*time) * N(d2);
    return call_price;
};

F32vec8 baw_vecopt(F32vec8 S,
		F32vec8 X, 
		F32vec8 r,
		F32vec8 b,
		F32vec8 sigma,
		F32vec8 time)
{
    F32vec8 sigma_sqr = sigma*sigma;
    F32vec8 time_sqrt = sqrt(time);
    F32vec8 nn_1 = twof*b*rcp_nr(sigma_sqr) - onef; 
    F32vec8 m = twof*r*rcp_nr(sigma_sqr);  
    F32vec8 K = onef-exp(-r*time); 
    F32vec8 rq2 = rcp_nr((-(nn_1)+sqrt(nn_1*nn_1+(fourf*m*rcp_nr(K))))*halff);
    F32vec8 rq2_inf = rcp_nr(halff * ( -(nn_1) + sqrt(nn_1*nn_1+fourf*m)));
    F32vec8 S_star_inf = X *rcp_nr(onef -rq2_inf );
    F32vec8 h2 = -(b*time+twof*sigma*time_sqrt)*(X*rcp_nr(S_star_inf-X));
    F32vec8 S_seed = X + (S_star_inf-X)*(onef-exp(h2));

    int no_iterations=0;
    F32vec8 Si=S_seed;         
    F32vec8 g=onef;
    F32vec8 gprime=onef;
    while ( no_iterations++<100) {
	F32vec8 c  = bs_call(Si,X,r,b,sigma,time);
	F32vec8 d1 = (log(Si/X)+(b+halff*sigma_sqr)*time)/(sigma*time_sqrt);
	g=(onef-rq2)*Si-X-c+rq2*Si*exp((b-r)*time)*N(d1);
	gprime=( onef-rq2)*(onef-exp((b-r)*time)*N(d1))
	    +rq2*exp((b-r)*time)*vn(d1)*(onef/(sigma*time_sqrt));
	Si=Si-(g/gprime); 
    };
    F32vec8 S_star = F32vec8(0.0f);
    __m256 solved = _mm256_cmp_ps(abs(g), F32vec8(ACCURACY), _CMP_GT_OQ);
    S_star = _mm256_blendv_ps(Si, S_seed, solved);
    F32vec8 C = F32vec8(0.0f);
    F32vec8 c  = bs_call(S,X,r,b,sigma,time);
    F32vec8 d1 = (log(S_star*rcp_nr(X))+(b+halff*sigma_sqr)*time)*rcp_nr(sigma*time_sqrt);
    F32vec8 A2 =  (onef-exp((b-r)*time)*N(d1))* (S_star*rq2); 
    C = _mm256_blendv_ps(c+A2*pow((S/S_star), onef/rq2), S-X, _mm256_cmp_ps(S, S_star, _CMP_GE_OQ));
    return simd_max(C,c); 
}

void main()
{
    unsigned long long start_cyc;
    unsigned long long end_cyc;

	float S = 100;   float X = 100;     float sigma = 0.20;
	float r = 0.08;  float b = -0.04;   float time = 0.25;
	F32vec8 res = baw_vecopt(F32vec8(S),F32vec8(X),F32vec8(r),F32vec8(b),F32vec8(sigma),F32vec8(time)); 
	cout << " Call price using Barone-Adesi Whaley approximation Optimized = " 
	 << std::setprecision(7) << res[0] << endl;
	int mem_size = sizeof(float) * OPT_N; 
	setlocale(LC_ALL,"");
	printf("Pricing American Options using BAW Approximation in SP, Number of options = %'d.\n", OPT_N);
	float *CallResult = (float *)scalable_aligned_malloc(mem_size, ALIGNMENT);
	float *CallResult2 = (float *)scalable_aligned_malloc(mem_size, ALIGNMENT);
	float *StockPrice = (float *)scalable_aligned_malloc(mem_size, ALIGNMENT);
	float *OptionStrike = (float *)scalable_aligned_malloc(mem_size, ALIGNMENT);
	float *OptionYears = (float *)scalable_aligned_malloc(mem_size, ALIGNMENT);
	float *CostofCarry = (float *)scalable_aligned_malloc(mem_size, ALIGNMENT);
	float *Volatility = (float *)scalable_aligned_malloc(mem_size, ALIGNMENT);

        unsigned int seed = 123;
	srand(seed);
	for(int i = OPT_N-1; i > -1; i--)
	{
		CallResult[i] = 0.0f;
		StockPrice[i]    = RandFloat(5.0f, 30.0f);
		OptionStrike[i]  = RandFloat(1.0f, 100.0f);
		OptionYears[i]   = RandFloat(0.25f, 2.0f);
		CostofCarry[i]   = RandFloat(0.02f, 0.05f);
		Volatility[i]   = RandFloat(0.10, 0.60);
	}
	sTime = second();
	for (int opt = 0; opt < OPT_N; opt+=8)
	{
		F32vec8 T = *(F32vec8*)&OptionYears[opt];
		F32vec8 S = *(F32vec8*)&StockPrice[opt];
		F32vec8 X = *(F32vec8*)&OptionStrike[opt];
		F32vec8 b = *(F32vec8*)&CostofCarry[opt];
		F32vec8 v = *(F32vec8*)&Volatility[opt];

		*(F32vec8 *)&(CallResult2[opt]) = baw_vecopt(S,X,R,b,v,T);
	}
	eTime = second();
       	printf("Completed pricing %7.5f million options in %7.5f seconds:\n", OPT_N/1e6, eTime-sTime);
       	printf("Vector version runs at %7.5f Thousand options per second.\n", OPT_N/(1e3*(eTime-sTime)));
	scalable_aligned_free(CallResult);
	scalable_aligned_free(CallResult2);
	scalable_aligned_free(StockPrice);
	scalable_aligned_free(OptionStrike);
	scalable_aligned_free(OptionYears);
	scalable_aligned_free(CostofCarry);
	scalable_aligned_free(Volatility);
};
