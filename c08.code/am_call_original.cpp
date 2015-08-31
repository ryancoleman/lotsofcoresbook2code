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
#include <cmath>
#include <algorithm> 
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <sys/time.h>
#include <time.h>

#if defined (__INTEL_COMPILER) 
#include <tbb/scalable_allocator.h>
#define aligned_malloc(memptr, alignment, size) (*memptr = _mm_malloc( size,  alignment))
#define aligned_free(memptr)  _mm_free(memptr) 
#else
#define aligned_malloc(memptr, alignment, size) (posix_memalign(memptr, alignment, size)) 
#define aligned_free(memptr)  free(memptr) 
#endif

using namespace std;
#ifndef PI 
#define PI 3.141592653589793238462643f
#endif

const float ACCURACY=1.0e-6;

const int ALIGNMENT = 1024;
//static int      OPT_N = 61*3*3*7*1024*32;
static int      OPT_N = 61*3*3*7*1024*4;
const float  RISKFREE = 0.02f;
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
///////////////////////////////////////////////////////////////////////////////
// Polynomial approximation of cumulative normal distribution function
///////////////////////////////////////////////////////////////////////////////
float cnd(float d){
    const float       A1 = 0.31938153;
    const float       A2 = -0.356563782;
    const float       A3 = 1.781477937;
    const float       A4 = -1.821255978;
    const float       A5 = 1.330274429;
    const float RSQRT2PI = 0.39894228040143267793994605993438;

    float
        K = 1.0 / (1.0 + 0.2316419 * fabs(d));

    float
        cnd = RSQRT2PI * exp(- 0.5 * d * d) *
        (K * (A1 + K * (A2 + K * (A3 + K * (A4 + K * A5)))));

    if(d > 0)
        cnd = 1.0 - cnd;

    return cnd;
}

float n(const float& z) {  // normal distribution function    
    return (1.0/sqrt(2.0*PI))*exp(-0.5*z*z);
};

float european_call_payout( const float& S,
		const float& X,
		const float& r,
		const float& q,
		const float& sigma,
		const float& time) {
    float sigma_sqr = sigma*sigma;
    float time_sqrt = sqrt(time);
    float d1 = (log(S/X) + (r-q + 0.5*sigma_sqr)*time)/(sigma*time_sqrt);
    float d2 = d1-(sigma*time_sqrt);
    float call_price = S * exp(-q*time)* cnd(d1) - X * exp(-r*time) * cnd(d2);
    return call_price;
};

float option_price_american_call_approximated_baw( const float& S,
						    const float& X, 
						    const float& r,
						    const float& b,
						    const float& sigma,
						    const float& time)
{
    float sigma_sqr = sigma*sigma;
    float time_sqrt = sqrt(time);
    float nn = 2.0*b/sigma_sqr; 
    float m = 2.0*r/sigma_sqr;  
    float K = 1.0-exp(-r*time); 
    float q2 = (-(nn-1)+sqrt(pow((nn-1),2.0)+(4*m/K)))*0.5;

    float q2_inf = 0.5 * ( -(nn-1) + sqrt(pow((nn-1),2.0)+4.0*m));
    float S_star_inf = X / (1.0 - 1.0/q2_inf);
    float h2 = -(b*time+2.0*sigma*time_sqrt)*(X/(S_star_inf-X));
    float S_seed = X + (S_star_inf-X)*(1.0-exp(h2));

    int no_iterations=0;
    float Si=S_seed;         
    float g=1;
    float gprime=1.0;
    while ( no_iterations++<100) {
	float c  = european_call_payout(Si,X,r,b,sigma,time);
	float d1 = (log(Si/X)+(b+0.5*sigma_sqr)*time)/(sigma*time_sqrt);
	g=(1.0-1.0/q2)*Si-X-c+(1.0/q2)*Si*exp((b-r)*time)*cnd(d1);
	gprime=( 1.0-1.0/q2)*(1.0-exp((b-r)*time)*cnd(d1))
	    +(1.0/q2)*exp((b-r)*time)*n(d1)*(1.0/(sigma*time_sqrt));
	Si=Si-(g/gprime); 
    };
    float S_star = 0;
    if (fabs(g)>ACCURACY) { S_star = S_seed; } // did not converge
    else { S_star = Si; };
    float C=0;
    float c  = european_call_payout(S,X,r,b,sigma,time);
    if (S>=S_star) {
	C=S-X;
    } 
    else {
	float d1 = (log(S_star/X)+(b+0.5*sigma_sqr)*time)/(sigma*time_sqrt);
	float A2 =  (1.0-exp((b-r)*time)*cnd(d1))* (S_star/q2); 
	C=c+A2*pow((S/S_star),q2);
    };
    return max(C,c);
};

int  main()
{
    unsigned long long start_cyc;
    unsigned long long end_cyc;

    float S = 100;   float X = 100;     float sigma = 0.20;
    float r = 0.08;  float b = -0.04;   float time = 0.25;
//    start_cyc = _rdtsc();
    float res =option_price_american_call_approximated_baw(S,X,r,b,sigma,time);
//    end_cyc = _rdtsc();
//    cout << " Call price using Barone-Adesi Whaley approximation Original = " 
//	 << std::setprecision(5) << res << endl << "  cycles consumed is " << end_cyc - start_cyc  << endl;
//    float  base_time = end_cyc -start_cyc;
//    start_cyc = _rdtsc();
//    float res = baw_scalaropt(S,X,r,b,sigma,time); 
//    end_cyc = _rdtsc();
    cout << " Call price using Barone-Adesi Whaley approximation Optimized = " 
	 << std::setprecision(7) << res << endl;
// << "  cycles consumed is " << end_cyc - start_cyc  << endl;
//    cout  << "Scalar improvement is " << std::setprecision(3) << base_time/float(end_cyc - start_cyc) << "X." << endl;
	int mem_size = sizeof(float) * OPT_N; 
	setlocale(LC_ALL,"");
	printf("Pricing American Options using BAW Approximation,  Number of options = %d.\n", OPT_N);
	int threadID = 0; 
	float *CallResult;
	float *CallResult2;
	float *StockPrice;
	float *OptionStrike;
	float *OptionYears;
	float *CostofCarry;
	float *Volatility;

	aligned_malloc((void **)&CallResult, ALIGNMENT, mem_size);
	aligned_malloc((void **)&CallResult2, ALIGNMENT, mem_size);
	aligned_malloc((void **)&StockPrice, ALIGNMENT, mem_size);
	aligned_malloc((void **)&OptionStrike, ALIGNMENT, mem_size);
	aligned_malloc((void **)&OptionYears, ALIGNMENT, mem_size);
	aligned_malloc((void **)&CostofCarry, ALIGNMENT, mem_size);
	aligned_malloc((void **)&Volatility, ALIGNMENT, mem_size);
        unsigned int seed = 123;
	unsigned int thread_seed = seed + threadID;
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
	for (int opt = 0; opt < OPT_N; opt++)
	{

		float T = OptionYears[opt];
		float S = StockPrice[opt];
		float X = OptionStrike[opt];
		float b = CostofCarry[opt];
		float v = Volatility[opt];

		CallResult[opt] = option_price_american_call_approximated_baw( S,X,RISKFREE,b, v, T);
	}
	eTime = second();
       	printf("Completed pricing %7.5f million options in %7.5f seconds:\n", OPT_N/1e6, eTime-sTime);
       	printf("Original version runs at %7.5f Thousand options per second.\n", OPT_N/(1e3*(eTime-sTime)));

	aligned_free(CallResult);
	aligned_free(CallResult2);
	aligned_free(StockPrice);
	aligned_free(OptionStrike);
	aligned_free(OptionYears);
	aligned_free(CostofCarry);
	aligned_free(Volatility);
};
