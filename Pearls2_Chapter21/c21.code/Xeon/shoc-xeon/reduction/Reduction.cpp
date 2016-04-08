//
// This example from a prerelease of the Scalable HeterOgeneous Computing
// (SHOC) Benchmark Suite Alpha v1.1.1i for Intel MIC architecture
// Contact: Kyle Spafford <kys@ornl.gov>
//         Rezaur Rahman <rezaur.rahman@intel.com>
//
// Copyright (c) 2011, UT-Battelle, LLC
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of Oak Ridge National Laboratory, nor UT-Battelle, LLC, nor
//    the names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// ==============================================================================

// ****************************************************************************
// File: Reduction.cpp
//
// Purpose:
//   Contains performance tests on a basic sum reduction with a few
//   optimizations
//
// Programmer:  Kyle Spafford
// Creation:    November 19, 2010
//
// ****************************************************************************


#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include "omp.h"

#include "offload.h"
#include "Timer.h"

#include "OptionParser.h"
#include "ResultDatabase.h"


#ifdef __MIC2__
//#include <lmmintrin.h>
#include <pthread.h>
#endif


using namespace std;

int micdev = 0;

__attribute__ ((target(mic))) float *indataf = NULL;
__attribute__ ((target(mic))) double *indatad = NULL;
__attribute__ ((target(mic))) double *outdata = NULL;

// Forward declaration
template <class T>
void RunTest(string, ResultDatabase &, OptionParser &);

//data do be shared between host and mic using no copy should be made global

// ****************************************************************************
// Function: reduceGold
//
// Purpose:
//   Simple cpu reduce routine to verify device results.  This could be
//   replaced with Kahan summation for better accuracy.
//
// Arguments:
//   data : the input data
//   size : size of the input data
//
// Returns:  sum of the data
//
// ****************************************************************************
template <class T>
T reduceGold(const T *data, int size)
{
   T sum = 0;
   for (int i = 0; i < size; i++)
   {
      sum += data[i];
   }
   return sum;
}


// CPU Tests
// Vanilla OpenMP
template <typename T>
T cpuOMP(T *data, size_t size)
{
    T ret = 0.0f;
    #pragma omp parallel for reduction(+:ret)
    for(int i = 0; i < size; i++)
    {
      ret += data[i];
    }
   return ret;
}

// Single-thread CEAN - C++ Extended Array Notation
template <typename T>
T cpuCEAN(T *data, size_t size)
{
   T ret = 0.0f;
   ret = __sec_reduce_add(data[0:size]);
   return ret;
}

// Host - OpenMP + CEAN
template <typename T>
T cpuOMP_CEAN(T *data, size_t size)
{
    T ret=0;

    int nThreads = omp_get_max_threads();
    int nPerThread = size / nThreads;

    #pragma omp parallel for reduction(+:ret)
    for (int i = 0; i < nThreads; i++)
    {
        ret = __sec_reduce_add(data[i*nPerThread:nPerThread]);
    }

    // Rest of the array
    for(int i = nThreads * nPerThread; i < size; i++)
    {
        ret += data[i];
    }

    return ret;
}

template <typename T>
__declspec(target(mic)) T KNF_OMP(T *ldata, size_t size)
{
    T ret = 0.0f;
    //#pragma offload target(mic:micdev) in(data:length(size) align(4*1024*1024))
    {
        #pragma omp parallel for reduction(+:ret)
        for (int i = 0; i < size; i++)
        {
            ret += ldata[i];
        }
    }
    return ret;
}

//template <typename T> __attribute((vector))
//myreduceadd

template <typename T>
__declspec(target(mic)) T KNF_OMP_CEAN(T *data, size_t size)
{
	T ret = 0.0;
	T intermed[512];

	int nThreads 	= omp_get_max_threads();
	int nPerThread 	= size / nThreads;

	/*
	#pragma omp parallel firstprivate(nPerThread, nThreads)
	{
		unsigned int i = omp_get_thread_num();

		T loopres = 0.0;

		unsigned int nVectors 	= nPerThread / 16;
		unsigned int nRemaining = nPerThread % 16;
		unsigned int nStart		= i * nPerThread;
		unsigned int nStop		= nStart + nVectors * 16;

		T tAccumulate[16] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

		for (int j = nStart; j < nStop; j += 16)
		{
			#pragma vector
			for (int k = 0; k < 16; k++)
				tAccumulate[k] += data[j + k];
		}

		for (int j = 0; j < 16; j++)
			loopres += tAccumulate[j];


		for (int j = nStop; j < nStart + nPerThread; j++)
			loopres += data[j];

		intermed[i] = loopres;
	}
	*/

	#pragma omp parallel for firstprivate(nThreads, nPerThread, data)
	for (unsigned int i = 0; i < nThreads; i++)
	{
		T loopres = 0.0;

		//#pragma ivdep
		for (unsigned int j = i * nPerThread; j < i * nPerThread + nPerThread; j++)
			loopres += data[j];

		intermed[i] = loopres;
	}

	for (int i = 0; i < nThreads; i++)
		ret += intermed[i];

	// Rest of the array
	for(int i = nThreads * nPerThread; i < size; i++)
		ret += data[i];

	return ret;
}

// Unlike the other functions, this function returns time
#if 0
template <typename T>
double KNF_OMP_CEAN_NT(T* ldata,  size_t size, double *avgTime, double *dataXferTime)
{
   #define REP (100)
   T ret = 0.0f;

   // Ideally, I'd like to be able to put these copy directives
   // in a different function, but when I tried that (compiler 0.0.17)
   // earlier, I ran into major problems (scrolling pages of errors
   // about buffer not being sent).

  double start = curr_second();

   //Note: In order to use persistency, one needs to use the global data
   #pragma offload target(mic:micdev) in(data:length(size) align(2*1024*1024) free_if(0))
   {}

   *dataXferTime = curr_second()-start;

   // Repeatedly execute the kernel to amortize PCIe
   // transfer time

   start = curr_second();
   #pragma offload target(mic:micdev) nocopy(data)
   for (int j = 0; j < REP; j++)
   {
       ret = 0;
       //#pragma offload target(mic:micdev) nocopy(data)
       {
           int nThreads = omp_get_max_threads();
           int nPerThread = size / nThreads;
           #pragma omp parallel for reduction(+:ret)
           for(int i = 0; i < nThreads;i++)
           {
               ret = __sec_reduce_add(data[i*nPerThread:nPerThread]);
           }
           //rest of the array
           for(int i = nThreads * nPerThread; i < size; i++)
           {
               ret += data[i];
           }
       }
   }
   double stop = curr_second();
   // Free data
   #pragma offload target(mic:micdev) in(data:length(size) alloc_if(0))
   {
   }


   *avgTime = (stop - start)/REP;

   return ret;
}
#endif

bool check(float result, float ref) {

   float diff = fabs(result - ref);

   float threshold = 1e-2;
   if (diff >= threshold * ref)
   {
      cout << "Pass FAILED\n";
      cout << "Diff: " << diff;
      exit(-1); // (don't report erroneous results)
   }
   else
   {
     cout<< "Passed\n";
   }
   return true;
}

void addBenchmarkSpecOptions(OptionParser& op) {
 op.addOption("iterations", OPT_INT, "256",
                 "specify reduction iterations");
   return; // none
}


void RunTestf(string testName, ResultDatabase& resultDB, OptionParser& op) {

   int probSizes[8] = { 4, 8, 32, 64, 128, 256, 512, 768 };
   int N = probSizes[op.getOptionInt("size")-1];
   N = (N * 1024 * 1024) / sizeof(float);

   indataf = (float*) _mm_malloc(N * sizeof(float), 64);
   if (!indataf) return;

   outdata = (double *)_mm_malloc(64 * sizeof(double), 64);
   if (!outdata) return;

   // Initialize host memory
   cout << "Initializing memory." << endl;
   for(int i = 0; i < N; i++)
   {
      indataf[i] = i % 3; // Fill with some pattern
   }

   float ref = reduceGold<float>(indataf, N);
   int passes =op.getOptionInt("passes");
   int iterations = op.getOptionInt("iterations");;
   micdev = op.getOptionInt("target");

   // Test attributes
   char atts[1024];
   sprintf(atts, "%d_items",N);

   cout<<"Running Benchmark\n";
   for (int k = 0; k < passes; k++)
   {
	   	float 	result;
		double 	start, stop;
		double 	avgTime;
		double 	transferTime=0;

		#pragma offload target(mic) in(outdata:length(64) alloc_if(1) free_if(0))
		{
		}

		//warm up
		#pragma offload target(mic) in(indataf:length(N) alloc_if(1) free_if(0))
		{
		}

		start=curr_second();

		#pragma offload target(mic) in(indataf:length(N) align(64) alloc_if(0) free_if(0))
		{

		}

      	stop = curr_second();
      	transferTime=stop-start;

		//#pragma offload target(mic:micdev) nocopy(data:length(N) align(4*1024*1024) alloc_if(0) free_if(0))
      	//#pragma offload target(mic:micdev) in(data:length(N) align(4*1024*1024))

      	start = curr_second();
      	#pragma offload target(mic:micdev) nocopy(indataf:length(N) align(64) alloc_if(0) free_if(0))
		// #pragma offload target(mic:micdev) nocopy(indataf:length(N) align(4*1024*1024) alloc_if(0) free_if(0))
      	{
         	for (int j=0; j<iterations; j++)
         	{
            	result = (float)KNF_OMP_CEAN<float>(indataf, N);
         	}
      	}
      	stop = curr_second();

      	avgTime = (stop - start) / (double)iterations;

		start=curr_second();
		#pragma offload target(mic:micdev) out(outdata:length(64) alloc_if(0) free_if(1) )
		{

		}
		stop = curr_second();
		transferTime += (stop - start);

		check(result, ref);

		// Free buffer on card
		#pragma offload target(mic:micdev) nocopy(indataf:length(N)  alloc_if(0) free_if(1))
		{

		}

		double gbytes = (double)(N * sizeof(float))/(1000.*1000.*1000.);
		resultDB.AddResult(testName, atts, "GB/s", gbytes / avgTime);
		resultDB.AddResult(testName+"_PCIe", atts, "GB/s", gbytes / (avgTime + transferTime));
		resultDB.AddResult(testName+"_Parity", atts, "N", transferTime / avgTime);
	}

	_mm_free( indataf);
	_mm_free( outdata);
}

void RunTestd(string testName, ResultDatabase& resultDB, OptionParser& op)
{
   int probSizes[8] = { 4, 8, 32, 64, 128, 256, 512, 768 };
   int N = probSizes[op.getOptionInt("size")-1];
   N = (N * 1024 * 1024) / sizeof(double);

   indatad = (double*) _mm_malloc(N * sizeof(double), 64);
   if (!indatad) return;

   outdata = (double *)_mm_malloc(64 * sizeof(double), 64);
   if (!outdata) return;

   // Initialize host memory
   cout << "Initializing memory." << endl;
   for(int i = 0; i < N; i++)
   {
      indatad[i] = i % 3; // Fill with some pattern
   }

   double ref = reduceGold<double>(indatad, N);
   int passes =op.getOptionInt("passes");
   int iterations = op.getOptionInt("iterations");;
   micdev = op.getOptionInt("target");

   // Test attributes
   char atts[1024];
   sprintf(atts, "%d_items",N);

   cout<<"Running Benchmark\n";
   for (int k = 0; k < passes; k++)
   {
      double result;
      double start, stop;
      double avgTime;
      double transferTime;

#pragma offload target(mic:micdev)\
     in(outdata:length(64)  align(64) alloc_if(1) free_if(0))
      {
      }
      //warm up
       #pragma offload target(mic:micdev)\
                in(indatad:length(N) alloc_if(1) free_if(0))
            {
            }

      start=curr_second();

      #pragma offload target(mic:micdev)\
                in(indatad:length(N) align(64) alloc_if(0) free_if(0))
      {
      }

      stop = curr_second();
      transferTime=stop-start;


//      #pragma offload target(mic:micdev) nocopy(data:length(N) align(4*1024*1024) alloc_if(0) free_if(0))
      //#pragma offload target(mic:micdev) in(data:length(N) align(64))
      start = curr_second();

      #pragma offload target(mic:micdev) nocopy(indatad:length(N) align(64) alloc_if(0) free_if(0))
      {
         for (int j=0; j<iterations; j++) {
            result = (double)KNF_OMP_CEAN<double>(indatad, N);
         }
      }

      stop = curr_second();

      avgTime = (stop - start) / (double)iterations;

      start=curr_second();
      #pragma offload target(mic:micdev) out(outdata:length(64) alloc_if(0) free_if(1))
      {
      }
      stop = curr_second();
      transferTime += (stop - start);

      check(result, ref);

      // Free buffer on card
      #pragma offload target(mic:micdev) nocopy(indatad:length(N) align(64) alloc_if(0) free_if(1))
      {
      }

      double gbytes = (double)(N*sizeof(double))/(1000.*1000.*1000.);
      resultDB.AddResult(testName, atts, "GB/s", gbytes / avgTime);
      resultDB.AddResult(testName+"_PCIe", atts, "GB/s", gbytes /
            (avgTime + transferTime));
      resultDB.AddResult(testName+"_Parity", atts, "N",
            transferTime / avgTime);
   }

   _mm_free( indatad);
   _mm_free( outdata);
}

/*
 * Best performace with:
 * setenv MIC_ENV_PREFIX MIC
 * setenv MIC_OMP_NUM_THREADS 90
 * setenv MIC_KMP_AFFINITY balanced,granularity=fine
 */
void RunBenchmark(OptionParser &op, ResultDatabase &resultDB)
{
    cout << "Running single precision test" << endl;
    RunTestf("Reduction", resultDB, op);

    cout << "Running double precision test" << endl;
    RunTestd("Reduction-DP", resultDB, op);
}
