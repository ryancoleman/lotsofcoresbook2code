// This example from an alpha release of the Scalable HeterOgeneous Computing
// (SHOC) Benchmark Suite Alpha v1.1.4a-mic for Intel MIC architecture
// Contact: Kyle Spafford <kys@ornl.gov>
//          Rezaur Rahman <rezaur.rahman@intel.com>
//
// Copyright (c) 2011, UT-Battelle, LLC
// Copyright (c) 2013, Intel Corporation
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
//  * Neither the name of Oak Ridge National Laboratory, nor UT-Battelle, LLC,
//    nor the names of its contributors may be used to endorse or promote
//    products derived from this software without specific prior written
//    permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
// OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
// THE POSSIBILITY OF SUCH DAMAGE.

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "omp.h"
#include "OptionParser.h"
#include "ResultDatabase.h"
#include "Timer.h"

#ifdef __MIC2__
#include <immintrin.h>
#endif

#include "Scan.h"

using namespace std;

#define ALIGN	4096
#define ERR 	1.0e-4
#define L1B		32768

#include "Scan_Kernel.h"

// ****************************************************************************
// Function: addBenchmarkSpecOptions
//
// Purpose:
//   Add benchmark specific options parsing
//
// Arguments:
//   op: the options parser / parameter database
//
// Returns:  nothing
//
// Programmer: Kyle Spafford
// Creation: August 13, 2009
//
// Modifications:
//
// ****************************************************************************
void addBenchmarkSpecOptions(OptionParser &op)
{
    op.addOption("iterations", OPT_INT, "256", "specify scan iterations");
}

// ****************************************************************************
// Function: RunBenchmark
//
// Purpose:
//   Executes the scan (parallel prefix sum) benchmark
//
// Arguments:
//   resultDB: results from the benchmark are stored in this db
//   op: the options parser / parameter database
//
// Returns:  nothing
//
// Programmer: Kyle Spafford
// Creation: August 13, 2009
//
// Modifications:
//
// ****************************************************************************
void
RunBenchmark(OptionParser &op, ResultDatabase &resultDB)
{
    cout << "Running single precision test" << endl;
    RunTest<float>("Scan", resultDB, op);

    // Test to see if this device supports double precision
    cout << "Running double precision test" << endl;
    RunTest<double>("Scan-DP", resultDB, op);
}

template <class T>
void RunTest(string testName, ResultDatabase &resultDB, OptionParser &op)
{
	int     pbIndex 		= op.getOptionInt("size") - 1;
	int 	passes 			= op.getOptionInt("passes");
    int 	iters 			= op.getOptionInt("iterations");
	int     micdev 			= op.getOptionInt("target");

	int		nThreads		= 240; // Default

	#pragma offload target(mic) inout(nThreads)
	{
		nThreads = sysconf(_SC_NPROCESSORS_ONLN) - 4; // Leave something for the OS
	}

	printf("Using %d available threads for MIC run.\n", nThreads);

	size_t	szOptimum		= L1B * nThreads;

	int 	pbSizesMB[] 	= { 1, 8, 32, 64, 128, 256, 512, 768 };
	size_t 	pbSizeBytes		= szOptimum * pbSizesMB[pbIndex] / 8;
	int 	pbSizeElements 	= pbSizeBytes / sizeof(T);

    // Allocate Host Memory
    __declspec(target(MIC)) static T* h_idata;
    __declspec(target(MIC)) static T* reference;
    __declspec(target(MIC)) static T* h_odata;

	h_idata 	= (T*)_mm_malloc(pbSizeBytes + ALIGN * sizeof(T), ALIGN);
	reference 	= (T*)_mm_malloc(pbSizeBytes + ALIGN * sizeof(T), ALIGN);
	h_odata 	= (T*)_mm_malloc(pbSizeBytes + ALIGN * sizeof(T), ALIGN);

	//Manually align memory
	h_idata += ALIGN - 1;
	h_odata += ALIGN - 1;

	srand(time(NULL));
    // Initialize host memory
    for (int i = 0; i < pbSizeElements; i++)
    {
        h_idata[i]    = rand() % 21 - 10; // Fill with some pattern
        h_odata[i]    = 0.0;
        reference[i]  = 0.0;
    }


    // Allocate data to mic
    #pragma offload target(mic:micdev) in(h_idata:length(pbSizeElements + 1) free_if(0)) \
        out(h_odata:length(pbSizeElements + 1) free_if(0))
    {
    }

    double start = curr_second();
    // Get data transfer time
    #pragma offload target(mic:micdev) in(h_idata:length(pbSizeElements + 1) alloc_if(0) \
            free_if(0)) out(h_odata:length(pbSizeElements + 1) alloc_if(0) free_if(0))
    {
    }

    float transferTime = curr_second()-start;



    cout << "Running benchmark with size " << pbSizeElements << endl;
    for (int k = 0; k < passes; k++)
    {

        double totalScanTime = 0.0f;
        start = curr_second();
        #pragma offload target(mic:micdev) nocopy(h_idata:length(pbSizeElements + 1) \
                alloc_if(0) free_if(0)) nocopy(h_odata:length(pbSizeElements + 1)    \
                alloc_if(0) free_if(0))
        {
			if (pbIndex > 0)
			{
				size_t elementsOptimum = szOptimum / sizeof(T);
				int  nChunks = pbSizesMB[pbIndex] / 8;
				T fOffset = 0;

				for (int iChunk = 0; iChunk < nChunks; iChunk++)
				{
					SCAN_KNC<T>(h_idata + iChunk * elementsOptimum, h_odata + iChunk * elementsOptimum, elementsOptimum, iters, fOffset);
					fOffset = (h_odata + iChunk * elementsOptimum)[elementsOptimum - 1];
				}
			}
			else
				SCAN_KNC<T>(h_idata, h_odata, szOptimum / (8 * sizeof(T)), iters, 0.0);
        }

        double stop = curr_second();
        totalScanTime = (stop-start);

        #pragma offload target(mic:micdev) out(h_odata:length(pbSizeElements + 1) \
                alloc_if(0) free_if(0))
        {
        }

        // If results aren't correct, don't report perf numbers
        if (! scanCPU<T>(h_idata, reference, h_odata, pbSizeElements))
        {
            return;
        }

        char atts[1024];
        double avgTime = (totalScanTime / (double) iters);
        sprintf(atts, "%d items", pbSizeElements);
        double gb = (double)(pbSizeElements * sizeof(T)) / (1000. * 1000. * 1000.);
        resultDB.AddResult(testName, atts, "GB/s", gb / avgTime);
        resultDB.AddResult(testName+"_PCIe", atts, "GB/s", gb / (avgTime + transferTime));
        resultDB.AddResult(testName+"_Parity", atts, "N", transferTime / avgTime);
    }

    // Clean up
    #pragma offload target(mic:micdev) in(h_idata:length(pbSizeElements + 1) alloc_if(0) ) \
                                out(h_odata:length(pbSizeElements + 1) alloc_if(0))
    {
    }
	_mm_free(h_idata - ALIGN + 1);
    _mm_free(h_odata - ALIGN + 1);
    _mm_free(reference);
}

// ****************************************************************************
// Function: scanCPU
//
// Purpose:
//   Simple cpu scan routine to verify device results
//
// Arguments:
//   data : the input data
//   reference : space for the cpu solution
//   dev_result : result from the device
//   size : number of elements
//
// Returns:  nothing, prints relevant info to stdout
//
// Modifications:
//
// ****************************************************************************
template <class T>
bool scanCPU(T *data, T* reference, T* dev_result, const size_t size)
{
    reference[0] = 0;
    bool passed = true;

    // NB: You cannot validate beyond a certain buffer size because
    // of rounding errors.
    if (size > 128)
    {
        // This is an inclusive scan while the OpenMP code is an exclusive scan
        reference[0] = data[0];
        for (unsigned int i = 1; i < size; ++i)
        {
            reference[i] = data[i] + reference[i - 1];
        }

        for (unsigned int i = 0; i < size; ++i)
        {
            if (abs(reference[i] - dev_result[i]) > ERR )
            {
#ifdef VERBOSE_OUTPUT
                cout << "Mismatch at i: " << i << " ref: " << reference[i - 1]
                     << " dev: " << dev_result[i] << endl;
#endif
                passed = false;
            }
        }
    }
    cout << "Test ";
    if (passed)
        cout << "Passed" << endl;
    else
        cout << "Failed" << endl;
    return passed;
}
