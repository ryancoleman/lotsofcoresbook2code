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

#include <string.h>

template <class T> __declspec(target(mic)) void SCAN_KNC(	T* 				pInput,
															T* 				pOutput,
															const 	size_t 	nElements,
															const	int		nIterations,
															T				fOffset)
{
	// Keep partial sums in it
	unsigned int nThreads	= sysconf(_SC_NPROCESSORS_ONLN) - 4;	// Leave something for the OS

	T*	pPartialSums    	= (T*)_mm_malloc((nThreads + 1) * sizeof(T), ALIGN);
	size_t	nThreadElements	= nElements / nThreads;

	__declspec(target(mic)) __declspec(align(64)) volatile int g_nThdIndex = -1;

	#pragma omp parallel // num_threads(nThreads) - Not really necessary
	{
		int i = _InterlockedIncrement((void *)&g_nThdIndex);

		for (int iteration = 0; iteration < nIterations; iteration++)
		{
			// Working pointers
			T*	pCrntInput	= pInput;
			T*	pCrntOutput	= pOutput;

			T fPartialSums = 0;

			#pragma simd reduction(+:fPartialSums) vectorlengthfor(T)
			for (int j = 0; j < nThreadElements; j++)
				fPartialSums += pCrntInput[j + i * nThreadElements];

			pPartialSums[i + 1] = fPartialSums;

			#pragma omp barrier
			;

			if(i == 0)
			{
				pPartialSums[0] = fOffset;
				for (int j = 1; j <= nThreads; j++)
				{
					pPartialSums[j] += pPartialSums[j-1];
				}
			}

			#pragma omp barrier
			;

			pCrntOutput[i * nThreadElements] = pCrntInput[i * nThreadElements] + pPartialSums[i];

			for (int j = i * nThreadElements + 1; j < (i+1) * nThreadElements; j += 32)
			{
				__declspec(align(64)) T	vTemp[32];

				// Don't use G/S
				#pragma novector
				for (int k = 0; k < 16; k++)
					vTemp[k * 2 + 1] = pCrntInput[j + k * 2] + pCrntInput[j + k * 2 + 1];

				// Don't use G/S
				#pragma novector
				for (int k = 1; k < 16; k ++)
					vTemp[k * 2 + 1] += vTemp[k * 2 - 1];

				vTemp[0] = pCrntInput[j];

				// Don't use G/S
				#pragma novector
				for (int k = 1; k < 16; k++)
					vTemp[k * 2] = vTemp[k * 2 - 1] + pCrntInput[j + k * 2];

				T* 	out = pCrntOutput + j;
				int inc = pCrntOutput[j - 1];

				__assume_aligned(out, 64);

				#pragma simd vectorlengthfor(T)
				for (int k = 0; k < 32; k++)
					out[k] = inc + vTemp[k];
			}
		}
	}
}
