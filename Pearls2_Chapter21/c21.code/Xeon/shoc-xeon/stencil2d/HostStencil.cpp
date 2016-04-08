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

#include <string.h> // for memcpy
#include "HostStencil.h"

#define LINESIZE 64

template<class T>
void
HostStencil<T>::operator()( Matrix2D<T>& mtx, unsigned int nIters )
{
    // we need a temp space buffer
    Matrix2D<T> tmpMtx( mtx.GetNumRows(), mtx.GetNumColumns() );

    // be able to access the matrices as 2D arrays
    typename Matrix2D<T>::DataPtr mtxData = mtx.GetData();
    typename Matrix2D<T>::DataPtr tmpMtxData = tmpMtx.GetData();

	// Compute only for target stencil
	unsigned int uHaloWidth = LINESIZE / sizeof(T);

    for( unsigned int iter = 0; iter < nIters; iter++ )
    {
        DoPreIterationWork( mtx, iter );

        /* copy the "real" data to the temp matrix */
        memcpy( tmpMtx.GetFlatData(), mtx.GetFlatData(), mtx.GetDataSize() );


        /* Apply the stencil operator */
        for( size_t i = uHaloWidth; i < mtx.GetNumRows() - uHaloWidth; i++ )
        {
            for( size_t j = uHaloWidth; j < mtx.GetNumColumns() - uHaloWidth; j++ )
            {
                T oldCenterValue = tmpMtxData[i][j];
                T oldNSEWValues = (tmpMtxData[i-1][j] +
                                        tmpMtxData[i+1][j] +
                                        tmpMtxData[i][j-1] +
                                        tmpMtxData[i][j+1]);
                T oldDiagonalValues = (tmpMtxData[i-1][j-1] +
                                            tmpMtxData[i+1][j-1] +
                                            tmpMtxData[i-1][j+1] +
                                            tmpMtxData[i+1][j+1]);

                mtxData[i][j] = this->wCenter * oldCenterValue +
                                this->wCardinal * oldNSEWValues +
                                this->wDiagonal * oldDiagonalValues;
            }
        }
    }
}


template<class T>
void
HostStencil<T>::DoPreIterationWork( Matrix2D<T>& mtx, unsigned int iter )
{
    // we have nothing to do
}