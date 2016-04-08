// This example from an alpha release of the Scalable HeterOgeneous Computing
// (SHOC) Benchmark Suite Alpha v1.1.4a-mic for Intel MIC architecture
// Contact: Kyle Spafford <kys@ornl.gov>
//          Rezaur Rahman <rezaur.rahman@intel.com>
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

#include <assert.h>
#include "StencilFactory.h"
#include "InvalidArgValue.h"


template<class T>
void
StencilFactory<T>::CheckOptions( const OptionParser& options ) const
{
    // number of iterations must be positive
    unsigned int nIters = (unsigned int)options.getOptionInt( "num-iters" );
    if( nIters == 0 )
    {
        throw InvalidArgValue( "number of iterations must be positive" );
    }

    // no restrictions on weight values, just that we have them
}

template<class T>
void
StencilFactory<T>::ExtractOptions( const OptionParser& options,
                                T& wCenter,
                                T& wCardinal,
                                T& wDiagonal )
{
    wCenter   = options.getOptionFloat( "weight-center" );
    wCardinal = options.getOptionFloat( "weight-cardinal" );
    wDiagonal = options.getOptionFloat( "weight-diagonal" );
}


template<class T>
std::vector<long long>
StencilFactory<T>::GetStandardProblemSize( int sizeClass )
{
    const int probSizes[4] = { 768, 1408, 2048, 4096 };
    if (!(sizeClass >= 0 && sizeClass < 5))
    {
        throw InvalidArgValue( "Size class must be between 1-4" );
    }

    std::vector<long long> ret( 2, probSizes[sizeClass - 1] );
    return ret;
}

