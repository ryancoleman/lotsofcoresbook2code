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

#ifndef STENCIL_H
#define STENCIL_H

#include <string>
#include <functional>
#include "Matrix2D.h"

// ****************************************************************************
// Class:  Stencil
//
// Purpose:
//   9-point stencil.
//
// Programmer:  Phil Roth
// Creation:    October 28, 2009
//
// ****************************************************************************
template<class T>
class Stencil : public std::binary_function<Matrix2D<T>&, unsigned int, void>
{
protected:
    T wCenter;
    T wCardinal;
    T wDiagonal;

protected:
    T GetCenterWeight( void ) const { return wCenter; }
    T GetCardinalWeight( void ) const { return wCardinal; }
    T GetDiagonalWeight( void ) const { return wDiagonal; }

public:
    Stencil( T _wCenter,
                T _wCardinal,
                T _wDiagonal )
      : wCenter( _wCenter ),
        wCardinal( _wCardinal ),
        wDiagonal( _wDiagonal )
    {
        // nothing else to do
    }


    /*
     * This is a 9-point stencil using three weights:
     *   wCenter is applied to the stencil 'center'
     *   wCardinal is applied to the sum of the stencil NSEW values
     *   wDiagonal is applied to the sum of the stencil diagonal values
     *
     * note two things:
     *   We use the overall boundary values but do not update them.
     *   We apply wCardinal and wDiagonal *only* to the sum of the NSEW and
     *     diagonal values. We don't do any other averaging, etc.
     */
    virtual void operator()( Matrix2D<T>& m, unsigned int nIters ) = 0;
};

#endif // STENCIL_H
