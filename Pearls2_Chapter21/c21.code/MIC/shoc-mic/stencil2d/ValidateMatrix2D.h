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

#ifndef VALIDATE_H
#define VALIDATE_H

#include <functional>
#include <vector>
#include "Matrix2D.h"


// ****************************************************************************
// Struct:  ValidationErrorInfo
//
// Purpose:
//   Stores information about validation errors originating in a 2D grid.
//
// Programmer:  Phil Roth
// Creation:    October 28, 2009
//
// ****************************************************************************
template<class T>
struct ValidationErrorInfo
{
    int i;
    int j;
    T val;
    T exp;
    double relErr;

    ValidationErrorInfo( int _i, int _j,
                            T _val,
                            T _exp,
                            double _relErr )
      : i( _i ),
        j( _j ),
        val( _val ),
        exp( _exp ),
        relErr( _relErr )
    {
        // nothing else to do
    }
};

// ****************************************************************************
// Class:  Validate
//
// Purpose:
//   Compares 2D matrices.
//
// Programmer:  Phil Roth
// Creation:    October 28, 2009
//
// ****************************************************************************
template<class T>
class Validate : public std::binary_function<const Matrix2D<T>&, const Matrix2D<T>&, std::vector<ValidationErrorInfo<T> > >
{
private:
    double relErrThreshold;

public:
    Validate( double _relErrThreshold )
      : relErrThreshold( _relErrThreshold )
    {
        // nothing else to do
    }

    std::vector<ValidationErrorInfo<T> > operator()( const Matrix2D<T>& s, const Matrix2D<T>& t );
};

#endif // VALIDATE_H
