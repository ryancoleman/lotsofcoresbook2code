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

#ifndef MATRIX2D_H
#define MATRIX2D_H

#include <stdio.h>
#include <iostream>
#include "config.h"

// ****************************************************************************
// Class:  Matrix2D
//
// Purpose:
//   Encapsulation of 2D matrices.
//
// Programmer:  Phil Roth
// Creation:    October 28, 2009
//
// ****************************************************************************
template<class T>
class Matrix2D
{
public:
    typedef T* restrict FlatDataPtr;
    typedef T* restrict* restrict DataPtr;
    typedef T* const restrict* const restrict ConstDataPtr;

private:
    size_t nRows;
    size_t nColumns;
    FlatDataPtr flatData;   // 1D array of data
    DataPtr data;           // data as 2D array (ptr to array of ptrs)

public:
    Matrix2D( size_t _nRows, size_t _nColumns )
      : nRows( _nRows ),
        nColumns( _nColumns ),
        flatData((T*) _mm_malloc(nRows*nColumns*sizeof(T), 64)),//new T[nRows * nColumns] ),
        data( new T*[nRows] )
    {
        for( size_t i = 0; i < nRows; i++ )
        {
            data[i] = &(flatData[i * nColumns]);
        }
    }

    ~Matrix2D( void )
    {
        delete[] data;
        data = NULL;

        _mm_free(flatData);
        flatData = NULL;
    }

    DataPtr GetData( void )
    {
        return data;
    }

    ConstDataPtr GetConstData( void ) const
    {
        return data;
    }

    FlatDataPtr GetFlatData( void )
    {
        return flatData;
    }

    size_t GetNumRows( void ) const { return nRows; }
    size_t GetNumColumns( void ) const { return nColumns; }

    size_t GetDataSize( void ) const { return nRows * nColumns * sizeof(T); }
};


template<class T>
std::ostream&
operator<<( std::ostream& s, const Matrix2D<T>& m )
{
    typename Matrix2D<T>::ConstDataPtr mdata = m.GetConstData();

    for( unsigned int i = 0; i < m.GetNumRows(); i++ )
    {
        for( unsigned int j = 0; j < m.GetNumColumns(); j++ )
        {
            if( j != 0 )
            {
                s << '\t';
            }
            s << mdata[i][j];
        }
        s << '\n';
    }
    return s;
}

#endif /* MATRIX2D_H */
