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

#ifndef STENCILFACTORY_H
#define STENCILFACTORY_H

#include <map>
#include "../common/OptionParser.h"
#include "Stencil.h"

// ****************************************************************************
// Class:  StencilFactory
//
// Purpose:
//   Class to generate stencils.
//
// Programmer:  Phil Roth
// Creation:    October 28, 2009
//
// ****************************************************************************
template<class T>
class StencilFactory
{
public:
    typedef std::map<std::string, StencilFactory*> FactoryMap;

private:
    // map of class name to a StencilFactory object
    // would be much easier if C++ classes were first class objects
    // so that we could programmatically construct a class name and
    // then create an instance of that class
    static FactoryMap* factoryMap;

    std::string sname;

protected:
    void ExtractOptions( const OptionParser& options,
                        T& wCenter,
                        T& wCardinal,
                        T& wDiagonal );

public:
    StencilFactory( std::string _sname )
      : sname( _sname )
    {
        // nothing else to do
    }
    virtual ~StencilFactory( void ) { }

    std::string GetStencilName( void ) { return sname; }
    
    virtual Stencil<T>* BuildStencil( const OptionParser& options ) = 0;
    virtual void CheckOptions( const OptionParser& options ) const = 0;

    static std::vector<long long> GetStandardProblemSize( int sizeClass );
};

#endif // STENCILFACTORY_H
