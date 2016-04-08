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

#ifndef UTILITY_H
#define UTILITY_H

#include <sstream>
#include <math.h>

// ****************************************************************************
// File:  Utility.h
//
// Purpose:
//   Various generic utility routines having to do with string and number
//   manipulation.
//
// Programmer:  Jeremy Meredith
// Creation:    September 18, 2009
// Modified:    Jan 2010, rothpc
//
// ****************************************************************************

inline std::string HumanReadable(long long value, long long *rounding=0)
{
    std::ostringstream vstr;
    long long pVal;
    if (value>10ll*1024*1024*1024)
    {
        pVal = (long long)round(value/(1024.0*1024*1024));
        if (rounding)
            *rounding = pVal*1024*1024*1024 - value;
        vstr << pVal << 'G';
    }
    else if (value>10ll*1024*1024)
    {
        pVal = (long long)round(value/(1024.0*1024));
        if (rounding)
            *rounding = pVal*1024*1024 - value;
        vstr << pVal << 'M';
    }
    else if (value>10ll*1024)
    {
        pVal = (long long)round(value/(1024.0));
        if (rounding)
            *rounding = pVal*1024 - value;
        vstr << pVal << 'k';
    }
    else
    {
        if (rounding)
            *rounding = 0;
        vstr << value;
    }
    return vstr.str();
}

inline vector<string> SplitValues(const std::string &buff, char delim)
{
    vector<std::string> output;   
    std::string tmp="";
    for (size_t i=0; i<buff.length(); i++)
    {
       if (buff[i] == delim)
       {
          if (!tmp.empty())
             output.push_back(tmp);
          tmp = "";  
       }
       else
       {
          tmp += buff[i];  
       }   
    } 
    if (!tmp.empty())
       output.push_back(tmp);

    return output;
}


#endif
