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

#ifndef RESULT_DATABASE_H
#define RESULT_DATABASE_H

#include <cfloat>
#include <iostream>
#include <string>
#include <vector>

using std::string;
using std::vector;
using std::ostream;

// ****************************************************************************
// Class:  ResultDatabase
//
// Purpose:
//   Track numerical results as they are generated.
//   Print statistics of raw results.
//
// Programmer:  Jeremy Meredith
// Creation:    June 12, 2009
//
// Modifications:
//    Jeremy Meredith, Wed Nov 10 14:20:47 EST 2010
//    Split timing reports into detailed and summary.  E.g. for serial code,
//    we might report all trial values, but skip them in parallel.
//
//    Jeremy Meredith, Thu Nov 11 11:40:18 EST 2010
//    Added check for missing value tag.
//
// ****************************************************************************
class ResultDatabase
{
    friend class ParallelResultDatabase; // TODO (JSM): this is a hack....
  protected:
    //
    // A performance result for a single SHOC benchmark run.
    //
    struct Result
    {
        string test;  // e.g. "readback"
        string atts;  // e.g. "pagelocked 4k^2"
        string unit;  // e.g. "MB/sec"
        vector<double> value; // e.g. "837.14"
        double GetMin();
        double GetMax();
        double GetMedian();
        double GetMean();
        double GetStdDev();

        bool operator<(const Result &rhs) const;

        bool HadAnyFLTMAXValues()
        {
            for (int i=0; i<value.size(); ++i)
            {
                if (value[i] >= FLT_MAX)
                    return true;
            }
            return false;
        }
    };

    vector<Result> results;

  public:
    void AddResult(const string &test,
                   const string &atts,
                   const string &unit,
                   double value);
    void AddResults(const string &test,
                    const string &atts,
                    const string &unit,
                    const vector<double> &values);
    void DumpDetailed(ostream&);
    void DumpSummary(ostream&);
};


#endif
