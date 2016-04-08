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

#include <algorithm>
#include <cfloat>
#include <cmath>
#include "ResultDatabase.h"

using namespace std;

bool ResultDatabase::Result::operator<(const Result &rhs) const
{
    if (test < rhs.test)
        return true;
    if (test > rhs.test)
        return false;
    if (atts < rhs.atts)
        return true;
    if (atts > rhs.atts)
        return false;
    return true; //?
}

double ResultDatabase::Result::GetMin()
{
    double r = FLT_MAX;
    for (int i=0; i<value.size(); i++)
    {
        r = min(r, value[i]);
    }
    return r;
}

double ResultDatabase::Result::GetMax()
{
    double r = -FLT_MAX;
    for (int i=0; i<value.size(); i++)
    {
        r = max(r, value[i]);
    }
    return r;
}

double ResultDatabase::Result::GetMedian()
{
    int n = value.size();
    if (n == 1)
        return value[0];
    if (n == 2)
        return (value[0] + value[1]) / 2.;

    vector<double> sorted = value;
    sort(sorted.begin(), sorted.end());

    if ((n % 2) == 1)
        return sorted[n/2];
    else
        return (sorted[n/2-1] + sorted[n/2]) / 2.;
}

double ResultDatabase::Result::GetMean()
{
    double r = 0;
    for (int i=0; i<value.size(); i++)
    {
        r += value[i];
    }
    return r / double(value.size());
}

double ResultDatabase::Result::GetStdDev()
{
    double r = 0;
    double u = GetMean();
    for (int i=0; i<value.size(); i++)
    {
        r += (value[i] - u) * (value[i] - u);
    }
    r = sqrt(r / value.size());
    return r;
}


void ResultDatabase::AddResults(const string &test,
                                const string &atts,
                                const string &unit,
                                const vector<double> &values)
{
    for (int i=0; i<values.size(); i++)
    {
        AddResult(test, atts, unit, values[i]);
    }
}

void ResultDatabase::AddResult(const string &test,
                               const string &atts,
                               const string &unit,
                               double value)
{
    int index;
    for (index = 0; index < results.size(); index++)
    {
        if (results[index].test == test &&
            results[index].atts == atts)
        {
            if (results[index].unit != unit)
                throw "Internal error: mixed units";

            break;
        }
    }

    if (index >= results.size())
    {
        Result r;
        r.test = test;
        r.atts = atts;
        r.unit = unit;
        results.push_back(r);
    }

    results[index].value.push_back(value);
}

// ****************************************************************************
//  Method:  ResultDatabase::DumpDetailed
//
//  Purpose:
//    Writes the full results, including all trials.
//
//  Arguments:
//    out        where to print
//
//  Programmer:  Jeremy Meredith
//  Creation:    August 14, 2009
//
//  Modifications:
//    Jeremy Meredith, Wed Nov 10 14:25:17 EST 2010
//    Renamed to DumpDetailed to make room for a DumpSummary.
//
//    Jeremy Meredith, Thu Nov 11 11:39:57 EST 2010
//    Added note about (*) missing value tag.
//
// ****************************************************************************
void ResultDatabase::DumpDetailed(ostream &out)
{
    vector<Result> sorted(results);

    sort(sorted.begin(), sorted.end());

    int maxtrials = 1;
    for (int i=0; i<sorted.size(); i++)
    {
        if (sorted[i].value.size() > maxtrials)
            maxtrials = sorted[i].value.size();
    }

    // TODO: in big parallel runs, the "trials" are the procs
    // and we really don't want to print them all out....
    out << "test\t"
        << "atts\t"
        << "units\t"
        << "median\t"
        << "mean\t"
        << "stddev\t"
        << "min\t"
        << "max\t";
    for (int i=0; i<maxtrials; i++)
        out << "trial"<<i<<"\t";
    out << endl;

    for (int i=0; i<sorted.size(); i++)
    {
        Result &r = sorted[i];
        out << r.test << "\t"
            << r.atts << "\t"
            << r.unit << "\t"
            << r.GetMedian() << "\t"
            << r.GetMean()   << "\t"
            << r.GetStdDev() << "\t"
            << r.GetMin()    << "\t"
            << r.GetMax()    << "\t";
        for (int j=0; j<r.value.size(); j++)        
        {
            out << r.value[j] << "\t";
        }
        
        out << endl;
    }
    out << endl
        << "Note: results marked with (*) had missing values such as" << endl
        << "might occur with a mixture of architectural capabilities." << endl;
}


// ****************************************************************************
//  Method:  ResultDatabase::DumpDetailed
//
//  Purpose:
//    Writes the summary results (min/max/stddev/med/mean), but not
//    every individual trial.
//
//  Arguments:
//    out        where to print
//
//  Programmer:  Jeremy Meredith
//  Creation:    November 10, 2010
//
//  Modifications:
//    Jeremy Meredith, Thu Nov 11 11:39:57 EST 2010
//    Added note about (*) missing value tag.
//
// ****************************************************************************
void ResultDatabase::DumpSummary(ostream &out)
{
    vector<Result> sorted(results);

    sort(sorted.begin(), sorted.end());

    // TODO: in big parallel runs, the "trials" are the procs
    // and we really don't want to print them all out....
    out << "test\t"
        << "atts\t"
        << "units\t"
        << "median\t"
        << "mean\t"
        << "stddev\t"
        << "min\t"
        << "max\t";
    out << endl;

    for (int i=0; i<sorted.size(); i++)
    {
        Result &r = sorted[i];
        out << r.test << "\t"
            << r.atts << "\t"
            << r.unit << "\t"
            << r.GetMedian() << "\t"
            << r.GetMean()   << "\t"
            << r.GetStdDev() << "\t"
            << r.GetMin()    << "\t"
            << r.GetMax()    << "\t";
        
        out << endl;
    }
    out << endl
        << "Note: results marked with (*) had missing values such as" << endl
        << "might occur with a mixture of architectural capabilities." << endl;
}

