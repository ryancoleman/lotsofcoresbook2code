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

#ifndef _PROGRESS_BAR_H_
#define _PROGRESS_BAR_H_

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>


// ****************************************************************************
// Class: ProgressBar
//
// Purpose:
//   Simple text progress bar class.
//
// Programmer: Gabriel Marin
// Creation:   October 12, 2009
//
// Modifications:
//
// ****************************************************************************
class ProgressBar
{
private:
    int itersDone;
    int totalIters;
    static const char barDone[81];
    double rTotal;
    double percDone;
    
public:
    //   Constructor
    //
    //   Arguments:
    //       _totalIters  total work amount to be tracked
    ProgressBar (int _totalIters = 0)
    {
        totalIters = _totalIters;
        itersDone = 0;
        if (totalIters)
        {
            rTotal = 100.0/totalIters;
        } else
        {
            rTotal = 0.0;
        }
        percDone = itersDone*rTotal;
    }
    
    //   Method: setTotalIters
    //
    //   Purpose: setter for the total work amount
    //
    //   Arguments:
    //       _totalIters  total work amount to be tracked
    void setTotalIters (int _totalIters)
    {
        totalIters = _totalIters;
        if (totalIters)
        {
            rTotal = 100.0/totalIters;
            percDone = itersDone*rTotal;
        }
    }
    
    //   Method: setItersDone
    //
    //   Purpose: setter for the completed work amount
    //
    //   Arguments:
    //       _itersDone  completed work amount
    void setItersDone (int _itersDone)
    {
        itersDone = _itersDone;
        percDone = itersDone*rTotal;
    }
    
    //   Method: addItersDone
    //
    //   Purpose: update amount of completed work
    //
    //   Arguments:
    //       _inc  amount of newly completed work
    void addItersDone (int _inc = 1)
    {
        itersDone += _inc;
        percDone = itersDone*rTotal;
    }
    
    //   Method: Show
    //
    //   Purpose: display progress bar
    //
    //   Arguments:
    //       fd  output file descriptor
    void Show (FILE *fd)
    {
        int lenDone = (int)(percDone/2.0 + 0.5);
        fprintf(fd, "\r|%.*s%*s| %5.1lf%%", lenDone, barDone, 50-lenDone, "", percDone);
        fflush(fd);
    }
};

#endif
