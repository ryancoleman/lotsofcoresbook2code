/*
  This file is provided under a BSD license.

  Copyright (c) 2015 Intel Corporation. All rights reserved.
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the
      distribution.
    * Neither the name of Intel Corporation nor the names of its
      contributors may be used to endorse or promote products derived
      from this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "get_time.h"

#if defined(__linux__) || defined(__CYGWIN__)
#include <time.h>

long long get_current_time()
{
	struct timespec curr_time;
	
	clock_gettime(CLOCK_REALTIME, &curr_time);
	
	return curr_time.tv_sec*1000000000l+curr_time.tv_nsec;
}

#else
#include <Windows.h>
long long get_current_time()
{
    LARGE_INTEGER           time;
    static double           freqToNano;
    static int             bInitialized = 0;
    if ( !bInitialized ) {
        static LARGE_INTEGER    freq;
        QueryPerformanceFrequency(&freq);
        freqToNano = (double)freq.QuadPart / 1000000000.;
        bInitialized = 1;
    }

    QueryPerformanceCounter(&time);
    return (long long) ((double)time.QuadPart / freqToNano);
}

#endif
